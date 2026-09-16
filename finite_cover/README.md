# Direct large-prime check for Verification 4.16

This directory contains a checker for Finite arithmetic verification 4.16.

The program works modulo

```text
p = 2013265921 = 15 * 2^27 + 1
```

with primitive root `31`.  Since `p > 50000`, every `0 <= n < 50000`
has `q=0` in the modular translation formula.  The program therefore checks
`D_10(g,n)` directly.

The full run tests precisely

```text
1 <= g <= 2000,
0 <= n < 50000,
3g + 2n >= 25,
(g,n) not in {(8,1),(12,0)}.
```

There are `99,999,950` such pairs.  The program treats this number as a hard
consistency check.

## Quick start

The convenience runner needs Python 3 and a C++17 compiler.  It takes no
arguments:

```bash
python3 check_4_16.py
```

It compiles the C++ source, runs the built-in tests, starts the full
calculation, displays progress immediately, and saves the complete transcript
as `check_4_16_large_prime.log`.

To compile and run directly:

```bash
c++ -O3 -DNDEBUG -std=c++17 \
  check_4_16_large_prime.cpp -o check_4_16_large_prime

./check_4_16_large_prime --self-test
./check_4_16_large_prime --smoke
./check_4_16_large_prime
```

The last command is the full computation.  The `--smoke` calculation uses
`g <= 20`, `n < 64`, degree `10`, and the smaller NTT prime `40961`; it is only
an implementation test, not part of Verification 4.16.

On a machine where the binary will only be used locally, adding
`-march=native` may improve performance.

## What is computed

The archived code first constructs the coefficients `R_{m,j}(w)` and
`H_h(w)`.  If `A_a(h,n)` denotes the coefficient of `w^a` in

```text
H_h(w) (h-1+w)_n / n!,
```

then the requested coefficient is the bivariate convolution

```text
D_10 = sum_{a+b <= 10} A_a * R_b.
```

The program performs this convolution directly in `F_p`.  Because `p` is
itself an NTT prime, the two auxiliary transforms and CRT reconstruction of
the old arbitrary-modulus program are unnecessary.

The `n`-range is divided into 13 blocks of width 4096.  Consecutive input
blocks overlap by 2000 values.  This is exact, rather than an approximation:
the second index of `R_{m,j}` satisfies `j <= m <= 1999`, so no coefficient
outside that overlap can contribute to the requested output block.

For each block, the program prints:

- the output and overlapping input ranges;
- the time needed to construct the `A` coefficients;
- progress after every forward transform;
- the block time and the number of zero residues;
- a running estimate for the remaining block time.

The construction of `R` also reports progress.  Thus the transcript gives a
useful runtime estimate before the full calculation finishes.

## Memory and output

For the production parameters, the main persistent array consists of eleven
transformed cumulative `R` series and occupies about 1.38 GiB.  The largest
compact `A` block is about 0.50 GiB, and the two transform work buffers use
about 0.25 GiB.  The `R` construction has separate temporary arrays.  A few
gigabytes of available RAM should therefore be sufficient, but the program
prints its principal memory estimates before starting.

All pairs with residue zero are written to

```text
check_4_16_large_prime_unresolved.tsv
```

If the final number of zero residues is zero, the output ends with

```text
RESULT check=4.16 modulus=2013265921 status=PASS unresolved=0
```

If zero residues remain, the result is labelled `INCOMPLETE`; a second prime
would then be required.  This is not treated as a program failure, since the
purpose of the first production run is precisely to determine whether one
large prime suffices.

## Validation performed during development

The built-in tests verify primality, the stated primitive roots, field
inversion, and NTT convolution.

In addition, the complete blocked computation was run at

```text
G=20, W=10, 0<=n<64, p=40961.
```

Its `D_10` checksum and its complete zero set agreed with the archived
`modular_low_truncation_table_large.cpp`.  In `n`-major, then `g`-major order,
both computations gave

```text
checksum = 570606223
zero residues = 50.
```

This comparison tests the construction of `R` and `H`, the rising-factorial
recurrence, the overlapping-block indexing, and the final cumulative
degree-ten convolution against the earlier implementation.
