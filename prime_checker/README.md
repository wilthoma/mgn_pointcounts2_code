# High-genus modular checker: planning package

This package scales the modular low-truncation calculation from the audited
`g <= 600` run to the proposed ranges `g <= 1000` and `g <= 2000`.

It has **not** been used in this sandbox to complete either large audit.  Only
a short compatibility test was run: the large-NTT generator exactly matches
the archived generator at `G=20, W=2, S=30, p=101`.

## Mathematical parameters

For weight `(11,0)` put `Gamma=10`; for weight `(15,0)` put `Gamma=14`.
If `n=q p+s`, then the modular translation uses `D_{Gamma-q}(g,s)`.
With the safe choice `p>G+2` and `n>=5000`, the smallest candidate primes give

| G | first p | floor(5000/p) | W for Gamma=10 | W for Gamma=14 |
|---:|---:|---:|---:|---:|
| 1000 | 1009 | 4 | 6 | 10 |
| 2000 | 2003 | 2 | 8 | 12 |

For direct evaluation below `n=5000`, choose `p>5000` and use `W=Gamma`.

## Build

```bash
g++ -O3 -std=c++17 -march=native \
  modular_low_truncation_table_large.cpp \
  -o modular_low_truncation_table_large

g++ -O3 -std=c++17 -march=native \
  check_modular_cover.cpp -o check_modular_cover
```

The generator uses the NTT primes `2013265921` and `469762049`, supporting
transform lengths through `2^26`.  The program checks the CRT coefficient
bound before each convolution family.

## Resource planning

```bash
python3 estimate_resources.py 2000 8 5002 --jobs 8
python3 estimate_resources.py 2000 12 5002 --jobs 8
python3 estimate_resources.py 2000 14 5002 --jobs 4
```

The estimates are calibrated against the audited `G=600,W=6,S=1600` run.
Benchmark one table on the target machine before launching a pool; NTT jobs are
memory-bandwidth sensitive.

## Pilot table pools for `5000 <= n <= 50000`

```bash
python3 make_jobs.py --G 2000 --Gamma 10 > jobs_a.txt
python3 make_jobs.py --G 2000 --Gamma 14 > jobs_b.txt
```

The printed pools use primes near `2000` and `5000`; their q-bands cover the
whole n-axis interval.  They are candidate pools, not an already verified
coefficient cover.  Run them in parallel, for example:

```bash
grep '^./' jobs_a.txt | parallel -j 8
```

Then check exact coverage:

```bash
./check_modular_cover 2000 10 5000 49999 uncovered_a.tsv tables/G2000_W8_p*.bin
./check_modular_cover 2000 14 5000 49999 uncovered_b.tsv tables/G2000_W12_p*.bin
```

If pairs remain, generate tables at nearby primes and rerun the checker.  The
checker reports the unresolved pairs and how many pairs each prime first
certifies.

## Direct finite rectangle below 5000

This is not ten million separate symbolic calculations.  One dense convolution
table produces all `G*5000` residues simultaneously.

```bash
python3 make_jobs.py --G 2000 --Gamma 10 --mode n-lt-5000 > direct_a.txt
python3 make_jobs.py --G 2000 --Gamma 14 --mode n-lt-5000 > direct_b.txt
```

After generating the tables:

```bash
./check_modular_cover 2000 10 11 4999 unresolved_a_small.tsv tables/G2000_W10_p*.bin
./check_modular_cover 2000 14 15 4999 unresolved_b_small.tsv tables/G2000_W14_p*.bin
```

A nonzero residue modulo any table prime proves exact nonvanishing.  Pairs that
remain zero modulo every chosen prime require another prime or an exact check;
they may also be genuine zero coefficients.

## Universal D1 tail

The exact tail checker requires SymPy for exact Bernoulli numerators:

```bash
python3 check_tail_prime_windows.py \
  --G 1000 --Gamma 10 --tail-start 50000 --scan-end 1000000 \
  --infinite-start 1000000 --workers 32

python3 check_tail_prime_windows.py \
  --G 2000 --Gamma 10 --tail-start 50000 --scan-end 2000000 \
  --infinite-start 2000000 --workers 32

python3 check_tail_prime_windows.py \
  --G 2000 --Gamma 14 --tail-start 60000 --scan-end 4000000 \
  --infinite-start 4000000 --workers 32
```

The finite part merges the admissible prime windows exactly.  The infinite
part uses exact rational logarithm enclosures and Dusart's explicit bounds for
`pi(x)`.  These full parameter runs were not executed in the sandbox.

## Parallel execution

The easiest parallelization is one single-thread generator process per
verification prime.  On a 64-core, 1 TB machine, start with 8 concurrent jobs,
then test 16.  Running all 64 simultaneously is likely to be limited by memory
bandwidth rather than RAM capacity.

## Simplest unified run when the range below 5000 is also needed

For `G=2000`, tables with `p` just above `5000` and `W=14` cover both families
without a separate low-n stage.  For such a table:

- `Gamma=10` is available for every `0 <= n < 10p` (about 50,000);
- `Gamma=14` is available for every `0 <= n < 14p` (about 70,000).

Thus the four `W=14` jobs printed by

```bash
python3 make_jobs.py --G 2000 --Gamma 14 --mode n-lt-5000
```

can be checked on the larger ranges

```bash
./check_modular_cover 2000 10 11 49999 unresolved_a.tsv tables/G2000_W14_p*.bin
./check_modular_cover 2000 14 15 59999 unresolved_b.tsv tables/G2000_W14_p*.bin
```

and then joined directly to the universal `D_1` tails.  This is computationally
slightly heavier per prime than the `W=8`/`W=12` split, but it gives the
simplest audit and manuscript description.
