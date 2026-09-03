# Referee verification scripts

This directory contains checkers for each numerical
verification box included in the manuscript.
For all but verification 4.16 the checkers are python programs that can be found in the root directory.
The Finite arithmetic verification 4.16 also uses some C++ parts and is located in a subfolder. To check 4.16, please run `check_4_16/check_4_16.py`. It is expected to run for around 10 minutes.

## Quick start

Python 3.11 or later is recommended.  From this directory, run:

```console
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 run_all.py
```

Each checker can instead be run on its own, for example:

```console
python3 check_5_5.py
```

No checker accepts or needs numerical input from the command line.  A successful
checker prints the quantities it reconstructed, the claimed threshold, useful
worst-case diagnostics, and a final `PASS`; a failed or inconclusive check exits
with a nonzero status.

On a representative Linux run, the complete suite took about six minutes (excluding 4.16).
Almost all of that time was the exhaustive 2,150,400-box replay in
`check_5_9.py`; the other five scripts together took about six seconds.

## Included checks

| Script | Manuscript box | What is reconstructed and checked |
|---|---:|---|
| `check_4_12.py` | 4.12 | The exact Bernoulli numerators \(A_g\), for \(1\leq g\leq2000\), and the bound \(\lvert A_g\rvert<200000^{782}\). |
| `check_4_14.py` | 4.14 | Deterministically constructs the prime lists \(P_g\), then independently rechecks primality, nondivisibility by \(A_g\), and complete interval coverage. |
| `check_5_5.py` | 5.5 | The finite error sum, the two geometric tail bounds, the two boundary values, and their combined \(<1/30\) estimate. |
| `check_5_9.py` | 5.9 | Every one of the 2,150,400 Case-A product boxes and the two compact Case-B boxes in the leading-term argument. |
| `check_6_2.py` | 6.2 | Reconstructs \(\mathscr L_q\), \(\ell_q\), and \(\widehat H_q\) through \(w^{14}\), then checks all three displayed norm bounds. |
| `check_6_6.py` | 6.6 | Reconstructs \(\log R\) through \(u^{10}\) modulo \(w^{15}\), recovers every displayed \(\lambda_m\), and checks the three scalar bounds. |

## Arithmetic and inspectability

The scripts do not load old certificate outputs or trust tables of final
answers.

- Integer, rational, Bernoulli, and formal-polynomial calculations use Python
  arbitrary-precision integers, `fractions.Fraction`, or exact FLINT arithmetic.
- `check_5_5.py` uses Arb real balls through `python-flint`; a strict comparison
  is accepted only when the complete certified ball lies below the threshold.
- `check_5_9.py` contains its interval operations and elementary-function
  enclosures in the same file.  Binary64 operations are rounded outward, while
  \(\pi\), logarithms, and exponentials are enclosed using finite series with
  explicit remainder bounds.  NumPy is used only to evaluate many independent
  boxes in parallel arrays.
- `check_6_2.py` constructs a rational interval for \(\pi\) from Machin's
  formula and alternating-series remainders.
- `check_6_6.py` reduces its final comparisons to rational inequalities using
  the elementary bounds \(\pi>3\) and \(\sqrt\pi>7/4\), as in the manuscript.

The dependency versions are pinned in `requirements.txt`.  The individual
docstrings and comments give the manuscript equation numbers and explain the
finite recurrences and interval conventions used in each check.
