# Referee verification scripts

This directory contains checkers for each numerical
verification box included in the manuscript.
For all but verification `finite_cover` the checkers are python programs that can be found in the root directory.
The Finite arithmetic verification finite_cover also uses some C++ parts and is located in a subfolder. To check it, please run `finite_cover/finite_cover.py`. It is expected to run for around 10 minutes.

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
python3 error_finite.py
```

No checker accepts or needs numerical input from the command line.  A successful
checker prints the quantities it reconstructed, the claimed threshold, useful
worst-case diagnostics, and a final `PASS`; a failed or inconclusive check exits
with a nonzero status.

On a representative Linux run, the complete suite took about six minutes (excluding `finite_cover`).
Almost all of that time was the exhaustive 2,150,400-box replay in
`interval_check_m.py`; the other five scripts together took about six seconds.

## Arithmetic and inspectability

The scripts do not load old certificate outputs or trust tables of final
answers.

- Integer, rational, Bernoulli, and formal-polynomial calculations use Python
  arbitrary-precision integers, `fractions.Fraction`, or exact FLINT arithmetic.
- `error_finite.py` uses Arb real balls through `python-flint`; a strict comparison
  is accepted only when the complete certified ball lies below the threshold.
- `interval_check_m.py` contains its interval operations and elementary-function
  enclosures in the same file.  Binary64 operations are rounded outward, while
  \(\pi\), logarithms, and exponentials are enclosed using finite series with
  explicit remainder bounds.  NumPy is used only to evaluate many independent
  boxes in parallel arrays.
- `llh_bounds.py` constructs a rational interval for \(\pi\) from Machin's
  formula and alternating-series remainders.
- `lambda_t_bounds.py` reduces its final comparisons to rational inequalities using
  the elementary bounds \(\pi>3\) and \(\sqrt\pi>7/4\).

The dependency versions are pinned in `requirements.txt`.  The individual
docstrings and comments give the manuscript equation numbers and explain the
finite recurrences and interval conventions used in each check.
