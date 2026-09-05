#!/usr/bin/env python3
"""Recheck Finite arithmetic verification 4.13.

The manuscript defines, for 3 <= g <= 2000,

    h = floor((g - 1)/2),
    (-1)^(g+1) B_(2h)/(2h) = A_g/C_g

in lowest terms with C_g > 0; it sets A_1 = A_2 = 1.  This script
constructs those rational numbers with FLINT exact arithmetic and verifies

    |A_g| < 200000^782                  (1 <= g <= 2000).

There are no floating-point calculations in the verification.  In
particular, ``fmpq.bernoulli`` returns the exact Bernoulli number and FLINT
keeps every ``fmpq`` in lowest terms with positive denominator.

Run from the repository directory with

    python3 check_4_12.py
"""

from __future__ import annotations

import sys
from time import perf_counter

try:
    from flint import fmpq, fmpz
except ImportError:
    print("ERROR: python-flint is required; install requirements.txt first.")
    raise SystemExit(2)


FIRST_GENUS = 1
LAST_GENUS = 2000
BOUND_BASE = 200000
BOUND_EXPONENT = 782


def numerator_A(g: int) -> fmpz:
    """Return the integer A_g from equation (4.21), using exact arithmetic."""

    if g in (1, 2):
        return fmpz(1)

    h = (g - 1) // 2
    value = fmpq.bernoulli(2 * h) / (2 * h)
    if (g + 1) % 2 == 1:  # multiply by (-1)^(g+1)
        value = -value

    # An fmpq is automatically reduced and has positive denominator, so its
    # numerator is precisely A_g (rather than merely a multiple of A_g).
    assert value.denominator > 0
    return value.numerator


def main() -> int:
    if len(sys.argv) != 1:
        print("This checker takes no command-line arguments.", file=sys.stderr)
        return 2

    started = perf_counter()
    bound = fmpz(BOUND_BASE) ** BOUND_EXPONENT

    print("Finite arithmetic verification 4.12")
    print("=" * 38)
    print("Reconstructing A_g from the exact Bernoulli numbers in (4.21).")
    print(f"Genus range: {FIRST_GENUS} <= g <= {LAST_GENUS}")
    print(f"Claim: |A_g| < {BOUND_BASE}^{BOUND_EXPONENT}")
    print("Arithmetic: exact FLINT integers and rational numbers (no floats)\n")

    largest_abs_A = fmpz(0)
    largest_genera: list[int] = []

    for g in range(FIRST_GENUS, LAST_GENUS + 1):
        abs_A = abs(numerator_A(g))
        if abs_A >= bound:
            print("FAIL")
            print(f"At g = {g}, the strict inequality does not hold.")
            print(f"|A_g| has {len(str(abs_A))} decimal digits.")
            return 1

        if abs_A > largest_abs_A:
            largest_abs_A = abs_A
            largest_genera = [g]
        elif abs_A == largest_abs_A:
            largest_genera.append(g)

    checked = LAST_GENUS - FIRST_GENUS + 1
    bound_digits = len(str(bound))
    numerator_digits = len(str(largest_abs_A))
    exact_quotient = bound // largest_abs_A
    elapsed = perf_counter() - started

    print(f"PASS: all {checked} values satisfy the strict inequality.")
    print("Worst case (largest |A_g|):")
    print(f"  g = {', '.join(map(str, largest_genera))}")
    print(f"  decimal digits in |A_g|: {numerator_digits}")
    print(f"  decimal digits in the bound: {bound_digits}")
    print(f"  exact floor(bound / |A_g|): {exact_quotient}")
    print(f"Elapsed time: {elapsed:.3f} seconds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
