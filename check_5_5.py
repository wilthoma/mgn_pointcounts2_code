#!/usr/bin/env python3
"""Certify the three numerical inequalities in Verification 5.5.

The manuscript first bounds the error term E(g,n) by a finite sum, an
infinite tail, and a boundary term.  Equations (5.19)--(5.21) reduce the
last two quantities to finite expressions.  This script evaluates exactly
those three expressions.

All non-rational quantities are evaluated with Arb ball arithmetic through
python-flint.  An Arb ball is a mathematically certified interval, with
outward rounding.  Thus ``value < threshold`` is accepted below only when
the entire computed ball lies strictly below the threshold ball.

The script takes no command-line arguments.
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction

try:
    from flint import arb, ctx
except ImportError as exc:  # Give a useful message if requirements were skipped.
    raise SystemExit(
        "This checker requires python-flint.  Install the pinned dependencies "
        "with:  python3 -m pip install -r requirements.txt"
    ) from exc


# A comfortable precision for strict inequalities whose margins are much larger
# than the radii of the resulting balls.  Arb measures precision in bits.
WORKING_PRECISION_BITS = 256
DISPLAY_DIGITS = 32


def show_ball(value: arb) -> str:
    """Return a readable midpoint-radius enclosure."""

    return value.str(DISPLAY_DIGITS, radius=True)


def certify_less(
    description: str, value: arb, threshold: arb, threshold_text: str
) -> bool:
    """Print and certify a strict inequality between two Arb balls."""

    print(f"\n{description}")
    print(f"  certified enclosure : {show_ball(value)}")
    print(f"  required threshold  : {threshold_text}")

    if value < threshold:
        print("  result              : PASS (the whole enclosure is below the threshold)")
        return True

    # Arb comparisons are deliberately conservative: overlapping balls are
    # neither < nor >=.  Distinguish that situation from a certified failure.
    if value >= threshold:
        print("  result              : FAIL (the enclosure is at or above the threshold)")
    else:
        print("  result              : INCONCLUSIVE (increase the working precision)")
    return False


def finite_part(pi: arb) -> arb:
    r"""Return the finite expression on the left of (5.22).

    The notation (1997)_{m-1} in the manuscript denotes the falling
    factorial

        (1997)_{m-1} = 1997 * 1996 * ... * (1999-m).

    We update this exact integer incrementally.  Only pi and Gamma require
    ball arithmetic; all integers in the formula remain exact.
    """

    two_pi = 2 * pi
    total = (
        arb(600)
        * (arb(55) / 12)
        * two_pi**2
        / (arb(1998) * 1997)
    )

    # For m=2 the falling factorial has one factor, namely 1997.  Before
    # evaluating the m-th summand we append the factor 1999-m, so that at
    # m=3 the exact product is 1997*1996.
    falling_factorial = 1997
    for m in range(3, 1998):
        falling_factorial *= 1999 - m
        gamma_factor = (arb(m) / 2 + 1).gamma()
        numerator = arb(600) * arb(6) ** m * gamma_factor * two_pi**m
        denominator = arb(1998) * falling_factorial
        total += numerator / denominator

    return total


def geometric_tail_bound(pi: arb) -> arb:
    r"""Return the two geometric-series bounds displayed in (5.23)."""

    pi_squared = pi**2
    two_pi = 2 * pi

    even_denominator_factor = 1 - arb(72) * pi_squared / 2001
    odd_denominator_factor = 1 - arb(72) * pi_squared / 2002

    # Positivity is part of what makes the geometric bounds meaningful.
    if not (even_denominator_factor > 0 and odd_denominator_factor > 0):
        raise ArithmeticError("Arb could not certify the geometric denominators positive")

    even_part = (
        arb(600)
        * arb(6) ** 1998
        * arb(1000).gamma()
        * two_pi**1998
        / (arb(math.factorial(1999)) * even_denominator_factor)
    )
    odd_part = (
        arb(600)
        * arb(6) ** 1999
        * (arb(2001) / 2).gamma()
        * two_pi**1999
        / (arb(math.factorial(2000)) * odd_denominator_factor)
    )
    return even_part + odd_part


def boundary_expression(h: int, pi: arb) -> arb:
    r"""Evaluate the last expression in (5.18) at the integer h."""

    first = arb(6) ** (h - 2) * (arb(h) / 2).gamma()
    second = 2 * arb(6) ** (h - 1) * (arb(h + 1) / 2).gamma()
    return (first + second) * (2 * pi) ** h / arb(math.factorial(h - 2))


def main() -> int:
    if len(sys.argv) != 1:
        print("This checker takes no command-line arguments.", file=sys.stderr)
        return 2

    ctx.prec = WORKING_PRECISION_BITS
    pi = arb.pi()

    print("Certified analytic verification 5.5: finite error-bound evaluation")
    print(f"Arithmetic: python-flint/Arb at {WORKING_PRECISION_BITS} bits")
    print("Every displayed enclosure includes all rounding and transcendental error.")

    print("\nComputing the finite part in (5.22):")
    print("  600*(55/12)*(2*pi)^2/(1998*1997)")
    print("  + sum_{m=3}^{1997} 600*6^m*Gamma(m/2+1)*(2*pi)^m")
    print("      / (1998*(1997)_{m-1}), with a falling factorial.")
    finite = finite_part(pi)
    finite_ok = certify_less(
        "Claim 1: finite part < 33/1000",
        finite,
        arb(33) / 1000,
        "33/1000 = 0.033",
    )

    print("\nComputing the two parity-tail bounds in (5.23), obtained from")
    print("the two-step ratio 72*pi^2/(m+3).")
    tail = geometric_tail_bound(pi)
    tail_ok = certify_less(
        "Claim 2: infinite-tail bound < 10^(-9)",
        tail,
        arb(1) / 10**9,
        "1/10^9 = 0.000000001",
    )

    print("\nComputing the boundary expression from the last line of (5.18):")
    print("  [6^(h-2)*Gamma(h/2) + 2*6^(h-1)*Gamma((h+1)/2)]")
    print("  * (2*pi)^h / (h-2)!")
    boundary_2000 = boundary_expression(2000, pi)
    boundary_2001 = boundary_expression(2001, pi)
    print(f"  h = 2000 enclosure: {show_ball(boundary_2000)}")
    print(f"  h = 2001 enclosure: {show_ball(boundary_2001)}")

    if boundary_2000 > boundary_2001:
        boundary_max = boundary_2000
        maximizing_h = 2000
    elif boundary_2001 > boundary_2000:
        boundary_max = boundary_2001
        maximizing_h = 2001
    else:
        print("  Could not order the two boundary balls; increase the precision.")
        return 1
    print(f"  certified larger value occurs at h = {maximizing_h}")

    boundary_ok = certify_less(
        "Claim 3: max(boundary(2000), boundary(2001)) < 377/10^13",
        boundary_max,
        arb(377) / 10**13,
        "377/10^13 = 0.0000000000377",
    )

    # The last implication in the box is exact rational arithmetic.  Keeping
    # it separate makes clear that no numerical approximation is involved.
    published_upper_bound = (
        Fraction(33, 1000) + Fraction(1, 10**9) + Fraction(377, 10**13)
    )
    target = Fraction(1, 30)
    combined_bounds_ok = published_upper_bound < target

    print("\nCombining the three published rational upper bounds exactly:")
    print("  33/1000 + 1/10^9 + 377/10^13")
    print(f"  exact value          : {published_upper_bound}")
    print(f"  required threshold   : 1/30")
    print(
        "  result               : "
        + ("PASS" if combined_bounds_ok else "FAIL")
        + " (exact integer comparison)"
    )

    # This extra direct check is redundant, but useful to a referee: it shows
    # the enclosure obtained before replacing the components by round bounds.
    actual_total = finite + tail + boundary_max
    actual_total_ok = certify_less(
        "Cross-check: sum of the three computed enclosures < 1/30",
        actual_total,
        arb(1) / 30,
        "1/30",
    )

    all_ok = all(
        (finite_ok, tail_ok, boundary_ok, combined_bounds_ok, actual_total_ok)
    )
    print("\n" + "=" * 72)
    if all_ok:
        print("PASS: all statements in Certified analytic verification 5.5 hold.")
        return 0

    print("FAIL OR INCONCLUSIVE: Verification 5.5 was not certified.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
