#!/usr/bin/env python3
"""Verify the three finite computations in Certified verification 6.2.

The script reconstructs every polynomial from the definitions in the paper.  It
does not read a table of precomputed norms.

All polynomial coefficients are computed with ``fractions.Fraction``.  Powers
of pi are bounded using a rational interval obtained from Machin's formula

    pi = 16 arctan(1/5) - 4 arctan(1/239),

with rigorous alternating-series remainders.  Thus every printed norm bound is
an upper bound, rather than a floating-point approximation.

Run from any directory with

    python3 check_6_2.py

No third-party Python package is required.
"""

from __future__ import annotations

from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from math import comb, factorial
import sys


# The norm in the manuscript keeps the coefficients of w^0,...,w^14.
W_DEGREE = 14
MAX_Q = 99
FACT = [factorial(n) for n in range(MAX_Q + 1)]

Polynomial = list[Fraction]


def zero_polynomial() -> Polynomial:
    """Return a zero polynomial represented through degree W_DEGREE."""

    return [Fraction(0) for _ in range(W_DEGREE + 1)]


def polynomial_add(left: Polynomial, right: Polynomial) -> Polynomial:
    """Add two truncated polynomials exactly."""

    return [a + b for a, b in zip(left, right)]


def polynomial_scale(polynomial: Polynomial, scalar: Fraction) -> Polynomial:
    """Multiply a polynomial by an exact rational scalar."""

    return [scalar * coefficient for coefficient in polynomial]


def polynomial_multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    """Multiply exactly, discarding powers w^15 and higher."""

    product = zero_polynomial()
    for left_degree, left_coefficient in enumerate(left):
        if left_coefficient == 0:
            continue
        for right_degree in range(W_DEGREE + 1 - left_degree):
            right_coefficient = right[right_degree]
            if right_coefficient:
                product[left_degree + right_degree] += (
                    left_coefficient * right_coefficient
                )
    return product


def rational_l1_norm(polynomial: Polynomial) -> Fraction:
    """Return the exact sum of the absolute rational coefficients."""

    return sum((abs(coefficient) for coefficient in polynomial), Fraction(0))


def bernoulli_numbers(last_index: int) -> list[Fraction]:
    """Compute B_0,...,B_last_index in the convention B_1=-1/2.

    This is the defining triangular recurrence obtained from

        sum_{k=0}^n binom(n+1,k) B_k = 0  (n >= 1).
    """

    numbers = [Fraction(1)]
    for n in range(1, last_index + 1):
        preceding_sum = sum(
            (Fraction(comb(n + 1, k)) * numbers[k] for k in range(n)),
            Fraction(0),
        )
        numbers.append(-preceding_sum / (n + 1))
    return numbers


BERNOULLI = bernoulli_numbers(MAX_Q + 1)


def bernoulli_at_one(index: int) -> Fraction:
    """Return B_index(1); only index 1 differs from B_index."""

    return Fraction(1, 2) if index == 1 else BERNOULLI[index]


def d_polynomial(q: int) -> Polynomial:
    r"""Construct d_q(w), truncated after w^14, from formula (6.6).

    Taylor expansion at 1 gives

      [w^j] B_{q+1}(1-w)
          = (-1)^j binom(q+1,j) B_{q+1-j}(1).

    Since q>=1, the constant term cancels in

      d_q(w) = -(B_{q+1}(1-w)-B_{q+1})/(q(q+1)).
    """

    degree = q + 1
    polynomial = zero_polynomial()
    for j in range(1, min(degree, W_DEGREE) + 1):
        polynomial[j] = (
            -Fraction((-1) ** j * comb(degree, j), q * degree)
            * bernoulli_at_one(degree - j)
        )
    return polynomial


def normalized_ell_core(q: int) -> Polynomial:
    r"""Return E_q with ell_q(w) = pi^(q+1) E_q(w).

    In the paper ell_q=d_q/S_{q+1}, where

      S_{q+1} = (q-1)!/(2*pi)^(q+1).

    Therefore E_q = 2^(q+1) d_q/(q-1)! has rational coefficients.
    """

    scalar = Fraction(2 ** (q + 1), FACT[q - 1])
    return polynomial_scale(d_polynomial(q), scalar)


def normalized_h_cores(ell_cores: dict[int, Polynomial]) -> dict[int, Polynomial]:
    r"""Reconstruct A_q with Hhat_q(w) = pi^(q+1) A_q(w).

    After taking the common power pi^(q+1) out of recurrence (6.7),

      A_q = E_q + (1/2) sum_{r=1}^{q-1} c_{q,r} E_{q-r} A_r,

      c_{q,r} = (r-1)! (q-r)! / q!.

    Products are truncated after w^14, exactly as in the stated norm.
    """

    h_cores: dict[int, Polynomial] = {}
    for q in range(1, MAX_Q + 1):
        h_q = ell_cores[q]
        for r in range(1, q):
            half_c_qr = Fraction(FACT[r - 1] * FACT[q - r], 2 * FACT[q])
            term = polynomial_multiply(ell_cores[q - r], h_cores[r])
            h_q = polynomial_add(h_q, polynomial_scale(term, half_c_qr))
        h_cores[q] = h_q
    return h_cores


def arctangent_bounds(reciprocal: int, terms: int) -> tuple[Fraction, Fraction]:
    r"""Enclose arctan(1/reciprocal) by an exact rational interval.

    The series sum (-1)^j x^(2j+1)/(2j+1) is alternating with decreasing
    terms for 0<x<=1.  Its value lies strictly between a partial sum and
    that partial sum plus the first omitted term.
    """

    partial_sum = sum(
        (
            Fraction((-1) ** j, (2 * j + 1) * reciprocal ** (2 * j + 1))
            for j in range(terms)
        ),
        Fraction(0),
    )
    first_omitted = Fraction(
        (-1) ** terms,
        (2 * terms + 1) * reciprocal ** (2 * terms + 1),
    )
    other_endpoint = partial_sum + first_omitted
    return min(partial_sum, other_endpoint), max(partial_sum, other_endpoint)


def certified_pi_interval() -> tuple[Fraction, Fraction]:
    """Return rational lower and upper bounds for pi via Machin's formula."""

    atan_fifth = arctangent_bounds(5, 40)
    atan_239th = arctangent_bounds(239, 20)
    lower = 16 * atan_fifth[0] - 4 * atan_239th[1]
    upper = 16 * atan_fifth[1] - 4 * atan_239th[0]
    assert lower < upper
    return lower, upper


def decimal_bound(value: Fraction, places: int, upper: bool) -> str:
    """Print a rational bound outward-rounded to a fixed number of places."""

    rounding = ROUND_CEILING if upper else ROUND_FLOOR
    quantum = Decimal(1).scaleb(-places)
    with localcontext() as context:
        context.prec = places + 40
        context.rounding = rounding
        decimal_value = Decimal(value.numerator) / Decimal(value.denominator)
        decimal_value = decimal_value.quantize(quantum, rounding=rounding)
    return f"{decimal_value:.{places}f}"


def l_polynomial(sigma: int) -> list[tuple[Fraction, int]]:
    r"""Construct L_sigma from (5.30).

    Each returned pair (c,k) represents the coefficient c*pi^k of w^k.
    Simplifying the powers of i in (5.30), only odd k survive for sigma=0
    and only even k survive for sigma=1.  In either case the surviving
    rational coefficient has absolute value 2^(k+1)/k!; its sign is kept
    here even though the l1 norm subsequently takes absolute values.
    """

    coefficients: list[tuple[Fraction, int]] = []
    wanted_parity = (sigma + 1) % 2
    for k in range(1, W_DEGREE + 1):
        if k % 2 != wanted_parity:
            coefficients.append((Fraction(0), k))
            continue
        sign_exponent = (k + 1) // 2 if sigma == 0 else k // 2 + 1
        coefficient = Fraction(((-1) ** sign_exponent) * 2 ** (k + 1), FACT[k])
        coefficients.append((coefficient, k))
    return coefficients


def l_norm_upper(sigma: int, pi_upper: Fraction) -> Fraction:
    """Return a rigorous upper bound for ||L_sigma||_{1,<=14}."""

    return sum(
        (abs(coefficient) * pi_upper**power for coefficient, power in l_polynomial(sigma)),
        Fraction(0),
    )


def common_pi_power_norm_upper(
    polynomial: Polynomial, power: int, pi_upper: Fraction
) -> Fraction:
    """Bound ||pi^power * polynomial|| using the certified upper bound."""

    return rational_l1_norm(polynomial) * pi_upper**power


def report_check(label: str, upper_bound: Fraction, threshold: int, q: int) -> bool:
    """Print one maximum and return whether it is strictly below threshold."""

    passed = upper_bound < threshold
    relation = "<" if passed else ">="
    print(
        f"  largest certified upper bound: {decimal_bound(upper_bound, 12, True)} "
        f"(at q = {q})"
    )
    print(f"  required comparison: {relation} {threshold}")
    print(f"  {label}: {'PASS' if passed else 'FAIL'}")
    return passed


def main() -> int:
    """Run all three finite checks in the box and return a shell status."""

    if len(sys.argv) != 1:
        print("This certificate takes no command-line parameters.", file=sys.stderr)
        return 2

    print("Certified analytic verification 6.2")
    print("===================================")
    print("All polynomial arithmetic is exact rational arithmetic.")

    pi_lower, pi_upper = certified_pi_interval()
    print("pi is enclosed using Machin's formula and alternating-series bounds:")
    print(f"  {decimal_bound(pi_lower, 50, False)} < pi")
    print(f"  pi < {decimal_bound(pi_upper, 50, True)}")

    all_passed = True

    print("\n(6.15) Two-mode polynomials L_q, reconstructed from (5.30)")
    l_bounds: list[tuple[Fraction, int]] = []
    for q in (0, 1):
        bound = l_norm_upper(q, pi_upper)
        l_bounds.append((bound, q))
        print(
            f"  q = {q}: ||L_{q}||_(1,<=14) "
            f"< {decimal_bound(bound, 12, True)}"
        )
    largest_l_bound, largest_l_q = max(l_bounds)
    all_passed &= report_check("claim (6.15)", largest_l_bound, 534, largest_l_q)

    print("\nReconstructing d_q, ell_q, and Hhat_q through degree w^14 ...")
    ell_cores = {q: normalized_ell_core(q) for q in range(1, MAX_Q + 1)}
    h_cores = normalized_h_cores(ell_cores)

    # The defining formulas imply zero constant terms; checking these is a
    # useful guard against a Bernoulli-convention or indexing mistake.
    assert all(ell_cores[q][0] == 0 for q in ell_cores)
    assert all(h_cores[q][0] == 0 for q in h_cores)

    print("\n(6.16) Normalized Bernoulli polynomials ell_q for 1 <= q <= 50")
    ell_bounds = [
        (common_pi_power_norm_upper(ell_cores[q], q + 1, pi_upper), q)
        for q in range(1, 51)
    ]
    largest_ell_bound, largest_ell_q = max(ell_bounds)
    all_passed &= report_check(
        "claim (6.16)", largest_ell_bound, 800, largest_ell_q
    )

    print("\n(6.17) Normalized factors Hhat_q for 1 <= q <= 99")
    h_bounds = [
        (common_pi_power_norm_upper(h_cores[q], q + 1, pi_upper), q)
        for q in range(1, 100)
    ]
    largest_h_bound, largest_h_q = max(h_bounds)
    all_passed &= report_check("claim (6.17)", largest_h_bound, 600, largest_h_q)

    print("\nOverall result:", "PASS" if all_passed else "FAIL")
    return 0 if all_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
