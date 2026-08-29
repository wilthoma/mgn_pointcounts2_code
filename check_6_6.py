#!/usr/bin/env python3
"""Recheck Certified analytic verification 6.6 (finite data for R).

The calculation is exact: every formal-series coefficient is a
``fractions.Fraction``.  We reconstruct

    log R(u,z,w) = sum_{m >= 1} Lambda_m(z,w) u^m  (mod w^15)

directly from equations (3.3), (3.11), (6.18), and (6.19) of the paper.
No precomputed Lambda polynomial or lambda value is used as input.

The final three comparisons contain powers of pi and, for odd indices,
sqrt(pi).  We replace them by the elementary strict lower bounds

    pi > 3,       sqrt(pi) > 7/4,

so those checks too reduce to exact rational arithmetic.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial
import sys
from typing import Dict, Iterable, Tuple


# We need u^1,...,u^10 and w^0,...,w^14.
U_MAX = 10
W_MAX = 14

# In the second sum in (3.3), a term with k has u-order at least k-1.
# Thus k <= 11 is sufficient.  Negative powers from (z/u)^j can shift an
# E^{-s} coefficient down by at most k-1 <= 10, so retaining its ordinary
# power-series part through degree U_MAX + 10 is sufficient.
K_MAX = U_MAX + 1
U_AUX_MAX = 2 * U_MAX

Rational = Fraction
UniPoly = Dict[int, Rational]
BiPoly = Dict[Tuple[int, int], Rational]  # exponents (z,w)


def add_to(poly: Dict, exponent, coefficient: Rational) -> None:
    """Add one term to a sparse polynomial, deleting exact zeroes."""

    if not coefficient:
        return
    value = poly.get(exponent, Fraction(0)) + coefficient
    if value:
        poly[exponent] = value
    elif exponent in poly:
        del poly[exponent]


def mul_univariate(a: UniPoly, b: UniPoly, degree_limit: int) -> UniPoly:
    """Multiply sparse univariate polynomials and truncate by degree."""

    result: UniPoly = {}
    for degree_a, coefficient_a in a.items():
        for degree_b, coefficient_b in b.items():
            degree = degree_a + degree_b
            if degree <= degree_limit:
                add_to(result, degree, coefficient_a * coefficient_b)
    return result


def power_univariate(base: UniPoly, exponent: int, degree_limit: int) -> UniPoly:
    """Return ``base**exponent`` truncated at ``degree_limit``."""

    result: UniPoly = {0: Fraction(1)}
    factor = dict(base)
    n = exponent
    while n:
        if n & 1:
            result = mul_univariate(result, factor, degree_limit)
        n //= 2
        if n:
            factor = mul_univariate(factor, factor, degree_limit)
    return result


def log_series(unit_series: UniPoly, degree_limit: int) -> UniPoly:
    """Compute log(unit_series), assuming its constant coefficient is one."""

    assert unit_series.get(0) == 1
    x = dict(unit_series)
    del x[0]
    result: UniPoly = {}
    power: UniPoly = {0: Fraction(1)}
    for n in range(1, degree_limit + 1):
        power = mul_univariate(power, x, degree_limit)
        if not power:
            break
        scalar = Fraction(1 if n % 2 else -1, n)
        for degree, coefficient in power.items():
            add_to(result, degree, scalar * coefficient)
    return result


def negative_power_series(
    unit_series: UniPoly, exponent: int, degree_limit: int
) -> UniPoly:
    """Compute ``unit_series**(-exponent)`` by the binomial series."""

    assert exponent >= 1
    assert unit_series.get(0) == 1
    x = dict(unit_series)
    del x[0]
    result: UniPoly = {0: Fraction(1)}
    power: UniPoly = {0: Fraction(1)}
    for n in range(1, degree_limit + 1):
        power = mul_univariate(power, x, degree_limit)
        if not power:
            break
        scalar = Fraction((-1) ** n * comb(exponent + n - 1, n))
        for degree, coefficient in power.items():
            add_to(result, degree, scalar * coefficient)
    return result


def mobius(n: int) -> int:
    """Return the Moebius function mu(n), by trial division."""

    value = 1
    prime = 2
    remaining = n
    while prime * prime <= remaining:
        if remaining % prime == 0:
            remaining //= prime
            value = -value
            if remaining % prime == 0:
                return 0
            while remaining % prime == 0:
                remaining //= prime
        prime += 1
    if remaining > 1:
        value = -value
    return value


def divisors(n: int) -> Iterable[int]:
    """Yield the positive divisors of n."""

    return (d for d in range(1, n + 1) if n % d == 0)


def bernoulli_numbers(n_max: int) -> list[Rational]:
    """Compute B_0,...,B_n_max from their standard exact recurrence."""

    numbers = [Fraction(1)]
    for n in range(1, n_max + 1):
        total = sum(
            (Fraction(comb(n + 1, k)) * numbers[k] for k in range(n)),
            Fraction(0),
        )
        numbers.append(-total / Fraction(n + 1))
    return numbers


def q_series(ell: int) -> UniPoly:
    """Return Q_ell(u)=sum_{d|ell} mu(ell/d) u^(ell-d)."""

    return {
        ell - d: Fraction(mobius(ell // d))
        for d in divisors(ell)
        if mobius(ell // d)
    }


def w_series(ell: int) -> UniPoly:
    """Return W_ell(w), truncated modulo w^15.

    For ell >= 2 the constant term vanishes, since sum_{d|ell}mu(ell/d)=0.
    Thus only the terms -mu(ell/d) w^d / ell need to be recorded.
    """

    result: UniPoly = {}
    for d in divisors(ell):
        if d <= W_MAX:
            add_to(result, d, Fraction(-mobius(ell // d), ell))
    return result


def e_inverse_power(ell: int, exponent: int) -> UniPoly:
    """Return E_ell^{-exponent} as a u-series to the auxiliary cutoff.

    Since E_ell = u^{-ell} Q_ell(u)/ell, this is
    ell^exponent u^(ell*exponent) Q_ell(u)^(-exponent).
    """

    inverse_q = negative_power_series(q_series(ell), exponent, U_AUX_MAX)
    shift = ell * exponent
    scale = Fraction(ell**exponent)
    return {
        shift + degree: scale * coefficient
        for degree, coefficient in inverse_q.items()
        if shift + degree <= U_AUX_MAX
    }


def delta_x_power(ell: int, exponent: int, a: UniPoly) -> list[tuple[int, int, int, Rational]]:
    """Expand (W_ell+c_ell*z/u)^exponent-(c_ell*z/u)^exponent.

    A returned tuple is (u exponent, z exponent, w exponent, coefficient).
    The u exponents can be negative at this intermediate stage.
    """

    c = Fraction(mobius(ell), ell)
    terms: list[tuple[int, int, int, Rational]] = []
    # j is the number of factors c*z/u.  The all-B term j=exponent is
    # exactly the term subtracted in the definition of R.
    for j in range(exponent):
        if j and not c:
            continue
        a_power = power_univariate(a, exponent - j, W_MAX)
        scalar = Fraction(comb(exponent, j)) * c**j
        for w_degree, coefficient in a_power.items():
            terms.append((-j, j, w_degree, scalar * coefficient))
    return terms


def add_factored_contribution(
    lambdas: list[BiPoly],
    x_difference: list[tuple[int, int, int, Rational]],
    u_factor: UniPoly,
    scalar: Rational,
) -> None:
    """Add scalar * x_difference * u_factor to the Lambda coefficients."""

    for u_shift, z_degree, w_degree, x_coefficient in x_difference:
        for u_degree_factor, u_coefficient in u_factor.items():
            u_degree = u_shift + u_degree_factor
            if 1 <= u_degree <= U_MAX:
                add_to(
                    lambdas[u_degree],
                    (z_degree, w_degree),
                    scalar * x_coefficient * u_coefficient,
                )


def reconstruct_lambdas() -> list[BiPoly]:
    """Reconstruct Lambda_1,...,Lambda_10 exactly from (3.3) and (3.11)."""

    lambdas: list[BiPoly] = [{} for _ in range(U_MAX + 1)]

    # The prefactor (1-u)^(1-w) in (3.11).
    for m in range(1, U_MAX + 1):
        add_to(lambdas[m], (0, 0), Fraction(-1, m))
        add_to(lambdas[m], (0, 1), Fraction(1, m))

    bernoulli = bernoulli_numbers(2 * U_MAX)

    # Lemma 3.3 implies ell <= 2*U_MAX.  Each summand below is the
    # numerator-minus-denominator evaluation in (3.11).
    for ell in range(2, 2 * U_MAX + 1):
        a = w_series(ell)

        # Linear term X log(lambda_ell E_ell).  In the difference of the
        # two evaluations, X is replaced simply by W_ell.
        q = q_series(ell)
        one_minus_u_ell: UniPoly = {0: Fraction(1)}
        if ell <= U_MAX:
            one_minus_u_ell[ell] = Fraction(-1)
        lambda_e = mul_univariate(q, one_minus_u_ell, U_MAX)
        log_lambda_e = log_series(lambda_e, U_MAX)
        linear_difference = [(0, 0, w_degree, coefficient) for w_degree, coefficient in a.items()]
        add_factored_contribution(
            lambdas, linear_difference, log_lambda_e, Fraction(1)
        )

        # First inverse-E sum: (1/2) sum_{k>=1} X^k E^{-k}/k.
        for k in range(1, K_MAX + 1):
            minimum_order = ell * k - (k - 1 if mobius(ell) else 0)
            if minimum_order > U_MAX:
                break
            add_factored_contribution(
                lambdas,
                delta_x_power(ell, k, a),
                e_inverse_power(ell, k),
                Fraction(1, 2 * k),
            )

        # Second inverse-E sum: -sum_{k>=2} X^k E^{-(k-1)}/(k(k-1)).
        for k in range(2, K_MAX + 1):
            minimum_order = ell * (k - 1) - (k - 1 if mobius(ell) else 0)
            if minimum_order > U_MAX:
                break
            add_factored_contribution(
                lambdas,
                delta_x_power(ell, k, a),
                e_inverse_power(ell, k - 1),
                Fraction(-1, k * (k - 1)),
            )

        # Bernoulli part.  The conservative finite loops include every
        # term that can reach u^10; the order test discards the rest.
        for r in range(2, 2 * U_MAX + 1, 2):
            for k in range(1, K_MAX + 1):
                e_exponent = r + k - 1
                minimum_order = ell * e_exponent - (k - 1 if mobius(ell) else 0)
                if minimum_order > U_MAX:
                    break
                scalar = (
                    -bernoulli[r]
                    * Fraction(comb(r + k - 2, k), r * (r - 1))
                )
                add_factored_contribution(
                    lambdas,
                    delta_x_power(ell, k, a),
                    e_inverse_power(ell, e_exponent),
                    scalar,
                )

    return lambdas


def coefficient_l1_norm(poly: BiPoly) -> Rational:
    """The norm ||.||_+ from (6.19)."""

    return sum((abs(coefficient) for coefficient in poly.values()), Fraction(0))


def gamma_factor_without_sqrt_pi(m: int) -> tuple[Rational, bool]:
    """Return (c, has_sqrt_pi) with Gamma(m/2+1)=c*sqrt(pi)^has_sqrt_pi."""

    if m % 2 == 0:
        return Fraction(factorial(m // 2)), False
    # m=2k+1: Gamma(k+3/2)=(2k+2)!/(4^(k+1)(k+1)!)*sqrt(pi).
    k = (m - 1) // 2
    coefficient = Fraction(factorial(2 * k + 2), 4 ** (k + 1) * factorial(k + 1))
    return coefficient, True


def scalar_sum_upper_bound(lambdas: list[Rational]) -> tuple[Rational, Rational, Rational]:
    """Return rational part, 1/sqrt(pi) coefficient, and a rational upper bound."""

    rational_part = Fraction(0)
    inverse_sqrt_pi_part = Fraction(0)
    for m in range(1, 11):
        gamma_coefficient, has_sqrt_pi = gamma_factor_without_sqrt_pi(m)
        term = Fraction(m) * lambdas[m] / (Fraction(6**m) * gamma_coefficient)
        if has_sqrt_pi:
            inverse_sqrt_pi_part += term
        else:
            rational_part += term
    # sqrt(pi)>7/4 implies 1/sqrt(pi)<4/7.
    upper = rational_part + Fraction(4, 7) * inverse_sqrt_pi_part
    return rational_part, inverse_sqrt_pi_part, upper


def t_upper_bound(m: int) -> Rational:
    """Exact rational upper bound for T_m from (6.32)-(6.33).

    We use pi>3 in every pi^{-b}; for odd m we additionally use
    1/sqrt(pi)<4/7 for the half-integral Gamma factor.
    """

    gamma_coefficient, has_sqrt_pi = gamma_factor_without_sqrt_pi(m)
    inverse_sqrt_bound = Fraction(4, 7) if has_sqrt_pi else Fraction(1)

    alpha_without_gamma = Fraction(1405 * m * 5**m, 264 * 6**m)
    total = alpha_without_gamma / gamma_coefficient * inverse_sqrt_bound

    for b in range(2, m // 2 + 1):
        beta_without_pi_gamma = Fraction(
            8 * m * 5 ** (m - 2 * b) * factorial(m - b - 1),
            6**m * factorial(m - 2 * b),
        )
        total += (
            beta_without_pi_gamma
            / (Fraction(3**b) * gamma_coefficient)
            * inverse_sqrt_bound
        )
    return total


def format_fraction(value: Rational) -> str:
    """Readable exact rational, followed by a decimal approximation."""

    exact = str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    return f"{exact}  (approximately {float(value):.12g})"


def main() -> int:
    if len(sys.argv) != 1:
        print("This checker takes no command-line arguments.", file=sys.stderr)
        return 2

    print("Certified analytic verification 6.6: finite data for R")
    print("---------------------------------------------------------")
    print("Reconstructing log R through u^10 and modulo w^15")
    print("from equations (3.3) and (3.11), using exact rational arithmetic.")

    polynomials = reconstruct_lambdas()
    # Equation (3.12) gives the support condition z-degree <= u-degree.
    # Checking it here is a useful independent guard on every z/u shift.
    assert all(
        z_degree <= m and w_degree <= W_MAX
        for m in range(1, U_MAX + 1)
        for z_degree, w_degree in polynomials[m]
    )
    computed = [Fraction(0)] + [coefficient_l1_norm(polynomials[m]) for m in range(1, 11)]
    expected = [
        Fraction(0),
        Fraction(3),
        Fraction(23, 6),
        Fraction(29, 6),
        Fraction(20, 3),
        Fraction(141, 10),
        Fraction(403, 18),
        Fraction(1985, 42),
        Fraction(10703, 120),
        Fraction(5273, 30),
        Fraction(256937, 660),
    ]

    print("\nExact coefficient norms lambda_m = ||Lambda_m||_+:")
    all_ok = True
    for m in range(1, 11):
        matches = computed[m] == expected[m]
        all_ok &= matches
        status = "matches" if matches else f"DOES NOT MATCH expected {expected[m]}"
        print(
            f"  m={m:2d}: lambda_m={str(computed[m]):>10s}; "
            f"{len(polynomials[m]):3d} nonzero coefficients; {status}"
        )

    rational_part, inverse_sqrt_part, finite_sum_upper = scalar_sum_upper_bound(computed)
    finite_sum_ok = finite_sum_upper < Fraction(17, 20)
    all_ok &= finite_sum_ok
    print("\nFirst scalar inequality in (6.34):")
    print(
        "  The exact sum has the form A + B/sqrt(pi), with\n"
        f"    A = {format_fraction(rational_part)}\n"
        f"    B = {format_fraction(inverse_sqrt_part)}"
    )
    print(
        "  Using sqrt(pi)>7/4 gives the strict rational upper bound\n"
        f"    A + 4B/7 = {format_fraction(finite_sum_upper)} < 17/20: "
        f"{'yes' if finite_sum_ok else 'NO'}"
    )

    t11_upper = t_upper_bound(11)
    t12_upper = t_upper_bound(12)
    t11_ok = t11_upper < Fraction(1, 32)
    t12_ok = t12_upper < Fraction(1, 96)
    all_ok &= t11_ok and t12_ok
    print("\nRemaining scalar inequalities in (6.34):")
    print("  T_11 and T_12 are evaluated from (6.32)-(6.33).")
    print(
        "  Replacing pi by 3 (and, for T_11, sqrt(pi) by 7/4) gives"
    )
    print(
        f"    T_11 < {format_fraction(t11_upper)} < 1/32: "
        f"{'yes' if t11_ok else 'NO'}"
    )
    print(
        f"    T_12 < {format_fraction(t12_upper)} < 1/96: "
        f"{'yes' if t12_ok else 'NO'}"
    )

    if all_ok:
        print("\nPASS: every statement in Certified analytic verification 6.6 holds.")
        return 0
    print("\nFAIL: at least one reconstructed value or inequality is incorrect.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
