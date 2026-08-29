#!/usr/bin/env python3
"""Recheck Certified analytic verification 5.9.

This script is a direct interval evaluation of equations (5.25)--(5.46) in
the manuscript.  It has no command-line options and uses no saved numerical
certificate values.  The only imported numerical package is NumPy, which is
used to evaluate many independent intervals in parallel.

The interval implementation below uses IEEE-754 binary64 arithmetic.  Every
basic arithmetic result is enlarged by one representable number in each
outward direction.  The only transcendental constants/functions needed by
the proof (pi, -log(1-y), and exp(-L)) are enclosed by finite power series
with explicit remainder bounds before the remaining outward arithmetic is
performed.

Exit status is zero only if every strict inequality claimed in the box has
been certified.
"""

from __future__ import annotations

import math
import sys
import time
from dataclasses import dataclass
from fractions import Fraction

import numpy as np


# Polynomial arithmetic is always truncated modulo w^15.
W_DEGREE = 14
H_RADIUS = 13_000
X_BOXES = 256

NEGATIVE_INFINITY = np.float64(-np.inf)
POSITIVE_INFINITY = np.float64(np.inf)


def _down(value):
    """Round a binary64 scalar/array one step toward minus infinity."""

    return np.nextafter(np.asarray(value, dtype=np.float64), NEGATIVE_INFINITY)


def _up(value):
    """Round a binary64 scalar/array one step toward plus infinity."""

    return np.nextafter(np.asarray(value, dtype=np.float64), POSITIVE_INFINITY)


@dataclass(frozen=True)
class Interval:
    """A closed interval, or an elementwise array of closed intervals.

    Point intervals should only be used for exactly representable binary64
    values (in this script: small integers and dyadic mesh coordinates).
    General rational constants are constructed with ``rational`` below.
    """

    lo: object
    hi: object

    def __post_init__(self) -> None:
        object.__setattr__(self, "lo", np.asarray(self.lo, dtype=np.float64))
        object.__setattr__(self, "hi", np.asarray(self.hi, dtype=np.float64))
        if np.any(self.lo > self.hi):
            raise ValueError("invalid interval")

    @staticmethod
    def point(value) -> "Interval":
        value = np.asarray(value, dtype=np.float64)
        return Interval(value, value)

    @staticmethod
    def hull(lo, hi) -> "Interval":
        return Interval(lo, hi)

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def _exact_scalar(self, value: float) -> bool:
        return self.lo.shape == () and self.lo == value and self.hi == value

    def __add__(self, other) -> "Interval":
        other = as_interval(other)
        if self._exact_scalar(0.0):
            return other
        if other._exact_scalar(0.0):
            return self
        return Interval(_down(self.lo + other.lo), _up(self.hi + other.hi))

    def __radd__(self, other) -> "Interval":
        return self + other

    def __sub__(self, other) -> "Interval":
        other = as_interval(other)
        if other._exact_scalar(0.0):
            return self
        if self._exact_scalar(0.0):
            return -other
        return Interval(_down(self.lo - other.hi), _up(self.hi - other.lo))

    def __rsub__(self, other) -> "Interval":
        return as_interval(other) - self

    def __mul__(self, other) -> "Interval":
        other = as_interval(other)
        shape = np.broadcast_shapes(self.lo.shape, other.lo.shape)
        if self._exact_scalar(0.0) or other._exact_scalar(0.0):
            return Interval.point(np.zeros(shape, dtype=np.float64))
        if self._exact_scalar(1.0):
            return other
        if other._exact_scalar(1.0):
            return self
        if self._exact_scalar(-1.0):
            return -other
        if other._exact_scalar(-1.0):
            return -self
        # Multiplication by an exact point needs only two endpoint products.
        if self.lo.shape == () and self.lo == self.hi:
            c = float(self.lo)
            if c >= 0.0:
                return Interval(_down(c * other.lo), _up(c * other.hi))
            return Interval(_down(c * other.hi), _up(c * other.lo))
        if other.lo.shape == () and other.lo == other.hi:
            return other * self
        products = np.stack(
            (
                self.lo * other.lo,
                self.lo * other.hi,
                self.hi * other.lo,
                self.hi * other.hi,
            )
        )
        return Interval(_down(np.min(products, axis=0)), _up(np.max(products, axis=0)))

    def __rmul__(self, other) -> "Interval":
        return self * other

    def reciprocal(self) -> "Interval":
        if np.any((self.lo <= 0.0) & (self.hi >= 0.0)):
            raise ZeroDivisionError("interval contains zero")
        return Interval(_down(1.0 / self.hi), _up(1.0 / self.lo))

    def __truediv__(self, other) -> "Interval":
        return self * as_interval(other).reciprocal()

    def __rtruediv__(self, other) -> "Interval":
        return as_interval(other) / self


def as_interval(value) -> Interval:
    if isinstance(value, Interval):
        return value
    # All bare numerical constants in the formulas are integers or dyadic
    # binary64 numbers, hence exactly representable.
    return Interval.point(value)


def rational(numerator: int, denominator: int = 1) -> Interval:
    """Return an outward binary64 enclosure of an exact rational number."""

    if denominator <= 0:
        raise ValueError("denominator must be positive")
    if numerator == 0:
        return Interval.point(0.0)
    nearest = np.float64(numerator / denominator)
    # Python's int/int conversion is correctly rounded; widening on both
    # sides also makes the construction independent of equality detection.
    return Interval(_down(nearest), _up(nearest))


# An interval, rather than the nearest binary64 literal, represents the exact
# threshold 1/5.  Whenever a scalar endpoint is used below, it is the upper
# endpoint, making the numerical test (very slightly) stronger than required.
ONE_FIFTH = rational(1, 5)


def interval_power(base: Interval, exponent: int) -> Interval:
    """Nonnegative integer power by repeated squaring."""

    if exponent < 0:
        raise ValueError("negative exponent")
    result = Interval.point(np.ones_like(base.lo))
    factor = base
    while exponent:
        if exponent & 1:
            result = result * factor
        exponent >>= 1
        if exponent:
            factor = factor * factor
    return result


def arctan_inverse(q: int, terms: int) -> Interval:
    """Enclose arctan(1/q) by its alternating series."""

    q_interval = Interval.point(float(q))
    q_squared = q_interval * q_interval
    term = Interval.point(1.0) / q_interval
    total = Interval.point(0.0)
    for k in range(terms):
        addend = term / Interval.point(float(2 * k + 1))
        total = total - addend if k & 1 else total + addend
        term = term / q_squared
    # The next alternating-series term bounds the remainder in magnitude.
    remainder = term / Interval.point(float(2 * terms + 1))
    radius = np.maximum(np.abs(remainder.lo), np.abs(remainder.hi))
    return total + Interval(-radius, radius)


def pi_interval() -> Interval:
    """Machin's formula: pi = 16 atan(1/5) - 4 atan(1/239)."""

    return 16 * arctan_inverse(5, 50) - 4 * arctan_inverse(239, 15)


def negative_log_one_minus(y: Interval) -> Interval:
    """Enclose -log(1-y) for the small nonnegative intervals used here.

    We sum through y^10/10.  The remaining positive tail is at most

        y^11 / (11 (1-y)).

    Throughout this script y < 1/400, so this remainder is far below the
    widths introduced by the input boxes.
    """

    # The arguments are products of nonnegative x and s.  A preceding
    # outward operation may have introduced the harmless value -5e-324 at
    # an exact zero endpoint, so impose this known domain information.
    y = Interval(np.maximum(0.0, y.lo), y.hi)
    if np.any(y.hi >= 1.0):
        raise ValueError("negative_log_one_minus called outside 0 <= y < 1")
    power = y
    total = Interval.point(np.zeros_like(y.lo))
    for k in range(1, 11):
        total = total + power / Interval.point(float(k))
        power = power * y
    remainder = power / (Interval.point(11.0) * (1 - y))
    # Every omitted term is nonnegative.
    return Interval(total.lo, (total + remainder).hi)


def exp_negative_points(y) -> Interval:
    """Enclose exp(-y) at nonnegative (dyadic) point arguments y <= 32.

    First enclose exp(y) using its Taylor polynomial through degree 130.
    Since 0 <= y <= 32 < 132, the positive remainder is bounded by

        next_term / (1 - y/132).

    Taking the reciprocal gives exp(-y).  The routine is vectorized because
    all endpoints of an L-mesh are handled at once.
    """

    y = np.asarray(y, dtype=np.float64)
    if np.any(y < 0.0) or np.any(y > 32.0):
        raise ValueError("exp_negative_points requires 0 <= y <= 32")
    argument = Interval.point(y)
    term = Interval.point(np.ones_like(y))
    total = Interval.point(np.ones_like(y))
    for k in range(1, 131):
        term = term * argument / Interval.point(float(k))
        total = total + term
    next_term = term * argument / Interval.point(131.0)
    remainder = next_term / (1 - argument / Interval.point(132.0))
    exp_positive = Interval(total.lo, (total + remainder).hi)
    return exp_positive.reciprocal()


def polynomial_zero(shape=()) -> list[Interval]:
    return [Interval.point(np.zeros(shape, dtype=np.float64)) for _ in range(W_DEGREE + 1)]


def polynomial_multiply(a: list[Interval], b: list[Interval]) -> list[Interval]:
    """Multiply interval polynomials modulo w^15."""

    shape = np.broadcast_shapes(a[0].lo.shape, b[0].lo.shape)
    out = polynomial_zero(shape)
    for i in range(W_DEGREE + 1):
        for j in range(W_DEGREE + 1 - i):
            out[i + j] = out[i + j] + a[i] * b[j]
    return out


def l_sigma(sigma: int, pi: Interval) -> list[Interval]:
    """Construct the real polynomial L_sigma from equation (5.30)."""

    sigma %= 4
    out = polynomial_zero()
    two_pi_power = Interval.point(1.0)
    factorial = 1
    for k in range(1, W_DEGREE + 1):
        two_pi_power = two_pi_power * (2 * pi)
        factorial *= k
        # The two complex powers cancel unless k == sigma+1 (mod 2).
        if (k - sigma - 1) % 2 == 0:
            sign_exponent = k + (k - sigma - 1) // 2
            coefficient = 2 * two_pi_power / Interval.point(float(factorial))
            out[k] = -coefficient if sign_exponent & 1 else coefficient
    return out


def h_polynomial(rho: int, m: int, x: Interval, pi: Interval) -> list[Interval]:
    """Construct H_{rho,m} in (5.31), including the exact-zero rule (5.33)."""

    sigma = (rho - 1 - m) % 4
    previous_sigma = (sigma - 1) % 4
    h = l_sigma(sigma, pi)
    previous = l_sigma(previous_sigma, pi)
    z = x / (1 - (m + 1) * x)

    w_minus_w2 = polynomial_zero()
    w_minus_w2[1] = pi * z
    w_minus_w2[2] = -(pi * z)
    correction = polynomial_multiply(w_minus_w2, previous)
    error_radius = H_RADIUS * z * z
    symmetric_error = Interval(-error_radius.hi, error_radius.hi)
    for k in range(W_DEGREE + 1):
        h[k] = h[k] + correction[k]
        if k >= 1:
            h[k] = h[k] + symmetric_error
    h[0] = Interval.point(0.0)

    # Here sigma is congruent to g-1-m modulo 4.  If it is odd, the
    # underlying index is > 1 and odd (g >= 2000), so (5.33) is exact.
    if sigma & 1:
        h[1] = Interval.point(0.0)
    return h


def nonnegative_endpoint(value: Interval) -> Interval:
    """Use the known nonnegative endpoint in (5.25)--(5.26).

    At an actual manuscript parameter, ``u=1/a`` and ``v=1/(a+n-j)``
    with ``n>=j``, so ``u-v`` (and every difference of their positive
    powers) is nonnegative.  An x-by-s product box also contains artificial
    independent pairs that need not satisfy this relation.  Intersecting
    with the analytically known half-line and retaining the computed upper
    endpoint gives an interval containing every actual endpoint.
    """

    return Interval(np.zeros_like(value.hi), np.maximum(0.0, value.hi))


def v_polynomial(m: int, j: int, x: Interval, s: Interval) -> list[Interval]:
    """Construct V_{m,j} from (5.25)--(5.28)."""

    a = m + j + 2
    b = m + 2
    u = x / (1 - b * x)
    v = x * s / (1 - a * x * s)

    beta: list[Interval | None] = [None] * (W_DEGREE + 1)
    minus_log_a = negative_log_one_minus(a * x * s)
    minus_log_b = negative_log_one_minus(b * x)
    beta[1] = minus_log_b - minus_log_a + nonnegative_endpoint(u - v)

    for q in range(2, W_DEGREE + 1):
        integral = (interval_power(u, q - 1) - interval_power(v, q - 1)) / (q - 1)
        endpoint = nonnegative_endpoint(interval_power(u, q) - interval_power(v, q))
        unsigned = (integral + endpoint) / q
        beta[q] = unsigned if q & 1 else -unsigned

    shape = np.broadcast_shapes(x.lo.shape, s.lo.shape)
    coefficients = polynomial_zero(shape)
    coefficients[0] = Interval.point(np.ones(shape, dtype=np.float64))
    for k in range(1, W_DEGREE + 1):
        value = Interval.point(np.zeros(shape, dtype=np.float64))
        for q in range(1, k + 1):
            assert beta[q] is not None
            value = value + q * beta[q] * coefficients[k - q]
        coefficients[k] = value / k
    return coefficients


def r_polynomial(m: int, j: int) -> list[Interval]:
    """Return R_{0,0}, R_{1,0}, or R_{1,1} from the manuscript."""

    out = polynomial_zero()
    if (m, j) == (0, 0):
        out[0] = Interval.point(1.0)
    elif (m, j) == (1, 0):
        out[0] = Interval.point(-1.0)
        out[1] = Interval.point(0.5)
        out[2] = Interval.point(0.5)
    elif (m, j) == (1, 1):
        out[1] = Interval.point(0.5)
        out[2] = Interval.point(-0.5)
    else:
        raise ValueError("unexpected (m,j)")
    return out


def omega(m: int, j: int, x: Interval, s: Interval, pi: Interval) -> Interval:
    """The prefactors in (5.34); omega_{0,0}=1."""

    if m == 0:
        return Interval.point(np.ones(np.broadcast_shapes(x.lo.shape, s.lo.shape)))
    value = 2 * pi * x * s * (1 - 3 * x)
    value = value / ((1 - 2 * x) * (1 - 3 * x * s))
    if j == 1:
        value = value * (1 - s) / (1 - 4 * x * s)
    return value


def m_base_polynomial(rho: int, x: Interval, s: Interval, pi: Interval) -> list[Interval]:
    """Return F with M_Gamma = T_<=Gamma(e^(Lw) F), cf. (5.34)."""

    shape = np.broadcast_shapes(x.lo.shape, s.lo.shape)
    total = polynomial_zero(shape)
    h_cache = {m: h_polynomial(rho, m, x, pi) for m in (0, 1)}
    for m, j in ((0, 0), (1, 0), (1, 1)):
        rh = polynomial_multiply(r_polynomial(m, j), h_cache[m])
        term = polynomial_multiply(rh, v_polynomial(m, j, x, s))
        prefactor = omega(m, j, x, s, pi)
        for k in range(W_DEGREE + 1):
            total[k] = total[k] + prefactor * term[k]
    return total


def m_coefficients(gamma: int, base: list[Interval]) -> list[Interval]:
    """Coefficients c_r of M_Gamma(L)=sum c_r L^r.

    Applying T_<=Gamma to e^(Lw) F(w) gives

        c_r = (F_0 + ... + F_{Gamma-r}) / r!.
    """

    prefix: list[Interval] = []
    running = 0 * base[0]
    for coefficient in base:
        running = running + coefficient
        prefix.append(running)
    result: list[Interval] = []
    factorial = 1
    for r in range(gamma + 1):
        if r:
            factorial *= r
        result.append(prefix[gamma - r] / factorial)
    return result


def evaluate_at_centered_interval(coefficients: list[Interval], variable: Interval) -> Interval:
    """Evaluate a polynomial after the exact change variable=center+t.

    Centering is essential: natural Horner evaluation on an L-box near 32
    would lose far more information than the narrow box warrants.
    """

    center_value = (variable.lo + variable.hi) * 0.5
    center = Interval.point(center_value)
    displacement = variable - center

    # Translate the coefficient list by synthetic Horner steps.  At each
    # stage ``shifted`` contains coefficients in powers of t.
    shifted = [coefficients[-1]]
    for original in reversed(coefficients[:-1]):
        new = [shifted[0] * center + original]
        for k in range(1, len(shifted)):
            new.append(shifted[k - 1] + shifted[k] * center)
        new.append(shifted[-1])
        shifted = new

    value = shifted[-1]
    for coefficient in reversed(shifted[:-1]):
        value = value * displacement + coefficient
    return value


def regularized_reciprocal_polynomial(coefficients: list[Interval], degree: int, theta: Interval) -> Interval:
    """Evaluate theta^degree P(1/theta), including theta=0, without division."""

    if len(coefficients) < degree + 1:
        raise ValueError("not enough coefficients")
    # Horner with the ascending input c_0,...,c_degree evaluates
    # c_0 theta^degree + ... + c_{degree-1} theta + c_degree.
    value = coefficients[0]
    for coefficient in coefficients[1 : degree + 1]:
        value = value * theta + coefficient
    return value


def l_base_mesh(parity: str) -> list[tuple[Fraction, Fraction]]:
    """Return the 526/524-leaf dyadic mesh used before the final 8-way split.

    Away from one narrow sign-change region it is the uniform mesh of width
    1/16.  The two explicit cut lists below replace that one cell.  Listing
    the cuts makes the numerical domain decomposition fully inspectable and
    avoids importing any data from an older certificate bundle.
    """

    if parity == "even":
        special_lo, special_hi = Fraction(1, 4), Fraction(5, 16)
        cuts = [
            Fraction(1, 4), Fraction(17, 64), Fraction(35, 128),
            Fraction(71, 256), Fraction(143, 512), Fraction(287, 1024),
            Fraction(575, 2048), Fraction(4601, 16384),
            Fraction(36809, 131072), Fraction(73619, 262144),
            Fraction(18405, 65536), Fraction(9203, 32768),
            Fraction(2301, 8192), Fraction(1151, 4096),
            Fraction(9, 32), Fraction(5, 16),
        ]
    elif parity == "odd":
        special_lo, special_hi = Fraction(19, 16), Fraction(5, 4)
        cuts = [
            Fraction(19, 16), Fraction(153, 128), Fraction(2449, 2048),
            Fraction(4899, 4096), Fraction(9799, 8192),
            Fraction(19599, 16384), Fraction(39199, 32768),
            Fraction(78399, 65536), Fraction(1225, 1024),
            Fraction(613, 512), Fraction(307, 256), Fraction(77, 64),
            Fraction(39, 32), Fraction(5, 4),
        ]
    else:
        raise ValueError("parity must be 'even' or 'odd'")

    mesh: list[tuple[Fraction, Fraction]] = []
    for k in range(32 * 16):
        lo, hi = Fraction(k, 16), Fraction(k + 1, 16)
        if (lo, hi) == (special_lo, special_hi):
            mesh.extend(zip(cuts[:-1], cuts[1:]))
        else:
            mesh.append((lo, hi))
    return mesh


def l_mesh(parity: str) -> tuple[np.ndarray, np.ndarray]:
    """Split every base leaf into eight equal dyadic intervals."""

    leaves: list[tuple[Fraction, Fraction]] = []
    for lo, hi in l_base_mesh(parity):
        step = (hi - lo) / 8
        leaves.extend((lo + k * step, lo + (k + 1) * step) for k in range(8))
    lo_array = np.array([float(lo) for lo, _ in leaves], dtype=np.float64)
    hi_array = np.array([float(hi) for _, hi in leaves], dtype=np.float64)
    # Every denominator is a power of two, so these conversions are exact.
    if not np.all(lo_array[1:] == hi_array[:-1]):
        raise AssertionError("L mesh is not contiguous")
    if lo_array[0] != 0.0 or hi_array[-1] != 32.0:
        raise AssertionError("L mesh does not cover [0,32]")
    return lo_array, hi_array


def finite_l_check(parity: str, rho: int, pi: Interval) -> dict[str, object]:
    """Verify (5.36) for one parity on all x-by-L product boxes."""

    l_lo, l_hi = l_mesh(parity)
    expected_l_boxes = 4208 if parity == "even" else 4192
    if len(l_lo) != expected_l_boxes:
        raise AssertionError(f"unexpected {parity} L-box count: {len(l_lo)}")
    l_interval = Interval(l_lo, l_hi)
    exp_at_lo = exp_negative_points(l_lo)
    exp_at_hi = exp_negative_points(l_hi)
    s_interval = Interval(exp_at_hi.lo, exp_at_lo.hi)

    minimum_margin = math.inf
    worst: tuple[int, int, int, int, float, float] | None = None
    counts = {(10, 1): 0, (10, -1): 0, (14, 1): 0, (14, -1): 0}

    for x_index in range(X_BOXES):
        x_left = rational(x_index, 512_000)
        x_right = rational(x_index + 1, 512_000)
        x = Interval(x_left.lo, x_right.hi)

        base = m_base_polynomial(rho, x, s_interval, pi)
        value10 = evaluate_at_centered_interval(m_coefficients(10, base), l_interval)
        value14 = evaluate_at_centered_interval(m_coefficients(14, base), l_interval)

        # A positive number below certifies epsilon*M - 1/5 > 0.
        margins = np.stack(
            (
                _down(value10.lo - ONE_FIFTH.hi),
                _down(-value10.hi - ONE_FIFTH.hi),
                _down(value14.lo - ONE_FIFTH.hi),
                _down(-value14.hi - ONE_FIFTH.hi),
            )
        )
        witnesses = np.argmax(margins, axis=0)
        best = np.max(margins, axis=0)
        if np.any(best <= 0.0):
            bad_l = int(np.flatnonzero(best <= 0.0)[0])
            raise RuntimeError(
                f"inconclusive {parity} product box x={x_index}, L={bad_l}: "
                f"M10=[{value10.lo[bad_l]:.17g},{value10.hi[bad_l]:.17g}], "
                f"M14=[{value14.lo[bad_l]:.17g},{value14.hi[bad_l]:.17g}]"
            )

        labels = ((10, 1), (10, -1), (14, 1), (14, -1))
        for witness_index, label in enumerate(labels):
            counts[label] += int(np.count_nonzero(witnesses == witness_index))

        local_l = int(np.argmin(best))
        local_margin = float(best[local_l])
        if local_margin < minimum_margin:
            minimum_margin = local_margin
            gamma, sign = labels[int(witnesses[local_l])]
            chosen = value10 if gamma == 10 else value14
            worst = (
                x_index, local_l, gamma, sign,
                float(chosen.lo[local_l]), float(chosen.hi[local_l]),
            )

        if (x_index + 1) % 64 == 0:
            print(
                f"    {parity}: checked {x_index + 1:3d}/{X_BOXES} x-intervals "
                f"({(x_index + 1) * expected_l_boxes:,} product boxes)",
                flush=True,
            )

    assert worst is not None
    return {
        "parity": parity,
        "l_boxes": expected_l_boxes,
        "product_boxes": X_BOXES * expected_l_boxes,
        "minimum_margin": minimum_margin,
        "worst": worst,
        "counts": counts,
    }


def tail_check(pi: Interval) -> dict[str, object]:
    """Verify the two regular compact inequalities (5.39) and (5.46)."""

    x = Interval(0.0, rational(1, 2000).hi)
    exp32 = exp_negative_points(np.float64(32.0))
    s = Interval(0.0, exp32.hi)
    theta = Interval(0.0, rational(1, 32).hi)

    # Odd g: theta^9 M(1/theta), evaluated from c_0,...,c_9 by Horner.
    odd_base = m_base_polynomial(1, x, s, pi)
    odd_coefficients = m_coefficients(10, odd_base)
    odd_regular = regularized_reciprocal_polynomial(
        odd_coefficients[:10], 9, theta
    )
    theta9 = interval_power(theta, 9)
    odd_signed = -odd_regular - ONE_FIFTH * theta9
    odd_margin = float(np.asarray(odd_signed.lo))
    if odd_margin <= 0.0:
        raise RuntimeError(
            "odd L>=32 compact box is inconclusive: "
            f"left side of (5.39) = [{odd_signed.lo},{odd_signed.hi}]"
        )

    # Even g: remove the s*J*L^9 term and regularize degrees 0,...,8.
    # This is the polynomial A in (5.42); no expression 1/theta is formed.
    even_base = m_base_polynomial(0, x, s, pi)
    even_coefficients = m_coefficients(10, even_base)
    even_regular = regularized_reciprocal_polynomial(
        even_coefficients[:9], 8, theta
    )

    h_even_1 = h_polynomial(0, 1, x, pi)
    j_value = -(2 * pi * x * (1 - 3 * x) * h_even_1[1])
    j_value = j_value / (
        Interval.point(float(math.factorial(9))) * (1 - 2 * x) * (1 - 3 * x * s)
    )
    kappa = Interval(0.0, (32 * exp32).hi)
    even_left = even_regular + kappa * j_value
    theta8 = interval_power(theta, 8)
    even_signed = -even_left - ONE_FIFTH * theta8
    even_margin = float(np.asarray(even_signed.lo))
    if even_margin <= 0.0:
        raise RuntimeError(
            "even L>=32 compact box is inconclusive: "
            f"left side of (5.46) = [{even_signed.lo},{even_signed.hi}]"
        )

    return {
        "odd_interval": (float(odd_signed.lo), float(odd_signed.hi)),
        "odd_margin": odd_margin,
        "even_interval": (float(even_signed.lo), float(even_signed.hi)),
        "even_margin": even_margin,
    }


def print_finite_summary(result: dict[str, object]) -> None:
    parity = str(result["parity"])
    x_index, l_index, gamma, sign, value_lo, value_hi = result["worst"]
    epsilon = "+1" if sign == 1 else "-1"
    print(f"  {parity.capitalize()} parity (rho={0 if parity == 'even' else 1}):")
    print(f"    L intervals:              {result['l_boxes']:,}")
    print(f"    x-by-L product boxes:     {result['product_boxes']:,}")
    print("    Chosen witnesses (Gamma, epsilon):")
    counts = result["counts"]
    for label in ((10, 1), (10, -1), (14, 1), (14, -1)):
        print(f"      ({label[0]:2d}, {label[1]:+d}): {counts[label]:,}")
    print(
        "    Smallest certified margin epsilon*M - 1/5: "
        f"{result['minimum_margin']:.9e}"
    )
    print(
        f"      attained for x-box {x_index}, L-box {l_index}, "
        f"Gamma={gamma}, epsilon={epsilon}; M in [{value_lo:.12g}, {value_hi:.12g}]"
    )


def main() -> int:
    if len(sys.argv) != 1:
        print("This certificate script takes no command-line arguments.", file=sys.stderr)
        return 2
    started = time.perf_counter()
    print("Certified analytic verification 5.9")
    print("====================================")
    print("Reconstructing equations (5.25)--(5.34) with outward interval arithmetic.")
    print("No saved function values or previously certified lower bounds are used.\n")

    try:
        pi = pi_interval()
        print(f"Certified pi interval: [{float(pi.lo):.17g}, {float(pi.hi):.17g}]\n")
        print("Case A: 0 <= L <= 32")
        print("Checking (5.36) with Gamma in {10,14} and epsilon in {-1,+1}.")
        even = finite_l_check("even", 0, pi)
        print_finite_summary(even)
        odd = finite_l_check("odd", 1, pi)
        print_finite_summary(odd)
        total = int(even["product_boxes"]) + int(odd["product_boxes"])
        if total != 2_150_400:
            raise AssertionError(f"total product-box count is {total}, expected 2,150,400")
        print(f"  Total Case-A product boxes: {total:,}\n")

        print("Case B: L >= 32")
        print("Checking the regularized compact polynomials, never evaluating 1/theta.")
        tail = tail_check(pi)
        print(
            "  Odd box, left side of (5.39):  "
            f"[{tail['odd_interval'][0]:.12g}, {tail['odd_interval'][1]:.12g}]"
        )
        print(
            "  Even box, left side of (5.46): "
            f"[{tail['even_interval'][0]:.12g}, {tail['even_interval'][1]:.12g}]"
        )

    except Exception as error:
        print(f"\nFAIL / INCONCLUSIVE: {error}", file=sys.stderr)
        return 1

    elapsed = time.perf_counter() - started
    print(f"\nPASS: every assertion in Certified analytic verification 5.9 was verified.")
    print(f"Elapsed time: {elapsed:.1f} seconds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
