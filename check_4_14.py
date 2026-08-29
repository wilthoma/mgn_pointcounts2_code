#!/usr/bin/env python3
"""Recheck Finite arithmetic verification 4.14.

For each 1 <= g <= 2000, the manuscript asks for a finite list P_g of
primes p such that

    p > 4012,                  p does not divide A_g,

and the closed intervals

    [9*p + g + 1, 10*p - d_g - 1]

cover every integer 50000 <= n < 2000000.  This script first *constructs*
such lists by a deterministic greedy search and then verifies the resulting
witnesses in a separate pass.

All arithmetic relevant to the certificate is exact.  Python integers are
arbitrary-precision, and python-flint supplies exact rational Bernoulli
numbers.  No command-line arguments or external data files are used.
"""

from bisect import bisect_right
from math import isqrt
import sys

try:
    from flint import fmpq
except ImportError as exc:
    raise SystemExit(
        "This checker requires python-flint. Install the pinned dependencies "
        "with: python3 -m pip install -r requirements.txt"
    ) from exc


G_MAX = 2_000
N_MIN = 50_000
N_MAX = 2_000_000 - 1       # The manuscript has the strict bound n < 2000000.
PRIME_LOWER_BOUND = 2 * G_MAX + 12


def require(condition: bool, message: str) -> None:
    """Raise a readable certificate failure instead of silently continuing."""

    if not condition:
        raise AssertionError(message)


def genus_data(g: int, bernoulli_cache: dict[int, fmpq]) -> tuple[int, int]:
    """Return the pair (A_g, d_g) defined in equation (4.21).

    For g >= 3, let h=floor((g-1)/2), d_g=2h-2, and reduce

        (-1)^(g+1) B_(2h)/(2h) = A_g/C_g,   C_g > 0.

    An ``fmpq`` is always stored in lowest terms with positive denominator,
    so its numerator is exactly A_g.  For g=1,2 the manuscript sets A_g=1
    and d_g=0.
    """

    if g <= 2:
        return 1, 0

    h = (g - 1) // 2
    if h not in bernoulli_cache:
        bernoulli_cache[h] = fmpq.bernoulli(2 * h)

    reduced_value = ((-1) ** (g + 1)) * bernoulli_cache[h] / (2 * h)
    return int(reduced_value.numerator), 2 * h - 2


def primes_up_to(limit: int) -> list[int]:
    """Generate prime candidates by the elementary sieve of Eratosthenes."""

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            first = p * p
            count = (limit - first) // p + 1
            sieve[first : limit + 1 : p] = b"\x00" * count
    return [p for p in range(PRIME_LOWER_BOUND + 1, limit + 1) if sieve[p]]


def is_prime_by_trial_division(n: int) -> bool:
    """Independently prove primality by deterministic trial division.

    The witnesses are at most about 222000, so this deliberately simple
    verifier is fast enough.  It is independent of the sieve used to find
    candidates: every possible divisor through floor(sqrt(n)) is tested.
    """

    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    divisor = 3
    while divisor * divisor <= n:
        if n % divisor == 0:
            return False
        divisor += 2
    return True


def construct_witnesses(
    genus_values: list[tuple[int, int]], prime_candidates: list[int]
) -> list[list[int]]:
    """Construct P_g by repeatedly extending the already covered interval.

    If ``next_n`` is the first integer not yet covered, a prime can start its
    interval no later than next_n exactly when

        p <= floor((next_n-g-1)/9).

    Among admissible candidates we take the largest such p.  This maximizes
    the right endpoint 10p-d_g-1 and is therefore the natural greedy choice.
    The success of the search is not itself used as a proof: the complete
    lists are checked afresh below.
    """

    all_lists: list[list[int]] = []
    for g, (a_g, d_g) in enumerate(genus_values, start=1):
        next_n = N_MIN
        witnesses: list[int] = []

        while next_n <= N_MAX:
            largest_allowed = (next_n - g - 1) // 9
            index = bisect_right(prime_candidates, largest_allowed) - 1

            # A_g has few large prime divisors, but skip any that occur.
            while index >= 0 and a_g % prime_candidates[index] == 0:
                index -= 1
            require(
                index >= 0,
                f"construction stopped at g={g}, n={next_n}: no candidate prime",
            )

            p = prime_candidates[index]
            left = 9 * p + g + 1
            right = 10 * p - d_g - 1
            require(
                left <= next_n <= right,
                "greedy interval does not contain the next uncovered integer "
                f"(g={g}, n={next_n}, p={p}, interval=[{left},{right}])",
            )
            require(
                not witnesses or p > witnesses[-1],
                f"constructed primes are not strictly increasing for g={g}",
            )

            witnesses.append(p)
            next_n = right + 1

        all_lists.append(witnesses)

    return all_lists


def verify_witnesses(
    genus_values: list[tuple[int, int]], all_lists: list[list[int]]
) -> dict[str, object]:
    """Verify conditions (i)--(iii) of the boxed statement from scratch."""

    require(len(all_lists) == G_MAX, "there must be exactly one list for each genus")
    distinct_primes = sorted({p for witnesses in all_lists for p in witnesses})

    # Condition (i): a transparent, deterministic primality proof for every
    # distinct recorded integer, followed by the explicit lower-bound check.
    for p in distinct_primes:
        require(p > PRIME_LOWER_BOUND, f"recorded p={p} does not exceed 4012")
        require(is_prime_by_trial_division(p), f"recorded p={p} is not prime")

    total_occurrences = 0
    list_sizes: list[int] = []
    tightest_join = None
    smallest_terminal_excess = None
    smallest_initial_excess = None

    for g, ((a_g, d_g), witnesses) in enumerate(
        zip(genus_values, all_lists), start=1
    ):
        require(witnesses, f"P_{g} is empty")
        require(
            len(witnesses) == len(set(witnesses)),
            f"P_{g} contains a duplicate prime",
        )

        # Condition (ii), using the exact numerator A_g computed above.
        for p in witnesses:
            require(a_g % p != 0, f"p={p} divides A_{g}")

        # Condition (iii).  For closed integer intervals, there is no gap
        # after [L,R] precisely when the next left endpoint is <= R+1.
        intervals = sorted(
            (9 * p + g + 1, 10 * p - d_g - 1, p) for p in witnesses
        )
        first_left = intervals[0][0]
        require(
            first_left <= N_MIN,
            f"P_{g} leaves the initial integer {N_MIN} uncovered",
        )
        initial_excess = N_MIN - first_left
        if smallest_initial_excess is None or initial_excess < smallest_initial_excess[0]:
            smallest_initial_excess = (initial_excess, g)

        covered_through = intervals[0][1]
        for left, right, p in intervals[1:]:
            join_slack = covered_through + 1 - left
            require(
                join_slack >= 0,
                f"P_{g} has a gap before interval for p={p}: "
                f"covered through {covered_through}, next left endpoint {left}",
            )
            if tightest_join is None or join_slack < tightest_join[0]:
                tightest_join = (join_slack, g, p)
            covered_through = max(covered_through, right)

        require(
            covered_through >= N_MAX,
            f"P_{g} covers only through {covered_through}, not through {N_MAX}",
        )
        terminal_excess = covered_through - N_MAX
        if (
            smallest_terminal_excess is None
            or terminal_excess < smallest_terminal_excess[0]
        ):
            smallest_terminal_excess = (terminal_excess, g)

        total_occurrences += len(witnesses)
        list_sizes.append(len(witnesses))

    largest_list_size = max(list_sizes)
    smallest_list_size = min(list_sizes)
    return {
        "distinct_primes": distinct_primes,
        "total_occurrences": total_occurrences,
        "largest_list_size": largest_list_size,
        "largest_list_genera": [
            g for g, size in enumerate(list_sizes, start=1) if size == largest_list_size
        ],
        "smallest_list_size": smallest_list_size,
        "smallest_list_genera": [
            g for g, size in enumerate(list_sizes, start=1) if size == smallest_list_size
        ],
        "tightest_join": tightest_join,
        "smallest_initial_excess": smallest_initial_excess,
        "smallest_terminal_excess": smallest_terminal_excess,
    }


def abbreviated_genus_list(genera: list[int]) -> str:
    """Keep diagnostic output readable when many genera attain an extremum."""

    if len(genera) <= 8:
        return ", ".join(map(str, genera))
    first = ", ".join(map(str, genera[:8]))
    return f"{first}, ... ({len(genera)} genera total)"


def main() -> None:
    if len(sys.argv) != 1:
        raise SystemExit("This checker takes no command-line arguments.")

    print("Finite arithmetic verification 4.14")
    print("=" * 39)
    print("Target: for every 1 <= g <= 2000, construct primes p > 4012")
    print("with p not dividing A_g whose closed intervals")
    print("    [9p + g + 1, 10p - d_g - 1]")
    print("cover every integer 50000 <= n < 2000000.\n")

    print("Construction phase")
    print("------------------")
    bernoulli_cache: dict[int, fmpq] = {}
    genus_values = [genus_data(g, bernoulli_cache) for g in range(1, G_MAX + 1)]
    require(
        all(0 <= d_g <= 1_996 for _, d_g in genus_values),
        "computed d_g lies outside the range stated in the manuscript",
    )
    print("Computed every A_g and d_g exactly from equation (4.21).")

    # No interval can need a prime larger than this during the greedy search.
    candidate_limit = (N_MAX - 1 - 1) // 9 + 1
    candidates = primes_up_to(candidate_limit)
    print(
        f"Generated {len(candidates):,} sieve candidates in "
        f"({PRIME_LOWER_BOUND}, {candidate_limit}]."
    )
    all_lists = construct_witnesses(genus_values, candidates)
    print("Constructed all 2,000 witness lists by the stated greedy rule.\n")

    print("Independent verification phase")
    print("------------------------------")
    summary = verify_witnesses(genus_values, all_lists)
    distinct_primes = summary["distinct_primes"]
    assert isinstance(distinct_primes, list)
    print(
        "(i) Primality and p > 4012: verified by trial division for "
        f"{len(distinct_primes):,} distinct recorded primes "
        f"({distinct_primes[0]} through {distinct_primes[-1]})."
    )
    print(
        "(ii) Nondivisibility: verified A_g mod p != 0 for all "
        f"{summary['total_occurrences']:,} (g,p) occurrences."
    )
    print(
        "(iii) Coverage: merged the closed intervals separately for each genus; "
        "all cover [50000, 1999999]."
    )

    print("\nDiagnostics")
    print("-----------")
    print(
        f"Smallest |P_g|: {summary['smallest_list_size']} "
        "(g = "
        f"{abbreviated_genus_list(summary['smallest_list_genera'])})"
    )
    print(
        f"Largest  |P_g|: {summary['largest_list_size']} "
        "(g = "
        f"{abbreviated_genus_list(summary['largest_list_genera'])})"
    )
    join_slack, join_g, join_p = summary["tightest_join"]
    print(
        "Tightest consecutive join: slack "
        f"{join_slack} integer(s) (g={join_g}, next p={join_p}); "
        "slack 0 means exactly adjacent intervals."
    )
    initial_excess, initial_g = summary["smallest_initial_excess"]
    terminal_excess, terminal_g = summary["smallest_terminal_excess"]
    print(
        f"Smallest initial excess: {initial_excess} integer(s) (g={initial_g})."
    )
    print(
        f"Smallest terminal excess: {terminal_excess} integer(s) (g={terminal_g})."
    )
    print("\nPASS: all three claims in Finite arithmetic verification 4.14 hold.")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"\nFAIL: {error}")
        raise
