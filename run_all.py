#!/usr/bin/env python3
"""Run the root-directory verifications in manuscript order.

Each checker remains an ordinary, independently runnable Python script.  This
small driver runs those six checkers in one command.  Computer Verification
4.17 is included in the finite_cover subdirectory and is run separately with
``python3 finite_cover/finite_cover.py``.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from time import perf_counter


CHECKERS = (
    "ag_bounds.py",
    "prime_cover.py",
    "error_finite.py",
    "interval_check_m.py",
    "llh_bounds.py",
    "lambda_t_bounds.py",
)


def main() -> int:
    repository = Path(__file__).resolve().parent
    started = perf_counter()

    print("Referee certificate suite")
    print("=========================")
    print("Running one no-argument checker per included verification box.\n")

    for index, checker in enumerate(CHECKERS, start=1):
        path = repository / checker
        if not path.is_file():
            print(f"ERROR: missing checker: {checker}", flush=True)
            return 2

        print(f"[{index}/{len(CHECKERS)}] {checker}", flush=True)
        completed = subprocess.run([sys.executable, str(path)], cwd=repository)
        if completed.returncode != 0:
            print(f"\nFAILED: {checker} exited with status {completed.returncode}.")
            return completed.returncode
        print()

    elapsed = perf_counter() - started
    print(f"ALL INCLUDED VERIFICATIONS PASSED ({elapsed:.1f} seconds total).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
