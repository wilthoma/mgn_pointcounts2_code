#!/usr/bin/env python3
"""Run every included verification in manuscript order.

Each checker remains an ordinary, independently runnable Python script.  This
small driver is only a convenience for a referee who wants to run the complete
suite in one command.  Finite arithmetic verification 4.16 is deliberately not
part of this repository.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from time import perf_counter


CHECKERS = (
    "check_4_12.py",
    "check_4_14.py",
    "check_5_5.py",
    "check_5_9.py",
    "check_6_2.py",
    "check_6_6.py",
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
