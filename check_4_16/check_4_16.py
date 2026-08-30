#!/usr/bin/env python3
"""Compile and run the experimental large-prime checker for box 4.16.

Running this file without arguments performs the full computation.  Compiler
and program output are shown immediately and also saved in
check_4_16_large_prime.log.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "check_4_16_large_prime.cpp"
BUILD = ROOT / "build"
BINARY = BUILD / "check_4_16_large_prime"
LOG = ROOT / "check_4_16_large_prime.log"


def find_compiler() -> str:
    requested = os.environ.get("CXX")
    if requested:
        path = shutil.which(requested)
        if path:
            return path
        raise SystemExit(f"CXX={requested!r}, but that compiler was not found")
    for candidate in ("c++", "g++", "clang++"):
        path = shutil.which(candidate)
        if path:
            return path
    raise SystemExit("No C++ compiler found (tried c++, g++, and clang++)")


def main() -> int:
    compiler = find_compiler()
    BUILD.mkdir(exist_ok=True)
    command = [
        compiler,
        "-O3",
        "-DNDEBUG",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        str(SOURCE),
        "-o",
        str(BINARY),
    ]
    print("Compiling the large-prime checker:", flush=True)
    print("  " + " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)

    print("\nRunning the built-in self-tests ...", flush=True)
    subprocess.run([str(BINARY), "--self-test"], cwd=ROOT, check=True)

    print("\nComparing a small blocked run with the archived generator ...", flush=True)
    subprocess.run([str(BINARY), "--smoke"], cwd=ROOT, check=True)

    print("\nStarting the full computation.", flush=True)
    print(f"A complete transcript will be saved to {LOG.name}.\n", flush=True)
    started = time.monotonic()
    with LOG.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            [str(BINARY)],
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
        status = process.wait()

    elapsed = time.monotonic() - started
    print(f"\nRunner elapsed time: {elapsed:.1f} seconds", flush=True)
    if status != 0:
        print(f"The checker exited with status {status}.", file=sys.stderr)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
