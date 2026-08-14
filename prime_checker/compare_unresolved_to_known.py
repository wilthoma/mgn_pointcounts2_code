#!/usr/bin/env python3
"""Compare check_modular_cover output with an allow-list of known zero pairs.

Both TSV files must have a header containing the two columns `g` and `n`.
The known-zero list is filtered to the requested rectangle before comparison.
Exit status is 0 iff the sets agree exactly.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def load_pairs(path: Path) -> set[tuple[int, int]]:
    with path.open() as f:
        header = next(f, "").split()
        if header != ["g", "n"]:
            raise ValueError(f"{path}: expected header 'g\\tn', got {header!r}")
        pairs: set[tuple[int, int]] = set()
        for lineno, line in enumerate(f, 2):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 2:
                raise ValueError(f"{path}:{lineno}: expected two integers")
            g, n = map(int, fields)
            pairs.add((g, n))
        return pairs


def write_pairs(path: Path, pairs: set[tuple[int, int]]) -> None:
    with path.open("w") as f:
        f.write("g\tn\n")
        for g, n in sorted(pairs, key=lambda x: (x[1], x[0])):
            f.write(f"{g}\t{n}\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("unresolved", type=Path)
    ap.add_argument("known_zeros", type=Path)
    ap.add_argument("--G", type=int, required=True)
    ap.add_argument("--n0", type=int, required=True)
    ap.add_argument("--n1", type=int, required=True)
    ap.add_argument("--unexpected-out", type=Path, default=Path("unexpected_unresolved.tsv"))
    ap.add_argument("--missing-out", type=Path, default=Path("known_zeros_not_unresolved.tsv"))
    args = ap.parse_args()

    unresolved = {
        pair for pair in load_pairs(args.unresolved)
        if 1 <= pair[0] <= args.G and args.n0 <= pair[1] <= args.n1
    }
    known = {
        pair for pair in load_pairs(args.known_zeros)
        if 1 <= pair[0] <= args.G and args.n0 <= pair[1] <= args.n1
    }

    unexpected = unresolved - known
    missing = known - unresolved
    write_pairs(args.unexpected_out, unexpected)
    write_pairs(args.missing_out, missing)

    print(f"unresolved in rectangle: {len(unresolved)}")
    print(f"known zeros in rectangle: {len(known)}")
    print(f"unexpected unresolved:    {len(unexpected)} -> {args.unexpected_out}")
    print(f"known zeros covered:      {len(missing)} -> {args.missing_out}")
    if unexpected or missing:
        return 1
    print("PASS: unresolved set agrees exactly with the known-zero list")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
