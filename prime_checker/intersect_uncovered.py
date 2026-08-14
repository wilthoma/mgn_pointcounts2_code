#!/usr/bin/env python3
"""Intersect two TSV files produced by check_modular_cover."""
from __future__ import annotations
import argparse

def load(path):
    with open(path) as f:
        head=next(f).split();assert head==['g','n']
        return {tuple(map(int,line.split())) for line in f if line.strip()}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('first');ap.add_argument('second');ap.add_argument('output');a=ap.parse_args()
    z=sorted(load(a.first)&load(a.second),key=lambda x:(x[1],x[0]))
    with open(a.output,'w') as f:
        f.write('g\tn\n');f.writelines(f'{g}\t{n}\n' for g,n in z)
    print(f'intersection size={len(z)}; output={a.output}')
if __name__=='__main__':main()
