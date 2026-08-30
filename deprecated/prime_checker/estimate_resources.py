#!/usr/bin/env python3
from __future__ import annotations
import argparse, math

def next_pow2_gt(x:int)->int:
    n=1
    while n<=x:n<<=1
    return n

def sizes(G:int,W:int,S:int):
    nr=next_pow2_gt(2*(G*(2*G+1)+G))
    base=S+G+1
    na=next_pow2_gt(((G-1)*base+G)+((G-1)*base+S))
    tri=(G+1)*(G+2)//2; out=G*(S+1)
    ints_ar=3*(W+1)*out+(W+1)*tri+(2*(W+1)+1)*na
    ints_r=(W+1)*nr+4*W*nr+4*nr+(W+1)*tri
    mem=4*max(ints_ar,ints_r)
    work=6*W*nr*int(math.log2(nr))+6*(W+1)*na*int(math.log2(na))
    return nr,na,mem,work,out
_,_,_,CAL_WORK,_=sizes(600,6,1600)
CAL_SEC=25.8

def hb(n):
    x=float(n)
    for u in ['B','KiB','MiB','GiB','TiB']:
        if x<1024 or u=='TiB':return f'{x:.2f} {u}'
        x/=1024

def main():
    ap=argparse.ArgumentParser();ap.add_argument('G',type=int);ap.add_argument('W',type=int);ap.add_argument('S',type=int)
    ap.add_argument('--jobs',type=int,default=1);ap.add_argument('--speed-factor',type=float,default=1.0)
    a=ap.parse_args();nr,na,mem,work,out=sizes(a.G,a.W,a.S)
    sec=CAL_SEC*(work/CAL_WORK)/a.speed_factor
    print(f'G={a.G} W={a.W} S={a.S}')
    print(f'R NTT:  {nr}=2^{int(math.log2(nr))}')
    print(f'AR NTT: {na}=2^{int(math.log2(na))}')
    print(f'peak/job estimate: {hb(mem)}')
    print(f'full table file:   {hb(16+2*(a.W+1)*out)}')
    print(f'time/job estimate on reference sandbox: {sec/60:.2f} min')
    print(f'{a.jobs} jobs aggregate memory model: {hb(mem*a.jobs)}')
    print('Benchmark one job; NTT workloads are memory-bandwidth sensitive.')
if __name__=='__main__':main()
