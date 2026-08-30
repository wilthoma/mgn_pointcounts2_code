#!/usr/bin/env python3
"""Exact D1 prime-window checker for the infinite n-tail.

SymPy is used only to obtain exact Bernoulli numerators.  All interval merging,
prime sieving, and the final logarithmic bounds use exact integer/Fraction
arithmetic.  The mathematical input is Dusart's explicit estimates
  x/(log x-1) <= pi(x) for x>5393,
  pi(x) <= x/(log x-1.1) for x>60184.
"""
from __future__ import annotations
import argparse
from bisect import bisect_left
from fractions import Fraction
from math import isqrt
from multiprocessing import Pool
import sympy as sp

def h_of_g(g:int)->int:
    if g<3: raise ValueError(g)
    return (g-1)//2 if g&1 else (g-2)//2

def rho(g:int)->int:
    return -2 if g<=2 else 2*h_of_g(g)-2

def sieve_primes(N:int)->list[int]:
    a=bytearray(b'\x01')*(N+1);a[0:2]=b'\x00\x00'
    for p in range(2,isqrt(N)+1):
        if a[p]:a[p*p:N+1:p]=b'\x00'*(((N-p*p)//p)+1)
    return [i for i,v in enumerate(a) if v]

def log2_interval(terms=70):
    q=Fraction(1,3)
    s=sum((2*q**(2*j+1)/Fraction(2*j+1) for j in range(terms+1)),Fraction())
    r=2*q**(2*terms+3)/Fraction(2*terms+3)/(1-q*q)
    return s-r,s+r
LOG2=log2_interval()

def log_interval(x:Fraction,terms=55):
    if x<=0:raise ValueError(x)
    k=x.numerator.bit_length()-x.denominator.bit_length()
    twok=Fraction(1<<k) if k>=0 else Fraction(1,1<<(-k))
    y=x/twok
    while y<1:k-=1;y*=2
    while y>=2:k+=1;y/=2
    q=(y-1)/(y+1)
    s=sum((2*q**(2*j+1)/Fraction(2*j+1) for j in range(terms+1)),Fraction())
    r=2*abs(q)**(2*terms+3)/Fraction(2*terms+3)/(1-q*q)
    if k>=0:return s-r+k*LOG2[0],s+r+k*LOG2[1]
    return s-r+k*LOG2[1],s+r+k*LOG2[0]

def genus_scan(task):
    g,Gamma,start,end,primes,num=task;cur=start;r=rho(g)
    i=bisect_left(primes,max(2*g+2,(cur+r+Gamma-1)//Gamma))
    for p in primes[i:]:
        if g>=3 and num%p==0:continue
        lo=(Gamma-1)*p+g+1;hi=Gamma*p-r-1
        if hi<cur:continue
        if lo>cur:return g,cur,lo-1
        cur=max(cur,hi+1)
        if cur>end:return g,None,None
    return g,cur,end

def possible_bad_count(nums:dict[int,int],base:int):
    best=(0,None)
    for g,num in nums.items():
        c=0;q=1
        while q*base<=num:q*=base;c+=1
        if c>best[0]:best=(c,g)
    return best

def infinite_check(G:int,Gamma:int,N:int,nums:dict[int,int]):
    rmax=max(rho(g) for g in range(1,G+1))
    U=Fraction(N-G,Gamma-1);L=Fraction(N+rmax,Gamma)
    assert U>L and U>60184 and L>60184
    lu=log_interval(U);ll=log_interval(L)
    count_lb=U/(lu[1]-1)-L/(ll[0]-Fraction(11,10))
    base=L.numerator//L.denominator+1
    bad,g=possible_bad_count(nums,base)
    assert count_lb>bad+2,(count_lb,bad,g)
    # Monotonicity of the count lower bound.  For Gamma=10 use
    # log(10/9)<53/500; for Gamma=14 use log(14/13)<3/40.
    if Gamma==10:
        P_at=2_500_000-14_870_000+31_954_190-19_529_599
        disc=(-29_740_000)**2-4*7_500_000*31_954_190
        assert P_at>0 and disc<0
    elif Gamma==14:
        # P(z)=3200z^3-20480z^2+46726z-29603 is increasing,
        # and P(11/10)>0.
        P_at=3200*1331-20480*1210+46726*1100-29603*1000
        disc=(-40_960)**2-4*9_600*46_726
        assert P_at>0 and disc<0
    else:
        raise ValueError(Gamma)
    return count_lb,bad,g,base,rmax

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--G',type=int,required=True);ap.add_argument('--Gamma',type=int,choices=(10,14),default=10);ap.add_argument('--tail-start',type=int,default=50000);ap.add_argument('--scan-end',type=int,default=2_000_000);ap.add_argument('--infinite-start',type=int,default=2_000_000);ap.add_argument('--workers',type=int,default=1);ap.add_argument('--finite-only',action='store_true');a=ap.parse_args()
    nums={g:abs(int(sp.numer(sp.bernoulli(2*h_of_g(g))))) for g in range(3,a.G+1)}
    pmax=(a.scan_end+max(rho(g) for g in range(1,a.G+1)))//(a.Gamma-1)+20
    primes=sieve_primes(pmax)
    tasks=[(g,a.Gamma,a.tail_start,a.scan_end,primes,nums.get(g,1)) for g in range(1,a.G+1)]
    if a.workers==1:results=list(map(genus_scan,tasks))
    else:
        with Pool(a.workers) as pool:results=list(pool.imap_unordered(genus_scan,tasks,chunksize=1))
    gaps=[x for x in results if x[1] is not None]
    if gaps:
        print('FAIL finite gaps',sorted(gaps)[:20]);raise SystemExit(1)
    print(f'finite PASS: Gamma={a.Gamma}, g<= {a.G}, n in [{a.tail_start},{a.scan_end}], pmax={pmax}, primes={len(primes)}')
    if not a.finite_only:
        lb,bad,g,base,rmax=infinite_check(a.G,a.Gamma,a.infinite_start,nums)
        print(f'infinite PASS from n={a.infinite_start}: rho_max={rmax}, common-window minimum prime={base}')
        print(f'exact prime-count lower bound={lb}')
        print(f'max possible Bernoulli-numerator exclusions={bad}, at g={g}')
if __name__=='__main__':main()
