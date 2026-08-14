#!/usr/bin/env python3
from __future__ import annotations
import argparse
POOLS={1000:[1009,1013,1019,1021,2003,2011,2017,2027,5003,5009,5011,5021],2000:[2003,2011,2017,2027,5003,5009,5011,5021]}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--G',type=int,choices=(1000,2000),required=True);ap.add_argument('--Gamma',type=int,choices=(10,14),required=True);ap.add_argument('--mode',choices=('n-ge-5000','n-lt-5000'),default='n-ge-5000');ap.add_argument('--exe',default='./modular_low_truncation_table_large');a=ap.parse_args()
 if a.mode=='n-ge-5000':p0=min(POOLS[a.G]);W=a.Gamma-5000//p0;ps=POOLS[a.G]
 else:W=a.Gamma;ps=[5003,5009,5011,5021]
 print(f'# G={a.G} Gamma={a.Gamma} W={W} mode={a.mode}\nmkdir -p tables logs')
 for p in ps:print(f'{a.exe} {a.G} {W} {p-1} {p} tables/G{a.G}_W{W}_p{p}.bin >logs/G{a.G}_W{W}_p{p}.out 2>logs/G{a.G}_W{W}_p{p}.err')
if __name__=='__main__':main()
