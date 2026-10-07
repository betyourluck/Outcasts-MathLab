import json
from collections import Counter
from math import gcd
from failure_atlas import jacobi, factor
rows=json.load(open("failure_atlas_rows.json"))["primes"]
def sqfree_divs(n):
    ps=[q for q in factor(n)]
    out=[1]
    for q in ps: out+= [d*q for d in out]
    return out
cnt=Counter(); unexplained=0; n=0
for r in rows:
    if r["ok"] or r["type"]!="QR" or r["psi_kernel"]: continue
    n+=1
    a,u,p=r["a"],r["u"],r["n"]; m0=4*a*u
    qs=[q for q in factor(r["M"]) if q%2==1 and gcd(q,m0)==1]
    ds=[-d for d in sqfree_divs(2*a*u)]          # negative d  <=> chi(-1) = -1
    hits=[d for d in ds if all(jacobi(d % q, q)==1 for q in qs)]
    if not hits: unexplained+=1
    for d in hits: cnt[d]+=1
print("QR-not-psi rows",n,"unexplained by Kronecker(d|.) with d|2au, d<0:",unexplained)
print("obstructing d (top):",cnt.most_common(12))
for p in (345601,670849): print(p, p%840)
