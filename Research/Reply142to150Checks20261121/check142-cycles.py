# Contract-free (2026-11-15): cycles of the map F on single-source vertices (NOT in >>142).
# U = non-residue primes q < p with 3 sigma(q) q > p (only tau = 1) whose x = (p + sigma q)/4 has a prime factor r > p/9;
# F(q) = r. A nonempty closed set made only of single-source vertices lies in U and contains an F-cycle (it need not be a
# union of cycles: in-trees can hang on a cycle -- corrected after >>144, failures #42); every vertex of U fails both branches.
# >>142 searched the F-cycles with x = r (cofactor s = 1) on every edge, all p, lengths 6..18.
# Inputs here: every hard prime 1009 <= p < 10^5 (288 - 19 of check142-quick part C) plus 2521, 12889 and the 17 primes.
import sys, time
from collections import Counter
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
HARD={1,121,169,289,361,529}
P17=[30241,61681,106681,111721,148201,153409,163249,206641,216841,280009,288649,307129,309241,325081,347041,352081,374929]
inputs=sorted(set(p for p in sympy.primerange(1009,100000) if p%840 in HARD)|{2521,12889}|set(P17))
def sig(r): return 1 if r%4==3 else 3
def dsq(x):
    out=[1]
    for r,e in sympy.factorint(x).items(): out=[a*r**k for a in out for k in range(2*e+1)]
    return out
tot=Counter(); edge_s=Counter(); maxchain=(0,None); cycles=[]; checked_fail=0
for p in inputs:
    F={}; S={}
    for q in sympy.primerange(p//9,p):
        if pow(q,(p-1)//2,p)!=p-1 or 3*sig(q)*q<p: continue
        x=(p+sig(q)*q)//4
        for s in range(1,9):
            if x%s==0 and sympy.isprime(x//s) and 9*(x//s)>p:
                F[q]=x//s; S[q]=s; break
        tot['single-source vertices']+=1
    tot['U']+=len(F)
    inU={q:F[q] for q in F if F[q] in F}           # edges that stay in U
    tot['edges staying in U']+=len(inU)
    for q in inU: edge_s[S[q]]+=1
    # longest F-path inside U and cycle detection (functional graph)
    state={}; depth={}
    for q0 in F:
        if q0 in state: continue
        path=[]; q=q0
        while q in F and q not in state:
            state[q]='open'; path.append(q); q=F[q]
        if q in state and state[q]=='open':             # new cycle
            cyc=path[path.index(q):]
            cycles.append((p,cyc))
            for v in cyc: depth[v]=len(cyc)
        base=depth.get(q,0) if q in F else 0
        for v in reversed(path):
            if v in depth: base=depth[v]; continue
            base+=1; depth[v]=base
        for v in path: state[v]='done'
    m=max(depth.values()) if depth else 0
    if m>maxchain[0]: maxchain=(m,p)
    # direct confirmation of the lemma on every vertex of U for the smaller inputs
    if p<20000:
        for q in F:
            c=sig(q)*q; x=(p+c)//4
            assert not any((d+p*x)%c==0 or (d+x)%c==0 for d in dsq(x)); checked_fail+=1
print('inputs:',len(inputs),'(hard primes 1009 <= p < 10^5 plus 2521, 12889 and the 17 primes)')
print('single-source vertices:',tot['single-source vertices'],'; in U (x has a prime factor > p/9):',tot['U'],'; F-edges with both ends in U:',tot['edges staying in U'])
print('cofactor s on the F-edges inside U:',sorted(edge_s.items()))
print('F-cycles inside U:',len(cycles),cycles[:5])
print('longest F-path inside U (number of vertices):',maxchain)
print('vertices of U checked directly over all divisors (p < 20000), all failing:',checked_fail)
print('%.1f s'%(time.time()-t0))
