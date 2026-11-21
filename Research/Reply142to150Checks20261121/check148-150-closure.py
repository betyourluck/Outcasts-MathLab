# Contract-free addition (2026-11-21), NOT a claim of >>148-150: on the 288 inputs, how the hypothesis and the conclusion of the
# target lemma of >>149 are distributed.  C_p = closure of the least non-residue prime l under all edges (all admissible sources).
# lower = 3c < p, upper = 3c > p; "succeeds" = Type I or Type II witness among the divisors of x^2.
import sys, time
from collections import Counter
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
HARD={1,121,169,289,361,529}
P17=[30241,61681,106681,111721,148201,153409,163249,206641,216841,280009,288649,307129,309241,325081,347041,352081,374929]
def A(q): return q if q%4==3 else 3*q
def nr(q,p): return pow(q,(p-1)//2,p)==p-1
def dsq(f):
    out=[1]
    for r,e in f.items(): out=[a*r**k for a in out for k in range(2*e+1)]
    return out
inputs=sorted(set(p for p in sympy.primerange(1009,100000) if p%840 in HARD)|{2521,12889}|set(P17))
cnt=Counter(); rows=[]; tot=Counter()
for p in inputs:
    l=3
    while not nr(l,p): l=sympy.nextprime(l)
    V={}
    for q in sympy.primerange(3,p):
        if not nr(q,p): continue
        At=A(q); t=1; kids=set(); lo=0; up=0; upsm=0
        while At*t*t<3*p:
            c=At*t*t; x=(p+c)//4; f=sympy.factorint(x)
            kids|={r for r in f if nr(r,p)}
            ok=any((d+p*x)%c==0 or (d+x)%c==0 for d in dsq(f))
            if ok:
                if 3*c<p: lo+=1
                else:
                    up+=1
                    if max(f)*l<=p: upsm+=1
            t+=2
        V[q]=(kids,lo,up,upsm)
    seen={l}; fr=[l]
    while fr:
        nx=[]
        for q in fr:
            for r in V[q][0]:
                if r not in seen: seen.add(r); nx.append(r)
        fr=nx
    lo=sum(V[q][1] for q in seen); up=sum(V[q][2] for q in seen); upsm=sum(V[q][3] for q in seen)
    tot['successful lower sources in C_p']+=lo; tot['successful upper sources in C_p']+=up
    tot['successful upper sources, whole graph']+=sum(v[2] for v in V.values()); tot['successful lower sources, whole graph']+=sum(v[1] for v in V.values())
    tot['vertices in C_p']+=len(seen); tot['vertices, whole graph']+=len(V)
    cnt['inputs']+=1
    cnt['C_p has a successful lower source']+=lo>0
    cnt['C_p has no successful lower source (hypothesis of the target lemma)']+=lo==0
    cnt['C_p has a successful upper source']+=up>0
    cnt['C_p has a successful upper source, all of them with every prime factor of x <= p/l']+=(up>0 and up==upsm)
    cnt['C_p has no successful upper source']+=up==0
    cnt['C_p has no successful source at all']+=(lo==0 and up==0)
    if lo==0: rows.append((p,l,len(seen),lo,up))
    if lo<=2: rows.append(('few lower',p,l,len(seen),lo,up))
for k in cnt: print('%-82s %d'%(k,cnt[k]))
for k in tot: print('%-82s %d'%(k,tot[k]))
print('inputs with at most 2 successful lower sources in C_p:',rows)
print('%.1f s'%(time.time()-t0))
