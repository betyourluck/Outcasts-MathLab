# Contract-free (2026-11-16): the j = 1 case by a Legendre symbol, and where successful sources live.
# Paper argument (any admissible tau): c = sigma q tau^2 with (q/p) = -1, so (c/p) = -1 (3 and tau^2 are squares mod p).
# j = 1 gives h c = alpha p + beta, so (h/p) (c/p) = (beta/p) = 1 (beta | s <= 8, and 2,3,5,7 are squares mod a hard p).
# Hence (h/p) = -1 and h has an odd power of a prime >= 11; with h <= 24 this means h >= 11.  But h <= 3 alpha <= 3 s and
# r > p/9 gives 4hs <= 9(h + alpha) <= 9(h + s), so s <= 9h/(4h - 9) < 3 for h >= 11, i.e. h <= 3s <= 6: contradiction.
# Part 1 checks the bookkeeping against the PR JSON; part 2 classifies every successful source of the 288 inputs.
import sys, json, time
from collections import Counter
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
R=json.load(open(sys.argv[1],encoding='utf-8'))
rows=R['rows']
cong=[x for x in rows if x['outcome']=='congruence']
big=[x for x in rows if any(l>=11 for l in sympy.primefactors(x['h']))]
print('PR tuples left for the congruence stage:',len(cong),'; their h values:',sorted(set(x['h'] for x in cong)),
      '; all with prime factors <= 7:',all(max(sympy.primefactors(x['h']) or [1])<=7 for x in cong))
print('PR tuples whose h has a prime factor >= 11:',len(big),'; outcomes:',dict(Counter(x['outcome'] for x in big)))
# exhaustive check of the inequality chain over all (s, alpha, beta, h) with alpha*beta | s, h <= 3 alpha, h >= 11
viol=[(s,a,b,h) for s in range(1,9) for a in range(1,9) for b in range(1,9) if s%(a*b)==0 and sympy.gcd(a,b)==1
      for h in range(11,3*a+1) if 4*h*s<=9*(h+a)]
print('(s, alpha, beta, h) with h >= 11, h <= 3 alpha and 4hs <= 9(h + alpha):',viol)
# ---- where successes live ----
HARD={1,121,169,289,361,529}
P17=[30241,61681,106681,111721,148201,153409,163249,206641,216841,280009,288649,307129,309241,325081,347041,352081,374929]
inputs=sorted(set(p for p in sympy.primerange(1009,100000) if p%840 in HARD)|{2521,12889}|set(P17))
def sig(r): return 1 if r%4==3 else 3
def dsq(x):
    out=[1]
    for r,e in sympy.factorint(x).items(): out=[a*r**k for a in out for k in range(2*e+1)]
    return out
src=Counter(); suc=Counter(); suc_branch=Counter()
for p in inputs:
    for q in sympy.primerange(3,p):
        if pow(q,(p-1)//2,p)!=p-1: continue
        A=sig(q)*q; t=1
        while A*t*t<3*p:
            c=A*t*t; x=(p+c)//4
            if 3*c<p: cls='lower (3c < p)'
            else:
                large=any(x%s==0 and sympy.isprime(x//s) and 9*(x//s)>p for s in range(1,9))
                cls='upper, x has a prime factor > p/9' if large else 'upper, x is p/9-smooth'
            src[cls]+=1
            ds=dsq(x); i1=any((d+p*x)%c==0 for d in ds); i2=any((d+x)%c==0 for d in ds)
            if i1 or i2:
                suc[cls]+=1; suc_branch[(cls,'I' if i1 else '-','II' if i2 else '-')]+=1
            t+=2
print('inputs:',len(inputs))
for k in sorted(src): print('   %-36s sources %8d   successful %6d'%(k,src[k],suc[k]))
print('   successful sources by branch:',sorted(suc_branch.items()))
print('%.1f s'%(time.time()-t0))
