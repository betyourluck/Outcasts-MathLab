# Contract-free quick check (2026-11-16) of >>144 / PR #4 7956611 (LargePrimeSource20261116).
# Written after reading >>144. The reduced j = 1 certificate is recomputed by a different route:
#   * tuples are generated from the divisors f | s^2 (not from coprime alpha, beta with alpha*beta | s) and reduced;
#   * the congruence stage is decided prime by prime (CRT): a tuple dies iff at some prime l in {2,3,5,7} or l | h*sigma
#     no residue of p modulo the local part of L survives; then compared with the full-period count in the PR JSON.
# The j = 2 paper arguments are checked symbolically on all (s, sigma, f) and numerically on random hard primes.
import sys, json, random, time
from math import gcd
from collections import Counter
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
HARD=(1,121,169,289,361,529)
nok=0; ndiff=0
def chk(name,got,exp):
    global nok,ndiff
    if got==exp: nok+=1; print('OK  ',name,got if len(str(got))<170 else '(long)')
    else: ndiff+=1; print('DIFF',name,'got',got,'expected',exp)

# ---- 1. tuples from f | s^2 ----
tuples=set()
for s in range(1,9):
    for f in sympy.divisors(s*s):
        g=gcd(f,s); a,b=f//g,s//g
        assert s%(a*b)==0 and gcd(a,b)==1          # >>144: alpha*beta | s
        for sigma in (1,3):
            for h in range(1,3*a+1): tuples.add((s,sigma,a,b,h))
# conversely every coprime (a, b) with a*b | s comes from f = s*a/b
conv=all(any(gcd(f,s)and f//gcd(f,s)==a and s//gcd(f,s)==b for f in sympy.divisors(s*s))
         for s in range(1,9) for a in range(1,9) for b in range(1,9) if s%(a*b)==0 and gcd(a,b)==1)
chk('tuples (s, sigma, alpha, beta, h) from f | s^2, and the converse map',(len(tuples),conv),(438,True))
def size_stage(s,sigma,a,b,h):
    if a>=sigma*h: return 'q_range'
    if 4*h*s>9*(h+a): return 'r_range'
    return None
stage=Counter(size_stage(*t) or 'periodic' for t in tuples)
chk('size stages (q range, r range, left for congruences)',(stage['q_range'],stage['r_range'],stage['periodic']),(85,286,67))

# ---- 2. congruence stage, prime by prime ----
def p_ok(l,pr):
    return {2:pr%8==1,3:pr%3==1,5:pr%5 in (1,4),7:pr%7 in (1,2,4)}.get(l,pr%l!=0)
def local_survives(s,sigma,a,b,h,l):
    # q = (a p + b)/(h sigma), r = ((h + a) p + b)/(4 h s); local conditions at l:
    # integrality, q = 3 or 1 mod 4 (l = 2), l does not divide q or r for l in {2,3,5,7}
    eq=sympy.multiplicity(l,h*sigma); er=sympy.multiplicity(l,4*h*s)
    need=max(eq,er)+(3 if l==2 else 1)
    mod=l**need
    for pr in range(mod):
        if not p_ok(l,pr%mod): continue
        qn=(a*pr+b)%mod; rn=((h+a)*pr+b)%mod
        if qn%l**eq or rn%l**er: continue
        qr=(qn//l**eq)%l**(need-eq); rr=(rn//l**er)%l**(need-er)
        if l==2:
            qv=qr*pow((h*sigma)//2**eq,-1,2**(need-eq))%2**(need-eq)
            if need-eq<2: raise SystemExit('precision')
            if qv%4!=(3 if sigma==1 else 1): continue
            if rr%2==0: continue
        elif l in (3,5,7):
            if qr%l==0 or rr%l==0: continue
        return True
    return False
dead_at=Counter(); per_tuple={}
for t in sorted(tuples):
    if size_stage(*t): continue
    s,sigma,a,b,h=t
    primes=sorted(set([2,3,5,7])|set(sympy.primefactors(h*sigma*s)))
    d=next((l for l in primes if not local_survives(s,sigma,a,b,h,l)),None)
    per_tuple[t]=d; dead_at[d]+=1
chk('congruence stage: every one of the 67 tuples dies at some prime (CRT, no full period)',dead_at.get(None,0),0)
print('     prime at which each tuple dies (first in the order 2,3,5,7,...):',sorted((k,v) for k,v in dead_at.items()))

# ---- 3. comparison with the PR results JSON ----
if len(sys.argv)>1:
    R=json.load(open(sys.argv[1],encoding='utf-8'))
    rows=R['rows']; summ=R['summary']
    theirs={(x['s'],x['sigma'],x['alpha'],x['beta'],x['h']):x for x in rows}
    chk('PR JSON: same 438 tuples',set(theirs)==tuples,True)
    chk('PR JSON: outcome of every tuple (q_range / r_range / congruence)',all(theirs[t]['outcome']==(size_stage(*t) or 'congruence') for t in tuples),True)
    chk('PR JSON: summary outcomes and unresolved',(summ['outcomes'],summ['unresolved']),({'q_range':85,'r_range':286,'congruence':67},[]))
    chk('PR JSON: periods L = 840 h s sigma and residues visited (6 per 840)',
        (all(theirs[t]['period']==840*t[4]*t[0]*t[1] for t in per_tuple),sum(theirs[t]['counts']['p_classes'] for t in per_tuple)),(True,sum(6*t[4]*t[0]*t[1] for t in per_tuple)))
    print('     residues visited in the full-period check:',sum(theirs[t]['counts']['p_classes'] for t in per_tuple))

# ---- 4. j = 2 paper arguments (>>144) ----
# sigma = 1: every prime l | s is a square mod q, so (-f/q) = -1 for every f | s^2.
# sigma = 3: 3 does not divide s and the witness needs f = 2 mod 3; s in {1,7}: impossible; s a power of 2: f an odd power of 2
#            and q = 5 mod 8 gives (-2/q) = -1; s = 5: f = 5 and (-5/q) = -1.
ok=True; rows4=[]
for s in range(1,9):
    for f in sympy.divisors(s*s):
        # sigma = 3 row: is f = 2 mod 3 possible with 3 not dividing s?
        if s%3 and (f+s*s)%3==0:
            rows4.append((s,f))
print('     sigma = 3, j = 2: (s, f) with 3 not dividing s and 3 | f + s^2:',rows4)
chk('sigma = 3, j = 2: the surviving (s, f) are the odd powers of 2 for s in {2,4,8} and f = 5 for s = 5',
    sorted(rows4),sorted([(2,2),(4,2),(4,8),(8,2),(8,8),(8,32),(5,5)]))
# numerical check of the two character facts on random hard primes and single-source q with x = s r, r > p/9 prime
random.seed(144)
hard=[p for p in sympy.primerange(10**6,3*10**6) if p%840 in HARD]
cnt=Counter(); bad=0
for p in random.sample(hard,60):
    for q in sympy.primerange(p//9,p):
        if pow(q,(p-1)//2,p)!=p-1: continue
        sg=1 if q%4==3 else 3
        if 3*sg*q<p: continue
        x=(p+sg*q)//4
        for s in range(1,9):
            if x%s==0 and sympy.isprime(x//s) and 9*(x//s)>p:
                for f in sympy.divisors(s*s):
                    if sg==1:
                        cnt['sigma1']+=1
                        if sympy.jacobi_symbol(-f%q,q)!=-1: bad+=1
                    elif (f+s*s)%3==0:
                        cnt['sigma3 f=2 mod 3']+=1
                        if sympy.jacobi_symbol(-f%q,q)!=-1: bad+=1
                break
print('     random sample: 60 hard primes in (10^6, 3*10^6); (vertex, f) pairs checked:',dict(cnt))
chk('j = 2: (-f/q) = -1 on every checked pair',bad,0)
print(nok,'OK /',ndiff,'DIFF,',round(time.time()-t0,1),'s')
