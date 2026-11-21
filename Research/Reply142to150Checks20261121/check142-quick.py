# Contract-free quick check (2026-11-15) of >>142 / PR #4 8399e66 (PrimeSourceCycles20261115).
# Written after reading >>142. The enumeration is redone by a different method (Duval's Lyndon-word generator and
# exact Fraction fixed points, P0 = lcm of the denominators); the PR results JSON (argv[1]) is read only for the
# comparison. Parts C and D are NOT in >>142.
# Units: "word" = primitive cyclic sigma word (one rotation); "vertex" = non-residue prime q < p.
import sys, json, math, time
from collections import Counter
from fractions import Fraction
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
HARD={1,121,169,289,361,529}
nok=0; ndiff=0
def chk(name,got,exp):
    global nok,ndiff
    if got==exp: nok+=1; print('OK  ',name,got if len(str(got))<170 else '(long)')
    else: ndiff+=1; print('DIFF',name,'got',got,'expected',exp)
def lyndon(k_alphabet,n):
    # Duval / FKM: all Lyndon words of length exactly n over range(k)
    w=[-1]; out=[]
    while w:
        w[-1]+=1
        m=len(w)
        if m==n: out.append(tuple(w))
        while len(w)<n: w.append(w[len(w)-m])
        while w and w[-1]==k_alphabet-1: w.pop()
    return out
def mobius_count(n,k=2):
    return sum(sympy.mobius(d)*k**(n//d) for d in sympy.divisors(n))//n
def cycle(steps):
    # steps: list of (sigma, s); y_{i+1} = (1 + sigma_i y_i)/(4 s_i); unique fixed point (strict contraction)
    A=Fraction(1); B=Fraction(0)
    for sg,s in steps:
        A,B=A*Fraction(sg,4*s),(B*sg+1)/(4*s)
    y0=B/(1-A); ys=[y0]
    for sg,s in steps[:-1]: ys.append((1+sg*ys[-1])/(4*s))
    assert (1+steps[-1][0]*ys[-1])/(4*steps[-1][1])==y0
    P0=math.lcm(*[y.denominator for y in ys])
    return P0,[int(P0*y) for y in ys]

# ---- A. the cofactor-1 family, lengths 6..18 ----
words={}; per_len={}
for n in range(6,19):
    ws=[tuple((1,3)[b] for b in w) for w in lyndon(2,n)]
    per_len[n]=len(ws); words[n]=ws
chk('A: Lyndon words per length (Duval) equal the Moebius formula, total',(all(per_len[n]==mobius_count(n) for n in per_len),sum(per_len.values())),(True,31028))
hard=[]; gstats=Counter(); n_int=0
for n,ws in words.items():
    for w in ws:
        P0,qs=cycle([(sg,1) for sg in w])
        assert all(P0+sg*q==4*qs[(i+1)%n] for i,(sg,q) in enumerate(zip(w,qs)))
        N=4**n-math.prod(w); g=N//P0; assert g*P0==N and math.gcd(*qs)==1
        n_int+=1
        gstats['g=1']+=g==1
        if P0%840 in HARD:
            hard.append((n,w,P0,g,qs))
            gstats['hard: g mod 8 in {5,7}']+=g%8 in (5,7)
            gstats['hard: g = -3^k mod 8']+=(g-(-pow(3,w.count(3),8)))%8==0
chk('A: integral reconstructions',n_int,31028)
chk('A: words in the six hard classes',len(hard),337)
dist=Counter(h[2] for h in hard)
chk('A: distinct p and number of words for each',sorted(dist.items()),[(9361,1),(80641,3),(598249,5),(2396641,15),(17602321,2),(330382081,22),(464316001,34),(2454266089,255)])
chk('A: every p composite (sympy), smallest prime factor of each',[(p,sympy.isprime(p),min(sympy.factorint(p))) for p in sorted(dist)],
    [(9361,False,11),(80641,False,11),(598249,False,281),(2396641,False,13),(17602321,False,11),(330382081,False,83),(464316001,False,21191),(2454266089,False,11)])
first=min(hard,key=lambda h:(h[0],h[2]))
chk('A: first hard word: length, word, p, g, factorisation',(first[0],first[1],first[2],first[3],sympy.factorint(first[2])),(8,(1,1,1,1,1,1,3,3),9361,7,{11:1,23:1,37:1}))
print('     A: words with g = 1:',gstats['g=1'],'of 31028 ; hard words with g mod 8 in {5,7}:',gstats['hard: g mod 8 in {5,7}'],'; hard words with g = -3^k mod 8:',gstats['hard: g = -3^k mod 8'],'of',len(hard))

# ---- B. comparison with the PR results JSON ----
if len(sys.argv)>1:
    R=json.load(open(sys.argv[1],encoding='utf-8'))
    theirs=sorted((c['length'],tuple(c['sigmas']),c['p'],c['g'],c['factor'],c['cofactor']) for c in R['hard_composite_candidates'])
    mine=sorted((n,w,p,g,min(sympy.factorint(p)),p//min(sympy.factorint(p))) for n,w,p,g,_ in hard)
    chk('B: the 337 hard words (length, word, p, g, least factor, cofactor) equal the PR JSON',theirs==mine,True)
    chk('B: per-length counts in the PR JSON',[(c['length'],c['primitive_necklaces'],c['integral'],c['hard_class'],c['prime_p']) for c in R['counts']],
        [(n,per_len[n],per_len[n],sum(1 for h in hard if h[0]==n),0) for n in range(6,19)])
    chk('B: witness / complete_through',(R['witness'],R['complete_through']),(None,18))

# ---- C. the local-failure lemma, checked numerically (NOT in >>142) ----
# q single-source: 3 sigma(q) q > p. x = (p + sigma q)/4. Cofactor family: x = s*r with r prime.
# >>142 (paper): s = 1 forces both branches to fail. Here: also s = 2 (paper argument in the reply draft).
P17=[30241,61681,106681,111721,148201,153409,163249,206641,216841,280009,288649,307129,309241,325081,347041,352081,374929]
hardp=[p for p in sympy.primerange(100,100000) if p%840 in HARD]
def sig(r): return 1 if r%4==3 else 3
def leg(a,p): return sympy.legendre_symbol(a%p,p)
def dsq(x):
    out=[1]
    for r,e in sympy.factorint(x).items(): out=[a*r**k for a in out for k in range(2*e+1)]
    return out
tab=Counter(); fails_s12=0; succ_by_s=Counter(); examples={}
inputs=sorted(set(hardp)|{2521,12889}|set(P17))
for p in inputs:
    for q in sympy.primerange(3,p):
        if leg(q,p)!=-1 or 3*sig(q)*q<p: continue
        c=sig(q)*q; x=(p+c)//4
        f=sympy.factorint(x); big=max(f)
        if f[big]!=1 or big*9<=p: tab['other']+=1; continue   # x has no prime factor > p/9 (or it is repeated)
        s=x//big
        ok=any((d+p*x)%c==0 or (d+x)%c==0 for d in dsq(x))
        tab[s]+=1
        if ok:
            succ_by_s[s]+=1
            examples.setdefault(s,(p,q,c,x,dict(f)))
print('     C: inputs',len(inputs),'(hard primes 100 < p < 10^5, plus 2521, 12889 and the 17 primes); single-source vertices by cofactor s = x / (largest prime factor of x, when that factor is > p/9):')
print('       ',sorted((k,v) for k,v in tab.items() if k!='other'),'; other',tab['other'])
print('        successes by s:',sorted(succ_by_s.items()))
for s in sorted(examples): print('        first success with s =',s,':',examples[s])
chk('C: single-source vertices with x prime (s = 1) or x = 2 * prime (s = 2): successes',(succ_by_s[1],succ_by_s[2]),(0,0))
chk('C: s = 1 and s = 2 both occur',(tab[1]>0,tab[2]>0),(True,True))

# ---- D. mod-3 condition for a prime cycle (NOT in >>142): sigma_i = 1 forces q_i = 1 mod 3 and q_(i+1) = 2 mod 3 ----
no11=lambda w: all(not(w[i]==1 and w[(i+1)%len(w)]==1) for i in range(len(w)))
allw=[w for n in words for w in words[n]]
print('     D: words with no cyclically adjacent "1 1":',sum(no11(w) for w in allw),'of',len(allw),'; among the 337 hard words:',sum(no11(h[1]) for h in hard))
print('     D: hard words whose reconstructed q_i include a multiple of 3:',sum(any(q%3==0 for q in h[4]) for h in hard),'; with every q_i prime:',sum(all(sympy.isprime(q) for q in h[4]) for h in hard))
# run 1 checked the window q_i, q_(i+1); the derivation gives q_i, q_(i+1), q_(i+2):
# q_i = 0 (mod 3), or q_i = 2 -> q_(i+1) = 0, or q_i = 1 -> q_(i+1) = 2 -> q_(i+2) = 0 (since p = 1 mod 3, 4 = 1 mod 3)
viol=0
for n,w,p,g,qs in hard:
    for i in range(n):
        if w[i]==1 and w[(i+1)%n]==1 and all(qs[(i+j)%n]%3!=0 for j in range(3)): viol+=1
chk('D: in the 337 hard words, every adjacent "1 1" at i, i+1 puts a multiple of 3 among q_i, q_(i+1), q_(i+2)',viol,0)
print(nok,'OK /',ndiff,'DIFF,',round(time.time()-t0,1),'s')
