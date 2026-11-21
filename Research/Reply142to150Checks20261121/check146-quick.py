# Contract-free quick check (2026-11-17) of >>146 / PR #4 20b7ac7 (UpperSourceCharacter20261117), written after reading >>146.
# Parts A-C check what >>146 states. Parts D-F are NOT in >>146: the threshold p/9 lowered to p/l, where l is the least
# non-residue prime of p (paper argument in the reply draft), checked on the 288 inputs of check142/check144.
# Units: "source" = one admissible (q, tau) with q a non-residue prime < p, tau odd, c = sigma(q) q tau^2 < 3p;
# "upper" = 3c > p; "succeeds" = some divisor d of x^2 has c | d + p x (Type I) or c | d + x (Type II).
import sys, json, time
from collections import Counter
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
HARD={1,121,169,289,361,529}
P17=[30241,61681,106681,111721,148201,153409,163249,206641,216841,280009,288649,307129,309241,325081,347041,352081,374929]
nok=0; ndiff=0
def chk(name,got,exp):
    global nok,ndiff
    if got==exp: nok+=1; print('OK  ',name,got if len(str(got))<170 else '(long)')
    else: ndiff+=1; print('DIFF',name,'got',got,'expected',exp)
def sig(r): return 1 if r%4==3 else 3
def dsq(x,f=None):
    out=[1]
    for r,e in (f or sympy.factorint(x)).items(): out=[a*r**k for a in out for k in range(2*e+1)]
    return out
def least_nr(p):
    l=2
    while pow(l,(p-1)//2,p)!=p-1: l=sympy.nextprime(l)
    return l

# ---- A. the coefficient bound of >>146 (h <= 9) by brute force over the integers it quantifies ----
tup=[(s,a,h) for s in range(1,9) for a in range(1,s+1) for h in range(1,3*a+1) if 4*h*s<=9*(h+a)]
chk('A: integer triples (s, a, h) with 1<=s<=8, a<=s, h<=3a, 4hs<=9(h+a): max h',max(h for s,a,h in tup),9)
print('     A: triples attaining h = 9:',[t for t in tup if t[2]==9],'; number of triples',len(tup))
# ---- B. the saved 438 rows of 2026-11-16 (PR JSON given as argv[1]) ----
if len(sys.argv)>1:
    R=json.load(open(sys.argv[1],encoding='utf-8'))['rows']
    cong=[x for x in R if x['outcome']=='congruence']
    big=[x for x in R if any(l>=11 for l in sympy.primefactors(x['h']))]
    chk('B: remaining tuples, their max h, tuples with a prime factor >= 11 in h, and the outcomes of those',
        (len(cong),max(x['h'] for x in cong),len(big),dict(Counter(x['outcome'] for x in big))),(67,9,34,{'r_range':34}))
# ---- C. small integers are residues for hard primes ----
hard6=[p for p in sympy.primerange(1000,10**6) if p%840 in HARD]
chk('C: hard primes below 10^6 for which some n in 1..9 is a non-residue',sum(any(pow(n,(p-1)//2,p)!=1 for n in range(1,10)) for p in hard6),0)
print('     C: hard primes below 10^6:',len(hard6),'; smallest:',hard6[0])

# ---- D. threshold p/l on the 288 inputs ----
inputs=sorted(set(p for p in sympy.primerange(1009,100000) if p%840 in HARD)|{2521,12889}|set(P17))
cnt=Counter(); by_s=Counter(); succ_rows=[]; ratio_max=(0,None); expo=Counter(); lvals=Counter()
for p in inputs:
    l=least_nr(p); lvals[l]+=1
    assert l*l<p and all(pow(n,(p-1)//2,p)==1 for n in range(1,l))
    for q in sympy.primerange(3,p):
        if pow(q,(p-1)//2,p)!=p-1: continue
        A=sig(q)*q; t=1
        while A*t*t<3*p:
            c=A*t*t; x=(p+c)//4; f=sympy.factorint(x); big=max(f)
            upper=3*c>p
            nr=[r for r in f if pow(r,(p-1)//2,p)==p-1]
            om=sum(f[r] for r in nr)
            ds=dsq(x,f)
            t1=[d for d in ds if (d+p*x)%c==0]; t2=[d for d in ds if (d+x)%c==0]
            ok=bool(t1 or t2)
            cnt['sources']+=1; cnt['successful sources']+=ok
            if upper:
                cnt['upper']+=1
                if big*l>p:
                    cnt['upper, prime factor > p/l']+=1; by_s[x//big]+=1
                    if big*9<=p: cnt['upper, largest prime factor in (p/l, p/9]']+=1
                    assert f[big]==1 and nr==[big]
                    if ok: succ_rows.append((p,l,q,t,c,x,dict(f),t1,t2))
                elif ok:
                    cnt['upper successes with every prime factor <= p/l']+=1
                    rr=big*l/p
                    if rr>ratio_max[0]: ratio_max=(rr,(p,l,q,t,c,x,dict(f)))
                if big*9>p: cnt['upper, prime factor > p/9']+=1
            # exponent of the unique non-residue prime in the witnesses (>>76 section 2, >>88 section 3)
            if om==1 and ok:
                r=nr[0]; cnt['successful sources with Omega_- = 1']+=1
                for d in t1: expo[('I',sympy.multiplicity(r,d))]+=1
                for d in t2: expo[('II',sympy.multiplicity(r,d))]+=1
            elif ok: cnt['successful sources with Omega_- >= 3']+=1
            t+=2
print('     D: inputs',len(inputs),'; least non-residue prime l:',sorted(lvals.items()))
for k in sorted(cnt): print('        %-52s %9d'%(k,cnt[k]))
chk('D: upper sources with a prime factor > p/l that succeed',len(succ_rows),0)
for row in succ_rows[:10]: print('        SUCCESS',row)
print('     D: cofactor s = x / r for those sources:',sorted(by_s.items()))
print('     D: exponent of the non-residue prime in the witnesses d (Omega_- = 1 sources):',sorted(expo.items()))
chk('D: Type I witnesses have exponent 1, Type II witnesses exponent 0 or 2',all((b=='I' and e==1) or (b=='II' and e in (0,2)) for b,e in expo),True)
print('     D: largest (largest prime factor of x) * l / p among upper successes: %.4f'%ratio_max[0],ratio_max[1])

# ---- E. l^2 < p for every hard prime below 10^8 ----
worst=(10**9,None); lmax=(0,None); n=0
for a in sorted(HARD):
    for p in range(a+840,10**8,840):
        if p%3==0 or p%5==0 or p%7==0 or p%11==0 and p>11 and False: pass
        if not sympy.isprime(p): continue
        n+=1; l=least_nr(p)
        if p/(l*l)<worst[0]: worst=(p/(l*l),(p,l))
        if l>lmax[0]: lmax=(l,p)
chk('E: hard primes below 10^8',n,179468)
print('     E: smallest p / l^2 = %.2f at (p, l) ='%worst[0],worst[1],'; largest l =',lmax[0],'at p =',lmax[1])
chk('E: l^2 < p for every hard prime below 10^8',worst[0]>1,True)
print(nok,'OK /',ndiff,'DIFF,',round(time.time()-t0,1),'s')
