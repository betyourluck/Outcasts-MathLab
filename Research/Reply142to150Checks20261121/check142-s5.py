# Contract-free (2026-11-15): the single case that survives check142-lemma.py.
# s = 5, sigma(q) = 3: x = (p + 3q)/4 = 5 r with r prime > p/9 (equivalently q > 11p/27), and Type I via e = 5 r^2,
# which holds iff 3q | 20 r^2 + 1 iff q | p^2 + 20 (since 20 (20 r^2 + 1) = (p + 3q)^2 + 20 = p^2 + 20 mod 3q).
# Search: every hard prime p < LIMIT; every prime divisor q of p^2 + 20 with 11p/27 < q < p, q = 1 mod 4,
# (q/p) = -1, 20 | p + 3q, r = (p + 3q)/20 prime. Each hit is then checked directly over all divisors of x^2.
import sys, time
import sympy
sys.stdout.reconfigure(encoding='utf-8')
LIMIT=int(sys.argv[1]) if len(sys.argv)>1 else 10**7
HARD={1,121,169,289,361,529}
t0=time.time(); n=0; hits=[]
def dsq(x):
    out=[1]
    for r,e in sympy.factorint(x).items(): out=[a*r**k for a in out for k in range(2*e+1)]
    return out
for p in sympy.primerange(1009,LIMIT):
    if p%840 not in HARD: continue
    n+=1
    for q in sympy.factorint(p*p+20):
        if not (27*q>11*p and q<p and q%4==1 and (p+3*q)%20==0): continue
        if sympy.legendre_symbol(q,p)!=-1: continue
        r=(p+3*q)//20
        if not sympy.isprime(r): continue
        c=3*q; x=(p+c)//4
        assert x==5*r and 9*r>p and 3*c>p
        ds=dsq(x)
        t1=sorted(d for d in ds if (d+p*x)%c==0); t2=sorted(d for d in ds if (d+x)%c==0)
        hits.append((p,q,r,x,t1,t2))
        print('HIT p=%d q=%d c=%d x=%d=5*%d  Type I d=%s  Type II d=%s'%(p,q,c,x,r,t1,t2),flush=True)
print('hard primes searched:',n,'(1009 <= p <',LIMIT,') ; hits:',len(hits),'; %.1f s'%(time.time()-t0))
