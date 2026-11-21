# Contract-free addition (2026-11-21), NOT a claim of >>148-150.
# For a few known inputs: among ALL upper denominators x (p/3 < c = 4x - p < 3p) that carry a Type I witness
# (e | x^2, c | 4e + 1), how many have c of the square-source form A(q) t^2 (q a non-residue prime; A = q if q = 3 mod 4, else 3q)?
# This is the base rate behind the "second existence guarantee" of >>150 (the constructed c must be an admissible source).
# Also: all multiplier-3 solutions (a, h) of ex1 (p = 289249), with the form of c.
import sys, time, math
from collections import Counter
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
def nr(q,p): return pow(q,(p-1)//2,p)==p-1
def core_form(c,p):
    """returns q if c = A(q) t^2 with q a non-residue prime of p, else None"""
    f=sympy.factorint(c); odd=sorted(r for r,e in f.items() if e%2)
    if len(odd)==1 and odd[0]!=3 and odd[0]%4==3 and nr(odd[0],p): return odd[0]
    if len(odd)==2 and odd[0]==3 and odd[1]%4==1 and nr(odd[1],p): return odd[1]
    return None
def dsq(f):
    out=[1]
    for r,e in f.items(): out=[a*r**k for a in out for k in range(2*e+1)]
    return out
for p in (2521,12889,30241,40609,289249):
    cnt=Counter(); kmin=Counter()
    for x in range(p//3+1,p):
        c=4*x-p
        if not (3*c>p and c<3*p): continue
        cnt['upper denominators x']+=1
        f=sympy.factorint(x)
        ks=[(4*e+1)//c for e in dsq(f) if (4*e+1)%c==0]
        q=core_form(c,p)
        if q: cnt['... with c = A(q) t^2']+=1
        if ks:
            cnt['x with a Type I witness']+=1
            if q: cnt['x with a Type I witness and c = A(q) t^2']+=1
            if 3 in ks:
                cnt['x with a multiplier-3 witness']+=1
                if q: cnt['x with a multiplier-3 witness and c = A(q) t^2']+=1
    print('p=%d: %s'%(p,dict(cnt)))
p=289249; N=(3*p+1)//4
print('ex1: N = (3p+1)/4 = %d = %s'%(N,sympy.factorint(N)))
D=sympy.divisors(N)
rows=[]
for a in D:
    for h in D:
        if math.gcd(a,h)==1 and N%(a*h)==0 and (a+h)%3==0:
            b=(a+h)//3; x=b*N//h; e=a*N//h; c=4*x-p
            assert (x*x)%e==0 and 4*e+1==3*c and h*c==a*p+b
            side='lower' if 3*c<p else ('upper' if c<3*p else 'c >= 3p')
            rows.append((a,h,b,c,side,core_form(c,p),dict(sympy.factorint(c))))
print('ex1: multiplier-3 pairs (a, h): %d'%len(rows))
for r in rows: print('   a=%d h=%d b=%d c=%d %s  vertex q=%s  c=%s'%r)
print('%.1f s'%(time.time()-t0))
