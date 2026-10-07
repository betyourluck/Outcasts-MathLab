from math import isqrt, gcd
def jacobi(a,n):
    assert n>0 and n%2==1
    a%=n; r=1
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5): r=-r
        a,n=n,a
        if a%4==3 and n%4==3: r=-r
        a%=n
    return r if n==1 else 0
def isprime(n):
    if n<2: return False
    for q in range(2,isqrt(n)+1):
        if n%q==0: return False
    return True
def divisors(n):
    out=[]
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            out.append(d)
            if d*d!=n: out.append(n//d)
    return out
# (1) every witness s with s|ap+u, 4au|s+1: Jacobi(p mod s / s) ?
from collections import Counter
cnt=Counter(); n=0
for p in range(3,20000):
    if not isprime(p): continue
    for a in range(1,9):
        for u in range(1,9):
            for s in divisors(a*p+u):
                if (s+1)%(4*a*u)==0:
                    n+=1; cnt[jacobi(p,s)]+=1
print("witnesses",n,"jacobi(p|s) distribution",dict(cnt))
# (2) residues of known hard primes mod 840
sq840=sorted({(x*x)%840 for x in range(840) if gcd(x,840)==1})
print("squares mod 840 coprime:",sq840)
for p in (409,577,2521,66529,1740481,5843041):
    print(p, p%840, (p%840) in sq840)
