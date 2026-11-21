# Contract-free (2026-11-17): the population behind the counts quoted in the draft reply to >>146, after the master-side
# audit could not reproduce "4,850 sources / 4,757 Type I witnesses" from the wording. Same 288 inputs.
# Definitions used here (and in check146-quick.py part D):
#   source   = admissible (q, tau): q non-residue prime < p, tau odd, c = sigma(q) q tau^2 < 3p, x = (p + c)/4
#   Omega_-  = number of non-residue prime factors of x counted with multiplicity
#   Type I witness = a divisor d of x^2 with c | d + p x     (one witness per divisor d)
#   Type II witness = a divisor d of x^2 with c | d + x
import sys, time
from fractions import Fraction
from collections import Counter
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
HARD={1,121,169,289,361,529}
P17=[30241,61681,106681,111721,148201,153409,163249,206641,216841,280009,288649,307129,309241,325081,347041,352081,374929]
inputs=sorted(set(p for p in sympy.primerange(1009,100000) if p%840 in HARD)|{2521,12889}|set(P17))
def sig(r): return 1 if r%4==3 else 3
C=Counter(); expo=Counter()
for p in inputs:
    for q in sympy.primerange(3,p):
        if pow(q,(p-1)//2,p)!=p-1: continue
        A=sig(q)*q; t=1
        while A*t*t<3*p:
            c=A*t*t; x=(p+c)//4; f=sympy.factorint(x)
            nr=[r for r in f if pow(r,(p-1)//2,p)==p-1]; om=sum(f[r] for r in nr)
            ds=[1]
            for r,e in f.items(): ds=[a*r**k for a in ds for k in range(2*e+1)]
            t1=[d for d in ds if (d+p*x)%c==0]; t2=[d for d in ds if (d+x)%c==0]
            key='Omega=1' if om==1 else 'Omega>=3'
            C[key,'sources']+=1
            C[key,'sources with a Type I witness']+=bool(t1)
            C[key,'sources with a Type II witness']+=bool(t2)
            C[key,'sources with either']+=bool(t1 or t2)
            C[key,'Type I witnesses']+=len(t1)
            C[key,'Type II witnesses']+=len(t2)
            if om==1:
                r=nr[0]
                for d in t1: expo['I',sympy.multiplicity(r,d)]+=1
                for d in t2: expo['II',sympy.multiplicity(r,d)]+=1
            t+=2
for k in sorted(C): print('%-10s %-34s %9d'%(k[0],k[1],C[k]))
print('total sources',C['Omega=1','sources']+C['Omega>=3','sources'],'; total Type I witnesses',C['Omega=1','Type I witnesses']+C['Omega>=3','Type I witnesses'],
      '; total Type II witnesses',C['Omega=1','Type II witnesses']+C['Omega>=3','Type II witnesses'])
print('exponent of r in the witnesses, Omega_- = 1 sources:',sorted(expo.items()))
# the boundary example: p/r versus 176/15
p,r=40609,3461
print('p/r =',Fraction(p,r),'; 176/15 - p/r =',Fraction(176,15)-Fraction(p,r),'; 1/(15*3461) =',Fraction(1,15*3461),'; 176*r - 15*p =',176*r-15*p)
print('%.1f s'%(time.time()-t0))
