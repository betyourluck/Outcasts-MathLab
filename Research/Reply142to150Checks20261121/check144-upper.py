# Contract-free exploration (2026-11-16, after >>144): does the large-prime failure lemma extend from tau = 1 to every
# "upper" source, i.e. any admissible source (q, tau) with c = sigma(q) q tau^2 > p/3 whose x = (p + c)/4 has a prime
# factor r > p/9?  (tau = 1 and c > p/3 is exactly the one-source case of >>143/>>144.)
# Inputs: every hard prime 1009 <= p < 10^5 plus 2521, 12889 and the 17 primes (the 288 inputs of check142).
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
cnt=Counter(); succ=[]; by_tau=Counter(); succ_tau=Counter()
for p in inputs:
    for q in sympy.primerange(3,p):
        if pow(q,(p-1)//2,p)!=p-1: continue
        A=sig(q)*q; t=1
        while A*t*t<3*p:
            c=A*t*t
            if 3*c>p:
                x=(p+c)//4
                for s in range(1,9):
                    if x%s==0 and sympy.isprime(x//s) and 9*(x//s)>p:
                        cnt['upper sources with a prime factor > p/9']+=1; by_tau[min(t,5)]+=1
                        ds=dsq(x)
                        t1=[d for d in ds if (d+p*x)%c==0]; t2=[d for d in ds if (d+x)%c==0]
                        if t1 or t2:
                            succ.append((p,q,t,c,x,s,t1,t2)); succ_tau[min(t,5)]+=1
                        break
                cnt['upper sources']+=1
            t+=2
print('inputs:',len(inputs))
print(dict(cnt))
print('by tau (5 means tau >= 5):',sorted(by_tau.items()),'; successes by tau:',sorted(succ_tau.items()))
print('successes:',len(succ))
for row in succ[:20]: print('   p=%d q=%d tau=%d c=%d x=%d s=%d TypeI=%s TypeII=%s'%row)
print('%.1f s'%(time.time()-t0))
