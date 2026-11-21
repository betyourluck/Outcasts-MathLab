# Contract-free (2026-11-15): a finite local check of the candidate lemma
#   (L)  p a hard prime (p mod 840 in {1,121,169,289,361,529}, so p >= 1009), q a non-residue prime with
#        3 sigma(q) q > p (only tau = 1 is admissible), c = sigma(q) q, x = (p + c)/4.  If x has a prime factor
#        r > p/9, then both branches fail at q.
# Reduction (paper, see the reply draft):
#   * x < p and r > p/9 give x = s r with s = x/r in 1..8, gcd(r, s) = 1, r prime > 7.
#   * Type II fails (size lemma of >>104 with 3c - p = 8 mod 12).
#   * Type I succeeds iff some e | x^2 has c | 4e + 1 (gcd(c, x) = 1). Write e = f r^j, f | s^2, j in {0,1,2}.
#     j = 0: c | 4f + 1 <= 257, impossible since c > p/3 > 336.
#     j = 1: 4fr + 1 = k c and 4 s r = p + c give  f p + s = m c  with m = k s - f, 1 <= m <= 3f (c > p/3, p > 24).
#            So c = (f p + s)/m is a linear function of p; every remaining condition is a congruence on p.
#     j = 2: c | 4 s^2 (4 f r^2 + 1) = f (p + c)^2 + 4 s^2, so c | f p^2 + 4 s^2; hence q | f p^2 + 4 s^2 and
#            (-f/q) = 1 (q > 64 >= f, q does not divide 2 s p); for sigma = 3 also 3 | f + s^2.
#   Local conditions used (all necessary): c = 3 mod 4; v3(c) = 1 if sigma = 3, 0 if sigma = 1; 5, 7 do not divide c;
#   v_l(x) = v_l(s) for l in {2,3,5,7} (r prime > 7, r not dividing s), with x = (p + c)/4;
#   p = 1 mod 8, p = 1 mod 3, p = +-1 mod 5, p in {1,2,4} mod 7, p not divisible by any other small prime;
#   for j = 1 and a prime l | m outside {2,3,5,7}: l does not divide c (q > p/3 > 336 >= 3f >= m >= l).
# A case dies if one prime has no admissible local residue (conditions at different primes are independent by CRT,
# except the Jacobi symbol in j = 2, which is handled by enumerating q mod 20160 jointly).
import sys, time
from math import gcd
import sympy
sys.stdout.reconfigure(encoding='utf-8')
t0=time.time()
def v(n,l):
    if n==0: return 99
    k=0
    while n%l==0: n//=l; k+=1
    return k
def p_ok(l,pr):
    # local class of a hard prime at l (pr is a residue modulo a power of l)
    if l==2: return pr%8==1
    if l==3: return pr%3==1
    if l==5: return pr%5 in (1,4)
    if l==7: return pr%7 in (1,2,4)
    return pr%l!=0
def need_c_prec(l,s):
    # precision (exponent of l) needed on c and on p + c
    if l==2: return max(2,v(s,2)+3)
    if l==3: return max(2,v(s,3)+1)
    if l in (5,7): return v(s,l)+1
    return 1
def local_c_ok(l,cr,pr,s,sigma,b):
    # cr = c mod l^b, pr = p mod l^b
    mod=l**b
    if l==2:
        if cr%4!=3: return False
        return v((pr+cr)%mod,2)==v(s,2)+2 if v(s,2)+2<b else (pr+cr)%(2**(v(s,2)+2))==0
    if l==3:
        if v(cr%mod,3)!=(1 if sigma==3 else 0): return False
        return v((pr+cr)%mod,3)==v(s,3)        # 4 is a unit at 3
    if l in (5,7):
        if cr%l==0: return False
        # exact valuation when the precision allows it; otherwise the necessary condition l^b | p + c
        # (run 2 returned True here, which dropped that condition in the j = 2 branch, where b = 1 for l = 5, 7)
        return v((pr+cr)%mod,l)==v(s,l) if v(s,l)<b else (pr+cr)%mod==0
    return cr%l!=0
survive=[]; dead_size=[]; ncase=0
for s in range(1,9):
    for sigma in (1,3):
        for f in sympy.divisors(s*s):
            # ---- j = 1 ----
            for m in range(1,3*f+1):
                if (m+f)%s: continue            # k = (m + f)/s must be an integer >= 1
                ncase+=1
                primes=set([2,3,5,7])|set(sympy.primefactors(m))
                dead_at=None
                for l in sorted(primes):
                    a=v(m,l); b=need_c_prec(l,s); mod=l**(a+b); found=False
                    for pr in range(mod):
                        if not p_ok(l,pr%l**max(1,a+b) if l!=2 else pr): continue
                        if l==2 and mod<8: raise SystemExit('precision')
                        num=(f*pr+s)%mod
                        if v(num,l)<a: continue
                        mm=m//l**a; cr=(num//l**a)*pow(mm,-1,l**b)%l**b
                        if local_c_ok(l,cr,pr%l**b,s,sigma,b): found=True; break
                    if not found: dead_at=l; break
                if dead_at is None:
                    # size: c = (f p + s)/m must satisfy c > p/3, c < p (sigma = 1) or c < 3p (sigma = 3), and
                    # r = (p + c)/(4 s) > p/9, i.e. 9 (p + c) > 4 s p. All are linear in p; test p >= 1009 exactly.
                    # each condition is  alpha p + beta > 0  with c = (f p + s)/m; intersect the solution sets with p >= 1009
                    from fractions import Fraction as Fr
                    conds=[(Fr(f,m)-Fr(1,3),Fr(s,m)),                       # c > p/3
                           ((1 if sigma==1 else 3)-Fr(f,m),Fr(-s,m)),        # c < p (sigma = 1) or c < 3p (sigma = 3)
                           (9+9*Fr(f,m)-4*s,9*Fr(s,m))]                      # r = (p + c)/(4s) > p/9
                    lo,hi=Fr(1009),None; empty=False
                    for al,be in conds:
                        if al>0: lo=max(lo,-be/al)
                        elif al<0: hi=-be/al if hi is None else min(hi,-be/al)
                        elif be<=0: empty=True
                    if empty or (hi is not None and hi<=lo):
                        dead_size.append(('j=1',s,sigma,f,m)); continue
                    survive.append(('j=1',s,sigma,f,m))
            # ---- j = 2 ----
            ncase+=1
            f0=1
            for l,e in sympy.factorint(f).items():
                if e%2: f0*=l
            if sigma==3 and (f+s*s)%3: continue      # dies at 3: 3 | f p^2 + 4 s^2 = f + s^2 mod 3
            M=2**6*9*5*7; ok=False
            for q in range(1,M,2):
                if q%4!=(3 if sigma==1 else 1) or q%3==0 or q%5==0 or q%7==0: continue
                if sympy.jacobi_symbol(-f0 % (8*f0*q), q)!=1:   # (-f0/q), periodic in q with period dividing 4 f0 | 20160
                    continue
                c=sigma*q
                good=True
                for l,b in ((2,6),(3,2),(5,1),(7,1)):
                    if not any(p_ok(l,pr) and local_c_ok(l,c%l**b,pr,s,sigma,b) for pr in range(l**b)):
                        good=False; break
                if good: ok=True; break
            if ok:
                # size: sigma = 1 means c = q < p, so x < p/2 and s = x/r < (p/2)/(p/9) = 4.5
                if sigma==1 and s>4: dead_size.append(('j=2',s,sigma,f,f0)); continue
                survive.append(('j=2',s,sigma,f,f0))
print('cases checked:',ncase,'; pass the local congruences but die by size:',len(dead_size),dead_size)
print('surviving cases:',len(survive))
for row in survive: print('   SURVIVES',row)
print('%.1f s'%(time.time()-t0))
