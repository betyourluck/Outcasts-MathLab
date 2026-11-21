# Written 2026-11-21 when the checks were placed on the fork. check146-example.log (2026-11-17) is the output of an inline
# command that was not saved; this script recomputes the same quantities for the boundary example quoted in >>147:
# p = 40609, q = 14767 (single source, tau = 1), x = 13844 = 2^2 * 3461, Type I witness d = 3461.
# Every equality is checked in integers or with Fraction (no rounded value is compared).
import sys
from fractions import Fraction
import sympy
sys.stdout.reconfigure(encoding='utf-8')
p, q = 40609, 14767
def nr(a): return pow(a, (p-1)//2, p) == p-1
l = 2
while not nr(l): l = sympy.nextprime(l)
sigma = 1 if q % 4 == 3 else 3
print('p prime', sympy.isprime(p), '; p mod 840 =', p % 840, '; least non-residue prime l =', l)
print('q prime', sympy.isprime(q), '; q non-residue', nr(q), '; sigma =', sigma)
c = sigma*q
x, rem = divmod(p + c, 4)
assert rem == 0
f = sympy.factorint(x)
print('c =', c, '; 3c > p:', 3*c > p, '; only tau = 1 admissible (9c >= 3p):', 9*c >= 3*p, '; x =', x, '=', dict(f))
r = max(f); s = x//r
print('r =', r, '; r non-residue', nr(r), '; s = x/r =', s, '; 9r > p:', 9*r > p, '; l*r > p:', l*r > p, '; 12r > p:', 12*r > p)
ds = [1]
for pr, e in f.items(): ds = [a*pr**k for a in ds for k in range(2*e+1)]
t1 = sorted(d for d in ds if (d + p*x) % c == 0)
t2 = sorted(d for d in ds if (d + x) % c == 0)
print('divisors of x^2:', len(ds), '; Type I witnesses d:', t1, '; Type II witnesses d:', t2)
d = t1[0]; N = p*x
y = (N + d)//c; z = (N + N*N//d)//c
assert (N + d) % c == 0 and (N*N) % d == 0 and (N + N*N//d) % c == 0
print('denominators', x, y, z, '; 4/p = 1/x + 1/y + 1/z exactly:', Fraction(1, x) + Fraction(1, y) + Fraction(1, z) == Fraction(4, p))
e = x*x//d
print('e = x^2/d =', e, '=', dict(sympy.factorint(e)), '; 4e + 1 =', 4*e + 1, '=', (4*e + 1)//c, '* c, remainder', (4*e + 1) % c)
h, alpha, beta = 11, 4, 1
print('h c = alpha p + beta with (h, alpha, beta) = (11, 4, 1):', h*c == alpha*p + beta, '; h non-residue', nr(h), '; h = l:', h == l)
bound = Fraction(4*h*s, h + alpha)
print('4hs/(h + alpha) =', bound, '; p/r =', Fraction(p, r), '; difference =', bound - Fraction(p, r),
      '; 1/(15*3461) =', Fraction(1, 15*3461), '; 176*r - 15*p =', 176*r - 15*p)
