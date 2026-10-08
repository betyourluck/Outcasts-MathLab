"""Exploratory checks on the ceiling-continued-fraction (FCT) condition of Ventas (arXiv:2605.04551, Thm 2.3)
as quoted in Outcasts >>41:  d | p+i,  d = 3 (mod 4),  4i | p+d.
(1) It is the Type I branch of PR #3 with c = d, x = (p+d)/4, t = i*x (ES identity checked exactly).
(2) If 4i | p-1, its divisor set equals the Type II witness set of Dahan's shift (a,u) = (1,i).
(3) At the fixed x of a doubly failing candidate (p=2521, c=11, x=633) it cannot help.
(4) First FCT depth i for the five primes of >>36/>>37 and over the hard primes below 10^6.
Exact integers; no randomness."""
import json, sys
from collections import Counter
from fractions import Fraction
from math import isqrt

SQ840 = {1, 121, 169, 289, 361, 529}

def isprime(n): return n > 1 and all(n % k for k in range(2, isqrt(n) + 1))

def divisors(n):
    s = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0: s |= {d, n // d}
    return sorted(s)

def fct_set(p, i): return {d for d in divisors(p + i) if d % 4 == 3 and (p + d) % (4 * i) == 0}
def dahan_set(p, i): return {d for d in divisors(p + i) if (d + 1) % (4 * i) == 0}

def main():
    out = {}
    # (1)
    checked = ok = 0
    for p in (q for q in range(5, 20000) if isprime(q) and q % 4 == 1):
        for i in range(1, 40):
            for d in fct_set(p, i):
                x = (p + d) // 4; t = i * x; q2 = x * x // t
                assert x % i == 0 and (x * x) % t == 0 and (t + p * x) % d == 0
                y = (t + p * x) // d; z = (p * p * q2 + p * x) // d
                checked += 1; ok += Fraction(1, x) + Fraction(1, y) + Fraction(1, z) == Fraction(4, p)
    out["1_type_I_branch"] = {"primes": "p = 1 (mod 4), p < 20000", "i": "1..39", "triples": checked, "identity_holds": ok}
    # (2)
    hard = [p for p in range(1009, 300000) if p % 840 in SQ840 and isprime(p)]
    eq = tot = 0; ne = tot2 = 0
    for p in hard:
        for i in range(1, 24):
            f, g = fct_set(p, i), dahan_set(p, i)
            if (p - 1) % (4 * i) == 0: tot += 1; eq += f == g
            else: tot2 += 1; ne += f != g
    out["2_equals_dahan_1_i"] = {"hard_primes_below": 300000, "count": len(hard), "i": "1..23",
                                 "4i|p-1: equal / cases": [eq, tot], "4i does not divide p-1: differ / cases": [ne, tot2]}
    # (3)
    p, c = 2521, 11; x = (p + c) // 4
    out["3_fixed_x_2521_c11"] = {"x": x, "i in divisors(x) with c | p+i": [i for i in divisors(x) if (p + i) % c == 0]}
    # (4)
    first = {}
    for p in (2521, 66529, 345601, 670849, 5843041):
        for i in range(1, 2001):
            s = sorted(fct_set(p, i))
            if s: first[p] = {"i": i, "d": s[:4], "x": [(p + d) // 4 for d in s[:2]]}; break
    out["4_first_depth_five_primes"] = first
    hard6 = [p for p in range(1009, 10**6) if p % 840 in SQ840 and isprime(p)]
    dep = Counter(); worst = (0, None)
    for p in hard6:
        for i in range(1, 5001):
            if fct_set(p, i): break
        dep[i if i <= 10 else ">10"] += 1
        if i > worst[0]: worst = (i, p)
    out["4_first_depth_hard_below_1e6"] = {"primes": len(hard6), "distribution": {str(k): v for k, v in dep.items()},
                                           "max_first_depth": worst[0], "at_p": worst[1]}
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
