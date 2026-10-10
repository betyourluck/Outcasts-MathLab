"""Independent check (Luna) of PR #4 fe2150e Research/Review20261009, items 2 and 3 of its contract.
Written from the README/contract only, before reading the author's code.

2. p = 3678481, B = ceil(log2 p) = 22, all coprime (a,u) with au <= B: witnesses s | ap+u, s = -1 mod 4au.
3. FCT (Ventas v2 Thm 2.3 as stated in the contract): d | p+i, d = 3 mod 4, 4i | p+d.
   Solution derived here: x = (p+d)/4, y = x(p+i)/d, z = p x (p+i)/(d i)
   (4/p - 1/x = d/(px); with N = px, c = ix | N^2 and d | N + c: d/N = 1/y + 1/z, y = (N+c)/d, z = N(N+c)/(dc))."""
import json, sys
from fractions import Fraction
from math import gcd, isqrt

def factor(n):
    f = {}; q = 2
    while q * q <= n:
        while n % q == 0: f[q] = f.get(q, 0) + 1; n //= q
        q += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def divisors(n):
    ds = [1]
    for q, e in factor(n).items(): ds = [d * q ** k for d in ds for k in range(e + 1)]
    return sorted(ds)

def divisors_scan(n):
    s = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0: s.add(d); s.add(n // d)
    return sorted(s)

def item2():
    p = 3678481; B = (p - 1).bit_length()
    rows = []
    for a in range(1, B + 1):
        for u in range(1, B // a + 1):
            if gcd(a, u) != 1: continue
            M = a * p + u; m = 4 * a * u
            d1, d2 = divisors(M), divisors_scan(M)
            assert d1 == d2
            wit = [s for s in d1 if (s + 1) % m == 0]
            rows.append({"a": a, "u": u, "M": M, "omega": len(factor(M)), "witnesses": wit})
    # identity for (1,1), s = 23: t = M/s, v = (s+1)/4
    M = p + 1; s = 23; t = M // s; v = (s + 1) // 4
    x, y, z = v * t, v * p, v * t * p
    return {"p": p, "B": B, "pairs": len(rows),
            "omega_ge4": sum(r["omega"] >= 4 for r in rows),
            "omega_ge4_with_witness": sum(r["omega"] >= 4 and bool(r["witnesses"]) for r in rows),
            "pairs_with_witness": sum(bool(r["witnesses"]) for r in rows),
            "pair_1_1": next(r for r in rows if (r["a"], r["u"]) == (1, 1)),
            "identity_1_1_23": {"x": x, "y": y, "z": z, "ok": Fraction(4, p) == Fraction(1, x) + Fraction(1, y) + Fraction(1, z)},
            "rows": rows}

def fct_first(p, imax=200):
    for i in range(1, imax + 1):
        for d in divisors(p + i):
            if d % 4 == 3 and (p + d) % (4 * i) == 0:
                x = (p + d) // 4; y = x * (p + i) // d; z = p * x * (p + i) // (d * i)
                assert x * (p + i) % d == 0 and p * x * (p + i) % (d * i) == 0
                ok = Fraction(4, p) == Fraction(1, x) + Fraction(1, y) + Fraction(1, z)
                typ = "I" if [v % p == 0 for v in (x, y, z)].count(True) == 1 else "other"
                return {"p": p, "i": i, "d": d, "x": x, "y": y, "z": z, "identity": ok, "positive": min(x, y, z) > 0,
                        "ordered_x_le_y_le_z": x <= y <= z, "type": typ}
    return {"p": p, "i": None}

def item3():
    firsts = [fct_first(p) for p in (2521, 66529, 345601, 670849, 5843041)]
    p, d = 2521, 11; x = (p + d) // 4
    cands = [i for i in divisors(x)]
    return {"first_i": firsts,
            "p2521_d11": {"x": x, "i_candidates_dividing_x": cands,
                          "d_divides_p_plus_i": [i for i in cands if (p + i) % d == 0]}}

def main():
    json.dump({"item2": item2(), "item3": item3()}, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
