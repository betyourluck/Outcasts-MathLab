"""Independent re-implementation (Luna) of the Kneser sufficient conditions in Outcasts >>46,
written from the post and the two contracts only, before reading the author's code.

m = 4au, M = ap + u, keep the prime factors q of M coprime to m with their exponents e_q.
G = <q mod m> in (Z/m)^*, D = prod {1, q, ..., q^e} mod m.  Target: -1 in D.
K (first contract): H0 = closure of saturating primes; -1 in H0, or (-1 in G and for every
    H0 <= H <= G with -1 notin H: sum_q min(e_q, ord_{G/H} q - 1) >= |G/H| - 1).
A (stronger theorem, admissible H only): -1 in G and for every admissible H (each q notin H has
    ord_{G/H} q > e_q + 1) with -1 notin H: sum_{q notin H} e_q >= |G/H| - 1.
"""
import json, sys, time
from math import gcd, isqrt
from fractions import Fraction
from functools import lru_cache

SQ840 = {1, 121, 169, 289, 361, 529}

def sieve(n):
    b = bytearray([1]) * (n + 1); b[0] = b[1] = 0
    for i in range(2, isqrt(n) + 1):
        if b[i]: b[i * i::i] = bytearray(len(b[i * i::i]))
    return [i for i in range(n + 1) if b[i]]

SMALL = sieve(60000)  # covers M = ap + u < 3.6e9 (p < 10^8, au <= 27)

def factor(n):
    f = {}
    for q in SMALL:
        if q * q > n: break
        while n % q == 0:
            f[q] = f.get(q, 0) + 1; n //= q
    if n > 1:
        assert n < SMALL[-1] ** 2
        f[n] = f.get(n, 0) + 1
    return f

def clog2(p):
    return (p - 1).bit_length()  # smallest k with 2^k >= p (p >= 2)

def closure(m, gens, start=frozenset([1])):
    S = set(start); frontier = list(S)
    gens = [g % m for g in gens]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = x * g % m
                if y not in S: S.add(y); new.append(y)
        frontier = new
    return frozenset(S)

@lru_cache(maxsize=None)
def all_subgroups(m):
    U = [x for x in range(1, m) if gcd(x, m) == 1] if m > 1 else [0]
    cyc = {closure(m, [x]) for x in U}
    subs = set(cyc) | {frozenset([1])}
    changed = True
    while changed:
        changed = False
        cur = list(subs)
        for A in cur:
            for C in cyc:
                J = closure(m, list(C), A)
                if J not in subs: subs.add(J); changed = True
    return tuple(sorted(subs, key=len))

def qord(m, r, H):
    """order of r in G/H (r in G, H subgroup)."""
    k, x = 1, r % m
    while x not in H:
        x = x * r % m; k += 1
    return k

def conditions(p, a, u):
    m = 4 * a * u; M = a * p + u
    fM = factor(M)
    prod = 1
    for q, e in fM.items(): prod *= q ** e
    assert prod == M
    F = {q: e for q, e in fM.items() if gcd(q, m) == 1}
    G = closure(m, list(F))
    minus = m - 1
    subs = [H for H in all_subgroups(m) if H <= G]
    hG = len(G)
    # K
    H0 = frozenset([1]); changed = True
    while changed:
        changed = False
        for q, e in F.items():
            if q % m not in H0 and e >= qord(m, q, H0) - 1:
                H0 = closure(m, [q], H0); changed = True
    if minus in H0:
        K = True; K_fail = []
    elif minus not in G:
        K = False; K_fail = ["-1 notin G"]
    else:
        K_fail = []
        for H in subs:
            if not (H0 <= H) or minus in H: continue
            L = sum(min(e, qord(m, q, H) - 1) for q, e in F.items())
            if L < hG // len(H) - 1: K_fail.append({"H_size": len(H), "L": L, "need": hG // len(H) - 1})
        K = not K_fail
    # A
    if minus not in G:
        A = False; A_fail = ["-1 notin G"]
    else:
        A_fail = []
        for H in subs:
            if minus in H: continue
            if any(q % m not in H and qord(m, q, H) <= e + 1 for q, e in F.items()): continue  # not admissible
            L = sum(e for q, e in F.items() if q % m not in H)
            if L < hG // len(H) - 1: A_fail.append({"H": sorted(H), "L": L, "need": hG // len(H) - 1})
        A = not A_fail
    # actual residue set D and witnesses (all divisors of M)
    D = {1}
    for q, e in F.items():
        D = {d * pow(q, k, m) % m for d in D for k in range(e + 1)}
    divs = [1]
    for q, e in fM.items():
        divs = [d * q ** k for d in divs for k in range(e + 1)]
    wit = sorted(s for s in divs if (s + 1) % m == 0)
    assert (minus in D) == bool(wit)
    return {"a": a, "u": u, "m": m, "M": M, "factors": {str(q): e for q, e in sorted(fM.items())},
            "G_size": hG, "minus1_in_G": minus in G, "H0_size": len(H0), "K": K, "A": A,
            "K_fail_count": len(K_fail), "A_fail": A_fail[:3], "D": sorted(D), "witnesses": wit}

def coprime_pairs(B):
    return [(a, u) for a in range(1, B + 1) for u in range(1, B // a + 1) if gcd(a, u) == 1]

def is_prime(n):
    return n > 1 and all(n % q for q in SMALL if q * q <= n)

def identity(p, a, u, s):
    M = a * p + u; t = M // s; v = (s + 1) // (4 * a * u)
    x, y, z = u * v * t, a * u * v * p, a * v * t * p
    return {"x": x, "y": y, "z": z, "ok": Fraction(4, p) == Fraction(1, x) + Fraction(1, y) + Fraction(1, z)}

def main():
    t0 = time.time(); out = {}
    p = 8615161
    out["p"] = {"p": p, "prime": is_prime(p), "mod840": p % 840, "clog2": clog2(p), "pow_check": [2 ** 23 < p, p <= 2 ** 24]}
    pairs = coprime_pairs(clog2(p))
    rows = [conditions(p, a, u) for a, u in pairs]
    out["pairs"] = len(pairs)
    out["counts"] = {"minus1_notin_G": sum(not r["minus1_in_G"] for r in rows),
                     "minus1_in_G_but_K_fails": sum(r["minus1_in_G"] and not r["K"] for r in rows),
                     "minus1_in_G_but_A_fails": sum(r["minus1_in_G"] and not r["A"] for r in rows),
                     "K_true": sum(r["K"] for r in rows), "A_true": sum(r["A"] for r in rows),
                     "pairs_with_witness": sum(bool(r["witnesses"]) for r in rows),
                     "total_witnesses": sum(len(r["witnesses"]) for r in rows)}
    out["witness_pairs"] = [{"a": r["a"], "u": r["u"], "witnesses": r["witnesses"]} for r in rows if r["witnesses"]]
    r41 = next(r for r in rows if (r["a"], r["u"]) == (4, 1))
    out["pair_4_1"] = {k: r41[k] for k in ("M", "factors", "G_size", "D", "witnesses", "A_fail", "K_fail_count")}
    out["pair_4_1"]["residues_mod16"] = {q: int(q) % 16 for q in r41["factors"]}
    out["pair_4_1"]["identity_s815"] = identity(p, 4, 1, 815)
    out["rows"] = rows
    # first-counterexample scan over hard primes below p
    t1 = time.time()
    hard = [q for q in sieve(p) if q % 840 in SQ840]
    out["scan"] = {"hard_primes_up_to_p": len(hard), "index_of_p": hard.index(p) + 1 if p in hard else None}
    first_fail = None; false_pos = []
    for q in hard:
        ok = False
        for a, u in coprime_pairs(clog2(q)):
            r = conditions(q, a, u)
            if (r["K"] or r["A"]) and not r["witnesses"]: false_pos.append((q, a, u))
            if r["K"]:
                ok = True; break
        if not ok:
            first_fail = q; break
    out["scan"].update({"first_prime_without_K_pair": first_fail, "false_positives": false_pos, "seconds": round(time.time() - t1, 1)})
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
