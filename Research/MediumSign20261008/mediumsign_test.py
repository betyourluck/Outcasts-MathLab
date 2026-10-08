"""Y test (contract: Research/MediumSign20261008/contract.md, pushed as 0e7ac68 before computing).
Is the sign bias of prime factors 29..997 of ap+u explained by 'q divides with probability 1/(q-1),
with the fixed value (d|q)'? Exact integers; no randomness."""
import gzip, json, sys, time
from collections import defaultdict
from math import gcd, isqrt, sqrt
import numpy as np

SQ840 = {1, 121, 169, 289, 361, 529}
SMALL = [q for q in range(2, 20000) if all(q % r for r in range(2, isqrt(q) + 1))]  # sqrt(23 * 6e6) = 11747
MED = [q for q in SMALL if 29 <= q <= 997]

def sieve(n):
    b = bytearray([1]) * (n + 1); b[0] = b[1] = 0
    for i in range(2, isqrt(n) + 1):
        if b[i]: b[i * i::i] = bytearray(len(b[i * i::i]))
    return b

def factor(n):
    f = {}
    for q in SMALL:
        if q * q > n: break
        while n % q == 0:
            f[q] = f.get(q, 0) + 1; n //= q
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def divisors_from(f):
    D = [1]
    for q, e in f.items():
        D = [d * q**k for d in D for k in range(e + 1)]
    return sorted(D)

def divisors_scan(n):
    r = np.arange(1, isqrt(n) + 1, dtype=np.int64)
    small = r[n % r == 0]
    return sorted(set(small.tolist()) | {n // int(d) for d in small})

def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0

def sqfree_divisors(n):
    ps = [q for q in range(2, n + 1) if n % q == 0 and all(q % r for r in range(2, isqrt(q) + 1))]
    ds = [1]
    for q in ps: ds += [d * q for d in ds]
    return ds

def summarize(obs, pred):
    s, n = obs; ps, w = pred
    S = s / n
    return {"n": n, "S_obs": round(S, 5), "SE": round(sqrt((1 - S * S) / n), 5), "S_pred": round(ps / w, 5),
            "diff_over_SE": round((S - ps / w) / sqrt((1 - S * S) / n), 2)}

def main():
    t0 = time.time()
    b = sieve(6 * 10**6)
    primes = [p for p in range(5 * 10**6 + 1, 6 * 10**6) if b[p] and p % 840 in SQ840]
    if len(primes) != 1985:
        raise ValueError(f"expected 1985 primes, got {len(primes)}")
    # precompute 1/(q-1)-weighted character sums per (m, d) lazily
    cache = {}
    def pred_sums(m, d, qs):
        key = (m, d, qs)
        if key not in cache:
            ps = sum(jacobi(d % q, q) / (q - 1) for q in MED if m % q != 0 and lo(qs) <= q <= hi(qs))
            w = sum(1 / (q - 1) for q in MED if m % q != 0 and lo(qs) <= q <= hi(qs))
            cache[key] = (ps, w)
        return cache[key]
    lo = lambda qs: {"all": 29, "29-100": 29, "101-997": 101}[qs]
    hi = lambda qs: {"all": 997, "29-100": 100, "101-997": 997}[qs]
    obs = defaultdict(lambda: [0, 0]); pred = defaultdict(lambda: [0.0, 0.0])
    Y3_checked = 0; Y3_viol = 0; rows = []
    for p in primes:
        B = (p - 1).bit_length()
        for a in range(1, B + 1):
            for u in range(1, B // a + 1):
                if gcd(a, u) != 1: continue
                M = a * p + u; m = 4 * a * u
                f = factor(M)
                D = divisors_from(f)
                if D != divisors_scan(M):
                    raise SystemExit(f"STOP: divisor methods disagree at {(p, a, u)}")
                witness = any((s + 1) % m == 0 for s in D)
                rows.append({"p": p, "a": a, "u": u, "B": B, "M": M, "factor": [[q, e] for q, e in sorted(f.items())],
                             "witness": witness})
                for q in f:
                    if gcd(q, m) == 1 and q != p:
                        Y3_checked += 1
                        if jacobi((-a * u) % q, q) != jacobi(p % q, q): Y3_viol += 1
                meds = [q for q in f if 29 <= q <= 997 and m % q != 0]
                for dp in sqfree_divisors(2 * a * u):
                    d = -dp
                    for qs in ("all", "29-100", "101-997"):
                        ps, w = pred_sums(m, d, qs)
                        for key in ((d, qs), ("ALL", qs)):
                            pred[key][0] += ps; pred[key][1] += w
                        for q in meds:
                            if lo(qs) <= q <= hi(qs):
                                sgn = jacobi(d % q, q)
                                for key in ((d, qs), ("ALL", qs)):
                                    obs[key][0] += sgn; obs[key][1] += 1
    if Y3_viol:
        raise SystemExit(f"STOP: Y3 violated {Y3_viol} times")
    out = {"contract": "Research/MediumSign20261008/contract.md @ 0e7ac68", "primes": len(primes), "rows": len(rows),
           "Y3_checked": Y3_checked, "Y3_violations": Y3_viol}
    pooled = summarize(obs[("ALL", "all")], pred[("ALL", "all")])
    out["Y1_pooled"] = pooled
    out["Y1_pass"] = abs(pooled["diff_over_SE"]) <= 3
    ds = sorted([k[0] for k in obs if k[1] == "all" and k[0] != "ALL"], key=lambda d: -obs[(d, "all")][1])
    top12 = {str(d): summarize(obs[(d, "all")], pred[(d, "all")]) for d in ds[:12]}
    out["Y2_top12"] = top12
    out["Y2_outside_3SE"] = sum(1 for v in top12.values() if abs(v["diff_over_SE"]) > 3)
    out["Y2_pass"] = out["Y2_outside_3SE"] <= 1
    out["Y4_top20"] = {str(d): summarize(obs[(d, "all")], pred[(d, "all")]) for d in ds[:20]}
    out["Y4_pooled_29-100"] = summarize(obs[("ALL", "29-100")], pred[("ALL", "29-100")])
    out["Y4_pooled_101-997"] = summarize(obs[("ALL", "101-997")], pred[("ALL", "101-997")])
    out["seconds"] = round(time.time() - t0, 1)
    with gzip.GzipFile("mediumsign_rows.json.gz", "wb", mtime=0) as fh:
        fh.write(json.dumps(rows).encode())
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
