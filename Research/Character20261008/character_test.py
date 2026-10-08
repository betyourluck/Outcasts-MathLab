"""X test (contract: Research/Character20261008/contract.md, pushed as 58fb041 before computing).
Quadratic-character obstructions of Type II witnesses versus a 'global sign + small-prime signs + fair coins'
model. Exact integers; no randomness."""
import gzip, json, sys, time
from collections import Counter, defaultdict
from math import gcd, isqrt, sqrt
import numpy as np

SQ840 = {1, 121, 169, 289, 361, 529}
SMALL = [q for q in range(2, 20000) if all(q % r for r in range(2, isqrt(q) + 1))]  # sqrt(23 * 5e6) = 10723; first run used 10000 (failed-run-1.log)

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

def coin_prob(signs, usable, G, T):
    """P(all usable signs +1) under the model with threshold T."""
    if G != 1: return 0.0
    if not all(signs[q] == 1 for q, _ in usable if q <= T): return 0.0
    large = [(q, e) for q, e in usable if q > T]
    if not large: return 1.0
    odd = sum(1 for _, e in large if e % 2); even = len(large) - odd
    return 0.5 ** even * (0.5 ** (odd - 1) if odd else 1.0)

def main():
    t0 = time.time()
    b = sieve(5 * 10**6)
    primes = [p for p in range(4 * 10**6 + 1, 5 * 10**6) if b[p] and p % 840 in SQ840]
    if len(primes) != 2001:
        raise ValueError(f"expected 2001 primes, got {len(primes)}")
    rows = []; obs_rows = []
    X1_viol = 0
    stats = {T: {"O": 0, "E": 0.0, "V": 0.0, "byL": defaultdict(lambda: [0, 0.0])} for T in (23, 100, 1000)}
    obs_large_le1 = 0; obs_total = 0
    per_d = defaultdict(lambda: [0, 0, 0])  # obstructions, G=+1 count, total
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
                usable = [(q, e) for q, e in sorted(f.items()) if gcd(q, m) == 1]
                rows.append({"p": p, "a": a, "u": u, "B": B, "M": M, "factor": [[q, e] for q, e in sorted(f.items())],
                             "witness": witness})
                for dp in sqfree_divisors(2 * a * u):
                    d = -dp
                    signs = {q: jacobi(d % q, q) for q, _ in usable}
                    G = 1
                    for q, e in usable: G *= signs[q] ** e
                    obs = all(signs[q] == 1 for q, _ in usable)
                    if obs and G != 1: X1_viol += 1
                    if obs:
                        obs_total += 1
                        obs_large_le1 += sum(1 for q, _ in usable if q > 1000) <= 1
                    pd = per_d[d]; pd[0] += obs; pd[1] += (G == 1); pd[2] += 1
                    for T, st in stats.items():
                        P = coin_prob(signs, usable, G, T)
                        st["O"] += obs; st["E"] += P; st["V"] += P * (1 - P)
                        k = min(sum(1 for q, _ in usable if q > T), 4)
                        st["byL"][k][0] += obs; st["byL"][k][1] += P
    if X1_viol:
        raise SystemExit(f"STOP: X1 violated {X1_viol} times")
    out = {"contract": "Research/Character20261008/contract.md @ 58fb041", "primes": len(primes), "rows": len(rows),
           "X1_violations": X1_viol}
    for T, st in stats.items():
        z = (st["O"] - st["E"]) / sqrt(st["V"]) if st["V"] > 0 else None
        out[f"T={T}"] = {"O": st["O"], "E": round(st["E"], 1), "O_over_E": round(st["O"] / st["E"], 4),
                         "z": round(z, 2) if z is not None else None,
                         "by_large_count": {k: [v[0], round(v[1], 1), round(v[0] / v[1], 4) if v[1] else None]
                                            for k, v in sorted(st["byL"].items())}}
    out["X2_pass"] = 0.98 <= out["T=1000"]["O_over_E"] <= 1.02 and abs(out["T=1000"]["z"]) < 3
    out["X3_pass"] = out["T=23"]["O_over_E"] < 1 and out["T=23"]["z"] < -3
    out["X5_share_obs_with_at_most_one_factor_above_1000"] = round(obs_large_le1 / obs_total, 4)
    out["X5_obstructions"] = obs_total
    top = sorted(per_d.items(), key=lambda kv: -kv[1][0])[:15]
    out["X6_top_d"] = {str(d): {"obstructions": v[0], "G_plus_share": round(v[1] / v[2], 4), "pairs": v[2]} for d, v in top}
    out["seconds"] = round(time.time() - t0, 1)
    with gzip.GzipFile("character_rows.json.gz", "wb", mtime=0) as fh:
        fh.write(json.dumps(rows).encode())
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
