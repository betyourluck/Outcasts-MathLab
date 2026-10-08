"""E test (contract: Research/Exponent20261008/contract.md, pushed as d9c4755 before computing).
Exponent shortfall versus a residue model: usable prime factors <= T fixed, larger ones (exponent 1) uniform in
(Z/4au)^x subject to their product; exact enumeration of the group. Exact integers; no randomness."""
import gzip, json, sys, time
from collections import Counter, defaultdict
from math import gcd, isqrt, log2, sqrt
import numpy as np

SQ840 = {1, 121, 169, 289, 361, 529}
SMALL = [q for q in range(2, 20000) if all(q % r for r in range(2, isqrt(q) + 1))]  # sqrt(23 * 7e6) = 12688
NAMES = ("success", "group", "exponent")

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

def closure(gens, m):
    S = {1 % m}; fr = [1 % m]
    while fr:
        nx = []
        for x in fr:
            for g in gens:
                y = x * g % m
                if y not in S: S.add(y); nx.append(y)
        fr = nx
    return S

def reachable(T_mult, m):
    S = {1 % m}
    for t, e in T_mult:
        pw = [1 % m]
        for _ in range(e): pw.append(pw[-1] * t % m)
        S = {s * x % m for s in S for x in pw}
    return S

_G = {}
def group(m):
    if m not in _G:
        G = [x for x in range(1, m) if gcd(x, m) == 1]
        _G[m] = (G, {x: pow(x, -1, m) for x in G})
    return _G[m]

def typ(T_mult, m):
    neg = m - 1
    if neg in reachable(T_mult, m): return "success"
    if neg not in closure(sorted({t for t, _ in T_mult}), m): return "group"
    return "exponent"

_cache = {}
def model_probs(m, fixed, n_large, r_large):
    key = (m, fixed, n_large, r_large)
    if key in _cache: return _cache[key]
    G, inv = group(m)
    one = lambda t: tuple(1.0 if n == t else 0.0 for n in NAMES)
    if n_large == 0:
        res = one(typ(list(fixed), m))
    elif n_large == 1:
        res = one(typ(list(fixed) + [(r_large, 1)], m))
    elif n_large == 2:
        c = Counter(typ(list(fixed) + [(t1, 1), (r_large * inv[t1] % m, 1)], m) for t1 in G)
        res = tuple(c[n] / len(G) for n in NAMES)
    elif n_large == 3:
        c = Counter(typ(list(fixed) + [(t1, 1), (t2, 1), (r_large * inv[t1 * t2 % m] % m, 1)], m) for t1 in G for t2 in G)
        res = tuple(c[n] / len(G) ** 2 for n in NAMES)
    else:
        raise ValueError(n_large)
    _cache[key] = res
    return res

def compare(rows_info, T):
    det_mis = 0; excluded = Counter(); agg = defaultdict(lambda: [0, 0.0, 0.0]); byL = defaultdict(lambda: defaultdict(lambda: [0, 0.0]))
    for info in rows_info:
        m, usable, obs = info["m"], info["usable"], info["type"]
        large = [(q, e) for q, e in usable if q > T]
        if any(e > 1 for _, e in large): excluded["squared_large_prime"] += 1; continue
        if len(large) >= 4: excluded["four_or_more_large"] += 1; continue
        fixed = tuple(sorted((q % m, e) for q, e in usable if q <= T))
        r = 1
        for q, _ in large: r = r * q % m
        P = model_probs(m, fixed, len(large), r)
        if len(large) <= 1:
            if P[NAMES.index(obs)] != 1.0: det_mis += 1
            continue
        for i, n in enumerate(NAMES):
            A = agg[n]; A[0] += (obs == n); A[1] += P[i]; A[2] += P[i] * (1 - P[i])
            B = byL[len(large)][n]; B[0] += (obs == n); B[1] += P[i]
    res = {"det_mismatch": det_mis, "excluded": dict(excluded),
           "rows_two_or_more_large": sum(agg[n][0] for n in NAMES)}
    for n in NAMES:
        O, E, V = agg[n]
        res[n] = {"O": O, "E": round(E, 1), "O_over_E": round(O / E, 4), "z": round((O - E) / sqrt(V), 2)}
    res["by_large_count"] = {L: {n: [v[0], round(v[1], 1), round(v[0] / v[1], 4) if v[1] else None] for n, v in d.items()}
                             for L, d in sorted(byL.items())}
    return res

def main():
    t0 = time.time()
    b = sieve(7 * 10**6)
    primes = [p for p in range(6 * 10**6 + 1, 7 * 10**6) if b[p] and p % 840 in SQ840]
    if len(primes) != 1977:
        raise ValueError(f"expected 1977 primes, got {len(primes)}")
    rows = []; info = []
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
                t = typ([(q % m, e) for q, e in usable], m)
                if (t == "success") != witness:
                    raise SystemExit(f"STOP: type and witness disagree at {(p, a, u)}")
                rows.append({"p": p, "a": a, "u": u, "B": B, "M": M, "factor": [[q, e] for q, e in sorted(f.items())],
                             "witness": witness, "type": t})
                info.append({"a": a, "u": u, "m": m, "usable": usable, "type": t})
    out = {"contract": "Research/Exponent20261008/contract.md @ d9c4755", "primes": len(primes), "rows": len(rows),
           "type_counts": dict(Counter(r["type"] for r in rows))}
    for T in (100, 1000):
        out[f"T={T}"] = compare(info, T)
    if out["T=100"]["det_mismatch"] or out["T=1000"]["det_mismatch"]:
        out["STOP"] = "E1 deterministic mismatch"
    t100 = out["T=100"]
    out["E1_pass"] = out["T=100"]["det_mismatch"] == 0 and out["T=1000"]["det_mismatch"] == 0
    out["E2_pass"] = t100["exponent"]["O_over_E"] > 1 and t100["exponent"]["z"] > 3
    out["E3_pass"] = t100["group"]["O_over_E"] < 1 and t100["group"]["z"] < -3
    # E4 and E6
    gp = [i for i in info if i["type"] != "group"]
    def share(sel):
        sel = list(sel); return [len(sel), round(sum(i["type"] == "exponent" for i in sel) / len(sel), 4) if sel else None]
    below = share(i for i in gp if len(i["usable"]) < log2(len(group(i["m"])[0])))
    above = share(i for i in gp if len(i["usable"]) >= log2(len(group(i["m"])[0])))
    out["E4"] = {"k<log2|G|": below, "k>=log2|G|": above, "ratio": round(below[1] / above[1], 2)}
    out["E4_pass"] = below[1] >= 3 * above[1]
    out["E6_by_k"] = {k: share(i for i in gp if len(i["usable"]) == k) for k in sorted({len(i["usable"]) for i in gp})}
    out["E6_by_G"] = {g: share(i for i in gp if len(group(i["m"])[0]) == g) for g in sorted({len(group(i["m"])[0]) for i in gp})}
    bins = defaultdict(list)
    for i in gp:
        d = len(i["usable"]) - log2(len(group(i["m"])[0]))
        bins[max(-4, min(3, int(np.floor(d))))].append(i)
    out["E6_by_k_minus_log2G_floor"] = {k: share(v) for k, v in sorted(bins.items())}
    byau = defaultdict(Counter)
    for i in info: byau[i["a"] * i["u"]][i["type"]] += 1
    out["E6_by_au"] = {au: {n: round(c[n] / sum(c.values()), 4) for n in NAMES} for au, c in sorted(byau.items())}
    out["seconds"] = round(time.time() - t0, 1)
    with gzip.GzipFile("exponent_rows.json.gz", "wb", mtime=0) as fh:
        fh.write(json.dumps(rows).encode())
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
