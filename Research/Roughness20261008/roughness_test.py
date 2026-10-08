"""R test (contract: Research/Roughness20261008/contract.md, pushed as c7338d0 before computing).
Is ap+u free of the primes 5..23 ("rough") more or less likely to carry a Type II witness?
Exact integers; two divisor methods; no randomness."""
import gzip, json, sys, time
from collections import Counter, defaultdict
from math import gcd, isqrt, erfc, sqrt
import numpy as np

SQ840 = {1, 121, 169, 289, 361, 529}

def sieve(n):
    b = bytearray([1]) * (n + 1); b[0] = b[1] = 0
    for i in range(2, isqrt(n) + 1):
        if b[i]: b[i * i::i] = bytearray(len(b[i * i::i]))
    return b

SMALL = [q for q in range(2, 10000) if all(q % r for r in range(2, isqrt(q) + 1))]

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

def cmh(rows, key_rough):
    strata = defaultdict(lambda: [0, 0, 0, 0])  # a: rough&wit, b: rough&no, c: not&wit, d: not&no
    for r in rows:
        t = strata[(r["a"], r["u"])]
        rough, w = r[key_rough], r["witness"]
        t[(0 if rough else 2) + (0 if w else 1)] += 1
    num = var = orn = ord_ = 0.0; used = 0
    for a, b, c, d in strata.values():
        n = a + b + c + d
        if n < 2 or (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0: continue
        used += 1
        num += a - (a + b) * (a + c) / n
        var += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
        orn += a * d / n; ord_ += b * c / n
    stat = num * num / var if var else float("nan")
    return {"strata": len(strata), "strata_used": used, "OR_MH": orn / ord_ if ord_ else float("inf"),
            "CMH_chi2": stat, "p_two_sided": erfc(sqrt(stat / 2)) if var else float("nan"),
            "sum_obs_minus_exp": num}

def main():
    t0 = time.time()
    b = sieve(3 * 10**6)
    primes = [p for p in range(2 * 10**6 + 1, 3 * 10**6) if b[p] and p % 840 in SQ840]
    assert len(primes) == 2109, len(primes)
    rows = []
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
                wit = [s for s in D if (s + 1) % m == 0]
                primes_M = sorted(f)
                rows.append({"p": p, "a": a, "u": u, "B": B, "M": M, "factor": [[q, e] for q, e in sorted(f.items())],
                             "witness": bool(wit), "n_witness": len(wit),
                             "rough23": not any(M % q == 0 for q in (5, 7, 11, 13, 17, 19, 23)),
                             "rough7": not any(M % q == 0 for q in (5, 7)),
                             "rough100": not any(M % q == 0 for q in SMALL if 5 <= q <= 97),
                             "L": sum(1 for q in primes_M if q > 23), "omega": len(primes_M),
                             "M_prime": f == {M: 1}})
    out = {"contract": "Research/Roughness20261008/contract.md @ c7338d0", "primes": len(primes), "rows": len(rows),
           "t_rows": round(time.time() - t0, 1)}
    out["R1_z23"] = cmh(rows, "rough23")
    out["sens_z7"] = cmh(rows, "rough7")
    out["sens_z100"] = cmh(rows, "rough100")
    rr = [r for r in rows if r["rough23"]]; nr = [r for r in rows if not r["rough23"]]
    out["pooled_rates"] = {"rough": [len(rr), round(sum(r["witness"] for r in rr) / len(rr), 4)],
                           "not_rough": [len(nr), round(sum(r["witness"] for r in nr) / len(nr), 4)]}
    # R2
    byL = defaultdict(lambda: [0, 0]); byW = defaultdict(lambda: [0, 0])
    for r in rows:
        k = min(r["L"], 4); byL[k][0] += 1; byL[k][1] += r["witness"]
        byW[r["omega"]][0] += 1; byW[r["omega"]][1] += r["witness"]
    out["R2_by_L(>=4 pooled)"] = {k: [v[0], round(v[1] / v[0], 4)] for k, v in sorted(byL.items())}
    out["R2_by_omega"] = {k: [v[0], round(v[1] / v[0], 4)] for k, v in sorted(byW.items())}
    out["R2_rows_L_le_1"] = sum(1 for r in rows if r["L"] <= 1)
    # R3
    rp = [r for r in rr if r["M_prime"]]
    out["R3"] = {"rough_rows": len(rr), "rough_and_M_prime": len(rp),
                 "frac_M_prime_among_rough": round(len(rp) / len(rr), 4),
                 "M_prime_with_a_dvd_u_plus_1": sum((r["u"] + 1) % r["a"] == 0 for r in rp),
                 "M_prime_with_witness": sum(r["witness"] for r in rp),
                 "M_prime_witness_without_a_dvd_u_plus_1": sum(r["witness"] and (r["u"] + 1) % r["a"] != 0 for r in rp)}
    # R4
    per = defaultdict(lambda: [0, 0, 0])
    for r in rows:
        t = per[r["p"]]; t[2] += r["witness"]
        if r["rough23"]: t[0] += 1; t[1] += r["witness"]
    Rp = [v[0] for v in per.values()]; Rpp = [v[1] for v in per.values()]; Np = [v[2] for v in per.values()]
    out["R4"] = {"R_p_min_median_max": [min(Rp), float(np.median(Rp)), max(Rp)],
                 "R_p_plus_min_median_max": [min(Rpp), float(np.median(Rpp)), max(Rpp)],
                 "N_p_plus_min_median_max": [min(Np), float(np.median(Np)), max(Np)],
                 "primes_R_p_plus_zero": sum(1 for x in Rpp if x == 0),
                 "primes_N_p_plus_zero (Hlog counterexamples)": sum(1 for x in Np if x == 0),
                 "R_p_plus_hist": dict(sorted(Counter(Rpp).items())),
                 "N_p_plus_hist": dict(sorted(Counter(Np).items()))}
    out["seconds"] = round(time.time() - t0, 1)
    with gzip.GzipFile("roughness_rows.json.gz", "wb", mtime=0) as fh:
        fh.write(json.dumps(rows).encode())
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
