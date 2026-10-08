"""M test (contract: Research/ManyFactor20261008/contract.md, pushed as c5f2b00 before computing).
For each hard prime p in (3*10^6, 4*10^6): how many pairs of the coprime log-budget box have
omega(ap+u) >= k, do the omega >= 4 pairs always keep a Type II witness, and how do three
candidate rules compare at equal size. Exact integers; two divisor methods; no randomness."""
import gzip, json, sys, time
from collections import Counter, defaultdict
from math import gcd, isqrt
import numpy as np

SQ840 = {1, 121, 169, 289, 361, 529}
SMALL = [q for q in range(2, 10000) if all(q % r for r in range(2, isqrt(q) + 1))]

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

def quantiles(xs):
    xs = sorted(xs)
    return {"min": xs[0], "p1": xs[len(xs) // 100], "median": xs[len(xs) // 2], "max": xs[-1]}

def select(rs, rule, m):
    if rule == "many_factor":
        return sorted(rs, key=lambda r: (-r["omega"], r["a"] * r["u"], r["a"]))[:m]
    if rule == "budget_order":
        return sorted(rs, key=lambda r: (r["a"] * r["u"], r["a"]))[:m]
    if rule == "rough":
        return sorted([r for r in rs if r["rough23"]], key=lambda r: (r["a"] * r["u"], r["a"]))[:m]
    raise ValueError(rule)

def main():
    t0 = time.time()
    b = sieve(4 * 10**6)
    primes = [p for p in range(3 * 10**6 + 1, 4 * 10**6) if b[p] and p % 840 in SQ840]
    if len(primes) != 2082:
        raise ValueError(f"expected 2082 primes, got {len(primes)}")
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
                nw = sum(1 for s in D if (s + 1) % m == 0)
                rows.append({"p": p, "a": a, "u": u, "B": B, "M": M, "factor": [[q, e] for q, e in sorted(f.items())],
                             "witness": nw > 0, "n_witness": nw, "omega": len(f),
                             "rough23": not any(M % q == 0 for q in (5, 7, 11, 13, 17, 19, 23))})
    per = defaultdict(list)
    for r in rows: per[r["p"]].append(r)
    out = {"contract": "Research/ManyFactor20261008/contract.md @ c5f2b00", "primes": len(primes), "rows": len(rows),
           "B_values": sorted({r["B"] for r in rows})}
    # M1
    zero4 = [p for p, rs in per.items() if not any(r["witness"] for r in rs if r["omega"] >= 4)]
    out["M1_Z4"] = len(zero4); out["M1_primes"] = zero4
    # M2 and M5
    for m in (10, 15, 20):
        z = {rule: sum(1 for rs in per.values() if not any(r["witness"] for r in select(rs, rule, m)))
             for rule in ("many_factor", "budget_order", "rough")}
        z["order_many<budget<rough"] = z["many_factor"] < z["budget_order"] < z["rough"]
        out[f"{'M2' if m == 10 else 'M5'}_m{m}"] = z
    # M3
    for k in (3, 4, 5):
        counts = {p: sum(1 for r in rs if r["omega"] >= k) for p, rs in per.items()}
        q = quantiles(list(counts.values()))
        q["argmin"] = sorted(p for p, c in counts.items() if c == q["min"])[:10]
        out[f"M3_M_p({k})"] = q
    # M4
    w4 = [sum(1 for r in rs if r["omega"] >= 4 and r["witness"]) for rs in per.values()]
    sel4 = [r for r in rows if r["omega"] >= 4]
    out["M4"] = {"witness_pairs_among_omega>=4": quantiles(w4), "hist": dict(sorted(Counter(w4).items())),
                 "rate_among_omega>=4": round(sum(r["witness"] for r in sel4) / len(sel4), 4), "rows_omega>=4": len(sel4)}
    out["all_box_min_witness_pairs"] = min(sum(r["witness"] for r in rs) for rs in per.values())
    out["seconds"] = round(time.time() - t0, 1)
    with gzip.GzipFile("manyfactor_rows.json.gz", "wb", mtime=0) as fh:
        fh.write(json.dumps(rows).encode())
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
