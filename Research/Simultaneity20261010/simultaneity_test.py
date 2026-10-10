"""S test (contract: Research/Simultaneity20261010/contract.md, pushed as 5d4fc3a before computing).
Rows for 7e6 < p < 1e7 (hard classes), coprime a*u <= ceil(log2 p); then the shared statistics."""
import gzip, json, random, sys, time
from math import gcd, isqrt
from pathlib import Path
import numpy as np
import simul_stats as ss

SQ840 = {1, 121, 169, 289, 361, 529}
LO, HI = 7 * 10**6, 10**7

def sieve(n):
    b = bytearray([1]) * (n + 1); b[0] = b[1] = 0
    for i in range(2, isqrt(n) + 1):
        if b[i]: b[i * i::i] = bytearray(len(b[i * i::i]))
    return [i for i in range(n + 1) if b[i]]

def factor_all(P, a, u, small):
    """Vectorised trial division of M = a*P + u by the primes in `small`; cofactor is 1 or prime."""
    M = a * P + u; R = M.copy(); facs = [dict() for _ in range(len(P))]
    for q in small:
        hit = np.nonzero(R % q == 0)[0]
        while hit.size:
            for i in hit: facs[i][q] = facs[i].get(q, 0) + 1
            R[hit] //= q
            hit = hit[R[hit] % q == 0]
    for i in np.nonzero(R > 1)[0]:
        facs[i][int(R[i])] = facs[i].get(int(R[i]), 0) + 1
    return M, facs

def judge(M, f, m):
    prod = 1
    for q, e in f.items(): prod *= q ** e
    if prod != M: raise SystemExit(f"STOP: factorisation product mismatch at M={M}")
    D = {1}
    for q, e in f.items():
        if gcd(q, m) != 1: continue
        r = q % m; new = set(D); pw = 1
        for _ in range(e):
            pw = pw * r % m; new |= {d * pw % m for d in D}
        D = new
    w1 = (m - 1) in D
    divs = [1]
    for q, e in f.items(): divs = [d * q ** k for d in divs for k in range(e + 1)]
    wit = sorted(s for s in divs if (s + 1) % m == 0)
    if w1 != bool(wit): raise SystemExit(f"STOP: two methods disagree at M={M}, m={m}")
    return w1, wit

def direct(M, m):
    for d in range(1, isqrt(M) + 1):
        if M % d == 0 and ((d + 1) % m == 0 or (M // d + 1) % m == 0): return True
    return False

def main():
    t0 = time.time(); out = {"contract": "Research/Simultaneity20261010/contract.md @ 5d4fc3a"}
    small = sieve(16000)
    assert small[-1] ** 2 > 24 * HI + 24 and 16001 ** 2 > 24 * HI + 24
    primes = [p for p in sieve(HI) if LO < p < HI and p % 840 in SQ840]
    P = np.array(primes, dtype=np.int64)
    budgets = sorted({(p - 1).bit_length() for p in primes})
    out["input"] = {"primes": len(primes), "first": primes[0], "last": primes[-1], "budgets": budgets}
    rows = []; rng = random.Random(20261020); direct_checked = 0
    for a in range(1, max(budgets) + 1):
        for u in range(1, max(budgets) // a + 1):
            if gcd(a, u) != 1: continue
            sel = np.array([(p - 1).bit_length() >= a * u for p in primes])
            if not sel.any(): continue
            Ps = P[sel]; m = 4 * a * u
            M, facs = factor_all(Ps, a, u, small)
            for p, Mi, f in zip(Ps.tolist(), M.tolist(), facs):
                w, wit = judge(Mi, f, m)
                if rng.random() < 0.01:
                    if direct(Mi, m) != w: raise SystemExit(f"STOP: direct scan disagrees at {(p, a, u)}")
                    direct_checked += 1
                rows.append({"p": p, "a": a, "u": u, "M": Mi, "factor": sorted(f.items()), "witness": w, "s_min": wit[0] if wit else None})
    out["rows"] = len(rows); out["direct_scan_checked_rows"] = direct_checked; out["t_rows"] = round(time.time() - t0, 1)
    with gzip.GzipFile("simultaneity_rows.json.gz", "wb", mtime=0) as fh: fh.write(json.dumps(rows).encode())
    print(f"[{round(time.time() - t0)}s] rows {len(rows)}", file=sys.stderr, flush=True)
    byp = ss.by_prime(rows)
    # S0: elementary check for p = 10, 8, 7 (mod 11)
    idx = {(r["p"], r["a"], r["u"]): r for r in rows}
    bad = []
    for p in byp:
        c = p % 11
        if c not in (10, 8, 7): continue
        s = {10: (1, 1), 8: (1, 3), 7: (3, 1)}[c]
        r = idx[(p,) + s]
        if not (r["M"] % 11 == 0 and (11 + 1) % (4 * s[0] * s[1]) == 0 and r["witness"]): bad.append(p)
        if all(not byp[p][t] for t in ss.SMALL): bad.append(p)
    out["S0"] = {"primes_checked": sum(1 for p in byp if p % 11 in (10, 8, 7)), "mismatches": bad}
    out["S0_pass"] = not bad
    if bad:
        json.dump(out, sys.stdout, indent=1); raise SystemExit("STOP: S0 failed")
    res, _ = ss.analyse(byp, 2000, 20261020)
    out["analysis"] = res
    out["S1_pass"] = res["mod11"]["delta"] > 0 and res["mod11"]["p_one_sided_delta_le_0"] < 0.01
    out["S2_pass"] = res["mod11"]["OR_ci95"][0] <= 1.0 <= res["mod11"]["OR_ci95"][1]
    out["S3_pass"] = res["specificity"]["delta11_minus_max_control"] > 0 and res["specificity"]["p_one_sided_le_0"] < 0.01
    out["S5_residue_table_mod11"] = ss.residue_table(byp, 11)
    out["S6_rescue"] = ss.rescue(byp)
    out["S7_single_exposure"] = ss.single_exposure(byp)
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, sys.stdout, indent=1, default=float)

if __name__ == "__main__":
    main()
