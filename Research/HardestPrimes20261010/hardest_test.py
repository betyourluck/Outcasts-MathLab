"""H test (contract: Research/HardestPrimes20261010/contract.md, pushed as 24990a2 before computing).
1.3e7 < p < 1.6e7 (hard classes), coprime a*u <= 24: n2(p) versus C(p), with the Proposition N check on every witness."""
import gzip, json, random, sys, time
from collections import defaultdict, Counter
from fractions import Fraction
from math import comb, gcd, sqrt, erfc
from pathlib import Path
import numpy as np
from explore import n2

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Simultaneity20261010"))
import simultaneity_test as st  # noqa: E402

SQ840 = {1, 121, 169, 289, 361, 529}
LO, HI = 13 * 10**6, 16 * 10**6

def avg_ranks(x):
    x = np.asarray(x, dtype=float); order = np.argsort(x, kind="mergesort"); r = np.empty(len(x)); i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]: j += 1
        r[order[i:j + 1]] = (i + j) / 2 + 1; i = j + 1
    return r

def fisher_one_sided_greater(a, b, c, d):
    """P(X >= a) for the hypergeometric with margins of [[a, b], [c, d]] (enrichment of the first cell)."""
    n1, n2_, k, N = a + b, c + d, a + c, a + b + c + d
    tot = comb(N, k); s = 0
    for x in range(a, min(n1, k) + 1): s += comb(n1, x) * comb(n2_, k - x)
    return float(Fraction(s, tot))

def S_of(n2v, primes):
    return sum(1 / q for q in primes if 11 <= q < n2v)

def main():
    t0 = time.time(); out = {"contract": "Research/HardestPrimes20261010/contract.md @ 24990a2"}
    small = st.sieve(20000); assert small[-1] ** 2 > 24 * HI + 24
    primes = [p for p in st.sieve(HI) if LO < p < HI and p % 840 in SQ840]
    budgets = sorted({(p - 1).bit_length() for p in primes}); B = max(budgets)
    P = np.array(primes, dtype=np.int64); rng = random.Random(20261022); checked = 0
    n2p = {p: n2(p) for p in primes}
    C = Counter(); wcount = 0; h0_bad = []; qstar = defaultdict(list); succ_au = defaultdict(Counter); rows = []
    for a in range(1, B + 1):
        for u in range(1, B // a + 1):
            if gcd(a, u) != 1: continue
            sel = np.array([(p - 1).bit_length() >= a * u for p in primes])
            Ps = P[sel]; m = 4 * a * u
            M, facs = st.factor_all(Ps, a, u, small)
            for p, Mi, f in zip(Ps.tolist(), M.tolist(), facs):
                w, wit = st.judge(Mi, f, m)
                if rng.random() < 0.01:
                    if st.direct(Mi, m) != w: raise SystemExit(f"STOP: direct scan disagrees at {(p, a, u)}")
                    checked += 1
                rows.append({"p": p, "a": a, "u": u, "M": Mi, "witness": w})
                if not w: continue
                C[p] += 1; succ_au[p][a * u] += 1
                for s in wit:  # H0: Proposition N on every witness
                    wcount += 1
                    nr = [q for q in f if s % q == 0 and q % 2 == 1 and pow(p % q, (q - 1) // 2, q) == q - 1]
                    if not nr or min(nr) < n2p[p]: h0_bad.append([p, a, u, s])
                s0 = wit[0]
                qstar[p].append(min(q for q in f if s0 % q == 0 and q % 2 == 1 and pow(p % q, (q - 1) // 2, q) == q - 1))
    out["input"] = {"primes": len(primes), "first": primes[0], "last": primes[-1], "budgets": budgets, "rows": len(rows), "direct_scan_checked_rows": checked}
    with gzip.GzipFile("hardest_rows.json.gz", "wb", mtime=0) as fh: fh.write(json.dumps(rows).encode())
    print(f"[{round(time.time() - t0)}s] rows {len(rows)}", file=sys.stderr, flush=True)
    out["H0"] = {"witnesses_checked": wcount, "successful_rows": sum(C.values()), "violations": h0_bad[:20], "violation_count": len(h0_bad)}
    out["H0_pass"] = not h0_bad
    if h0_bad:
        json.dump(out, sys.stdout, indent=1); raise SystemExit("STOP: H0 failed")
    Cv = np.array([C[p] for p in primes]); Nv = np.array([n2p[p] for p in primes])
    # H1
    r = float(np.corrcoef(avg_ranks(Cv), avg_ranks(Nv))[0, 1]); z = r * sqrt(len(primes) - 1)
    out["H1"] = {"spearman": r, "z": z, "p_one_sided": 0.5 * erfc(-z / sqrt(2))}
    out["H1_pass"] = r < 0 and out["H1"]["p_one_sided"] < 0.001
    # H2
    thr = float(np.quantile(Cv, 0.05)); low = Cv <= thr; hi = Nv >= 17
    a_, b_, c_, d_ = int((hi & low).sum()), int((~hi & low).sum()), int((hi & ~low).sum()), int((~hi & ~low).sum())
    OR = (a_ * d_) / (b_ * c_) if b_ * c_ else float("inf")
    out["H2"] = {"C_threshold": thr, "table": [[a_, b_], [c_, d_]], "OR": OR, "p_one_sided": fisher_one_sided_greater(a_, b_, c_, d_)}
    out["H2_pass"] = OR >= 2 and out["H2"]["p_one_sided"] < 0.001
    # H3
    sp = st.sieve(400); by = defaultdict(list)
    for p in primes: by[n2p[p]].append(C[p])
    xs = np.array([S_of(k, sp) for k in sorted(by)]); ys = np.array([np.mean(by[k]) for k in sorted(by)]); ws = np.array([len(by[k]) for k in sorted(by)])
    A = np.vstack([np.ones_like(xs), xs]).T; W = np.diag(ws); coef = np.linalg.solve(A.T @ W @ A, A.T @ W @ ys)
    out["H3"] = {"by_n2": {int(k): {"primes": int(w), "mean_C": round(float(y), 3), "pred_explore": round(14.61 - 11.65 * x, 2)} for k, x, y, w in zip(sorted(by), xs, ys, ws)},
                 "refit_a": float(coef[0]), "refit_kappa": float(-coef[1])}
    # H4
    hiN = [p for p in primes if n2p[p] >= 29]
    qs_hi = [q for p in hiN for q in qstar[p]]; qs_all = [q for p in primes for q in qstar[p]]
    out["H4"] = {"qstar_quantiles_all": [float(np.quantile(qs_all, x)) for x in (0.1, 0.25, 0.5, 0.75, 0.9)],
                 "qstar_quantiles_n2_ge_29": [float(np.quantile(qs_hi, x)) for x in (0.1, 0.25, 0.5, 0.75, 0.9)] if qs_hi else None,
                 "share_qstar_gt_1000_all": round(float(np.mean(np.array(qs_all) > 1000)), 4),
                 "share_qstar_gt_1000_n2_ge_29": round(float(np.mean(np.array(qs_hi) > 1000)), 4) if qs_hi else None,
                 "primes_n2_ge_29": len(hiN)}
    # H5
    out["H5_hardest"] = [[p, int(C[p]), int(n2p[p])] for p in sorted(primes, key=lambda p: (C[p], -n2p[p]))[:15]]
    # H6
    def au_dist(ps):
        tot = Counter()
        for p in ps: tot.update(succ_au[p])
        n = sum(tot.values()); return {"successes": n, "share_au_le_3": round(sum(v for k, v in tot.items() if k <= 3) / n, 4),
                                       "share_au_4_12": round(sum(v for k, v in tot.items() if 4 <= k <= 12) / n, 4),
                                       "share_au_13_24": round(sum(v for k, v in tot.items() if k >= 13) / n, 4)}
    out["H6"] = {"n2_ge_29": au_dist(hiN), "n2_eq_11": au_dist([p for p in primes if n2p[p] == 11])}
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, sys.stdout, indent=1, default=float)

if __name__ == "__main__":
    main()
