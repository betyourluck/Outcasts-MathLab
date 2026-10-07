"""H10' (contract: briefs/p1-h10prime-contract.md): small-prime-conditioned null model."""
import gzip, json, random, sys
from collections import defaultdict
from math import gcd, sqrt
from failure_atlas import factor, group
from exponent_atlas import inv_table, load

SEED = 20261005
MC = 300

def has_neg(ts, m0):
    neg = m0 - 1; S = {1}
    for t in ts:
        S = S | {s * t % m0 for s in S}
        if neg in S: return True
    return False

def model(m0, r, F, L, rng):
    """P(no subset product == -1) with fixed residues F and L large residues, uniform s.t. product == r."""
    G, _ = group(m0); inv = inv_table(m0); neg = m0 - 1
    pf = 1
    for t in F: pf = pf * t % m0
    rp = r * inv[pf] % m0
    if L == 0: return 0.0 if has_neg(F, m0) else 1.0
    if L == 1: return 0.0 if has_neg(F + (rp,), m0) else 1.0
    if L == 2:
        cnt = 0
        for t1 in G:
            if not has_neg(F + (t1, rp * inv[t1] % m0), m0): cnt += 1
        return cnt / len(G)
    fails = 0
    for _ in range(MC):
        ts = [rng.choice(G) for _ in range(L - 1)]
        prod = 1
        for t in ts: prod = prod * t % m0
        ts.append(rp * inv[prod] % m0)
        if not has_neg(F + tuple(ts), m0): fails += 1
    return fails / MC

def main():
    rows = load()
    sq = []
    for r in rows:
        a, u, M = r["a"], r["u"], r["M"]; m0 = 4 * a * u
        odd = {q: e for q, e in factor(M).items() if q % 2 == 1}
        if any(e > 1 for e in odd.values()): continue
        assert all(gcd(q, m0) == 1 for q in odd)
        Mp = 1
        for q in odd: Mp *= q
        r["m0"] = m0; r["odd"] = sorted(odd); r["r"] = Mp % m0; r["omega"] = len(odd)
        sq.append(r)
    S = {"sqfree_rows": len(sq)}
    for B in (0, 7, 23, 100):   # B = 0 is the original H10 (all residues uniform)
        rng = random.Random(SEED); cache = {}
        aggL = defaultdict(lambda: [0, 0.0, 0, 0.0])   # L -> rows, sum p, observed, sum p(1-p)
        aggW = defaultdict(lambda: [0, 0.0, 0, 0.0])   # omega (L>=2 only)
        strat = defaultdict(lambda: [0, 0.0, 0])        # (omega, has_small<=23) for B==0
        det_mismatch = 0
        for r in sq:
            F = tuple(q % r["m0"] for q in r["odd"] if q <= B)
            L = sum(q > B for q in r["odd"])
            key = (r["m0"], r["r"], F, L)
            if key not in cache: cache[key] = model(r["m0"], r["r"], F, L, rng)
            p = cache[key]; fail = int(not r["ok"])
            if L <= 1 and p != float(fail): det_mismatch += 1
            A = aggL[L]; A[0] += 1; A[1] += p; A[2] += fail; A[3] += p * (1 - p)
            if L >= 2:
                W = aggW[r["omega"]]; W[0] += 1; W[1] += p; W[2] += fail; W[3] += p * (1 - p)
            if B == 0 and r["omega"] in (2, 3):
                T = strat[(r["omega"], any(q <= 23 for q in r["odd"]))]; T[0] += 1; T[1] += p; T[2] += fail
        def rep(A):
            z = (A[2] - A[1]) / sqrt(A[3]) if A[3] > 0 else None
            return {"rows": A[0], "model": round(A[1], 1), "observed": A[2],
                    "ratio": round(A[2] / A[1], 4) if A[1] else None, "z_indep": round(z, 2) if z is not None else None}
        tot = [0, 0.0, 0, 0.0]
        for L, A in aggL.items():
            if L >= 2:
                for i in range(4): tot[i] += A[i]
        S[f"B={B}"] = {"cells": len(cache), "deterministic_mismatch(L<=1)": det_mismatch,
                       "by_L": {L: rep(A) for L, A in sorted(aggL.items())},
                       "L>=2_by_omega": {w: rep(A) for w, A in sorted(aggW.items())},
                       "L>=2_total": rep(tot)}
        if B == 0:
            S["H10prime_b_strata(B=0 model, omega 2-3)"] = {
                f"omega={w},has_small={h}": {"rows": T[0], "model": round(T[1], 1), "observed": T[2],
                                             "ratio": round(T[2] / T[1], 4)} for (w, h), T in sorted(strat.items())}
    json.dump(S, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
