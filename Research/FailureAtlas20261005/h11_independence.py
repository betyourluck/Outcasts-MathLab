"""H11 (contract: briefs/p1-h11-contract.md): are failures of different pairs (a,u) at the same p independent?"""
import json, random, sys
from collections import defaultdict, Counter
from math import gcd, sqrt
from failure_atlas import factor, group
from exponent_atlas import inv_table, load

SEED = 20261005
MC = 1000

def has_neg(F, m0):
    """F: list of (residue, multiplicity). Any sub-multiset product == -1?"""
    neg = m0 - 1; S = {1}
    for t, v in F:
        pw, x = [1], 1
        for _ in range(v): x = x * t % m0; pw.append(x)
        S = {s * q % m0 for s in S for q in pw}
        if neg in S: return True
    return False

def model(m0, rp, F, L, rng):
    """P(fail) with fixed (residue, mult) list F and L large squarefree residues with product == rp."""
    G, _ = group(m0); inv = inv_table(m0)
    if L == 0: return 0.0 if has_neg(F, m0) else 1.0
    if L == 1: return 0.0 if has_neg(F + ((rp, 1),), m0) else 1.0
    if L == 2:
        cnt = sum(1 for t1 in G if not has_neg(F + ((t1, 1), (rp * inv[t1] % m0, 1)), m0))
        return cnt / len(G)
    fails = 0
    for _ in range(MC):
        ts = [rng.choice(G) for _ in range(L - 1)]
        prod = 1
        for t in ts: prod = prod * t % m0
        ts.append(rp * inv[prod] % m0)
        if not has_neg(F + tuple((t, 1) for t in ts), m0): fails += 1
    return fails / MC

def run(rows, B):
    rng = random.Random(SEED); cache = {}
    by_p = defaultdict(list); excluded = 0
    for r in rows:
        a, u, M, p = r["a"], r["u"], r["M"], r["n"]; m0 = 4 * a * u
        odd = {q: e for q, e in factor(M).items() if q % 2 == 1}
        assert all(gcd(q, m0) == 1 for q in odd)
        if any(e > 1 for q, e in odd.items() if q > B): excluded += 1; continue
        F = tuple((q % m0, e) for q, e in sorted(odd.items()) if q <= B)
        large = [q for q in odd if q > B]; L = len(large)
        rp = 1
        for q in large: rp = rp * q % m0
        key = (m0, rp, F, L)
        if key not in cache: cache[key] = model(m0, rp, F, L, rng)
        by_p[p].append((cache[key], int(r["ok"]), large))
    out = {"B": B, "cells": len(cache), "excluded_pairs": excluded, "primes": len(by_p)}
    S_sum = E_sum = V_sum = dev_sum = Z_sum = 0.0
    hist_obs = Counter(); hist_model = defaultdict(float); zeros = []; shared = 0; mismatch_det = 0
    for p, lst in by_p.items():
        S = sum(ok for _, ok, _ in lst); E = sum(1 - P for P, _, _ in lst); V = sum(P * (1 - P) for P, _, _ in lst)
        Z = 1.0
        for P, _, _ in lst: Z *= P
        for P, ok, _ in lst:
            if P in (0.0, 1.0) and (1 - P) != ok: mismatch_det += 1
        S_sum += S; E_sum += E; V_sum += V; dev_sum += (S - E) ** 2; Z_sum += Z
        hist_obs[S] += 1
        pmf = [1.0]
        for P, _, _ in lst:
            new = [0.0] * (len(pmf) + 1)
            for k, w in enumerate(pmf): new[k] += w * P; new[k + 1] += w * (1 - P)
            pmf = new
        for k, w in enumerate(pmf): hist_model[k] += w
        if S == 0: zeros.append((p, len(lst), round(Z, 4)))
        seen = Counter(q for _, _, large in lst for q in large)
        if any(c > 1 for c in seen.values()): shared += 1
    out["H11a_ratio_S_over_E"] = round(S_sum / E_sum, 4)
    out["det_mismatch"] = mismatch_det
    out["H11b_variance_ratio"] = round(dev_sum / V_sum, 4)
    out["H11c_zero_primes_observed"] = len(zeros)
    out["H11c_zero_primes_model"] = round(Z_sum, 2)
    out["H11c_zero_primes_model_pm2sd"] = [round(Z_sum - 2 * sqrt(Z_sum), 1), round(Z_sum + 2 * sqrt(Z_sum), 1)]
    out["H11c_zero_list(p, pairs, Z_p)"] = sorted(zeros)
    out["H11c_zero_with_Z_below_0.01"] = [z for z in zeros if z[2] < 0.01]
    out["H11d_hist_observed"] = {k: hist_obs[k] for k in sorted(hist_obs)}
    out["H11d_hist_model"] = {k: round(hist_model[k], 1) for k in sorted(hist_model) if hist_model[k] >= 0.05}
    out["H11e_primes_with_shared_large_factor"] = shared
    out["pairs_per_prime"] = dict(sorted(Counter(len(l) for l in by_p.values()).items()))
    return out

def main():
    rows = load()
    res = {"rows": len(rows)}
    for B in (7, 23): res[f"B={B}"] = run(rows, B)
    json.dump(res, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
