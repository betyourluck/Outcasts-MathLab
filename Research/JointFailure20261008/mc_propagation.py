"""Post-hoc (OUTSIDE the >>38 contract): propagate the Monte Carlo estimation error of the cell
probabilities P_i, shared by all rows in the same cell, into the null distribution of T_k.

Null with plug-in truth P = P_hat. Each replicate: f* ~ Bernoulli(P_hat) independently per pair;
for every Monte Carlo cell (L >= 3) draw P_hat*_c = Binomial(N_MC, P_hat_c) / N_MC once and give it to
all rows of that cell (shared error); exact/deterministic cells keep P_hat. Statistic T_k(f*, P_hat*),
compared with the observed T_k(f, P_hat). Also reports how much of the data sits in MC cells."""
import gzip, json, random, sys
from collections import defaultdict
from math import gcd, sqrt
import numpy as np
from joint_test import esym_from_power_sums, fa, h11  # joint_test puts ../FailureAtlas20261005 on sys.path

B, MC, SEED, REPS, BOOT_SEED = 23, 10000, 20261008, 2000, 20261011

def table(rows):
    h11.MC = MC
    rng = random.Random(SEED); cache = {}; keys = {}
    recs = []
    for r in rows:
        a, u, M, p = r["a"], r["u"], r["M"], r["n"]; m0 = 4 * a * u
        odd = {q: e for q, e in fa.factor(M).items() if q % 2 == 1}
        if any(e > 1 for q, e in odd.items() if q > B): continue
        F = tuple((q % m0, e) for q, e in sorted(odd.items()) if q <= B)
        large = [q for q in odd if q > B]; L = len(large); rp = 1
        for q in large: rp = rp * q % m0
        key = (m0, rp, F, L)
        if key not in cache:
            cache[key] = h11.model(m0, rp, F, L, rng); keys[key] = len(keys)
        recs.append((p, cache[key], int(not r["ok"]), keys[key], L >= 3))
    return recs, len(keys)

def T_stats(gid, G, P, F):
    d = F - P
    def gs(x):
        if x.ndim == 1: return np.bincount(gid, weights=x, minlength=G)
        return np.stack([np.bincount(gid, weights=row, minlength=G) for row in x])
    e2, e3, e4 = esym_from_power_sums(gs(d), gs(d**2), gs(d**3), gs(d**4))
    return {"T2": e2.sum(-1), "T3": e3.sum(-1), "T4": e4.sum(-1)}

def main():
    rows = json.load(gzip.open("confirmatory_rows.json.gz", "rt"))
    recs, ncells = table(rows)
    ps = sorted({r[0] for r in recs}); idx = {p: i for i, p in enumerate(ps)}
    gid = np.array([idx[r[0]] for r in recs]); G = len(ps)
    P = np.array([r[1] for r in recs]); F = np.array([r[2] for r in recs], float)
    cell = np.array([r[3] for r in recs]); is_mc = np.array([r[4] for r in recs])
    obs = {k: float(v) for k, v in T_stats(gid, G, P, F).items()}
    mc_cells = np.unique(cell[is_mc])
    Pc = np.zeros(ncells); Pc[cell] = P
    se = np.sqrt(Pc[mc_cells] * (1 - Pc[mc_cells]) / MC)
    rng = np.random.default_rng(BOOT_SEED)
    out = {"role": "post-hoc, outside the >>38 contract", "B": B, "mc": MC, "reps": REPS, "boot_seed": BOOT_SEED,
           "pairs": int(len(P)), "rows_in_mc_cells": int(is_mc.sum()), "mc_cells": int(len(mc_cells)), "cells": ncells,
           "max_cell_se": float(se.max()), "mean_cell_se_weighted_by_rows": float(np.sqrt(P[is_mc] * (1 - P[is_mc]) / MC).mean())}
    bootA = {k: [] for k in obs}; bootB = {k: [] for k in obs}
    for start in range(0, REPS, 100):
        r = min(100, REPS - start)
        Fb = (rng.random((r, len(P))) < P).astype(float)
        # A: P fixed (as in the contract)
        for k, v in T_stats(gid, G, P, Fb).items(): bootA[k].append(v)
        # B: shared MC error propagated
        Ts = {k: [] for k in obs}
        for j in range(r):
            Pstar_c = Pc.copy()
            Pstar_c[mc_cells] = rng.binomial(MC, Pc[mc_cells]) / MC
            Pstar = Pstar_c[cell]
            for k, v in T_stats(gid, G, Pstar, Fb[j]).items(): Ts[k].append(v)
        for k in obs: bootB[k].append(np.array(Ts[k]))
    for name, boot in (("P_fixed", bootA), ("MC_error_propagated", bootB)):
        res = {}
        for k in obs:
            b = np.concatenate(boot[k]); t = obs[k]
            lo = (1 + int((b <= t).sum())) / (REPS + 1); hi = (1 + int((b >= t).sum())) / (REPS + 1)
            res[k] = {"obs": t, "boot_mean": float(b.mean()), "boot_sd": float(b.std()), "p_two_sided": min(1.0, 2 * min(lo, hi))}
        out[name] = res
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
