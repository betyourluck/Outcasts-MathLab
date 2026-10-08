"""J test (contract: Outcasts >>38, briefs/p1-joint-contract.md): joint failures of pairs at the same p
versus the independence model with H11's marginal model. Exact integers for rows; MC only inside the model."""
import gzip, json, random, sys, time
from collections import defaultdict, Counter
from math import gcd, sqrt
from pathlib import Path
import numpy as np

ATLAS = Path(__file__).resolve().parent.parent / "FailureAtlas20261005"
sys.path.insert(0, str(ATLAS))
import failure_atlas as fa  # noqa: E402  (published with the atlas, unchanged)
import h11_independence as h11  # noqa: E402


def load_exploratory():
    """Same rows and order as exponent_atlas.load(): P1 primes, then H2'."""
    rows = json.load(gzip.open(ATLAS / "failure_atlas_rows.json.gz", "rt"))["primes"]
    rows += json.load(gzip.open(ATLAS / "h2prime_rows.json.gz", "rt"))
    return rows

SQ840 = {1, 121, 169, 289, 361, 529}
PAIRS_P1 = [(a, u) for a in range(1, 7) for u in range(1, 7) if gcd(a, u) == 1]
PAIRS_H2 = [(a, u) for u in (11, 13, 17, 19, 23) for a in range(1, 7)]
BOX = PAIRS_P1 + PAIRS_H2
assert len(BOX) == 53

def confirmatory_rows():
    b = fa.sieve(2 * 10**6)
    primes = [p for p in b if 10**6 < p < 2 * 10**6 and p % 840 in SQ840]
    assert len(primes) == 2149, len(primes)
    rows = []
    for p in primes:
        for a, u in BOX:
            r = fa.analyze(p, a, u)
            if not r["agree"]:
                raise SystemExit(f"STOP: two methods disagree at {(p, a, u)}")
            rows.append(r)
    return rows

def model_table(rows, B, mc, seed):
    """Per p: list of (P_fail, fail) using H11's model with B, mc samples/cell, seed. Returns (by_p, info)."""
    h11.MC = mc
    rng = random.Random(seed); cache = {}
    by_p = defaultdict(list); excluded = 0; det_mismatch = 0
    for r in rows:
        a, u, M, p = r["a"], r["u"], r["M"], r["n"]; m0 = 4 * a * u
        odd = {q: e for q, e in fa.factor(M).items() if q % 2 == 1}
        assert all(gcd(q, m0) == 1 for q in odd)
        if any(e > 1 for q, e in odd.items() if q > B): excluded += 1; continue
        F = tuple((q % m0, e) for q, e in sorted(odd.items()) if q <= B)
        large = [q for q in odd if q > B]; L = len(large)
        rp = 1
        for q in large: rp = rp * q % m0
        key = (m0, rp, F, L)
        if key not in cache: cache[key] = h11.model(m0, rp, F, L, rng)
        P = cache[key]; f = int(not r["ok"])
        if P in (0.0, 1.0) and P != float(f): det_mismatch += 1
        by_p[p].append((P, f))
    return by_p, {"B": B, "mc": mc, "seed": seed, "cells": len(cache), "excluded_pairs": excluded,
                  "det_mismatch": det_mismatch, "primes": len(by_p), "pairs": sum(len(v) for v in by_p.values())}

def esym_from_power_sums(p1, p2, p3, p4):
    e1 = p1
    e2 = (e1 * p1 - p2) / 2
    e3 = (e2 * p1 - e1 * p2 + p3) / 3
    e4 = (e3 * p1 - e2 * p2 + e1 * p3 - p4) / 4
    return e2, e3, e4

def stats_arrays(gid, ngroups, P, F, n_p):
    """F may be 2-D (reps, N). Returns dict of T2,T3,T4 and N0..N3 (arrays over reps)."""
    d = F - P
    def gsum(x):
        if x.ndim == 1: return np.bincount(gid, weights=x, minlength=ngroups)
        return np.stack([np.bincount(gid, weights=row, minlength=ngroups) for row in x])
    p1, p2, p3, p4 = gsum(d), gsum(d**2), gsum(d**3), gsum(d**4)
    e2, e3, e4 = esym_from_power_sums(p1, p2, p3, p4)
    fails = gsum(F.astype(float))
    S = n_p - fails
    out = {"T2": e2.sum(axis=-1), "T3": e3.sum(axis=-1), "T4": e4.sum(axis=-1)}
    for K in range(4): out[f"N{K}"] = (S <= K + 1e-9).sum(axis=-1)
    return out

def analyze_sample(by_p, reps, seed):
    ps = sorted(by_p)
    gid = np.concatenate([np.full(len(by_p[p]), i) for i, p in enumerate(ps)])
    P = np.array([P for p in ps for P, _ in by_p[p]])
    F = np.array([f for p in ps for _, f in by_p[p]], dtype=float)
    n_p = np.array([len(by_p[p]) for p in ps], dtype=float)
    G = len(ps)
    obs = stats_arrays(gid, G, P, F, n_p)
    # variances under independence: e_k of v per p
    v = P * (1 - P)
    q1 = np.bincount(gid, weights=v, minlength=G); q2 = np.bincount(gid, weights=v**2, minlength=G)
    q3 = np.bincount(gid, weights=v**3, minlength=G); q4 = np.bincount(gid, weights=v**4, minlength=G)
    ve2, ve3, ve4 = esym_from_power_sums(q1, q2, q3, q4)
    var = {"T2": ve2.sum(), "T3": ve3.sum(), "T4": ve4.sum()}
    # calibration, all-fail expectation, Poisson-binomial expectations for N_K
    S_sum = float((n_p - np.bincount(gid, weights=F, minlength=G)).sum()); E_sum = float((1 - P).sum())
    Z = np.ones(G)
    np.multiply.at(Z, gid, P)
    pmf_exp = np.zeros(int(n_p.max()) + 1)
    for i, p in enumerate(ps):
        pmf = np.array([1.0])
        for Pi, _ in by_p[p]:
            new = np.zeros(len(pmf) + 1); new[:-1] += pmf * Pi; new[1:] += pmf * (1 - Pi); pmf = new
        pmf_exp[:len(pmf)] += pmf
    # bootstrap
    rng = np.random.default_rng(seed)
    boot = defaultdict(list)
    chunk = 100
    for start in range(0, reps, chunk):
        r = min(chunk, reps - start)
        Fb = (rng.random((r, len(P))) < P).astype(float)
        st = stats_arrays(gid, G, P, Fb, n_p)
        for k, val in st.items(): boot[k].append(np.asarray(val))
    boot = {k: np.concatenate(v) for k, v in boot.items()}
    R = reps
    res = {"primes": G, "pairs": int(len(P)), "S_over_E": round(S_sum / E_sum, 4), "S_sum": S_sum, "E_sum": round(E_sum, 2),
           "sum_Z": round(float(Z.sum()), 3)}
    for k in ("T2", "T3", "T4"):
        t = float(obs[k]); b = boot[k]
        lo = (1 + int((b <= t).sum())) / (R + 1); hi = (1 + int((b >= t).sum())) / (R + 1)
        res[k] = {"obs": t, "Z": t / sqrt(var[k]), "sd_formula": sqrt(var[k]),
                  "boot_mean": float(b.mean()), "boot_sd": float(b.std()),
                  "p_two_sided": min(1.0, 2 * min(lo, hi)), "p_low": lo, "p_high": hi}
    for K in range(4):
        n = int(obs[f"N{K}"]); b = boot[f"N{K}"]
        res[f"N{K}"] = {"obs": n, "model_expect": round(float(pmf_exp[:K + 1].sum()), 2),
                        "boot_mean": float(b.mean()),
                        "p_low": (1 + int((b <= n).sum())) / (R + 1), "p_high": (1 + int((b >= n).sum())) / (R + 1)}
    res["all_fail_primes"] = [ps[i] for i in range(G) if (n_p[i] - np.bincount(gid, weights=F, minlength=G)[i]) == 0]
    S_obs = n_p - np.bincount(gid, weights=F, minlength=G)
    res["hist_obs"] = {int(k): int(v) for k, v in sorted(Counter(S_obs.astype(int)).items())}
    res["hist_model"] = {k: round(float(w), 1) for k, w in enumerate(pmf_exp) if w >= 0.05}
    return res

def main():
    t0 = time.time(); out = {"contract": "Outcasts >>38 (2026-10-08T00:43:34Z)", "reps": 2000, "boot_seed": 20261009}
    conf = confirmatory_rows()
    out["confirmatory_rows"] = len(conf); out["t_rows"] = round(time.time() - t0, 1)
    with gzip.GzipFile("confirmatory_rows.json.gz", "wb", mtime=0) as fh: fh.write(json.dumps(conf).encode())
    runs = [("confirmatory", conf, 23, 10000, 20261008, "MAIN"),
            ("confirmatory", conf, 7, 10000, 20261008, "sensitivity"),
            ("confirmatory", conf, 100, 10000, 20261008, "sensitivity"),
            ("confirmatory", conf, 23, 10000, 20261010, "sensitivity (model seed)"),
            ("exploratory", load_exploratory(), 23, 10000, 20261008, "measurement")]
    for name, rows, B, mc, seed, role in runs:
        t1 = time.time()
        by_p, info = model_table(rows, B, mc, seed)
        if info["det_mismatch"]:
            out[f"{name}_B{B}_seed{seed}"] = {"role": role, "info": info, "STOP": "deterministic mismatch"}
            json.dump(out, sys.stdout, indent=1); return
        res = analyze_sample(by_p, 2000, 20261009)
        res["info"] = info; res["role"] = role; res["seconds"] = round(time.time() - t1, 1)
        out[f"{name}_B{B}_seed{seed}"] = res
        print(f"[{round(time.time() - t0)}s] done {name} B={B} seed={seed}", file=sys.stderr, flush=True)
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
