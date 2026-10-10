"""Shared statistics for the S test (small-shift wipe-out vs failures of the wider box).

Rows: (p, a, u, witness).  SMALL = the five shifts with a*u <= 3.  F_S(p) = all five fail.
Exposure = F_S(p); outcome = a wider pair (a*u >= 4) fails.  Mantel-Haenszel odds ratio over
strata = pair (raw) or (pair, p mod q) (stratified).  Uncertainty: bootstrap over primes
(the rows of one p are not independent), the same resample for every statistic (paired)."""
from collections import defaultdict
from math import log
import numpy as np

SMALL = ((1, 1), (1, 2), (2, 1), (1, 3), (3, 1))
CONTROL_MODULI = (13, 17, 19, 23)

def by_prime(rows):
    byp = defaultdict(dict)
    for r in rows: byp[r["p"]][(r["a"], r["u"])] = bool(r["witness"])
    return byp

def design(byp):
    """Arrays over the wider-pair rows: prime index, pair index, exposure, fail."""
    ps = sorted(byp); pid = {p: i for i, p in enumerate(ps)}
    FS = np.array([all(not byp[p][s] for s in SMALL) for p in ps])
    pairs = sorted({au for p in ps for au in byp[p] if au not in SMALL})
    kid = {au: i for i, au in enumerate(pairs)}
    P, K, Y = [], [], []
    for p in ps:
        for au, w in byp[p].items():
            if au in SMALL: continue
            P.append(pid[p]); K.append(kid[au]); Y.append(not w)
    P = np.array(P); K = np.array(K); Y = np.array(Y)
    return {"ps": np.array(ps), "FS": FS, "P": P, "K": K, "Y": Y, "E": FS[P], "npairs": len(pairs)}

def mh_or(d, strata, weights=None):
    """MH odds ratio with prime weights (bootstrap multiplicities).  strata: int array over rows."""
    w = np.ones(len(d["ps"])) if weights is None else weights
    rw = w[d["P"]]
    S = int(strata.max()) + 1
    E, Y = d["E"], d["Y"]
    a = np.bincount(strata, rw * (E & Y), S); b = np.bincount(strata, rw * (E & ~Y), S)
    c = np.bincount(strata, rw * (~E & Y), S); e = np.bincount(strata, rw * (~E & ~Y), S)
    n = a + b + c + e; ok = n > 0
    num = (a[ok] * e[ok] / n[ok]).sum(); den = (b[ok] * c[ok] / n[ok]).sum()
    return num / den if den > 0 else float("nan")

def strata_raw(d):
    return d["K"]

def strata_mod(d, q):
    return d["K"] * q + (d["ps"][d["P"]] % q)

def analyse(byp, reps, seed, moduli=(11,) + CONTROL_MODULI):
    d = design(byp)
    st = {"raw": strata_raw(d)}
    for q in moduli: st[q] = strata_mod(d, q)
    est = {k: mh_or(d, s) for k, s in st.items()}
    rng = np.random.default_rng(seed); G = len(d["ps"])
    boot = {k: np.empty(reps) for k in st}
    for b in range(reps):
        w = np.bincount(rng.integers(0, G, G), minlength=G).astype(float)
        for k, s in st.items(): boot[k][b] = mh_or(d, s, w)
    lr = np.log(boot["raw"])
    out = {"primes": G, "FS_primes": int(d["FS"].sum()), "wider_rows": int(len(d["Y"])), "reps": reps, "seed": seed,
           "OR_raw": est["raw"], "OR_raw_ci95": [float(np.quantile(boot["raw"], 0.025)), float(np.quantile(boot["raw"], 0.975))],
           "p_one_sided_raw_le_1": (1 + int((lr <= 0).sum())) / (reps + 1)}
    for q in moduli:
        delta = log(est["raw"]) - log(est[q]); bd = lr - np.log(boot[q])
        out[f"mod{q}"] = {"OR": est[q], "OR_ci95": [float(np.quantile(boot[q], 0.025)), float(np.quantile(boot[q], 0.975))],
                          "delta": delta, "delta_ci95": [float(np.quantile(bd, 0.025)), float(np.quantile(bd, 0.975))],
                          "p_one_sided_delta_le_0": (1 + int((bd <= 0).sum())) / (reps + 1)}
    # specificity: delta_11 - max over control moduli
    if 11 in moduli and all(q in moduli for q in CONTROL_MODULI):
        b11 = lr - np.log(boot[11]); bmax = np.max([lr - np.log(boot[q]) for q in CONTROL_MODULI], axis=0)
        diff = b11 - bmax
        out["specificity"] = {"delta11_minus_max_control": out["mod11"]["delta"] - max(out[f"mod{q}"]["delta"] for q in CONTROL_MODULI),
                              "ci95": [float(np.quantile(diff, 0.025)), float(np.quantile(diff, 0.975))],
                              "p_one_sided_le_0": (1 + int((diff <= 0).sum())) / (reps + 1)}
    return out, d

def residue_table(byp, q=11):
    rows = {}
    for c in range(q):
        ps = [p for p in byp if p % q == c]
        if not ps: continue
        fs = sum(all(not byp[p][s] for s in SMALL) for p in ps)
        wf = sum(1 - w for p in ps for au, w in byp[p].items() if au not in SMALL)
        wn = sum(1 for p in ps for au in byp[p] if au not in SMALL)
        rows[c] = {"primes": len(ps), "FS_rate": round(fs / len(ps), 4), "wider_fail_rate": round(wf / wn, 4),
                   "small_pairs_divisible_by_q": [list(s) for s in SMALL if (s[0] * c + s[1]) % q == 0]}
    return rows

def rescue(byp):
    fs = [p for p in byp if all(not byp[p][s] for s in SMALL)]
    succ = sorted(sum(1 for au, w in byp[p].items() if au not in SMALL and w) for p in fs)
    other = [sum(1 for au, w in byp[p].items() if au not in SMALL and w) for p in byp if p not in set(fs)]
    return {"FS_primes": len(fs), "min_wider_successes": succ[0] if succ else None,
            "mean_wider_successes_FS": round(sum(succ) / len(succ), 3) if succ else None,
            "mean_wider_successes_other": round(sum(other) / len(other), 3),
            "primes_without_any_witness": sum(1 for p in byp if not any(byp[p].values()))}

def single_exposure(byp, seed=0):
    """MH OR (strata = pair) with exposure = one small pair fails, raw and stratified by p mod 11."""
    out = {}
    for s in SMALL:
        ps = sorted(byp); pid = {p: i for i, p in enumerate(ps)}
        d = design(byp)
        d["E"] = np.array([not byp[p][s] for p in ps])[d["P"]]
        out[f"{s[0]},{s[1]}"] = {"raw": round(mh_or(d, strata_raw(d)), 4), "mod11": round(mh_or(d, strata_mod(d, 11)), 4)}
    return out
