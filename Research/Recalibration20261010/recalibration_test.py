"""D test (contract: Research/Recalibration20261010/contract.md, pushed as 440518e before computing).
J1 re-evaluated with the recalibration step inside the null distribution (Outcasts >>45)."""
import gzip, json, sys, time
from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "Stratified20261008"))
import stratified_test as st  # noqa: E402  (model_rows, pbin, by_p_from unchanged; also loads joint_test)
jt = st.jt

C3_OBS = {"T2": -161.705, "T3": -49.255, "T4": -5.963}
CHUNK = 50

def cell_C3(r):
    if r["method"] == "deterministic": return (r["method"], "det")
    return (r["method"], st.pbin(r["P"]) or f"P={int(r['P'])}")

def cell_box(r):
    return cell_C3(r) + (r["box"],)

def ratios(rows, cell):
    cal = defaultdict(lambda: [0, 0.0])
    for r in rows:
        k = cell(r); cal[k][0] += r["f"]; cal[k][1] += r["P"]
    return {k: (v[0] / v[1] if v[1] else 1.0) for k, v in cal.items()}

def apply(rows, ratio, cell, name):
    for r in rows:
        r[name] = r["P"] if r["method"] == "deterministic" else min(1.0, max(0.0, r["P"] * ratio.get(cell(r), 1.0)))

def arrays(rows, cell):
    order = sorted(range(len(rows)), key=lambda i: rows[i]["p"])  # stable: same grouping as analyze_sample
    rs = [rows[i] for i in order]
    ps = sorted({r["p"] for r in rs}); pid = {p: i for i, p in enumerate(ps)}
    gid = np.array([pid[r["p"]] for r in rs])
    keys = sorted({cell(r) for r in rs}, key=str); kid = {k: i for i, k in enumerate(keys)}
    cid = np.array([kid[cell(r)] for r in rs])
    P = np.array([r["P"] for r in rs]); F = np.array([r["f"] for r in rs], dtype=float)
    det = np.array([r["method"] == "deterministic" for r in rs])
    n_p = np.bincount(gid, minlength=len(ps)).astype(float)
    return rs, gid, len(ps), cid, len(keys), P, F, det, n_p

def refit(P, F, cid, C, det):
    """Procedure R applied to fails F (1-D or 2-D reps x N). Returns P'' with the same shape as F."""
    E = np.bincount(cid, weights=P, minlength=C)
    if F.ndim == 1:
        O = np.bincount(cid, weights=F, minlength=C)
        rat = np.where(E > 0, O / np.where(E > 0, E, 1.0), 1.0)
        Pp = np.clip(P * rat[cid], 0.0, 1.0)
    else:
        O = np.stack([np.bincount(cid, weights=row, minlength=C) for row in F])
        rat = np.where(E > 0, O / np.where(E > 0, E, 1.0), 1.0)
        Pp = np.clip(P[None, :] * rat[:, cid], 0.0, 1.0)
    return np.where(det, P, Pp)

def refit_bootstrap(rows, cell, reps, seed):
    rs, gid, G, cid, C, P, F, det, n_p = arrays(rows, cell)
    P1 = refit(P, F, cid, C, det)                     # P' fitted on the observed fails
    obs = jt.stats_arrays(gid, G, P1, F, n_p)
    rng = np.random.default_rng(seed); boot = defaultdict(list)
    for start in range(0, reps, CHUNK):
        r = min(CHUNK, reps - start)
        Fb = (rng.random((r, len(P1))) < P1).astype(float)   # null: independent rows with probabilities P'
        P2 = refit(P, Fb, cid, C, det)                      # redo procedure R on each resample
        s = jt.stats_arrays(gid, G, P2, Fb, n_p)
        for k in ("T2", "T3", "T4"): boot[k].append(np.asarray(s[k]))
    out = {"cells": C, "S_over_E_fitted": round(float((n_p.sum() - F.sum()) / (1 - P1).sum()), 4)}
    for k in ("T2", "T3", "T4"):
        b = np.concatenate(boot[k]); t = float(obs[k])
        lo = (1 + int((b <= t).sum())) / (reps + 1); hi = (1 + int((b >= t).sum())) / (reps + 1)
        out[k] = {"obs": round(t, 3), "boot_mean": round(float(b.mean()), 3), "boot_sd": round(float(b.std()), 3),
                  "p_two_sided": round(min(1.0, 2 * min(lo, hi)), 4), "p_low": round(lo, 4), "p_high": round(hi, 4)}
    return out, P1, rs

def toy_example():
    """>>45: two independent Bernoulli(1/2) rows recentred by the same-sample mean."""
    half = Fraction(1, 2); e_fit = Fraction(0); e_true = Fraction(0)
    for f1, f2 in product((0, 1), repeat=2):
        Ph = Fraction(f1 + f2, 2)
        e_fit += Fraction(1, 4) * (f1 - Ph) * (f2 - Ph)
        e_true += Fraction(1, 4) * (f1 - half) * (f2 - half)
    return {"E_refit": str(e_fit), "E_true_P": str(e_true)}

def main():
    t0 = time.time()
    out = {"contract": "Research/Recalibration20261010/contract.md @ 440518e"}
    conf = json.loads(gzip.open(jt_path("confirmatory_rows.json.gz")).read())
    rows = st.model_rows(conf, 23, 10000, 20261008)
    # D0: reproduce C3
    apply(rows, ratios(rows, cell_C3), cell_C3, "P_recal")
    rep = jt.analyze_sample(st.by_p_from(rows, "P_recal"), 2000, 20261012)  # identical to C3 (also D3's fixed-P' null)
    got = {k: round(rep[k]["obs"], 3) for k in C3_OBS}
    out["D0_reproduce"] = {"obs": got, "S_over_E": rep["S_over_E"], "expected": {**C3_OBS, "S_over_E": 1.0}}
    out["D0_pass"] = got == C3_OBS and rep["S_over_E"] == 1.0
    if not out["D0_pass"]:
        json.dump(out, sys.stdout, indent=1); raise SystemExit("STOP: D0 failed")
    print(f"[{round(time.time() - t0)}s] D0 ok", file=sys.stderr, flush=True)
    # D1: refit inside the bootstrap
    d1, P1, rs = refit_bootstrap(rows, cell_C3, 2000, 20261013)
    assert np.allclose(P1, [r["P_recal"] for r in rs])
    out["D1_refit_bootstrap"] = d1
    out["D1_pass"] = d1["T3"]["p_two_sided"] >= 0.01 and d1["T4"]["p_two_sided"] >= 0.01
    print(f"[{round(time.time() - t0)}s] D1 done", file=sys.stderr, flush=True)
    # D2: ratios from the exploratory sample, applied to the confirmatory rows
    expl = st.model_rows(jt.load_exploratory(), 23, 10000, 20261008)
    rex = ratios(expl, cell_C3)
    apply(rows, rex, cell_C3, "P_ex")
    d2 = jt.analyze_sample(st.by_p_from(rows, "P_ex"), 2000, 20261014)
    out["D2_split_calibration"] = {k: {"obs": round(d2[k]["obs"], 3), "Z": round(d2[k]["Z"], 2), "boot_mean": round(d2[k]["boot_mean"], 3),
                                       "boot_sd": round(d2[k]["boot_sd"], 3), "p_two_sided": round(d2[k]["p_two_sided"], 4)}
                                   for k in ("T2", "T3", "T4")}
    out["D2_split_calibration"]["S_over_E"] = d2["S_over_E"]
    out["D2_ratios_exploratory"] = {"|".join(k): round(v, 4) for k, v in sorted(rex.items(), key=str)}
    out["D2_ratios_confirmatory"] = {"|".join(k): round(v, 4) for k, v in sorted(ratios(rows, cell_C3).items(), key=str)}
    out["D2_pass"] = d2["T3"]["p_two_sided"] >= 0.01 and d2["T4"]["p_two_sided"] >= 0.01
    print(f"[{round(time.time() - t0)}s] D2 done", file=sys.stderr, flush=True)
    # D3: size of the refit effect
    out["D3_fixed_vs_refit_null"] = {k: {"fixed_P_prime": {"boot_mean": round(rep[k]["boot_mean"], 3), "boot_sd": round(rep[k]["boot_sd"], 3),
                                                           "p_two_sided": round(rep[k]["p_two_sided"], 4)},
                                         "refit": {"boot_mean": d1[k]["boot_mean"], "boot_sd": d1[k]["boot_sd"], "p_two_sided": d1[k]["p_two_sided"]}}
                                     for k in ("T2", "T3", "T4")}
    out["D3_toy_example"] = toy_example()
    # D4: finer cells (method x P-bin x box), refit bootstrap
    d4, _, _ = refit_bootstrap(rows, cell_box, 2000, 20261015)
    out["D4_box_cells_refit_bootstrap"] = d4
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, sys.stdout, indent=1)

def jt_path(name):
    return Path(jt.__file__).resolve().parent / name

if __name__ == "__main__":
    main()
