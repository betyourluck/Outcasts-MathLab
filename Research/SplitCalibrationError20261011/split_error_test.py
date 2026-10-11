"""V test (contract: Research/SplitCalibrationError20261011/contract.md, pushed before computing).
Estimation error of the split-sample calibration D2 (Outcasts >>52 section 3): first-order size of the
shift of E[T2] and a two-stage parametric bootstrap that re-estimates the cell ratios on each resample."""
import gzip, json, sys, time
from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "Recalibration20261010"))
import recalibration_test as rt  # noqa: E402  (cell_C3, ratios, apply unchanged; loads stratified_test and joint_test)
st, jt = rt.st, rt.jt

D2_OBS = {"T2": -160.186, "T3": -49.696, "T4": -3.575}
D2_NULL = {"T2": (-1.35, 86.886), "T3": (1.594, 78.514), "T4": (-0.702, 58.546)}
D2_S_OVER_E = 0.9995
REPS, SEED, CHUNK = 2000, 20261016, 50

def arrays(rows, kid, cell):
    """Rows sorted by p (same grouping as analyze_sample); cells indexed by the shared table kid."""
    order = sorted(range(len(rows)), key=lambda i: rows[i]["p"])
    rs = [rows[i] for i in order]
    ps = sorted({r["p"] for r in rs}); pid = {p: i for i, p in enumerate(ps)}
    gid = np.array([pid[r["p"]] for r in rs]); cid = np.array([kid[cell(r)] for r in rs])
    P = np.array([r["P"] for r in rs]); F = np.array([r["f"] for r in rs], dtype=float)
    det = np.array([r["method"] == "deterministic" for r in rs])
    n_p = np.bincount(gid, minlength=len(ps)).astype(float)
    return rs, gid, len(ps), cid, P, F, det, n_p

def ratio_from(P, F, cid, C):
    """Procedure R: r_c = sum_c f / sum_c P, 1 when the denominator is 0. F is 1-D or 2-D (reps x N)."""
    E = np.bincount(cid, weights=P, minlength=C)
    if F.ndim == 1: O = np.bincount(cid, weights=F, minlength=C)
    else: O = np.stack([np.bincount(cid, weights=row, minlength=C) for row in F])
    return np.where(E > 0, O / np.where(E > 0, E, 1.0), 1.0)

def scaled(P, rat, cid, det):
    """P^ex = clip(P * r_c) on non-deterministic rows, P on deterministic rows. rat is (C,) or (reps, C)."""
    if rat.ndim == 1: return np.where(det, P, np.clip(P * rat[cid], 0.0, 1.0))
    return np.where(det[None, :], P[None, :], np.clip(P[None, :] * rat[:, cid], 0.0, 1.0))

def toy(n_train):
    """>>52: n training rows and 2 test rows, all independent Bernoulli(1/2); h = training mean."""
    e_split = Fraction(0); e_plug = Fraction(0)
    w = Fraction(1, 2 ** (n_train + 2))
    for ys in product((0, 1), repeat=n_train):
        h = Fraction(sum(ys), n_train)
        for f1, f2 in product((0, 1), repeat=2):
            e_split += w * (f1 - h) * (f2 - h)
        # plug-in: given h, F* ~ Bernoulli(h) independent, so E[(F1-h)(F2-h) | h] = 0 exactly
        e_plug += Fraction(1, 2 ** n_train) * (h - h) * (h - h)
    return {"n_train": n_train, "E_split": str(e_split), "predicted_1_over_4n": str(Fraction(1, 4 * n_train)), "E_plug_in": str(e_plug)}

def toy_same_sample():
    """>>45: two rows recentred by their own mean (D3 reproduced this as -1/8)."""
    e = Fraction(0)
    for f1, f2 in product((0, 1), repeat=2):
        h = Fraction(f1 + f2, 2); e += Fraction(1, 4) * (f1 - h) * (f2 - h)
    return str(e)

def first_order_shift(keys, tr, te):
    """V2: Delta E[T2] = sum_c Var(r_c) * sum over same-p same-cell pairs (i<j, non-deterministic) of P_i P_j."""
    rs_tr, _, _, cid_tr, P_tr, _, det_tr, _ = tr
    rs_te, _, _, cid_te, P_te, _, det_te, _ = te
    C = len(keys)
    rat = ratio_from(P_tr, tr[5], cid_tr, C); Pex_tr = scaled(P_tr, rat, cid_tr, det_tr)
    E = np.bincount(cid_tr, weights=P_tr, minlength=C)
    V = np.bincount(cid_tr, weights=np.where(det_tr, 0.0, Pex_tr * (1 - Pex_tr)), minlength=C)
    var_r = np.where(E > 0, V / np.where(E > 0, E, 1.0) ** 2, 0.0)
    groups = defaultdict(list)
    for r, c, d in zip(rs_te, cid_te, det_te):
        if not d: groups[(r["p"], int(c))].append(r["P"])
    S2 = np.zeros(C); npairs = np.zeros(C, dtype=int)
    for (p, c), Ps in groups.items():
        s1 = sum(Ps); s2 = sum(x * x for x in Ps); S2[c] += (s1 * s1 - s2) / 2; npairs[c] += len(Ps) * (len(Ps) - 1) // 2
    contrib = var_r * S2
    cells = {"|".join(keys[c]): {"train_rows": int((cid_tr == c).sum()), "E_c": round(float(E[c]), 1), "var_ratio": float(var_r[c]),
                                 "same_p_pairs_conf": int(npairs[c]), "sum_PiPj": round(float(S2[c]), 2), "contribution": round(float(contrib[c]), 4)}
             for c in range(C) if npairs[c] > 0}
    return float(contrib.sum()), cells, rat

def two_stage_bootstrap(keys, tr, te, Pex_te, reps, seed):
    """V3: resample the training fails, refit the ratios, rescale the test probabilities, resample the test fails."""
    _, _, _, cid_tr, P_tr, F_tr, det_tr, _ = tr
    _, gid_te, G, cid_te, P_te, F_te, det_te, n_p = te
    C = len(keys)
    rat = ratio_from(P_tr, F_tr, cid_tr, C); Pex_tr = scaled(P_tr, rat, cid_tr, det_tr)
    obs = jt.stats_arrays(gid_te, G, Pex_te, F_te, n_p)
    rng = np.random.default_rng(seed); boot = defaultdict(list); shift_rat = []
    for start in range(0, reps, CHUNK):
        r = min(CHUNK, reps - start)
        Ftr = (rng.random((r, len(P_tr))) < Pex_tr).astype(float)   # (a) training fails under the fitted null
        rat_b = ratio_from(P_tr, Ftr, cid_tr, C)                      # (b) procedure R redone on the resample
        Pex_b = scaled(P_te, rat_b, cid_te, det_te)
        Fte = (rng.random((r, len(P_te))) < Pex_te).astype(float)    # (c) test fails: the truth under the null is the fitted P^ex
        s = jt.stats_arrays(gid_te, G, Pex_b, Fte, n_p)               # (d) centred at the re-estimated P^ex*
        for k in ("T2", "T3", "T4"): boot[k].append(np.asarray(s[k]))
        shift_rat.append(rat_b - rat[None, :])
    out = {"reps": reps, "seed": SEED}
    for k in ("T2", "T3", "T4"):
        b = np.concatenate(boot[k]); t = float(obs[k])
        lo = (1 + int((b <= t).sum())) / (reps + 1); hi = (1 + int((b >= t).sum())) / (reps + 1)
        out[k] = {"obs": round(t, 3), "boot_mean": round(float(b.mean()), 3), "boot_sd": round(float(b.std()), 3),
                  "p_two_sided": round(min(1.0, 2 * min(lo, hi)), 4), "p_low": round(lo, 4), "p_high": round(hi, 4)}
    sr = np.concatenate(shift_rat)
    out["ratio_resample_sd"] = {"|".join(keys[c]): round(float(sr[:, c].std()), 5) for c in range(C) if sr[:, c].std() > 0}
    return out

def main():
    t0 = time.time()
    out = {"contract": "Research/SplitCalibrationError20261011/contract.md (pushed before computing)"}
    conf = json.loads(gzip.open(rt.jt_path("confirmatory_rows.json.gz")).read())
    rows = st.model_rows(conf, 23, 10000, 20261008)
    expl = st.model_rows(jt.load_exploratory(), 23, 10000, 20261008)
    print(f"[{round(time.time() - t0)}s] model rows: conf {len(rows)}, expl {len(expl)}", file=sys.stderr, flush=True)
    # V0: reproduce D2 exactly as recalibration_test did it
    rex = rt.ratios(expl, rt.cell_C3)
    rt.apply(rows, rex, rt.cell_C3, "P_ex")
    d2 = jt.analyze_sample(st.by_p_from(rows, "P_ex"), 2000, 20261014)
    got = {k: round(d2[k]["obs"], 3) for k in D2_OBS}
    null = {k: (round(d2[k]["boot_mean"], 3), round(d2[k]["boot_sd"], 3)) for k in D2_OBS}
    out["V0_reproduce_D2"] = {"obs": got, "null_mean_sd": null, "S_over_E": d2["S_over_E"],
                              "expected": {"obs": D2_OBS, "null_mean_sd": D2_NULL, "S_over_E": D2_S_OVER_E}}
    out["V0_pass"] = got == D2_OBS and null == D2_NULL and d2["S_over_E"] == D2_S_OVER_E
    if not out["V0_pass"]:
        json.dump(out, sys.stdout, indent=1); raise SystemExit("STOP: V0 failed")
    print(f"[{round(time.time() - t0)}s] V0 ok", file=sys.stderr, flush=True)
    # shared cell table and arrays; check the vectorised P^ex against recalibration_test.apply
    keys = sorted({rt.cell_C3(r) for r in expl} | {rt.cell_C3(r) for r in rows}, key=str); kid = {k: i for i, k in enumerate(keys)}
    tr = arrays(expl, kid, rt.cell_C3); te = arrays(rows, kid, rt.cell_C3)
    rat = ratio_from(tr[4], tr[5], tr[3], len(keys)); Pex_te = scaled(te[4], rat, te[3], te[6])
    assert np.allclose(Pex_te, [r["P_ex"] for r in te[0]])
    assert all(abs(rat[kid[k]] - v) < 1e-12 for k, v in rex.items())
    # V1: >>52's small example and its n-row generalisation
    out["V1_toy"] = {"split": [toy(n) for n in (2, 4, 8)], "same_sample_two_rows": toy_same_sample()}
    # V2: first-order shift of E[T2]
    shift, cells, _ = first_order_shift(keys, tr, te)
    sd2 = D2_NULL["T2"][1]
    out["V2_first_order_shift_T2"] = {"delta_E_T2": round(shift, 3), "delta_over_plug_in_sd": round(shift / sd2, 4),
                                      "Z_plug_in": round((D2_OBS["T2"] - D2_NULL["T2"][0]) / sd2, 3),
                                      "Z_shifted": round((D2_OBS["T2"] - (D2_NULL["T2"][0] + shift)) / sd2, 3),
                                      "cells": cells, "cells_in_table": len(keys)}
    print(f"[{round(time.time() - t0)}s] V2 done: delta {shift:.3f}", file=sys.stderr, flush=True)
    # V3: two-stage bootstrap
    v3 = two_stage_bootstrap(keys, tr, te, Pex_te, REPS, SEED)
    for k in ("T2", "T3", "T4"):
        v3[k]["mean_minus_plug_in"] = round(v3[k]["boot_mean"] - D2_NULL[k][0], 3)
        v3[k]["plug_in"] = {"boot_mean": D2_NULL[k][0], "boot_sd": D2_NULL[k][1]}
    out["V3_two_stage_bootstrap"] = v3
    out["V3_pass"] = v3["T3"]["p_two_sided"] >= 0.01 and v3["T4"]["p_two_sided"] >= 0.01
    print(f"[{round(time.time() - t0)}s] V3 done", file=sys.stderr, flush=True)
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
