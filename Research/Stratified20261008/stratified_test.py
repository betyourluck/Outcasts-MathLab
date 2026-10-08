"""C test (contract: Research/Stratified20261008/contract.md, pushed as 83f3bed before computing).
Stratified calibration of the J-test model on the J confirmatory rows (Outcasts >>42)."""
import gzip, json, random, sys, time
from collections import defaultdict
from math import erfc, gcd, sqrt
from pathlib import Path

JF = Path(__file__).resolve().parent.parent / "JointFailure20261008"
sys.path.insert(0, str(JF))
import joint_test as jt  # noqa: E402  (puts ../FailureAtlas20261005 on sys.path; model and statistics unchanged)
fa, h11 = jt.fa, jt.h11

P1 = {(a, u) for a in range(1, 7) for u in range(1, 7) if gcd(a, u) == 1}
P_BINS = [(0.0, 0.1), (0.1, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.0)]

def model_rows(rows, B, mc, seed):
    """Same loop as joint_test.model_table, but keeps per-row strata."""
    h11.MC = mc
    rng = random.Random(seed); cache = {}; out = []
    for r in rows:
        a, u, M, p = r["a"], r["u"], r["M"], r["n"]; m0 = 4 * a * u
        odd = {q: e for q, e in fa.factor(M).items() if q % 2 == 1}
        if any(e > 1 for q, e in odd.items() if q > B): continue
        F = tuple((q % m0, e) for q, e in sorted(odd.items()) if q <= B)
        large = [q for q in odd if q > B]; L = len(large)
        rp = 1
        for q in large: rp = rp * q % m0
        key = (m0, rp, F, L)
        if key not in cache: cache[key] = h11.model(m0, rp, F, L, rng)
        P = cache[key]
        out.append({"p": p, "a": a, "u": u, "P": P, "f": int(not r["ok"]), "L": L,
                    "method": "deterministic" if L <= 1 else ("exact" if L == 2 else "monte_carlo"),
                    "omega": len(fa.factor(M)), "box": "P1" if (a, u) in P1 else "H2prime"})
    return out

def pbin(P):
    for lo, hi in P_BINS:
        if lo < P <= hi or (hi == 1.0 and lo < P < 1.0): return f"({lo},{hi}]"
    return None

def tally(rows, key):
    agg = defaultdict(lambda: [0, 0.0, 0.0, 0])
    for r in rows:
        k = key(r)
        if k is None: continue
        A = agg[k]; A[0] += r["f"]; A[1] += r["P"]; A[2] += r["P"] * (1 - r["P"]); A[3] += 1
    return {k: {"rows": v[3], "O": v[0], "E": round(v[1], 1), "O_over_E": round(v[0] / v[1], 4) if v[1] else None,
                "z": round((v[0] - v[1]) / sqrt(v[2]), 2) if v[2] > 0 else None} for k, v in sorted(agg.items(), key=lambda kv: str(kv[0]))}

def z_star(alpha_two_sided):
    lo, hi = 0.0, 10.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if erfc(mid / sqrt(2)) > alpha_two_sided: lo = mid
        else: hi = mid
    return hi

def strata(rows):
    nd = [r for r in rows if r["method"] != "deterministic"]
    det = [r for r in rows if r["method"] == "deterministic"]
    res = {"rows": len(rows), "deterministic_rows": len(det),
           "deterministic_mismatches": sum(1 for r in det if r["P"] != float(r["f"]))}
    res["by_method"] = tally(rows, lambda r: r["method"])
    res["by_P_bin"] = tally(nd, lambda r: pbin(r["P"]))
    res["by_omega"] = tally(nd, lambda r: min(r["omega"], 5))
    res["by_box"] = tally(nd, lambda r: r["box"])
    return res

def by_p_from(rows, Pkey="P"):
    by_p = defaultdict(list)
    for r in rows: by_p[r["p"]].append((r[Pkey], r["f"]))
    return by_p

def main():
    t0 = time.time()
    conf = json.loads(gzip.open(JF / "confirmatory_rows.json.gz").read())
    rows = model_rows(conf, 23, 10000, 20261008)
    # C0: reproduce J
    rep = jt.analyze_sample(by_p_from(rows), 200, 20261009)
    out = {"contract": "Research/Stratified20261008/contract.md @ 83f3bed",
           "C0_reproduce": {"S_over_E": rep["S_over_E"], "T2_obs": round(rep["T2"]["obs"], 2),
                            "expected": {"S_over_E": 1.0115, "T2_obs": -146.44}}}
    out["C0_pass"] = rep["S_over_E"] == 1.0115 and round(rep["T2"]["obs"], 2) == -146.44
    if not out["C0_pass"]:
        json.dump(out, sys.stdout, indent=1); raise SystemExit("STOP: C0 failed")
    st = strata(rows)
    ext = [r for r in rows if r["method"] != "deterministic" and r["P"] in (0.0, 1.0)]
    out["nondeterministic_rows_with_P_0_or_1"] = {"rows": len(ext), "mismatches": sum(1 for r in ext if r["P"] != float(r["f"]))}
    out["C1_deterministic_mismatches"] = st["deterministic_mismatches"]
    out["C1_pass"] = st["deterministic_mismatches"] == 0
    out["strata"] = st
    cells = [v for name in ("by_P_bin", "by_omega", "by_box") for v in st[name].values() if v["z"] is not None]
    K = len(cells); zs = z_star(0.01 / K)
    flagged = [(name, k, v) for name in ("by_P_bin", "by_omega", "by_box") for k, v in st[name].items()
               if v["z"] is not None and abs(v["z"]) > zs]
    out["C2"] = {"K": K, "z_star": round(zs, 3), "flagged": [{"layer": n, "stratum": str(k), **v} for n, k, v in flagged]}
    out["C2_no_flag"] = not flagged
    # C3: recalibrated J1 (method x P-bin)
    def ckey(r):
        if r["method"] == "deterministic": return (r["method"], "det")
        return (r["method"], pbin(r["P"]) or f"P={int(r['P'])}")
    cal = defaultdict(lambda: [0, 0.0])
    for r in rows:
        k = ckey(r)
        cal[k][0] += r["f"]; cal[k][1] += r["P"]
    ratio = {k: (v[0] / v[1] if v[1] else 1.0) for k, v in cal.items()}
    for r in rows:
        k = ckey(r)
        r["P_recal"] = min(1.0, max(0.0, r["P"] * ratio[k])) if r["method"] != "deterministic" else r["P"]
    rec = jt.analyze_sample(by_p_from(rows, "P_recal"), 2000, 20261012)
    out["C3_recalibrated_J1"] = {k: {"obs": round(rec[k]["obs"], 3), "Z": round(rec[k]["Z"], 2), "p_two_sided": round(rec[k]["p_two_sided"], 4)}
                                 for k in ("T2", "T3", "T4")}
    out["C3_recalibrated_J1"]["S_over_E"] = rec["S_over_E"]
    out["C3_ratios"] = {f"{k[0]}|{k[1]}": round(v, 4) for k, v in sorted(ratio.items())}
    # C4: exploratory sample and B = 100
    out["C4_exploratory_B23"] = strata(model_rows(jt.load_exploratory(), 23, 10000, 20261008))
    out["C4_confirmatory_B100"] = strata(model_rows(conf, 100, 10000, 20261008))
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
