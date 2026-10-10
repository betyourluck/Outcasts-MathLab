"""Q test (contract: Research/QR11Bit20261010/contract.md, pushed as 9935f4a before computing).
Rows for 1e7 < p < 1.3e7 (hard classes), coprime a*u <= ceil(log2 p) = 24; then the pre-registered statistics."""
import gzip, json, random, sys, time
from collections import defaultdict
from math import gcd
from pathlib import Path
import numpy as np
import qr_stats as qs
from explore import FILES

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Simultaneity20261010"))
import simultaneity_test as st  # noqa: E402  (sieve, factor_all, judge, direct unchanged)

SQ840 = {1, 121, 169, 289, 361, 529}
LO, HI = 10**7, 13 * 10**6
MODULI = (11, 13, 17, 19, 23)

def make_rows():
    small = st.sieve(18000)
    assert small[-1] ** 2 > 24 * HI + 24
    primes = [p for p in st.sieve(HI) if LO < p < HI and p % 840 in SQ840]
    budgets = sorted({(p - 1).bit_length() for p in primes})
    P = np.array(primes, dtype=np.int64); rows = []; rng = random.Random(20261021); checked = 0
    for a in range(1, max(budgets) + 1):
        for u in range(1, max(budgets) // a + 1):
            if gcd(a, u) != 1: continue
            sel = np.array([(p - 1).bit_length() >= a * u for p in primes])
            if not sel.any(): continue
            Ps = P[sel]; m = 4 * a * u
            M, facs = st.factor_all(Ps, a, u, small)
            for p, Mi, f in zip(Ps.tolist(), M.tolist(), facs):
                w, wit = st.judge(Mi, f, m)
                if rng.random() < 0.01:
                    if st.direct(Mi, m) != w: raise SystemExit(f"STOP: direct scan disagrees at {(p, a, u)}")
                    checked += 1
                rows.append({"p": p, "a": a, "u": u, "M": Mi, "witness": w})
    return primes, budgets, rows, checked

def residue_dev_range(rows, q, drop_q_au):
    acc = defaultdict(lambda: [0, 0])
    for r in rows:
        au = (r["a"], r["u"])
        if au in qs.SMALL or r["M"] % q == 0: continue
        if drop_q_au and (r["a"] * r["u"]) % q == 0: continue
        A = acc[(au, r["p"] % q)]; A[0] += (not r["witness"]); A[1] += 1
    dev = defaultdict(list)
    for au in {a for a, _ in acc}:
        rates = {c: acc[(au, c)][0] / acc[(au, c)][1] for c in range(1, q) if acc[(au, c)][1] > 20}
        if len(rates) < q - 2: continue
        mu = np.mean(list(rates.values()))
        for c, v in rates.items(): dev[c].append(v - mu)
    d = {c: float(np.mean(v)) for c, v in dev.items()}
    lo = min(d, key=d.get); hi = max(d, key=d.get)
    return {"min": [lo, round(d[lo], 4)], "max": [hi, round(d[hi], 4)]}

def main():
    t0 = time.time(); out = {"contract": "Research/QR11Bit20261010/contract.md @ 9935f4a"}
    primes, budgets, rows, checked = make_rows()
    out["input"] = {"primes": len(primes), "first": primes[0], "last": primes[-1], "budgets": budgets, "rows": len(rows), "direct_scan_checked_rows": checked}
    with gzip.GzipFile("qr11_rows.json.gz", "wb", mtime=0) as fh: fh.write(json.dumps(rows).encode())
    print(f"[{round(time.time() - t0)}s] rows {len(rows)}", file=sys.stderr, flush=True)
    # Q0
    idx = {(r["p"], r["a"], r["u"]): r for r in rows}
    bad = [p for p in primes if p % 11 == 10 and not ((p + 23) % 88 == 0 and idx[(p, 1, 22)]["witness"])]
    out["Q0"] = {"primes_checked": sum(1 for p in primes if p % 11 == 10), "mismatches": bad}
    out["Q0_pass"] = not bad
    if bad:
        json.dump(out, sys.stdout, indent=1); raise SystemExit("STOP: Q0 failed")
    # Q1, Q2
    out["Q1_Q2"] = {}
    for q in MODULI:
        r = qs.analyse_pair(rows, q, 1000, 20261021)
        out["Q1_Q2"][f"q{q}"] = r
    out["Q1_pass"] = all(v["log_ratio_neg_over_pos"]["est"] < 0 and v["log_ratio_neg_over_pos"]["p_one_sided_ge_0"] < 0.002 for v in out["Q1_Q2"].values())
    out["Q2_pass"] = all(v["psi_neg"]["OR_div_vs_nodiv"] <= 0.35 and v["psi_pos"]["OR_div_vs_nodiv"] >= 0.55 for v in out["Q1_Q2"].values())
    print(f"[{round(time.time() - t0)}s] Q1/Q2 done", file=sys.stderr, flush=True)
    # Q3
    fit = [r for k in "XYES" for r in json.load(gzip.open(FILES[k])) if (r["a"] * r["u"]) % 11 != 0]
    model = qs.fit_pair_model(fit, 11)
    test = [r for r in rows if (r["a"] * r["u"]) % 11 != 0]
    out["Q3"] = qs.residue_calibration(test, model, 11)
    out["Q3_pass"] = out["Q3"]["p"] >= 0.01
    # Q4-Q6 measurements
    out["Q4_one_bit_all_wider"] = qs.prime_anova(rows, 11)
    out["Q4_one_bit_drop_11au"] = qs.prime_anova(test, 11)
    acc = defaultdict(lambda: [0, 0])
    for r in rows:
        if (r["a"], r["u"]) in ((1, 22), (1, 11)):
            A = acc[(f"{r['a']},{r['u']}", r["p"] % 11)]; A[0] += (not r["witness"]); A[1] += 1
    out["Q5_pair_fail_by_residue"] = {k: {c: round(acc[(k, c)][0] / acc[(k, c)][1], 4) for c in range(1, 11) if acc[(k, c)][1]} for k in ("1,22", "1,11")}
    out["Q6_special_residues"] = {f"q{q}": {"with_q_au": residue_dev_range(rows, q, False), "without_q_au": residue_dev_range(rows, q, True)} for q in MODULI}
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, sys.stdout, indent=1, default=float)

if __name__ == "__main__":
    main()
