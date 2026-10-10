"""Exploration, third pass (the checks run in the terminal before writing the contract): where the residue-10 effect
sits (pairs with 11 | au), the 1-bit model after dropping those pairs, and the cross-range pair-composition model."""
import gzip, json, sys
from collections import defaultdict
import numpy as np
import qr_stats as qs
from explore import FILES

def drop(rows, q): return [r for r in rows if (r["a"] * r["u"]) % q != 0]

def residue10_by_pair(rows):
    acc = defaultdict(lambda: [0, 0])
    for r in rows:
        au = (r["a"], r["u"])
        if au in qs.SMALL or r["M"] % 11 == 0: continue
        A = acc[(au, r["p"] % 11 == 10)]; A[0] += (not r["witness"]); A[1] += 1
    eff = []
    for au in {a for a, _ in acc}:
        if acc[(au, True)][1] > 100 and acc[(au, False)][1] > 100:
            eff.append([round(acc[(au, True)][0] / acc[(au, True)][1] - acc[(au, False)][0] / acc[(au, False)][1], 4), list(au)])
    return sorted(eff)[:6]

def main():
    data = {k: json.load(gzip.open(FILES[k])) for k in "XYES"}
    allr = [r for v in data.values() for r in v]
    out = {"residue10_effect_most_negative_pairs": residue10_by_pair(allr)}
    out["anova_drop_11au"] = {k: qs.prime_anova(drop(v, 11), 11) for k, v in [("all", allr)] + list(data.items())}
    fit = drop([r for k in "XYE" for r in data[k]], 11); test = drop(data["S"], 11)
    out["pair_model_XYE_to_S"] = qs.residue_calibration(test, qs.fit_pair_model(fit, 11), 11)
    out["pair_model_S_to_XYE"] = qs.residue_calibration(fit, qs.fit_pair_model(test, 11), 11)
    out["pair_model_with_11au_XYE_to_S"] = qs.residue_calibration(data["S"], qs.fit_pair_model([r for k in "XYE" for r in data[k]], 11), 11)
    json.dump(out, sys.stdout, indent=1, default=float)

if __name__ == "__main__":
    main()
