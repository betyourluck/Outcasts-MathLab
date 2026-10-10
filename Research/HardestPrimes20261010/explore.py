"""Exploration on published rows (X, Y, E, S, Q tests; 4e6 < p < 1.3e7): the least quadratic non-residue n2(p) of
hard primes versus the number C(p) of successful coprime pairs in the box a*u <= ceil(log2 p). Not a judgement."""
import gzip, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

R = Path(__file__).resolve().parent.parent
FILES = {"X": R / "Character20261008/character_rows.json.gz", "Y": R / "MediumSign20261008/mediumsign_rows.json.gz",
         "E": R / "Exponent20261008/exponent_rows.json.gz", "S": R / "Simultaneity20261010/simultaneity_rows.json.gz",
         "Q": R / "QR11Bit20261010/qr11_rows.json.gz"}

def small_primes(n):
    return [q for q in range(2, n) if all(q % r for r in range(2, int(q ** .5) + 1))]
PR = small_primes(400)

def n2(p):
    """Least prime quadratic non-residue mod p. Hard primes have p = 1 (mod 8), so 2 is a residue; for odd q and
    p = 1 (mod 4), (q|p) = (p|q) by reciprocity, computed by Euler's criterion mod q."""
    assert p % 8 == 1
    for q in PR[1:]:
        if pow(p % q, (q - 1) // 2, q) == q - 1: return q
    return None

def load():
    byp = defaultdict(dict); rng = {}
    for k, f in FILES.items():
        for r in json.load(gzip.open(f)):
            byp[r["p"]][(r["a"], r["u"])] = bool(r["witness"]); rng[r["p"]] = k
    return byp, rng

def main():
    byp, rng = load()
    rows = [{"p": p, "range": rng[p], "n2": n2(p), "C": sum(d.values()), "pairs": len(d)} for p, d in byp.items()]
    out = {"primes": len(rows)}
    by = defaultdict(list)
    for r in rows: by[r["n2"]].append(r["C"])
    out["C_by_n2"] = {int(k): {"primes": len(v), "mean_C": round(float(np.mean(v)), 2), "min_C": int(min(v))} for k, v in sorted(by.items())}
    Cs = np.array([r["C"] for r in rows]); n2s = np.array([r["n2"] for r in rows])
    out["spearman_like"] = float(np.corrcoef(np.argsort(np.argsort(Cs)), np.argsort(np.argsort(n2s)))[0, 1])
    for frac in (0.01, 0.05):
        thr = np.quantile(Cs, frac); low = n2s[Cs <= thr]
        out[f"bottom_{frac}"] = {"C_threshold": float(thr), "primes": int(len(low)),
                                 "share_n2_ge_17": round(float(np.mean(low >= 17)), 3), "share_n2_ge_17_all": round(float(np.mean(n2s >= 17)), 3),
                                 "share_n2_ge_29": round(float(np.mean(low >= 29)), 3), "share_n2_ge_29_all": round(float(np.mean(n2s >= 29)), 3)}
    hardest = sorted(rows, key=lambda r: (r["C"], -r["n2"]))[:15]
    out["hardest_15"] = hardest
    out["n2_of_old_hard_primes"] = {p: n2(p) for p in (345601, 670849, 2521, 66529, 1740481, 5843041, 3678481, 8615161, 53722321)}
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
