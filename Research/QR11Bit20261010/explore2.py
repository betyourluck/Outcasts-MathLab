"""Exploration (second pass, after the first pass showed the within-pair confound): within-pair contrast split by
psi = (-au|q), and the per-prime nested F-test (residue mod 11 vs the QR bit). Same published rows."""
import gzip, json, sys
from pathlib import Path
import qr_stats as qs
from explore import FILES, R
def main():
    rows = []
    for f in FILES.values(): rows += json.load(gzip.open(f))
    out = {"rows": len(rows), "anova_q11": qs.prime_anova(rows, 11)}
    for q in (11, 13, 17, 19, 23): out[f"pair_q{q}"] = qs.analyse_pair(rows, q, 1000, 20261010)
    for k in ("X", "Y", "E", "S"):
        out[f"anova_q11_{k}"] = qs.prime_anova(json.load(gzip.open(FILES[k])), 11)
    json.dump(out, sys.stdout, indent=1, default=float)
if __name__ == "__main__":
    main()
