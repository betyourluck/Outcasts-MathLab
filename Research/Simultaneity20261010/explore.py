"""Exploratory run of the S test statistics on rows already published (X, Y, E tests, 4e6 < p < 7e6).
These are the data the hypotheses were formed on; nothing here is a judgement."""
import gzip, json, sys
from pathlib import Path
import simul_stats as ss

R = Path(__file__).resolve().parent.parent
FILES = {"X": R / "Character20261008/character_rows.json.gz", "Y": R / "MediumSign20261008/mediumsign_rows.json.gz",
         "E": R / "Exponent20261008/exponent_rows.json.gz"}

def main():
    rows = []
    for f in FILES.values(): rows += json.load(gzip.open(f))
    byp = ss.by_prime(rows)
    res, _ = ss.analyse(byp, 2000, 20261010)
    out = {"inputs": {k: str(v.relative_to(R)) for k, v in FILES.items()}, "pooled": res,
           "residue_table_mod11": ss.residue_table(byp, 11), "rescue": ss.rescue(byp), "single_exposure": ss.single_exposure(byp)}
    for k, f in FILES.items():
        b = ss.by_prime(json.load(gzip.open(f)))
        r, _ = ss.analyse(b, 500, 20261010, moduli=(11,))
        out[f"range_{k}"] = {"primes": r["primes"], "OR_raw": r["OR_raw"], "OR_raw_ci95": r["OR_raw_ci95"], "mod11": r["mod11"]}
    json.dump(out, sys.stdout, indent=1, default=float)

if __name__ == "__main__":
    main()
