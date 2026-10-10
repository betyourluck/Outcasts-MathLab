"""Exploration on rows already published (X, Y, E tests 4e6-7e6; S test 7e6-1e7).  Nothing here is a judgement."""
import gzip, json, sys
from pathlib import Path
import qr_stats as qs
R = Path(__file__).resolve().parent.parent
FILES = {"X": R / "Character20261008/character_rows.json.gz", "Y": R / "MediumSign20261008/mediumsign_rows.json.gz",
         "E": R / "Exponent20261008/exponent_rows.json.gz", "S": R / "Simultaneity20261010/simultaneity_rows.json.gz"}
def main():
    rows = []
    for f in FILES.values(): rows += json.load(gzip.open(f))
    out = {"inputs": {k: str(v.relative_to(R)) for k, v in FILES.items()}, "rows": len(rows)}
    for q in (11, 13, 17, 19, 23): out[f"q{q}"] = qs.analyse_q(rows, q, 1000, 20261010)
    json.dump(out, sys.stdout, indent=1, default=float)
if __name__ == "__main__":
    main()
