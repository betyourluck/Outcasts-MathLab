"""Post-hoc (OUTSIDE the M contract): count only the prime factors of ap+u that are coprime to 4au
("usable": a witness s must be coprime to 4au, so 2 and the primes dividing au can never appear in s).
Same rows as manyfactor_test.py; nothing new is computed about divisors."""
import gzip, json, sys
from collections import Counter, defaultdict

rows = json.loads(gzip.open("manyfactor_rows.json.gz").read())
for r in rows:
    m = 4 * r["a"] * r["u"]
    r["omega_usable"] = sum(1 for q, _ in r["factor"] if m % q != 0)
per = defaultdict(list)
for r in rows: per[r["p"]].append(r)
out = {"role": "post-hoc, outside the M contract", "rows": len(rows)}
for key in ("omega", "omega_usable"):
    by = defaultdict(lambda: [0, 0])
    for r in rows:
        by[r[key]][0] += 1; by[r[key]][1] += r["witness"]
    out[f"rate_by_{key}"] = {k: [v[0], round(v[1] / v[0], 4)] for k, v in sorted(by.items())}
    for k in (2, 3, 4):
        counts = [sum(1 for r in rs if r[key] >= k) for rs in per.values()]
        zeros = [p for p, rs in per.items() if not any(r["witness"] for r in rs if r[key] >= k)]
        out[f"{key}>={k}"] = {"M_p_min": min(counts), "M_p_median": sorted(counts)[len(counts) // 2],
                             "primes_without_witness": len(zeros), "examples": zeros[:10]}
p = 3678481
out["p3678481"] = [(r["a"], r["u"], r["omega"], r["omega_usable"], r["witness"]) for r in per[p]
                   if r["omega"] >= 4 or r["witness"]]
json.dump(out, sys.stdout, indent=1)
