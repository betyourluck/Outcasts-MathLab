"""H11 post-hoc (OUTSIDE the preregistered contract; labelled as such in the results).
R' = sum_p (S_p-E_p)^2 / sum_rows (ok-(1-P))^2 isolates within-p correlation from marginal miscalibration.
Also: rank of the two all-fail primes by Z_p, and B=100 for the calibration trend."""
import json, random, sys
from collections import defaultdict
from math import gcd
from failure_atlas import factor
from exponent_atlas import load
from h11_independence import model, SEED

def stats(rows, B):
    rng = random.Random(SEED); cache = {}; by_p = defaultdict(list)
    for r in rows:
        a, u, M, p = r["a"], r["u"], r["M"], r["n"]; m0 = 4 * a * u
        odd = {q: e for q, e in factor(M).items() if q % 2 == 1}
        if any(e > 1 for q, e in odd.items() if q > B): continue
        F = tuple((q % m0, e) for q, e in sorted(odd.items()) if q <= B)
        large = [q for q in odd if q > B]; L = len(large); rp = 1
        for q in large: rp = rp * q % m0
        key = (m0, rp, F, L)
        if key not in cache: cache[key] = model(m0, rp, F, L, rng)
        by_p[p].append((cache[key], int(r["ok"])))
    dev = res2 = V = S = E = 0.0; Z = {}
    for p, lst in by_p.items():
        Sp = sum(ok for _, ok in lst); Ep = sum(1 - P for P, _ in lst)
        dev += (Sp - Ep) ** 2; res2 += sum((ok - (1 - P)) ** 2 for P, ok in lst)
        V += sum(P * (1 - P) for P, _ in lst); S += Sp; E += Ep
        z = 1.0
        for P, _ in lst: z *= P
        Z[p] = z
    order = sorted(Z, key=lambda p: -Z[p])
    return {"S_over_E": round(S / E, 4), "R_model_variance": round(dev / V, 4),
            "R_prime_corr_only": round(dev / res2, 4), "sum_Z": round(sum(Z.values()), 2),
            "rank_by_Z(1=hardest)": {p: order.index(p) + 1 for p in (345601, 670849)},
            "Z": {p: round(Z[p], 5) for p in (345601, 670849)},
            "top5_Z": [(p, round(Z[p], 4)) for p in order[:5]]}

rows = load()
json.dump({f"B={B}": stats(rows, B) for B in (7, 23, 100)}, sys.stdout, indent=1)
