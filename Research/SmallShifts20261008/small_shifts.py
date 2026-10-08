"""Small shifts a*u <= 3 on the hard primes (exploratory; reads the published E-test rows, no new ap+u).
For each pair and each odd quadratic character chi mod 4au, the global sign chi(usable part of ap+u) is
compared with the elementary prediction from p = 1 (mod 24) and p mod 16."""
import gzip, json, sys
from collections import Counter, defaultdict
from math import gcd
from pathlib import Path

ROWS = Path(__file__).resolve().parent.parent / "Exponent20261008" / "exponent_rows.json.gz"
CHARS = {4: [-1], 8: [-1, -2], 12: [-1, -3]}   # odd quadratic characters mod m, as Kronecker (d | .)
PAIRS = [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1)]

def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0

def predicted(pair, d, p):
    """Elementary derivation from p = 1 (mod 24); see README."""
    if pair == (1, 1): return {-1: 1}[d]
    if pair in ((1, 2), (2, 1)): return {-1: -1, -2: 1}[d]
    if pair in ((1, 3), (3, 1)): return {-3: 1, -1: 1 if p % 16 == 1 else -1}[d]
    raise ValueError(pair)

def main():
    rows = json.loads(gzip.open(ROWS).read())
    per = defaultdict(dict)
    for r in rows:
        if (r["a"], r["u"]) in PAIRS: per[r["p"]][(r["a"], r["u"])] = r
    out = {"primes": len(per), "all_p_mod_24_is_1": all(p % 24 == 1 for p in per)}
    mism = 0; guaranteed = 0; table = defaultdict(Counter); allfail = 0; nfail = Counter()
    for p, d in per.items():
        g = False; fails = 0
        for pair, r in d.items():
            m = 4 * pair[0] * pair[1]
            usable = [(q, e) for q, e in r["factor"] if gcd(q, m) == 1]
            signs = {}
            for dd in CHARS[m]:
                G = 1
                for q, e in usable: G *= jacobi(dd % q, q) ** e
                signs[dd] = G
                if G != predicted(pair, dd, p): mism += 1
                table[f"{pair} d={dd} p%16={p % 16}"][G] += 1
            if all(v == -1 for v in signs.values()): g = True
            fails += not r["witness"]; nfail[pair] += not r["witness"]
        guaranteed += g; allfail += (fails == len(d))
    out["prediction_mismatches"] = mism
    out["primes_with_a_pair_guaranteed_by_global_signs"] = guaranteed
    out["global_sign_table"] = {k: dict(v) for k, v in sorted(table.items())}
    out["failures_per_pair"] = {str(k): v for k, v in nfail.items()}
    out["primes_failing_all_five"] = allfail
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
