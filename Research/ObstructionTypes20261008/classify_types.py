"""Exploratory classification of Type II failures (plan.md). Reads the published M-test rows;
no new ap+u is computed. Definitions follow FailureAtlas20261005/failure_atlas.py (analyze)."""
import gzip, json, sys
from collections import Counter, defaultdict
from math import gcd
from pathlib import Path

ROWS = Path(__file__).resolve().parent.parent / "ManyFactor20261008" / "manyfactor_rows.json.gz"
FOUR = [3613009, 3678481, 3900649, 3973729]

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

def kronecker(d, q):
    """(d | q) for odd prime q."""
    return jacobi(d % q, q)

def closure(gens, m):
    S = {1 % m}; fr = [1 % m]
    while fr:
        nx = []
        for x in fr:
            for g in gens:
                y = x * g % m
                if y not in S: S.add(y); nx.append(y)
        fr = nx
    return S

def reachable(T_mult, m, k=None):
    """Residues of products q1^f1 ... with 0 <= fi <= (k if k is not None else ei)."""
    S = {1 % m}
    for t, e in T_mult:
        top = e if k is None else k
        pw = [1 % m]
        for _ in range(top): pw.append(pw[-1] * t % m)
        S = {s * x % m for s in S for x in pw}
    return S

def sqfree_divisors(n):
    ds = [1]
    for q in [q for q in range(2, n + 1) if n % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]:
        ds += [d * q for d in ds]
    return ds

def classify(r):
    a, u = r["a"], r["u"]; m = 4 * a * u; neg = m - 1
    usable = [(q % m, e) for q, e in r["factor"] if gcd(q, m) == 1]
    T = sorted({t for t, _ in usable})
    G = [x for x in range(1, m) if gcd(x, m) == 1]
    G2 = {x * x % m for x in G}
    H = closure(T, m)
    HG2 = closure(T + sorted(G2), m)
    out = {"omega_usable": len(usable)}
    if neg not in HG2:
        qs = [q for q, _ in r["factor"] if gcd(q, m) == 1]
        ds = [-d for d in sqfree_divisors(2 * a * u)]
        out["type"] = "QR"
        out["d"] = [d for d in ds if all(kronecker(d, q) == 1 for q in qs)]
    elif neg not in H:
        out["type"] = "subgroup"
    else:
        out["type"] = "exponent"
        assert neg not in reachable(usable, m)
        k = 1
        while neg not in reachable(usable, m, k): k += 1
        out["k_min"] = k
        out["max_exponent"] = max((e for _, e in usable), default=0)
    return out

def main():
    rows = json.loads(gzip.open(ROWS).read())
    fails = [r for r in rows if not r["witness"]]
    # consistency: a row with a witness must not be classified as failing
    res = {"rows": len(rows), "failing_rows": len(fails)}
    cls = []
    for r in fails:
        c = classify(r); c.update({"p": r["p"], "a": r["a"], "u": r["u"], "factor": r["factor"]}); cls.append(c)
    # 1. p = 3678481
    res["p3678481"] = [{k: c[k] for k in ("a", "u", "type", "omega_usable", "d", "k_min", "factor") if k in c}
                       for c in cls if c["p"] == 3678481]
    # 2. type share by omega_usable
    tab = defaultdict(Counter)
    for c in cls: tab[min(c["omega_usable"], 6)][c["type"]] += 1
    res["type_by_omega_usable"] = {k: dict(v) for k, v in sorted(tab.items())}
    res["type_total"] = dict(Counter(c["type"] for c in cls))
    # 3. the four primes whose usable>=4 pairs all fail
    res["four_primes_usable>=4"] = {p: [{k: c[k] for k in ("a", "u", "type", "omega_usable", "d", "k_min") if k in c}
                                        for c in cls if c["p"] == p and c["omega_usable"] >= 4] for p in FOUR}
    # 4. among usable>=4 failures: d distribution and k_min distribution
    hi = [c for c in cls if c["omega_usable"] >= 4]
    res["usable>=4_failures"] = len(hi)
    res["usable>=4_QR_d"] = dict(Counter(tuple(c["d"]) for c in hi if c["type"] == "QR").most_common(15))
    res["usable>=4_QR_d"] = {str(k): v for k, v in res["usable>=4_QR_d"].items()}
    res["usable>=4_exponent_kmin_minus_maxexp"] = dict(sorted(Counter(c["k_min"] - c["max_exponent"] for c in hi if c["type"] == "exponent").items()))
    res["all_QR_rows_explained_by_some_d"] = all(c["d"] for c in cls if c["type"] == "QR")
    json.dump(res, sys.stdout, indent=1, default=str)

if __name__ == "__main__":
    main()
