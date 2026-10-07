"""Exponent-type classification (contract: briefs/p1-exponent-contract.md). Exact integers; MC only in H10."""
import gzip, json, random, sys
from collections import Counter, defaultdict
from math import gcd
from failure_atlas import factor, closure, group

SEED = 20261005
MC_SAMPLES = 300

def load():
    rows = json.load(gzip.open("failure_atlas_rows.json.gz", "rt"))["primes"]
    rows += json.load(gzip.open("h2prime_rows.json.gz", "rt"))
    return rows

_inv = {}
def inv_table(m):
    if m not in _inv:
        G, _ = group(m)
        _inv[m] = {x: pow(x, -1, m) for x in G}
    return _inv[m]

def bounded_products(T_mult, m, k):
    S = {1}
    for t in T_mult:
        pw = [1]
        for _ in range(k): pw.append(pw[-1] * t % m)
        S = {s * q % m for s in S for q in pw}
    return S

def enrich(r):
    a, u, M = r["a"], r["u"], r["M"]; m0 = 4 * a * u
    f = factor(M)
    odd = {q: e for q, e in f.items() if q % 2 == 1}
    Mp = 1
    for q, e in odd.items(): Mp *= q ** e
    assert all(gcd(q, m0) == 1 for q in odd), (M, m0)
    r["omega"] = len(odd); r["Omega"] = sum(odd.values())
    r["sqfree"] = all(e == 1 for e in odd.values())
    r["r"] = Mp % m0
    if r["ok"]: return
    G, G2 = group(m0); inv = inv_table(m0)
    T = sorted({q % m0 for q in odd}); r["Tsize"] = len(T); r["phi"] = len(G)
    r["orbit_ok"] = all(((-inv[x]) % m0) not in T for x in T)
    H = closure(T, m0); r["Hsize"] = len(H); r["H_full"] = len(H) == len(G)
    # exact divisor residues (odd part), exponent of each prime bounded by its multiplicity
    D = {1}
    for q, e in odd.items():
        D = {d * x % m0 for d in D for x in bounded_products([q % m0], m0, e)}
    r["Dsize"] = len(D); r["prod_v1"] = 1
    for e in odd.values(): r["prod_v1"] *= e + 1
    if r["type"] == "exponent":
        T_mult = [q % m0 for q in odd]
        k = 1
        while (m0 - 1) not in bounded_products(T_mult, m0, k):
            k += 1
            assert k <= len(G)
        r["k_min"] = k; r["max_v"] = max(odd.values())
        if len(odd) == 1:
            (q, v), = odd.items(); t = q % m0
            o, x = 1, t
            while x != 1: x = x * t % m0; o += 1
            r["ord"] = o; r["v"] = v
            r["half_is_neg1"] = (o % 2 == 0) and pow(t, o // 2, m0) == m0 - 1

def model_fail_prob(m0, omega, r, rng):
    """P(no subset product == -1) with t_1..t_{w-1} uniform in G, t_w = r / prod."""
    G, _ = group(m0); inv = inv_table(m0); neg = m0 - 1
    if r == neg: return 0.0
    if omega == 1: return 1.0
    if omega == 2:
        bad = {neg, (neg * inv[r]) % m0}  # t1 = -1 or t2 = r/t1 = -1
        return 1 - len(bad) / len(G)
    if omega == 3:
        cnt = 0
        for t1 in G:
            if t1 == neg or (r * inv[t1]) % m0 == neg: continue
            # t2 must avoid: -1, -r (=> t3=-1), -1/t1 (t1 t2 = -1), -r/t1 (t2 t3 = r/t1 => -1 when t2 = -r/t1)
            bad = {neg, (neg * r) % m0, (neg * inv[t1]) % m0, (neg * r * inv[t1]) % m0}
            cnt += len(G) - len(bad)
        return cnt / (len(G) ** 2)
    fails = 0
    for _ in range(MC_SAMPLES):
        ts = [rng.choice(G) for _ in range(omega - 1)]
        prod = 1
        for t in ts: prod = prod * t % m0
        ts.append(r * inv[prod] % m0)
        S = {1}
        hit = False
        for t in ts:
            S = S | {s * t % m0 for s in S}
            if neg in S: hit = True; break
        if not hit: fails += 1
    return fails / MC_SAMPLES

def main():
    rows = load()
    for r in rows: enrich(r)
    fails = [r for r in rows if not r["ok"]]
    expo = [r for r in fails if r["type"] == "exponent"]
    S = {"rows": len(rows), "fail_rows": len(fails), "exponent_rows": len(expo)}
    # H6
    S["H6_orbit_violations"] = sum(not r["orbit_ok"] for r in fails)
    S["H6_Tsize_over_half_phi"] = sum(2 * r["Tsize"] > r["phi"] for r in fails)
    # H7
    S["H7_kmin_lt_2"] = sum(r["k_min"] < 2 for r in expo)
    w1 = [r for r in expo if r["omega"] == 1]
    S["H7_omega1_rows"] = len(w1)
    S["H7_omega1_violations"] = sum(not (r["half_is_neg1"] and r["v"] < r["ord"] // 2) for r in w1)
    # H8
    def dist(rs, key): return dict(sorted(Counter(r[key] for r in rs).items()))
    S["H8_omega"] = {
        "all": dist(rows, "omega"), "ok": dist([r for r in rows if r["ok"]], "omega"),
        "QR": dist([r for r in fails if r["type"] == "QR"], "omega"),
        "subgroup": dist([r for r in fails if r["type"] == "subgroup"], "omega"),
        "exponent": dist(expo, "omega")}
    S["H8_sqfree_frac"] = {k: round(sum(r["sqfree"] for r in rs) / len(rs), 4) for k, rs in
                           (("all", rows), ("exponent", expo), ("QR", [r for r in fails if r["type"] == "QR"]))}
    S["H8_exponent_Omega"] = dist(expo, "Omega")
    # H9
    S["H9_kmin"] = dist(expo, "k_min")
    S["H9_kmin_by_omega"] = {w: dist([r for r in expo if r["omega"] == w], "k_min") for w in sorted({r["omega"] for r in expo})}
    S["H9_H_full_frac"] = round(sum(r["H_full"] for r in expo) / len(expo), 4)
    S["H9_D_lt_prodv1"] = sum(r["Dsize"] < r["prod_v1"] for r in expo)
    S["H9_omega1_ord"] = dist(w1, "ord")
    S["H9_omega1_gap"] = dict(sorted(Counter(r["ord"] // 2 - r["v"] for r in w1).items()))
    S["H9_Tsize_frac_of_half_phi_exponent"] = dict(sorted(Counter(round(2 * r["Tsize"] / r["phi"], 1) for r in expo).items()))
    # failure type by omega (all rows)
    tab = defaultdict(Counter)
    for r in rows: tab[r["omega"]][r["type"] if not r["ok"] else "ok"] += 1
    S["type_by_omega"] = {w: dict(tab[w]) for w in sorted(tab)}
    # H10
    rng = random.Random(SEED); cache = {}
    agg = defaultdict(lambda: [0, 0.0, 0])  # omega -> [rows, sum model, observed fails]
    agg_expo = defaultdict(lambda: [0, 0.0])
    for r in rows:
        if not r["sqfree"]: continue
        key = (4 * r["a"] * r["u"], r["omega"], r["r"])
        if key not in cache: cache[key] = model_fail_prob(*key, rng)
        pm = cache[key]
        A = agg[r["omega"]]; A[0] += 1; A[1] += pm; A[2] += (not r["ok"])
    S["H10_cells"] = len(cache)
    S["H10_by_omega"] = {w: {"rows": A[0], "model_fails": round(A[1], 1), "observed_fails": A[2],
                             "ratio_obs_over_model": round(A[2] / A[1], 3) if A[1] else None}
                         for w, A in sorted(agg.items())}
    tot = [sum(A[0] for A in agg.values()), sum(A[1] for A in agg.values()), sum(A[2] for A in agg.values())]
    S["H10_total"] = {"rows": tot[0], "model_fails": round(tot[1], 1), "observed_fails": tot[2],
                      "ratio": round(tot[2] / tot[1], 4)}
    json.dump(S, sys.stdout, indent=1)
    with open("exponent_rows.json", "w") as fh:
        json.dump([{k: v for k, v in r.items() if k not in ("agree",)} for r in fails], fh)

if __name__ == "__main__":
    main()
