"""P1 failure atlas (contract: briefs/p1-failure-atlas-contract.md). Exact integer arithmetic only."""
import json, sys
from math import gcd, isqrt
from collections import Counter, defaultdict

SQ840 = {1, 121, 169, 289, 361, 529}
J = 6
PAIRS = [(a, u) for a in range(1, J + 1) for u in range(1, J + 1) if gcd(a, u) == 1]

def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0

def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def sieve(n):
    b = bytearray([1]) * (n + 1); b[0] = b[1] = 0
    for i in range(2, isqrt(n) + 1):
        if b[i]: b[i*i::i] = bytearray(len(b[i*i::i]))
    return [i for i in range(n + 1) if b[i]]

_group = {}
def group(m):
    if m not in _group:
        G = [x for x in range(1, m) if gcd(x, m) == 1]
        G2 = {x * x % m for x in G}
        _group[m] = (G, G2)
    return _group[m]

def closure(gens, m):
    S = {1}; frontier = [1]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = x * g % m
                if y not in S: S.add(y); nxt.append(y)
        frontier = nxt
    return S

def analyze(n, a, u):
    M = a * n + u; m0 = 4 * a * u
    f = factor(M)
    odd = {q: e for q, e in f.items() if q % 2 == 1}
    T = sorted({q % m0 for q in odd if gcd(q, m0) == 1})
    # method 1: divisor residues from factorization (exponents included)
    D = {1}
    for q, e in odd.items():
        r = q % m0; new = set(D); pw = 1
        for _ in range(e):
            pw = pw * r % m0
            new |= {d * pw % m0 for d in D}
        D = new
    ok1 = (m0 - 1) in D
    # method 2: direct integer traversal of all divisors
    ok2 = False
    for d in range(1, isqrt(M) + 1):
        if M % d == 0 and ((d + 1) % m0 == 0 or (M // d + 1) % m0 == 0):
            ok2 = True; break
    rec = {"n": n, "a": a, "u": u, "M": M, "ok": ok1, "agree": ok1 == ok2}
    if not ok1:
        G, G2 = group(m0)
        H = closure(T, m0)
        HG2 = closure(list(T) + list(G2), m0)
        if (m0 - 1) not in HG2: typ = "QR"
        elif (m0 - 1) not in H: typ = "subgroup"
        else: typ = "exponent"
        rec["type"] = typ
        rec["T"] = T
        rec["psi_kernel"] = all(jacobi(-a * u, q) == 1 for q in odd if gcd(q, m0) == 1)
    return rec

def main():
    primes = [p for p in sieve(10**6) if p % 840 in SQ840]
    squares = [m * m for m in range(3, 1000, 2)]
    out = {"primes": [], "squares": []}
    for kind, ns in (("primes", primes), ("squares", squares)):
        for n in ns:
            for a, u in PAIRS:
                out[kind].append(analyze(n, a, u))
    # summaries
    S = {}
    for kind in out:
        rows = out[kind]
        fails = [r for r in rows if not r["ok"]]
        S[kind] = {
            "inputs": len({r["n"] for r in rows}), "rows": len(rows),
            "disagree": sum(not r["agree"] for r in rows),
            "fail_rows": len(fails),
            "types": dict(Counter(r["type"] for r in fails)),
            "QR_with_psi": sum(r["type"] == "QR" and r["psi_kernel"] for r in fails),
            "psi_kernel_any_type": sum(r["psi_kernel"] for r in fails),
        }
    # H1
    sq = out["squares"]
    S["H1_all_fail"] = all(not r["ok"] for r in sq)
    S["H1_all_QR_psi"] = all((not r["ok"]) and r["type"] == "QR" and r["psi_kernel"] for r in sq)
    # H2
    pr = out["primes"]
    h2 = [r for r in pr if not r["ok"] and r["psi_kernel"] and r["n"] % 8 == 1]
    S["H2_rows"] = len(h2)
    S["H2_violations"] = [ (r["n"], r["a"], r["u"]) for r in h2 if jacobi(r["u"], r["n"]) != 1 ][:20]
    # H4: QR failures not explained by psi
    qr_not_psi = [r for r in pr if not r["ok"] and r["type"] == "QR" and not r["psi_kernel"]]
    S["H4_QR_not_psi"] = len(qr_not_psi)
    S["H4_examples"] = [(r["n"], r["a"], r["u"], r["T"]) for r in qr_not_psi[:10]]
    # per-pair failure counts and primes failing all pairs
    by_pair = Counter((r["a"], r["u"]) for r in pr if not r["ok"])
    S["prime_fail_by_pair"] = {f"{a},{u}": by_pair[(a, u)] for a, u in PAIRS}
    allfail = defaultdict(int)
    for r in pr:
        if not r["ok"]: allfail[r["n"]] += 1
    S["primes_failing_all_pairs"] = sorted(p for p, c in allfail.items() if c == len(PAIRS))
    # u-residuosity vs failure (exploratory, for H2 context)
    tab = Counter()
    for r in pr:
        if r["n"] % 8 == 1:
            tab[(jacobi(r["u"], r["n"]), r["ok"])] += 1
    S["u_residue_vs_ok"] = {f"J(u|p)={k[0]},ok={k[1]}": v for k, v in sorted(tab.items())}
    json.dump(S, sys.stdout, indent=1, ensure_ascii=False)
    with open("failure_atlas_rows.json", "w") as fh:
        json.dump(out, fh)

if __name__ == "__main__":
    main()
