"""H2' (contract: briefs/p1-h2prime-contract.md)."""
import json, sys
from collections import Counter
from failure_atlas import analyze, jacobi, sieve, SQ840
US = [11, 13, 17, 19, 23]; AS = range(1, 7)
primes = [p for p in sieve(10**6) if p % 840 in SQ840]
rows = []
for p in primes:
    for u in US:
        for a in AS:
            r = analyze(p, a, u); r["Jup"] = jacobi(u, p); rows.append(r)
S = {"inputs": len(primes), "rows": len(rows), "disagree": sum(not r["agree"] for r in rows)}
for sgn in (-1, 1):
    sub = [r for r in rows if r["Jup"] == sgn]
    fails = [r for r in sub if not r["ok"]]
    S[f"J(u|p)={sgn}"] = {
        "rows": len(sub), "fail_rows": len(fails),
        "fail_rate": round(len(fails) / len(sub), 4) if sub else None,
        "types": dict(Counter(r["type"] for r in fails)),
        "psi_kernel_fails": sum(r["psi_kernel"] for r in fails),
    }
S["H2prime_violations"] = [(r["n"], r["a"], r["u"]) for r in rows
                           if r["Jup"] == -1 and not r["ok"] and r["psi_kernel"]][:20]
S["by_u"] = {u: {sgn: (lambda sub: [len(sub), sum(not r["ok"] for r in sub)])(
    [r for r in rows if r["u"] == u and r["Jup"] == sgn]) for sgn in (-1, 1)} for u in US}
S["primes_with_all_rows_failing"] = sorted({r["n"] for r in rows} - {r["n"] for r in rows if r["ok"]})[:20]
json.dump(S, sys.stdout, indent=1)
with open("h2prime_rows.json", "w") as fh: json.dump(rows, fh)
