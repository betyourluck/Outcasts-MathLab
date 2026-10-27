#!/usr/bin/env python3
"""Independent check of >>96 (PR #4 16ff5ac, Research/GraphDepth20261027) and of the
>>92 checksums (PR #4 c6a0617, Research/Reply20261025).

Contract: contract.md in this directory, committed before this script.
Standard library only. The author's two scripts are neither read nor imported; the
author's result JSON files are read from git blobs only after our own computation
is finished, and only for the comparison of item G11.
"""
import hashlib
import json
import subprocess
import sys
import time
from collections import deque
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

P = 2521
SIX = {1, 121, 169, 289, 361, 529}
GD_COMMIT = "16ff5ac97aefcab04b97483e7a1ffa9d3123a083"
SR_COMMIT = "c6a0617ebad31ee9543b9cf2caff34e3af39435a"
GD_DIR = "Research/GraphDepth20261027"
SR_DIR = "Research/Reply20261025"

items = []


def record(key, claim, ours, ok, note=""):
    items.append({"item": key, "claim": claim, "ours": ours,
                  "result": "OK" if ok else "DIFF", "note": note})
    print(f"[{'OK' if ok else 'DIFF'}] {key}: {note or ours}")


# ---------------------------------------------------------------- arithmetic

def sieve(n):
    """Primes < n by the sieve of Eratosthenes."""
    flags = bytearray([1]) * n
    flags[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if flags[i]:
            flags[i * i::i] = bytearray(len(flags[i * i::i]))
    return [i for i in range(n) if flags[i]]


def is_prime_trial(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def factor(n):
    """Trial-division factorization {prime: exponent}."""
    out, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def chi(a, p):
    v = pow(a, (p - 1) // 2, p)
    assert v in (1, p - 1), (a, p, v)
    return 1 if v == 1 else -1


def divisors_from_factors(fac, mult):
    """All divisors of prod q^(mult*e)."""
    divs = [1]
    for q, e in fac.items():
        divs = [d * q ** f for d in divs for f in range(mult * e + 1)]
    return sorted(divs)


def divisors_by_scan(x):
    """Divisors of x^2 without factoring: scan d <= x and take the cofactor."""
    n = x * x
    out = set()
    for d in range(1, x + 1):
        if n % d == 0:
            out.add(d)
            out.add(n // d)
    return sorted(out)


def sigma(q):
    return 1 if q % 4 == 3 else 3


def taus(p, q):
    A = sigma(q) * q
    t, out = 1, []
    while A * t * t < 3 * p:
        out.append(t)
        t += 2
    return out


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def eval_source(p, q, tau):
    """One square source: three independent success tests must agree."""
    c = sigma(q) * q * tau * tau
    assert (p + c) % 4 == 0, (p, q, tau)
    x = (p + c) // 4
    assert 0 < x < p, (p, q, tau, x)
    assert gcd(c, p) == 1 and gcd(c, x) == 1, (p, q, tau)
    fac = factor(x)
    divs = divisors_from_factors(fac, 2)                      # (i)
    assert divs == divisors_by_scan(x), (p, q, tau)          # (ii) divisor set
    hits_II = [d for d in divs if (d + x) % c == 0]
    hits_I = [d for d in divs if (d + p * x) % c == 0]
    scan = divisors_by_scan(x)
    scan_II = [d for d in scan if (d + x) % c == 0]
    scan_I = [d for d in scan if (d + p * x) % c == 0]
    assert (scan_II, scan_I) == (hits_II, hits_I), (p, q, tau)  # (ii) hits
    N = p * x                                                 # (iii) original
    fac_N = factor(N)
    orig = any((D + N) % c == 0 for D in divisors_from_factors(fac_N, 2))
    success = bool(hits_I or hits_II)
    assert orig == success, (p, q, tau, orig, success)
    children = sorted(r for r in fac if chi(r, p) == -1)
    assert all(r < p and r != q for r in children), (p, q, tau)
    return {"q": q, "tau": tau, "c": c, "x": x,
            "factors": {str(r): e for r, e in sorted(fac.items())},
            "n_divisors_x_squared": len(divs),
            "divisors_x_squared": divs,
            "hits": {"I": hits_I, "II": hits_II},
            "success": success, "children": children}


def build_graph(p, vertices=None):
    primes_sieve = sieve(p)
    assert primes_sieve == [n for n in range(p) if is_prime_trial(n)]
    V = [q for q in primes_sieve if chi(q, p) == -1]
    if vertices is not None:
        V = [q for q in V if q in vertices]
    sources = {q: [eval_source(p, q, t) for t in taus(p, q)] for q in V}
    edges = {q: sorted({r for s in sources[q] for r in s["children"]}) for q in V}
    succ = sorted(q for q in V if any(s["success"] for s in sources[q]))
    return V, sources, edges, succ


def distances_reverse(V, edges, succ):
    rev = {q: [] for q in V}
    for q in V:
        for r in edges[q]:
            rev[r].append(q)
    dist = {q: 0 for q in succ}
    dq = deque(succ)
    while dq:
        r = dq.popleft()
        for q in rev[r]:
            if q not in dist:
                dist[q] = dist[r] + 1
                dq.append(q)
    return dist


def distances_forward(V, edges, succ):
    S = set(succ)
    dist = {}
    for v in V:
        seen, frontier, level, found = {v}, [v], 0, None
        while frontier:
            if any(u in S for u in frontier):
                found = level
                break
            nxt = []
            for u in frontier:
                for w in edges[u]:
                    if w not in seen:
                        seen.add(w)
                        nxt.append(w)
            frontier, level = nxt, level + 1
        if found is not None:
            dist[v] = found
    return dist


# ---------------------------------------------------------------- git blobs

def blob(commit, path):
    return subprocess.run(["git", "-C", str(REPO), "cat-file", "blob", f"{commit}:{path}"],
                          check=True, capture_output=True).stdout


def check_sums(key, commit, d):
    sums = blob(commit, f"{d}/SHA256SUMS").decode("utf-8").splitlines()
    rows = []
    for line in sums:
        if not line.strip():
            continue
        h, name = line.split(None, 1)
        name = name.lstrip("*").strip()
        ours = hashlib.sha256(blob(commit, f"{d}/{name}")).hexdigest()
        rows.append({"file": name, "listed": h, "blob": ours, "match": ours == h})
    n_ok = sum(r["match"] for r in rows)
    record(key, f"{len(rows)}/{len(rows)}", rows, n_ok == len(rows),
           f"SHA256SUMS {n_ok}/{len(rows)} match the blobs at {commit[:7]}")
    return rows


# ---------------------------------------------------------------- main

def main():
    t0 = time.time()
    p = P
    # G1
    ell = next(q for q in sieve(p) if chi(q, p) == -1)
    g1 = {"prime": is_prime_trial(p), "mod840": p % 840, "ell": ell}
    record("G1", "prime, class 1 mod 840, ell = 11", g1,
           g1["prime"] and g1["mod840"] in SIX and g1["mod840"] == 1 and ell == 11)

    V, sources, edges, succ = build_graph(p)
    n_src = sum(len(v) for v in sources.values())
    n_div = sum(s["n_divisors_x_squared"] for v in sources.values() for s in v)
    record("G2", "|V| = 187, sources 310, divisors 7926",
           {"V": len(V), "sources": n_src, "divisors": n_div},
           (len(V), n_src, n_div) == (187, 310, 7926))
    record("G3", "S = {23, 29, 37, 1583}", succ, succ == [23, 29, 37, 1583])

    dist_a = distances_reverse(V, edges, succ)
    dist_b = distances_forward(V, edges, succ)
    assert dist_a == dist_b
    unreachable = sorted(q for q in V if q not in dist_a)
    hist = {}
    for q in V:
        if q in dist_a:
            hist[dist_a[q]] = hist.get(dist_a[q], 0) + 1
    hist = dict(sorted(hist.items()))
    six = sorted(q for q in V if dist_a.get(q) == 6)
    record("G4", "histogram 0..6 = 4,33,73,51,18,6,2; D=6 {1249,2153}; none unreachable",
           {"histogram": hist, "distance_six": six, "unreachable": unreachable},
           list(hist.values()) == [4, 33, 73, 51, 18, 6, 2] and list(hist) == list(range(7))
           and six == [1249, 2153] and not unreachable,
           f"histogram {hist}, D=6 {six}, unreachable {unreachable} (reverse BFS = forward BFS)")

    # G5
    claim17 = [(1, 51, 643, {643: 1}, [643]), (3, 459, 745, {5: 1, 149: 1}, [149]),
               (5, 1275, 949, {13: 1, 73: 1}, [73]), (7, 2499, 1255, {5: 1, 251: 1}, [251]),
               (9, 4131, 1663, {1663: 1}, [1663]), (11, 6171, 2173, {41: 1, 53: 1}, [53])]
    ours17 = [(s["tau"], s["c"], s["x"], {int(k): v for k, v in s["factors"].items()},
               s["children"]) for s in sources[17]]
    fail17 = all(not s["success"] for s in sources[17])
    record("G5", "q=17 table, all six sources fail", ours17, ours17 == claim17 and fail17,
           f"{len(ours17)} sources, table equal = {ours17 == claim17}, all fail = {fail17}")

    # G6
    ch17 = edges[17]
    counts = [len(sources[r]) for r in ch17]
    allfail = all(not s["success"] for r in ch17 for s in sources[r])
    total = len(sources[17]) + sum(counts)
    record("G6", "children {53,73,149,251,643,1663}, counts 3,3,2,3,2,1, all fail, 20 sources",
           {"children": ch17, "counts": counts, "all_fail": allfail, "total": total},
           ch17 == [53, 73, 149, 251, 643, 1663] and counts == [3, 3, 2, 3, 2, 1]
           and allfail and total == 20)

    # G7
    s1 = next(s for s in sources[17] if s["tau"] == 11)
    s2 = next(s for s in sources[53] if s["tau"] == 5)
    s3 = next(s for s in sources[29] if s["tau"] == 1)
    d = 1304
    c, x = s3["c"], s3["x"]
    N = p * x
    y = (N + d) // c
    z = (N + N * N // d) // c
    exact = Fraction(4, p) == Fraction(1, x) + Fraction(1, y) + Fraction(1, z)
    ident = 4 * x * y * z == p * (x * y + x * z + y * z)
    g7 = {"x1": s1["x"], "f1": s1["factors"], "x2": s2["x"], "f2": s2["factors"],
          "c": c, "x": x, "d_in_hits_I": d in s3["hits"]["I"], "d_divides_x2": (x * x) % d == 0,
          "denominators": [x, y, z], "fraction": exact, "identity": ident,
          "D17": dist_a.get(17), "path_edges": [53 in edges[17], 29 in edges[53]]}
    record("G7", "17 -11-> 53 -5-> 29 -1-> Type I; d=1304; (652,18908,23833534); D(17)=2", g7,
           s1["x"] == 2173 and s1["factors"] == {"41": 1, "53": 1} and 53 in s1["children"]
           and s2["x"] == 1624 and s2["factors"] == {"2": 3, "7": 1, "29": 1} and 29 in s2["children"]
           and (c, x) == (87, 652) and g7["d_in_hits_I"] and g7["d_divides_x2"]
           and [x, y, z] == [652, 18908, 23833534] and exact and ident and dist_a.get(17) == 2)

    # G8
    record("G8", "D(11) = 1", dist_a.get(11), dist_a.get(11) == 1)

    # G9 and the U1 set
    U1 = sorted(q for q in V if q not in succ and not any(r in succ for r in edges[q]))
    via_dist = sorted(q for q in V if dist_a.get(q, 10 ** 9) >= 2)
    assert U1 == via_dist
    record("G9", "smallest U1-refuting q = 17", {"min": U1[0], "count": len(U1)},
           U1[0] == 17, f"smallest {U1[0]}, U1 counterexample vertices {len(U1)} "
           f"(= vertices with D >= 2), completely failed vertices {len(V) - len(succ)}")

    record("G10", "three success tests agree on all 310 sources", n_src, True,
           "factorization / divisor scan / (px)^2 criterion agree on all "
           f"{n_src} sources (asserted per source)")

    # G12: our own earlier results
    root_fail = all(not s["success"] for s in sources[11])
    bfs = json.loads((REPO / "Research/Bfs12Inputs20261024/bfs12_results.json")
                     .read_text(encoding="utf-8"))
    row = next(r for r in bfs["inputs"] if r["p"] == p)
    prev = row["sweep"]["child_status"]
    ours_status = {str(r): ("success" if r in succ else "fail") for r in edges[11]}
    prev_status = {k: v["status"] for k, v in prev.items()}
    prev_ntau = {k: v["n_tau"] for k, v in prev.items() if "n_tau" in v}
    ours_ntau = {str(r): len(sources[r]) for r in edges[11] if r not in succ}
    path_prev = row["bfs"]["path"]
    ok12 = (root_fail and len(sources[11]) == 13
            and edges[11] == [17, 19, 29, 73, 107, 131, 139, 211, 233, 541, 853, 1249]
            and ours_status == prev_status and ours_ntau == prev_ntau
            and path_prev == [{"parent": 11, "tau": 25, "child": 29, "v_child": 1}]
            and 29 in next(s for s in sources[11] if s["tau"] == 25)["children"])
    record("G12", "11: 13 sources fail, 12 children, matches Bfs12Inputs20261024",
           {"children": edges[11], "status": ours_status}, ok12)

    # C1: cited known pair at p = 12889 (only the two vertices)
    pc = 12889
    Vc, src_c, edges_c, succ_c = build_graph(pc, vertices={43, 61})
    c1 = {"prime": is_prime_trial(pc), "mod840": pc % 840, "in_V": Vc}
    if Vc == [43, 61]:
        c1.update({"sources": sum(len(src_c[q]) for q in Vc),
                   "failed": sum(1 for q in Vc for s in src_c[q] if not s["success"]),
                   "taus": {str(q): [s["tau"] for s in src_c[q]] for q in Vc},
                   "43->61": 61 in edges_c[43], "61->43": 43 in edges_c[61],
                   "edges43": edges_c[43], "edges61": edges_c[61]})
    record("C1", "12889: 43 <-> 61, all 22 sources fail (cited, not the author's claim)", c1,
           c1["prime"] and c1["mod840"] in SIX and Vc == [43, 61] and c1["sources"] == 22
           and c1["failed"] == 22 and c1["43->61"] and c1["61->43"])

    t_compute = time.time() - t0
    ours_out = {"p": p, "vertices": V, "edges": {str(q): edges[q] for q in V},
                "successful_vertices": succ, "distances": {str(q): dist_a[q] for q in V if q in dist_a},
                "histogram": {str(k): v for k, v in hist.items()}, "U1": U1,
                "sources": {str(q): [{k: v for k, v in s.items() if k != "divisors_x_squared"}
                                     for s in sources[q]] for q in V}}

    # ------------------------------------------------ G11: compare with the author's JSON
    theirs = json.loads(blob(GD_COMMIT, f"{GD_DIR}/graph-depth-results.json"))
    run = theirs["runs"][0]
    diffs = []
    if run["p"] != p or run["vertex_count"] != len(V) or run["source_count"] != n_src:
        diffs.append("header")
    if sorted(int(k) for k in run["sources"]) != V:
        diffs.append("vertex set")
    n_rows = 0
    for qs, rows in run["sources"].items():
        q = int(qs)
        mine = {s["tau"]: s for s in sources[q]}
        if sorted(r["tau"] for r in rows) != sorted(mine):
            diffs.append(f"taus of {q}")
            continue
        for r in rows:
            m = mine[r["tau"]]
            n_rows += 1
            for key in ("q", "c", "x", "children"):
                if r[key] != m[key]:
                    diffs.append(f"{q} tau={r['tau']} {key}")
            if {str(k): v for k, v in r["factors"].items()} != m["factors"]:
                diffs.append(f"{q} tau={r['tau']} factors")
            if sorted(r["divisors_x_squared"]) != m["divisors_x_squared"]:
                diffs.append(f"{q} tau={r['tau']} divisors")
            if sorted(r["hits"]["I"]) != m["hits"]["I"] or sorted(r["hits"]["II"]) != m["hits"]["II"]:
                diffs.append(f"{q} tau={r['tau']} hits")
    if {int(k): sorted(v) for k, v in run["edges"].items()} != edges:
        diffs.append("edges")
    if sorted(run["successful_vertices"]) != succ:
        diffs.append("successful_vertices")
    if {int(k): v for k, v in run["distances"].items()} != dist_a:
        diffs.append("distances")
    if {int(k): v for k, v in run["distance_histogram"].items()} != hist:
        diffs.append("distance_histogram")
    if sorted(run["unreachable"]) != unreachable:
        diffs.append("unreachable")
    if sorted(run["U1_counterexample_vertices"]) != U1:
        diffs.append("U1_counterexample_vertices")
    sel = run["selected"]
    replay = json.loads(blob(GD_COMMIT, f"{GD_DIR}/independent-replay.json"))
    rc = replay["checks"][0]
    fs_ours = {str(q): len(sources[q]) for q in [17] + ch17}
    if rc["all_graph_sources"] != n_src or rc["divisors_tested"] != n_div:
        diffs.append("replay counts")
    if rc["failure_sources"] != fs_ours:
        diffs.append("replay failure_sources")
    if sorted(rc["distance_six_vertices"]) != six:
        diffs.append("replay distance_six_vertices")
    q17_ours = [{"tau": s["tau"], "c": s["c"], "x": s["x"], "factors": s["factors"],
                 "children": s["children"]} for s in sources[17]]
    if rc["q17_sources"] != q17_ours:
        diffs.append("replay q17_sources")
    record("G11", "author's JSON (310 rows, edges, distances, U1 list, replay)",
           {"rows_compared": n_rows, "diffs": diffs, "their_U1_count": len(run["U1_counterexample_vertices"]),
            "selected": sel}, not diffs,
           f"{n_rows} source rows x 7 fields + graph fields + replay: {len(diffs)} differences {diffs}")

    # ------------------------------------------------ G15, G16: checksums
    check_sums("G15", GD_COMMIT, GD_DIR)
    check_sums("G16", SR_COMMIT, SR_DIR)
    tree_a = subprocess.run(["git", "-C", str(REPO), "rev-parse", f"{SR_COMMIT}:{SR_DIR}"],
                            check=True, capture_output=True, text=True).stdout.strip()
    tree_b = subprocess.run(["git", "-C", str(REPO), "rev-parse", f"{GD_COMMIT}:{SR_DIR}"],
                            check=True, capture_output=True, text=True).stdout.strip()
    record("G16b", "Reply20261025 unchanged between c6a0617 and 16ff5ac", [tree_a, tree_b],
           tree_a == tree_b, f"tree {tree_a[:12]} at both commits = {tree_a == tree_b}")

    n_ok = sum(i["result"] == "OK" for i in items)
    n_diff = len(items) - n_ok
    out = {"contract": "contract.md", "script": "check96.py",
           "head": subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], check=True,
                                  capture_output=True, text=True).stdout.strip(),
           "python": sys.version.split()[0], "seconds_compute": round(t_compute, 3),
           "seconds_total": round(time.time() - t0, 3),
           "summary": {"OK": n_ok, "DIFF": n_diff}, "items": items, "graph": ours_out}
    (HERE / "check96_results.json").write_text(json.dumps(out, ensure_ascii=False, indent=1),
                                               encoding="utf-8", newline="\n")
    print(f"OK {n_ok} / DIFF {n_diff}; compute {t_compute:.2f} s, total {time.time() - t0:.2f} s")
    return 0 if n_diff == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
