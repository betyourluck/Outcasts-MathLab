"""Compare check_results.json with the author's files in PR #4 (fe2150e), read through git.
Usage: python compare.py <git ref of the PR #4 branch, e.g. origin/codex/witness-budget-falsification>
Needs a clone that has the PR #4 branch and the PR #7 commit 4b5a002 (fork betyourluck/Outcasts-MathLab)."""
import hashlib, json, subprocess, sys
from math import prod
from pathlib import Path

HERE = Path(__file__).resolve().parent

def show(ref, path):
    return subprocess.run(["git", "show", f"{ref}:{path}"], capture_output=True, check=True).stdout

def main():
    ref = sys.argv[1]; base = "Research/Review20261009"
    mine = json.loads((HERE / "check_results.json").read_text(encoding="utf-8"))
    out = {}
    # 61 rows
    A = json.loads(show(ref, f"{base}/p3678481-all-coprime-rows.json"))
    M = {(r["a"], r["u"]): r for r in mine["item2"]["rows"]}
    bad = [(r["a"], r["u"]) for r in A
           if not (r["M"] == M[(r["a"], r["u"])]["M"] and r["omega"] == M[(r["a"], r["u"])]["omega"]
                   and sorted(r["witnesses"]) == M[(r["a"], r["u"])]["witnesses"]
                   and prod(int(q) ** e for q, e in r["factorization"].items()) == r["M"])]
    out["rows"] = {"author": len(A), "same_pairs": set(M) == {(r["a"], r["u"]) for r in A}, "mismatches": bad}
    # FCT
    F = json.loads(show(ref, f"{base}/fixed-review-results.json"))
    keys = ("p", "i", "d", "x", "y", "z", "type")
    out["fct_equal"] = [{k: r[k] for k in keys} for r in F["FCT_fixed_five"]] == [{k: r[k] for k in keys} for r in mine["item3"]["first_i"]]
    out["p2521_d11_equal"] = F["FCT_p2521_d11_possible_i"] == mine["item3"]["p2521_d11"]["i_candidates_dividing_x"] and F["FCT_p2521_d11_successful_i"] == []
    out["small_shift_solution_equal"] = all(F["small_shift_solution"][k] == mine["item2"]["identity_1_1_23"][k] for k in "xyz")
    # source manifest
    S = json.loads(show(ref, f"{base}/source-manifest.json"))
    res = []
    for f in S["files"]:
        data = show(f["commit"], f["path"])
        blob = subprocess.run(["git", "rev-parse", f"{f['commit']}:{f['path']}"], capture_output=True, text=True).stdout.strip()
        res.append(hashlib.sha256(data).hexdigest() == f["sha256"] and blob == f["git_blob_sha1"] and len(data) == f["bytes"])
    out["manifest"] = {"files": len(res), "all_match": all(res)}
    # SHA256SUMS of the two author directories
    ok = total = 0
    for d in ("Research/Review20261009", "Research/KneserBudget20261010"):
        for line in show(ref, f"{d}/SHA256SUMS").decode().replace("\r", "").splitlines():
            h, f = line.split(maxsplit=1); f = f.lstrip("*")
            total += 1; ok += hashlib.sha256(show(ref, f"{d}/{f}")).hexdigest() == h
    out["sha256sums"] = {"entries": total, "match": ok}
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
