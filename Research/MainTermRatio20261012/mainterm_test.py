"""K test (contract: Research/MainTermRatio20261012/contract.md, pushed before computing).
1.6e7 < p < 1.9e7 (hard classes), coprime a*u <= bit_length(p-1): witness count W(p) versus the trivial-character
term M_t and the two-term main term M_t2 (trivial character + the odd real character chi_{-au}, whose values on the
prime factors of ap+u are the Legendre symbols (p|q))."""
import gzip, json, random, sys, time
from collections import defaultdict
from math import gcd
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
# Load ./explore.py under a distinct module name: it does `from explore import n2` itself (HardestPrimes20261010/explore.py),
# so importing it as `explore` here made that line hit the half-initialised module (failed-run-1.log).
import importlib.util  # noqa: E402
_spec = importlib.util.spec_from_file_location("mainterm_explore", HERE / "explore.py")
ex = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(ex)  # witness_count, phi, quant; ex.st = simultaneity_test, ex.n2
st, n2 = ex.st, ex.n2

SQ840 = {1, 121, 169, 289, 361, 529}
LO, HI = 16 * 10**6, 19 * 10**6
K1, K2, K3 = (0.90, 0.98), (-0.05, 0.20), (0.95, 1.05)
CHECK_RATE, CHECK_SEED = 0.01, 20261023

def avg_ranks(x):
    x = np.asarray(x, dtype=float); order = np.argsort(x, kind="mergesort"); r = np.empty(len(x)); i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]: j += 1
        r[order[i:j + 1]] = (i + j) / 2 + 1; i = j + 1
    return r

def spearman(a, b):
    return float(np.corrcoef(avg_ranks(a), avg_ranks(b))[0, 1])

def main():
    t0 = time.time()
    out = {"contract": "Research/MainTermRatio20261012/contract.md (pushed before computing)", "range": [LO, HI]}
    small = st.sieve(22000)
    primes = [p for p in st.sieve(HI) if LO < p < HI and p % 840 in SQ840]
    P = np.array(primes, dtype=np.int64); budget = np.array([(p - 1).bit_length() for p in primes]); B = int(budget.max())
    assert small[-1] ** 2 > B * HI + B
    N = len(primes); out["primes"] = N; out["budgets"] = {str(b): int((budget == b).sum()) for b in sorted(set(budget.tolist()))}
    print(f"[{round(time.time() - t0)}s] {N} primes, budgets {out['budgets']}", file=sys.stderr, flush=True)
    rng = random.Random(CHECK_SEED); checked = 0; rows_done = 0
    W = np.zeros(N); Mt = np.zeros(N); Mt2 = np.zeros(N); Mt3 = np.zeros(N); S = np.zeros(N, dtype=int); npairs = np.zeros(N, dtype=int)
    pooled = defaultdict(lambda: [0.0, 0.0, 0])      # phi -> [sum W, sum Mt2, rows]
    acc1 = {}; acc2 = {}                              # (a,u) -> [n, sum e, sum e^2] for e = W - Mt and W - Mt2
    pair_stats = {}
    for a in range(1, B + 1):
        for u in range(1, B // a + 1):
            if gcd(a, u) != 1: continue
            m = 4 * a * u; ph = ex.phi(m)
            sel = np.nonzero(budget >= a * u)[0]
            Ms, facs = st.factor_all(P[sel], a, u, small)
            w = np.zeros(len(sel)); mt = np.zeros(len(sel)); mt2 = np.zeros(len(sel)); mt3 = np.zeros(len(sel))
            for j, (idx, Mi, f) in enumerate(zip(sel.tolist(), Ms.tolist(), facs)):
                prod = 1
                for q, e in f.items(): prod *= q ** e
                if prod != Mi: raise SystemExit(f"STOP: factorisation mismatch at {(primes[idx], a, u)}")
                cnt, tau, chi_sum, fair = ex.witness_count(f, m, primes[idx])
                w[j] = cnt; mt[j] = tau / ph; mt2[j] = (tau - chi_sum) / ph; mt3[j] = (tau - fair) / ph
                if m == 4 and cnt != mt2[j]: raise SystemExit(f"STOP: (1,1) identity W = (tau - chi_sum)/2 fails at p={primes[idx]}")
                if rng.random() < CHECK_RATE:
                    if st.direct(Mi, m) != (cnt > 0): raise SystemExit(f"STOP: direct scan disagrees at {(primes[idx], a, u)}")
                    checked += 1
            W[sel] += w; Mt[sel] += mt; Mt2[sel] += mt2; Mt3[sel] += mt3; S[sel] += (w > 0); npairs[sel] += 1
            g = pooled[ph]; g[0] += float(w.sum()); g[1] += float(mt2.sum()); g[2] += len(sel)
            e1 = w - mt; e2 = w - mt2
            acc1[(a, u)] = [len(sel), float(e1.sum()), float((e1 ** 2).sum())]; acc2[(a, u)] = [len(sel), float(e2.sum()), float((e2 ** 2).sum())]
            pair_stats[f"{a},{u}"] = {"rows": int(len(sel)), "phi_4au": ph, "mean_W": round(float(w.mean()), 4), "mean_Mt": round(float(mt.mean()), 4),
                                      "mean_Mt2": round(float(mt2.mean()), 4), "mean_Mt3": round(float(mt3.mean()), 4),
                                      "W_over_Mt2": round(float(w.sum() / mt2.sum()), 4) if mt2.sum() else None, "success_rate": round(float((w > 0).mean()), 4)}
            rows_done += len(sel)
        print(f"[{round(time.time() - t0)}s] a={a} done, rows {rows_done}, direct checks {checked}", file=sys.stderr, flush=True)
    out["rows"] = rows_done; out["direct_checks"] = checked; out["W_zero_primes"] = int((W == 0).sum())
    bad = Mt2 <= 0; out["Mt2_nonpositive_primes"] = int(bad.sum()); ok = ~bad
    n2p = np.array([n2(p) for p in primes]); R1 = W[ok] / Mt[ok]; R2 = W[ok] / Mt2[ok]; R3 = W[ok] / Mt3[ok]
    # judgements
    k1 = float(R2.mean()); rho2 = spearman(R2, n2p[ok]); rho1 = spearman(R1, n2p[ok]); rho_m = spearman(Mt2[ok], n2p[ok])
    small_phi = [ph for ph in pooled if ph <= 8]; k3 = sum(pooled[ph][0] for ph in small_phi) / sum(pooled[ph][1] for ph in small_phi)
    out["K1_mean_R2"] = {"value": round(k1, 4), "window": K1, "pass": K1[0] <= k1 <= K1[1]}
    out["K2_spearman_R2_n2"] = {"value": round(rho2, 4), "window": K2, "pass": K2[0] < rho2 < K2[1],
                                "rho_R1_n2": round(rho1, 4), "rho_Mt2_n2": round(rho_m, 4), "rho_W_n2": round(spearman(W[ok], n2p[ok]), 4),
                                "rho_Mt_n2": round(spearman(Mt[ok], n2p[ok]), 4)}
    out["K3_pooled_small_phi"] = {"value": round(k3, 4), "window": K3, "pass": K3[0] <= k3 <= K3[1], "phi_values": sorted(small_phi),
                                  "rows": int(sum(pooled[ph][2] for ph in small_phi))}
    out["all_pass"] = out["K1_mean_R2"]["pass"] and out["K2_spearman_R2_n2"]["pass"] and out["K3_pooled_small_phi"]["pass"]
    # measurements (K4)
    out["K4"] = {"mean_W": round(float(W.mean()), 3), "mean_Mt": round(float(Mt.mean()), 3), "mean_Mt2": round(float(Mt2.mean()), 3), "mean_Mt3": round(float(Mt3.mean()), 3),
                 "mean_S": round(float(S.mean()), 3), "mean_pairs": round(float(npairs.mean()), 2),
                 "mean_R1": round(float(R1.mean()), 4), "R1_quantiles": ex.quant(R1), "min_R1": round(float(R1.min()), 4),
                 "R2_quantiles": ex.quant(R2), "min_R2": round(float(R2.min()), 4), "mean_R3": round(float(R3.mean()), 4),
                 "mean_W_minus_Mt2": round(float((W[ok] - Mt2[ok]).mean()), 3), "sd_W_minus_Mt2": round(float((W[ok] - Mt2[ok]).std()), 3),
                 "mean_W_minus_Mt": round(float((W - Mt).mean()), 3), "sd_W_minus_Mt": round(float((W - Mt).std()), 3),
                 "corr_W_Mt": round(float(np.corrcoef(W, Mt)[0, 1]), 4), "corr_W_Mt2": round(float(np.corrcoef(W, Mt2)[0, 1]), 4),
                 "pooled_by_phi": {str(ph): {"rows": v[2], "W_over_Mt2": round(v[0] / v[1], 4)} for ph, v in sorted(pooled.items())}}
    order = np.argsort(W / np.where(Mt2 > 0, Mt2, np.inf))[:10]
    out["K4"]["lowest_R2_primes"] = [{"p": int(P[i]), "B": int(budget[i]), "W": int(W[i]), "Mt": round(float(Mt[i]), 3), "Mt2": round(float(Mt2[i]), 3),
                                      "R2": round(float(W[i] / Mt2[i]), 4), "S": int(S[i]), "pairs": int(npairs[i]), "n2": int(n2p[i])} for i in order]
    by_n2 = defaultdict(list)
    for i in range(N): by_n2[int(n2p[i])].append(i)
    out["K4"]["by_n2"] = {str(k): {"primes": len(v), "mean_W": round(float(W[v].mean()), 3), "mean_Mt2": round(float(Mt2[v].mean()), 3),
                                   "mean_R1": round(float((W[v] / Mt[v]).mean()), 4), "mean_R2": round(float((W[v] / Mt2[v]).mean()), 4)} for k, v in sorted(by_n2.items())}
    out["K4"]["by_budget"] = {str(b): {"primes": int((budget == b).sum()), "mean_W": round(float(W[budget == b].mean()), 3),
                                       "mean_Mt2": round(float(Mt2[budget == b].mean()), 3), "mean_R2": round(float((W / Mt2)[budget == b].mean()), 4)}
                              for b in sorted(set(budget.tolist()))}
    def decomposition(acc, d):
        s = 0.0
        for (a, u), (n, s1, s2) in acc.items():
            if n > 1: s += (s2 / n - (s1 / n) ** 2) * n / N
        return {"var_of_sum": round(float(d.var()), 3), "sum_of_pair_vars": round(s, 3), "ratio": round(float(d.var()) / s, 4) if s else None}
    out["K4"]["alignment_Mt"] = decomposition(acc1, W - Mt); out["K4"]["alignment_Mt2"] = decomposition(acc2, W - Mt2)
    out["pairs"] = pair_stats
    out["seconds"] = round(time.time() - t0, 1)
    with gzip.GzipFile("mainterm_test_rows.json.gz", "wb", mtime=0) as fh:
        fh.write(json.dumps([{"p": int(P[i]), "B": int(budget[i]), "pairs": int(npairs[i]), "W": int(W[i]), "Mt": round(float(Mt[i]), 6),
                              "Mt2": round(float(Mt2[i]), 6), "Mt3": round(float(Mt3[i]), 6), "S": int(S[i]), "n2": int(n2p[i])} for i in range(N)]).encode())
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
