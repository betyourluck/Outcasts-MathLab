"""Exploration on published rows (X, Y, E, S, Q, H tests; 4e6 < p < 1.6e7, hard classes), not a judgement.
For each prime p and each coprime pair (a, u) with a*u <= ceil(log2 p), the witness count
    W(p; a, u) = #{ s | ap+u : s = -1 (mod 4au) }
is compared with the trivial-character term of its character expansion ("crude main term")
    Mt(p; a, u) = tau(odd part of ap+u) / phi(4au),
i.e. the count expected if the odd divisors of ap+u were equidistributed in (Z/4au)^x.  Question (Outcasts >>53):
how is the ratio R(p) = W(p)/Mt(p) distributed over p, in particular its lower tail, and do the per-pair errors
W - Mt line up across the pairs of one p (variance of the sum versus the sum of the variances)?"""
import gzip, json, sys, time
from collections import defaultdict
from math import gcd
from pathlib import Path
import numpy as np

R = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(R / "Simultaneity20261010")); sys.path.insert(0, str(R / "HardestPrimes20261010"))
import simultaneity_test as st  # noqa: E402  (sieve, factor_all: vectorised trial division, cofactor 1 or prime)
from explore import n2  # noqa: E402  (least prime quadratic non-residue, HardestPrimes20261010/explore.py)

FILES = {"X": "Character20261008/character_rows.json.gz", "Y": "MediumSign20261008/mediumsign_rows.json.gz",
         "E": "Exponent20261008/exponent_rows.json.gz", "S": "Simultaneity20261010/simultaneity_rows.json.gz",
         "Q": "QR11Bit20261010/qr11_rows.json.gz", "H": "HardestPrimes20261010/hardest_rows.json.gz"}

def phi(m):
    r = m; q = 2; n = m
    while q * q <= n:
        if n % q == 0:
            r -= r // q
            while n % q == 0: n //= q
        q += 1
    if n > 1: r -= r // n
    return r

def witness_count(f, m, p):
    """Number of odd divisors s of M (prime factorisation f) with s = -1 (mod m); tau of the odd part;
    the divisor sum of the odd real character chi(q) = (p|q) (= (-au|q) for q | ap+u, a character mod 4au with
    chi(-1) = -1), prod_q (1 + chi(q) + ... + chi(q)^e); and its fair-coin expectation prod_q ((e+1) + [e even]) / 2."""
    odd = {q: e for q, e in f.items() if q % 2 == 1}
    tau = 1; chi_sum = 1; fair = 1.0
    for q, e in odd.items():
        tau *= e + 1
        leg = 1 if pow(p % q, (q - 1) // 2, q) == 1 else -1
        chi_sum *= (e + 1) if leg == 1 else (1 if e % 2 == 0 else 0)
        fair *= ((e + 1) + (1 if e % 2 == 0 else 0)) / 2
    res = {1: 1}  # residue mod m -> number of odd divisors with that residue
    for q, e in odd.items():
        new = defaultdict(int); r = q % m
        for d, c in res.items():
            pw = 1
            for _ in range(e + 1):
                new[d * pw % m] += c; pw = pw * r % m
        res = new
    return res.get(m - 1, 0), tau, chi_sum, fair

def load_published():
    key_all = []; key_wit = []; rng = {}
    for k, f in FILES.items():
        rows = json.load(gzip.open(R / f, "rt"))
        for r in rows:
            kk = r["p"] * 10000 + r["a"] * 100 + r["u"]; key_all.append(kk)
            if r["witness"]: key_wit.append(kk)
            rng[r["p"]] = k
        del rows
    return np.array(sorted(key_all), dtype=np.int64), np.array(sorted(key_wit), dtype=np.int64), rng

def quant(x, qs=(0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)):
    x = np.asarray(x, dtype=float)
    return {str(q): round(float(np.quantile(x, q)), 4) for q in qs}

def main():
    t0 = time.time()
    key_all, key_wit, rng = load_published()
    primes = sorted(rng); P = np.array(primes, dtype=np.int64)
    budget = np.array([(p - 1).bit_length() for p in primes]); B = int(budget.max())
    small = st.sieve(20000); assert small[-1] ** 2 > B * primes[-1] + B
    print(f"[{round(time.time() - t0)}s] {len(primes)} primes, budgets {sorted(set(budget.tolist()))}, published rows {len(key_all)}", file=sys.stderr, flush=True)
    W = np.zeros(len(primes)); Mt = np.zeros(len(primes)); S = np.zeros(len(primes), dtype=int); npairs = np.zeros(len(primes), dtype=int)
    Mt2 = np.zeros(len(primes)); Mt3 = np.zeros(len(primes)); triv = np.zeros(len(primes), dtype=int)
    pair_stats = {}; mismatches = 0; rows_done = 0
    # per-range per-pair accumulators of e = W - Mt: n, sum e, sum e^2 (for the variance decomposition)
    acc = defaultdict(lambda: [0, 0.0, 0.0])
    labels = np.array([rng[p] for p in primes])
    for a in range(1, B + 1):
        for u in range(1, B // a + 1):
            if gcd(a, u) != 1: continue
            m = 4 * a * u; ph = phi(m)
            sel = np.nonzero(budget >= a * u)[0]
            Ms, facs = st.factor_all(P[sel], a, u, small)
            w = np.zeros(len(sel)); mt = np.zeros(len(sel)); mt2 = np.zeros(len(sel)); mt3 = np.zeros(len(sel)); tv = np.zeros(len(sel), dtype=int)
            for j, (idx, Mi, f) in enumerate(zip(sel.tolist(), Ms.tolist(), facs)):
                prod = 1
                for q, e in f.items(): prod *= q ** e
                if prod != Mi: raise SystemExit(f"STOP: factorisation mismatch at {(primes[idx], a, u)}")
                cnt, tau, chi_sum, fair = witness_count(f, m, primes[idx])
                w[j] = cnt; mt[j] = tau / ph; mt2[j] = (tau - chi_sum) / ph; mt3[j] = (tau - fair) / ph; tv[j] = (Mi + 1) % m == 0
                if m == 4 and cnt != mt2[j]: raise SystemExit(f"STOP: (1,1) identity W = (tau - chi_sum)/2 fails at p={primes[idx]}")
            keys = P[sel] * 10000 + a * 100 + u
            pub = np.isin(keys, key_wit); present = np.isin(keys, key_all)
            if not present.all(): raise SystemExit(f"STOP: pair {(a, u)} missing from the published rows for some p")
            mismatches += int((pub != (w > 0)).sum())
            W[sel] += w; Mt[sel] += mt; Mt2[sel] += mt2; Mt3[sel] += mt3; S[sel] += (w > 0); npairs[sel] += 1; triv[sel] += tv
            e = w - mt
            for lab in set(labels[sel].tolist()):
                k = labels[sel] == lab; A = acc[(lab, a, u)]; A[0] += int(k.sum()); A[1] += float(e[k].sum()); A[2] += float((e[k] ** 2).sum())
            pair_stats[f"{a},{u}"] = {"rows": int(len(sel)), "phi_4au": ph, "mean_W": round(float(w.mean()), 4), "mean_Mt": round(float(mt.mean()), 4),
                                      "mean_Mt2_with_odd_real_character": round(float(mt2.mean()), 4), "mean_Mt3_fair_coin": round(float(mt3.mean()), 4),
                                      "mean_W_minus_Mt": round(float(e.mean()), 4), "sd_W_minus_Mt": round(float(e.std()), 4),
                                      "success_rate": round(float((w > 0).mean()), 4), "trivial_rate": round(float(tv.mean()), 4),
                                      "mean_W_given_trivial": round(float(w[tv == 1].mean()), 4) if tv.any() else None,
                                      "mean_W_given_nontrivial": round(float(w[tv == 0].mean()), 4)}
            rows_done += len(sel)
        print(f"[{round(time.time() - t0)}s] a={a} done, rows {rows_done}", file=sys.stderr, flush=True)
    n2p = np.array([n2(p) for p in primes])
    Rr = W / Mt; R2 = W / Mt2; R3 = W / Mt3
    out = {"primes": len(primes), "rows": rows_done, "published_witness_mismatches": mismatches, "budgets": sorted(set(budget.tolist())),
           "W_zero_primes": int((W == 0).sum()), "ranges": {}}
    out["overall"] = {"mean_W": round(float(W.mean()), 3), "mean_Mt": round(float(Mt.mean()), 3), "mean_S": round(float(S.mean()), 3),
                      "mean_pairs": round(float(npairs.mean()), 2), "mean_ratio": round(float(Rr.mean()), 4), "ratio_quantiles": quant(Rr),
                      "min_ratio": round(float(Rr.min()), 4), "frac_ratio_lt_0.5": round(float((Rr < 0.5).mean()), 4),
                      "frac_ratio_lt_0.25": round(float((Rr < 0.25).mean()), 4), "mean_trivial_witness_pairs": round(float(triv.mean()), 3),
                      "corr_W_Mt": round(float(np.corrcoef(W, Mt)[0, 1]), 4), "spearman_ratio_n2": None,
                      "mean_Mt2_with_odd_real_character": round(float(Mt2.mean()), 3), "mean_ratio_W_over_Mt2": round(float(R2.mean()), 4),
                      "ratio2_quantiles": quant(R2), "min_ratio2": round(float(R2.min()), 4), "corr_W_Mt2": round(float(np.corrcoef(W, Mt2)[0, 1]), 4),
                      "sd_W_minus_Mt2": round(float((W - Mt2).std()), 3), "mean_W_minus_Mt2": round(float((W - Mt2).mean()), 3),
                      "mean_Mt3_fair_coin": round(float(Mt3.mean()), 3), "mean_ratio_W_over_Mt3": round(float(R3.mean()), 4)}
    rk = lambda x: np.argsort(np.argsort(x, kind="mergesort"), kind="mergesort").astype(float)
    out["overall"]["spearman_ratio_n2"] = round(float(np.corrcoef(rk(Rr), rk(n2p))[0, 1]), 4)
    out["overall"]["spearman_W_n2"] = round(float(np.corrcoef(rk(W), rk(n2p))[0, 1]), 4)
    out["overall"]["spearman_Mt_n2"] = round(float(np.corrcoef(rk(Mt), rk(n2p))[0, 1]), 4)
    out["overall"]["spearman_Mt2_n2"] = round(float(np.corrcoef(rk(Mt2), rk(n2p))[0, 1]), 4)
    out["overall"]["spearman_ratio2_n2"] = round(float(np.corrcoef(rk(R2), rk(n2p))[0, 1]), 4)
    for lab in sorted(set(labels.tolist()), key=lambda l: int(P[labels == l].min())):
        k = labels == lab; d = W[k] - Mt[k]
        var_sum = float(d.var()); sum_var = 0.0
        for (l2, a, u), (n, s1, s2) in acc.items():
            if l2 == lab and n > 1: sum_var += (s2 / n - (s1 / n) ** 2) * n / len(d)  # weight by the fraction of primes that have the pair
        out["ranges"][lab] = {"p_min": int(P[k].min()), "p_max": int(P[k].max()), "primes": int(k.sum()), "budget": sorted(set(budget[k].tolist())),
                              "mean_W": round(float(W[k].mean()), 3), "mean_Mt": round(float(Mt[k].mean()), 3), "mean_S": round(float(S[k].mean()), 3),
                              "mean_ratio": round(float(Rr[k].mean()), 4), "ratio_quantiles": quant(Rr[k]), "min_ratio": round(float(Rr[k].min()), 4),
                              "frac_ratio_lt_0.5": round(float((Rr[k] < 0.5).mean()), 4), "frac_ratio_lt_0.25": round(float((Rr[k] < 0.25).mean()), 4),
                              "var_of_sum_W_minus_Mt": round(var_sum, 3), "sum_of_pair_vars": round(sum_var, 3),
                              "var_ratio": round(var_sum / sum_var, 4) if sum_var else None,
                              "mean_W_minus_Mt": round(float(d.mean()), 3), "sd_W_minus_Mt": round(float(d.std()), 3),
                              "mean_Mt2": round(float(Mt2[k].mean()), 3), "mean_ratio2": round(float(R2[k].mean()), 4), "ratio2_quantiles": quant(R2[k]),
                              "min_ratio2": round(float(R2[k].min()), 4), "sd_W_minus_Mt2": round(float((W[k] - Mt2[k]).std()), 3),
                              "mean_Mt3": round(float(Mt3[k].mean()), 3), "mean_ratio3": round(float(R3[k].mean()), 4)}
    order = np.argsort(Rr)[:15]
    out["lowest_ratio_primes"] = [{"p": int(P[i]), "range": labels[i], "W": int(W[i]), "Mt": round(float(Mt[i]), 3), "ratio": round(float(Rr[i]), 4),
                                   "Mt2": round(float(Mt2[i]), 3), "ratio2": round(float(R2[i]), 4),
                                   "S": int(S[i]), "pairs": int(npairs[i]), "n2": int(n2p[i])} for i in order]
    out["by_n2_ratio2"] = {}
    by_n2 = defaultdict(list)
    for i in range(len(primes)): by_n2[int(n2p[i])].append(i)
    out["by_n2"] = {str(k): {"primes": len(v), "mean_W": round(float(W[v].mean()), 3), "mean_Mt": round(float(Mt[v].mean()), 3),
                             "mean_ratio": round(float(Rr[v].mean()), 4), "min_ratio": round(float(Rr[v].min()), 4),
                             "mean_Mt2": round(float(Mt2[v].mean()), 3), "mean_ratio2": round(float(R2[v].mean()), 4)} for k, v in sorted(by_n2.items())}
    out["pairs"] = pair_stats
    out["seconds"] = round(time.time() - t0, 1)
    with gzip.GzipFile("mainterm_rows.json.gz", "wb", mtime=0) as fh:
        fh.write(json.dumps([{"p": int(P[i]), "range": labels[i], "B": int(budget[i]), "pairs": int(npairs[i]), "W": int(W[i]),
                              "Mt": round(float(Mt[i]), 6), "Mt2": round(float(Mt2[i]), 6), "Mt3": round(float(Mt3[i]), 6),
                              "S": int(S[i]), "trivial": int(triv[i]), "n2": int(n2p[i])} for i in range(len(primes))]).encode())
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
