"""Measurement requested in Outcasts >>54 (not a judgement): for each prime p, in >>54's notation,
    M = sum U_{a,u} (U = tau_odd(ap+u)/phi(4au)),  E_{a,u} = W_{a,u} - U_{a,u},  E+ = sum max(E,0),  E- = sum max(-E,0),  R = E+ + E-,
so W_tot = M + E+ - E- and M - R = W_tot - 2 E+.  Habslog (>>54) is "M > R".  Also splits E- into the part carried by the
odd real character chi_{-au} (chi(q) = (p|q) on the prime factors), T_chi = prod_q (1 + chi(q) + ... + chi(q)^e) / phi(4au) >= 0,
and the rest.  Runs on the K-test range (1.6e7 < p < 1.9e7) and on the primes of the published rows (4e6 < p < 1.6e7)."""
import gzip, json, sys, time
from math import gcd
from pathlib import Path
import importlib.util
import numpy as np

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("mainterm_explore", HERE / "explore.py")
ex = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(ex)
st = ex.st
SQ840 = {1, 121, 169, 289, 361, 529}

def run(primes, label):
    t0 = time.time()
    P = np.array(primes, dtype=np.int64); budget = np.array([(p - 1).bit_length() for p in primes]); B = int(budget.max()); N = len(primes)
    small = st.sieve(22000); assert small[-1] ** 2 > B * int(P.max()) + B
    W = np.zeros(N); M = np.zeros(N); Ep = np.zeros(N); Em = np.zeros(N); Tchi = np.zeros(N); Em_chi = np.zeros(N); Em_rest = np.zeros(N)
    for a in range(1, B + 1):
        for u in range(1, B // a + 1):
            if gcd(a, u) != 1: continue
            m = 4 * a * u; ph = ex.phi(m); sel = np.nonzero(budget >= a * u)[0]
            Ms, facs = st.factor_all(P[sel], a, u, small)
            for idx, Mi, f in zip(sel.tolist(), Ms.tolist(), facs):
                cnt, tau, chi_sum, _ = ex.witness_count(f, m, primes[idx])
                U = tau / ph; e = cnt - U; t = chi_sum / ph
                W[idx] += cnt; M[idx] += U; Tchi[idx] += t
                if e > 0: Ep[idx] += e
                else:
                    Em[idx] += -e
                    # negative error split: the chi_{-au} term explains min(t, -e) of it; the rest is other characters
                    Em_chi[idx] += min(t, -e); Em_rest[idx] += max(0.0, -e - t)
    R = Ep + Em; MR = M - R
    q = lambda x: {k: round(float(np.quantile(x, float(k))), 3) for k in ("0.01", "0.05", "0.25", "0.5", "0.75", "0.95", "0.99")}
    out = {"label": label, "primes": N, "p_min": int(P.min()), "p_max": int(P.max()), "seconds": round(time.time() - t0, 1),
           "mean": {"W_tot": round(float(W.mean()), 3), "M": round(float(M.mean()), 3), "E_plus": round(float(Ep.mean()), 3), "E_minus": round(float(Em.mean()), 3),
                    "R": round(float(R.mean()), 3), "M_minus_R": round(float(MR.mean()), 3), "T_chi": round(float(Tchi.mean()), 3),
                    "E_minus_from_chi": round(float(Em_chi.mean()), 3), "E_minus_rest": round(float(Em_rest.mean()), 3)},
           "habslog_M_gt_R": {"primes": int((MR > 0).sum()), "fraction": round(float((MR > 0).mean()), 4)},
           "M_minus_R_quantiles": q(MR), "E_plus_quantiles": q(Ep), "E_minus_quantiles": q(Em), "E_minus_over_M_quantiles": q(Em / M),
           "W_minus_M_plus_Tchi_quantiles": q(W - M + Tchi), "max_M_minus_R": round(float(MR.max()), 3),
           "primes_with_M_gt_R": [int(P[i]) for i in np.nonzero(MR > 0)[0][:20]],
           "example_lowest_W": [{"p": int(P[i]), "W_tot": int(W[i]), "M": round(float(M[i]), 3), "E_plus": round(float(Ep[i]), 3), "E_minus": round(float(Em[i]), 3),
                                 "M_minus_R": round(float(MR[i]), 3), "T_chi": round(float(Tchi[i]), 3)} for i in np.argsort(W)[:5]]}
    rows = [{"p": int(P[i]), "W_tot": int(W[i]), "M": round(float(M[i]), 6), "E_plus": round(float(Ep[i]), 6), "E_minus": round(float(Em[i]), 6),
             "T_chi": round(float(Tchi[i]), 6)} for i in range(N)]
    return out, rows

def main():
    out = {}
    k_primes = [p for p in st.sieve(19 * 10**6) if 16 * 10**6 < p < 19 * 10**6 and p % 840 in SQ840]
    out["K_range"], rows_k = run(k_primes, "K range 1.6e7-1.9e7")
    print(f"K range done {out['K_range']['seconds']}s", file=sys.stderr, flush=True)
    pub = sorted({r["p"] for r in json.load(gzip.open(HERE / "mainterm_rows.json.gz", "rt"))})
    out["published_range"], rows_e = run(pub, "published rows 4e6-1.6e7")
    print(f"published range done {out['published_range']['seconds']}s", file=sys.stderr, flush=True)
    with gzip.GzipFile(HERE / "habslog_rows.json.gz", "wb", mtime=0) as fh: fh.write(json.dumps(rows_e + rows_k).encode())
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
