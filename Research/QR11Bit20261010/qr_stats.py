"""Shared statistics for the QR-bit test (does the p mod q effect reduce to the bit (p|q)?).
Rows: (p, a, u, M, witness). For a prime q: rows with q | M and rows with q not dividing M; exposure = p is a
quadratic residue mod q ((p|q) = +1); outcome = the pair fails. MH odds ratio over strata = pair;
bootstrap over primes (rows of one p are dependent), the same resample for every statistic."""
from collections import defaultdict
import numpy as np

SMALL = ((1, 1), (1, 2), (2, 1), (1, 3), (3, 1))

def legendre(p, q):
    r = pow(p % q, (q - 1) // 2, q)
    return 1 if r == 1 else (-1 if r == q - 1 else 0)

def design(rows, q, wider_only=True):
    ps = sorted({r["p"] for r in rows}); pid = {p: i for i, p in enumerate(ps)}
    pairs = sorted({(r["a"], r["u"]) for r in rows}); kid = {au: i for i, au in enumerate(pairs)}
    P, K, D, E, Y = [], [], [], [], []
    for r in rows:
        if wider_only and (r["a"], r["u"]) in SMALL: continue
        P.append(pid[r["p"]]); K.append(kid[(r["a"], r["u"])])
        D.append(r["M"] % q == 0); E.append(legendre(r["p"], q) == 1); Y.append(not r["witness"])
    return {"ps": np.array(ps), "P": np.array(P), "K": np.array(K), "D": np.array(D), "E": np.array(E), "Y": np.array(Y),
            "nK": len(pairs)}

def mh(d, mask, w):
    rw = w[d["P"]] * mask
    K = d["K"]; S = d["nK"]; E = d["E"]; Y = d["Y"]
    a = np.bincount(K, rw * (E & Y), S); b = np.bincount(K, rw * (E & ~Y), S)
    c = np.bincount(K, rw * (~E & Y), S); e = np.bincount(K, rw * (~E & ~Y), S)
    n = a + b + c + e; ok = n > 0
    num = (a[ok] * e[ok] / n[ok]).sum(); den = (b[ok] * c[ok] / n[ok]).sum()
    return num / den if den > 0 else float("nan")

def analyse_q(rows, q, reps, seed):
    d = design(rows, q); G = len(d["ps"])
    masks = {"div": d["D"].astype(float), "nodiv": (~d["D"]).astype(float)}
    est = {k: mh(d, m, np.ones(G)) for k, m in masks.items()}
    rng = np.random.default_rng(seed); boot = {k: np.empty(reps) for k in masks}
    for b in range(reps):
        w = np.bincount(rng.integers(0, G, G), minlength=G).astype(float)
        for k, m in masks.items(): boot[k][b] = mh(d, m, w)
    out = {"q": q, "primes": G, "rows_div": int(d["D"].sum()), "rows_nodiv": int((~d["D"]).sum())}
    for k in masks:
        lb = np.log(boot[k])
        out[k] = {"OR": est[k], "ci95": [float(np.quantile(boot[k], .025)), float(np.quantile(boot[k], .975))],
                  "ci99": [float(np.quantile(boot[k], .005)), float(np.quantile(boot[k], .995))],
                  "p_one_sided_le_1": (1 + int((lb <= 0).sum())) / (reps + 1)}
    # fail rates by (divisible, QR) cell
    for dv in (True, False):
        for qr in (True, False):
            m = (d["D"] == dv) & (d["E"] == qr)
            out[f"fail_rate_{'div' if dv else 'nodiv'}_{'QR' if qr else 'NR'}"] = round(float(d["Y"][m].mean()), 4) if m.any() else None
    return out

def design_pair(rows, q, wider_only=True):
    """Within-pair contrast: exposure = q | M, strata = pair, pairs split by psi = (-a u | q)."""
    ps = sorted({r["p"] for r in rows}); pid = {p: i for i, p in enumerate(ps)}
    pairs = sorted({(r["a"], r["u"]) for r in rows}); kid = {au: i for i, au in enumerate(pairs)}
    psi = np.array([legendre((-a * u) % q, q) for a, u in pairs])
    P, K, D, Y = [], [], [], []
    for r in rows:
        if wider_only and (r["a"], r["u"]) in SMALL: continue
        if (r["a"] * r["u"]) % q == 0: continue  # q | au: q never divides M; drop
        P.append(pid[r["p"]]); K.append(kid[(r["a"], r["u"])]); D.append(r["M"] % q == 0); Y.append(not r["witness"])
    P = np.array(P); K = np.array(K)
    return {"ps": np.array(ps), "P": P, "K": K, "E": np.array(D), "Y": np.array(Y), "nK": len(pairs), "psi_row": psi[K]}

def analyse_pair(rows, q, reps, seed):
    d = design_pair(rows, q); G = len(d["ps"])
    masks = {"psi_neg": (d["psi_row"] == -1).astype(float), "psi_pos": (d["psi_row"] == 1).astype(float)}
    est = {k: mh(d, m, np.ones(G)) for k, m in masks.items()}
    rng = np.random.default_rng(seed); boot = {k: np.empty(reps) for k in masks}
    for b in range(reps):
        w = np.bincount(rng.integers(0, G, G), minlength=G).astype(float)
        for k, m in masks.items(): boot[k][b] = mh(d, m, w)
    ratio = np.log(boot["psi_neg"]) - np.log(boot["psi_pos"])
    out = {"q": q, "primes": G}
    for k in masks:
        out[k] = {"OR_div_vs_nodiv": est[k], "ci95": [float(np.quantile(boot[k], .025)), float(np.quantile(boot[k], .975))]}
    out["log_ratio_neg_over_pos"] = {"est": float(np.log(est["psi_neg"]) - np.log(est["psi_pos"])),
                                     "ci95": [float(np.quantile(ratio, .025)), float(np.quantile(ratio, .975))],
                                     "p_one_sided_ge_0": (1 + int((ratio >= 0).sum())) / (reps + 1)}
    return out

def prime_anova(rows, q=11, wider_only=True):
    """Per-prime wider fail fraction f_p; nested F-test: residue (q-1 classes) vs QR bit (2 classes)."""
    from collections import defaultdict
    acc = defaultdict(lambda: [0, 0])
    for r in rows:
        if wider_only and (r["a"], r["u"]) in SMALL: continue
        acc[r["p"]][0] += (not r["witness"]); acc[r["p"]][1] += 1
    ps = sorted(acc); f = np.array([acc[p][0] / acc[p][1] for p in ps]); res = np.array([p % q for p in ps])
    bit = np.array([legendre(p, q) for p in ps])
    n = len(f); grand = f.mean()
    def rss(groups):
        return sum(((f[groups == g] - f[groups == g].mean()) ** 2).sum() for g in np.unique(groups))
    rss0 = ((f - grand) ** 2).sum(); rss_bit = rss(bit); rss_res = rss(res)
    k_res = len(np.unique(res)); k_bit = len(np.unique(bit))
    F_within = ((rss_bit - rss_res) / (k_res - k_bit)) / (rss_res / (n - k_res))
    F_bit = ((rss0 - rss_bit) / (k_bit - 1)) / (rss_bit / (n - k_bit))
    from math import lgamma, exp, log
    def f_sf(F, d1, d2):  # survival function of F via regularized incomplete beta (continued fraction)
        x = d2 / (d2 + d1 * F); a, b = d2 / 2, d1 / 2
        return betainc_reg(a, b, x)
    means = {int(c): round(float(f[res == c].mean()), 4) for c in np.unique(res)}
    return {"primes": n, "F_bit": F_bit, "df_bit": [k_bit - 1, n - k_bit], "p_bit": f_sf(F_bit, k_bit - 1, n - k_bit),
            "F_within_class": F_within, "df_within": [k_res - k_bit, n - k_res], "p_within": f_sf(F_within, k_res - k_bit, n - k_res),
            "mean_f_by_residue": means, "mean_f_QR": round(float(f[bit == 1].mean()), 4), "mean_f_NR": round(float(f[bit == -1].mean()), 4)}

def betainc_reg(a, b, x):
    """Regularized incomplete beta I_x(a, b) (Numerical Recipes continued fraction)."""
    from math import lgamma, exp, log
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    lbeta = lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log(1 - x)
    def cf(a, b, x):
        qab, qap, qam = a + b, a + 1, a - 1; c, d = 1.0, 1 - qab * x / qap
        d = 1 / (d if abs(d) > 1e-300 else 1e-300); h = d
        for m in range(1, 1000):
            m2 = 2 * m; aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300); c = 1 + aa / c if abs(c) > 1e-300 else 1e300; h *= d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300); c = 1 + aa / c if abs(c) > 1e-300 else 1e300
            de = d * c; h *= de
            if abs(de - 1) < 1e-14: break
        return h
    if x < (a + 1) / (a + b + 2): return exp(lbeta) * cf(a, b, x) / a
    return 1 - exp(lbeta) * cf(b, a, 1 - x) / b

def fit_pair_model(rows, q=11, wider_only=True):
    """Per pair k: fail probability when q | M and when q does not divide M (pooled over the fitting rows)."""
    from collections import defaultdict
    acc = defaultdict(lambda: [0, 0, 0, 0])  # fail_div, n_div, fail_nodiv, n_nodiv
    for r in rows:
        au = (r["a"], r["u"])
        if wider_only and au in SMALL: continue
        A = acc[au]; dv = r["M"] % q == 0
        if dv: A[0] += (not r["witness"]); A[1] += 1
        else: A[2] += (not r["witness"]); A[3] += 1
    return {au: ((A[0] / A[1]) if A[1] else None, (A[2] / A[3]) if A[3] else None) for au, A in acc.items()}

def residue_calibration(rows, model, q=11, wider_only=True):
    """Per prime: observed wider fail fraction minus model prediction (pairs in its box, divisibility from p mod q).
    Nested F-test of the residuals: residue classes vs one common offset."""
    from collections import defaultdict
    acc = defaultdict(lambda: [0.0, 0.0, 0])
    for r in rows:
        au = (r["a"], r["u"])
        if wider_only and au in SMALL: continue
        pd, pn = model.get(au, (None, None))
        dv = r["M"] % q == 0; pr = pd if dv else pn
        if pr is None: continue
        A = acc[r["p"]]; A[0] += (not r["witness"]); A[1] += pr; A[2] += 1
    ps = sorted(acc); resid = np.array([(acc[p][0] - acc[p][1]) / acc[p][2] for p in ps]); res = np.array([p % q for p in ps])
    n = len(resid); rss0 = ((resid - resid.mean()) ** 2).sum()
    rss1 = sum(((resid[res == c] - resid[res == c].mean()) ** 2).sum() for c in np.unique(res)); k = len(np.unique(res))
    F = ((rss0 - rss1) / (k - 1)) / (rss1 / (n - k))
    p_value = betainc_reg((n - k) / 2, (k - 1) / 2, (n - k) / ((n - k) + (k - 1) * F))
    return {"primes": n, "offset": float(resid.mean()), "F_residue": F, "df": [k - 1, n - k], "p": p_value,
            "mean_resid_by_residue": {int(c): round(float(resid[res == c].mean()), 4) for c in np.unique(res)}}
