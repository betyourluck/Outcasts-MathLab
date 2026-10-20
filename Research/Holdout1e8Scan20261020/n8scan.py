"""N8 scan (Outcasts 80c02e64, contract.md in this directory): preregistered scan of the unused range
10^8 < p < 1.2*10^8 (six classes mod 840, 33,663 primes, B = 27 throughout) for
  H_Elog / H_E exceptions and the minima of E_log / E_all          (>>67, >>75 section 4, >>80 section 2)
  W_B = 0 (Hlog at budget 27: 75 pairs, actual exponents)          (>>80 section 2)
  the L = 8 bound c_mag of >>69 combined with c_old = (W - eta)_+  (>>75 section 4, >>80 section 3)
Own implementation on top of Research/CertificateRate20261016/gtest.py + chars20261015.py (same author, 2026-10-15/16)
and a copy of magnitude() from Research/L8Magnitude69Check20261018/check69_l8.py (same author, 2026-10-18).
Nothing of the other participants' code is imported or executed.

Usage (from anywhere):
  python -I n8scan.py selftest   # outside the range: 11 fixed inputs vs check69_rows.json.gz, the 38296729 box,
                                 # 300 development rows, timing on 50 primes >= 1.2e8 (writes selftest_results.json)
  python -I n8scan.py cohort     # the 33,663 primes; refuses to run unless the prime-list hash and the git state match
Outputs (cohort): n8_results.json, n8_primes.jsonl.gz, n8_rows.csv.gz next to this file.
"""
import sys, os, json, gzip, math, cmath, random, subprocess, time, hashlib
from collections import Counter
from fractions import Fraction as Fr
from sympy import isprime, totient

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
GT = os.path.join(ROOT, 'Research', 'CertificateRate20261016')
sys.path.insert(0, GT)
import gtest  # noqa: E402
from chars20261015 import characters, F_poly, reduce_mod_cyclo, ser  # noqa: E402

CLASSES = (1, 121, 169, 289, 361, 529)
LO, HI = 10 ** 8, 12 * 10 ** 7
NPRIMES = 33663
PRIME_SHA = 'e2cf39689d2a6d1b05b19b42b551533995060e16509769235ae8f94bc0440c21'
BUDGET = 27
SEED = 20261020
LIMIT = 4 * 3600
BAND_A = (0.04, 0.11)          # N5(a): band for the row rescue rate (development data 355/4709 = 0.075, per pair 5-9%)
FIXED11 = [66529, 345601, 670849, 1740481, 3678481, 5843041, 8615161, 18490201, 53722321, 84525841, 38296729]
SQ2 = math.sqrt(2)
T0 = time.time()


# ---------------------------------------------------------------- primes, n2, Legendre

def sieve_range(lo, hi):
    r = math.isqrt(hi) + 1
    small = bytearray([1]) * (r + 1)
    small[0] = small[1] = 0
    for i in range(2, math.isqrt(r) + 1):
        if small[i]:
            small[i * i::i] = bytes(len(range(i * i, r + 1, i)))
    seg = bytearray([1]) * (hi - lo)
    for q in (i for i in range(2, r + 1) if small[i]):
        s = max(q * q, (lo + q - 1) // q * q)
        seg[s - lo::q] = bytes(len(range(s, hi, q)))
    return [lo + i for i in range(hi - lo) if seg[i]]


def legendre(n, p):
    v = pow(n % p, (p - 1) // 2, p)
    return -1 if v == p - 1 else v


def n2(p):
    """least quadratic non-residue of p (odd primes in increasing order, Euler's criterion; p = 1 mod 4 so (q|p) = (p|q))."""
    q = 3
    while True:
        if isprime(q) and pow(p % q, (q - 1) // 2, q) == q - 1:
            return q
        q += 2


# ---------------------------------------------------------------- c_mag (copy of check69_l8.py, 2026-10-18)

def ceil_sqrt_1024(U):
    """smallest t with (t/1024)^2 >= U, U a non-negative Fraction; integer comparisons only."""
    assert U >= 0
    n, d = U.numerator * 1024 ** 2, U.denominator
    t = math.isqrt(n // d)
    while t * t * d < n:
        t += 1
    assert t * t * d >= n and (t == 0 or (t - 1) * (t - 1) * d < n)
    return Fr(t, 1024)


def mul_dict(a, b, L):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            k = (i + j) % L
            out[k] = out.get(k, 0) + x * y
    return out


def pad4(cs):
    cs = list(cs) + [0] * (4 - len(cs))
    assert len(cs) == 4
    return cs


def local_norm_cos(h, e):
    """|sum_{j=0}^{e} zeta8^{hj}|^2 = (e+1) + 2 sum_{j=1}^{e} (e+1-j) cos(pi h j / 4) as (g, d) with value g + d*sqrt2."""
    g, d = e + 1, 0
    for j in range(1, e + 1):
        w = e + 1 - j
        t = (h * j) % 8
        if t == 0:
            g += 2 * w
        elif t == 4:
            g -= 2 * w
        elif t in (2, 6):
            pass
        elif t in (1, 7):
            d += w
        else:
            d -= w
    return g, d


def mul_sqrt2(x, y):
    return x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def magnitude(p, a, u, odd, m, W, R):
    """c_mag of >>69 from the factorization and the modulus only (plus R, which comes from the real characters)."""
    L, units, chars, bykey = characters(m)
    assert L == 8
    phi = int(totient(m))
    zeta = cmath.exp(1j * math.pi / 4)
    seen = set()
    T4 = 0
    terms = []
    sum_re8 = 0
    sum_irr8 = 0
    for c in chars:
        if c['real'] or c['param'] in seen:
            continue
        conj = bykey[tuple((-c['phase'][r]) % L for r in units)]
        assert conj['param'] != c['param'] and conj['sign'] == c['sign']
        seen.update((c['param'], conj['param']))
        d = L // math.gcd(L, *c['phase'].values())
        F = F_poly(c['phase'], odd, m, L)
        red = pad4(reduce_mod_cyclo(F, L))
        Ff = 1
        for q, e in odd:
            z = zeta ** c['phase'][q % m]
            Ff *= sum(z ** k for k in range(e + 1))
        assert abs(Ff - (red[0] + red[1] * zeta + red[2] * zeta ** 2 + red[3] * zeta ** 3)) < 1e-6
        if d == 4:
            assert red[1] == 0 and red[3] == 0
            T4 += c['sign'] * 2 * red[0]
        else:
            assert d == 8
            Fbar = F_poly(conj['phase'], odd, m, L)
            prod = pad4(reduce_mod_cyclo(mul_dict(F, Fbar, L), L))
            G, D = prod[0], prod[1]
            assert prod == [G, D, 0, -D], (p, a, u, prod)
            gd = (1, 0)
            for q, e in odd:
                gd = mul_sqrt2(gd, local_norm_cos(c['phase'][q % m], e))
            assert gd == (G, D), (p, a, u, gd, (G, D))
            assert abs(abs(Ff) ** 2 - (G + D * SQ2)) < 1e-6
            U = Fr(G) + (Fr(3, 2) if D >= 0 else Fr(7, 5)) * D
            assert U >= G + D * SQ2 - 1e-9
            rho = ceil_sqrt_1024(U)
            assert rho * rho >= U
            assert float(rho) >= abs(Ff) - 1e-9
            assert c['sign'] * 2 * Ff.real >= -2 * float(rho) - 1e-9
            sum_re8 += c['sign'] * 2 * red[0]
            sum_irr8 += c['sign'] * (red[1] - red[3])
            terms.append({'gamma': G, 'delta': D, 'U': U, 'rho': rho})
    assert len(seen) == sum(1 for c in chars if not c['real'])
    assert sum_irr8 == 0, (p, a, u, sum_irr8)
    assert phi * W == phi * R + T4 + sum_re8, (p, a, u, phi * W, phi * R, T4, sum_re8)
    raw = R + Fr(T4, phi) - Fr(2, phi) * sum((t['rho'] for t in terms), Fr(0))
    c_mag = max(Fr(0), raw)
    assert 0 <= c_mag <= W
    return {'T4': T4, 'terms': terms, 'raw_lower': raw, 'c_mag': c_mag}


# ---------------------------------------------------------------- one prime

def l8(p, a, u, r):
    """c_mag fields for an L = 8 box row r (output of gtest.row); re-factors N (6 rows per prime, cheap)."""
    N, m = a * p + u, 4 * a * u
    f, odd, divs, tau, W = gtest.divisor_data(N, m)
    assert W == r['W'] and r['L'] == 8 and r['m'] == m
    R = Fr(r['W_G2'], r['G2size'])
    nonres = sum(e for q, e in f.items() if legendre(q, p) == -1)
    assert p % 840 in CLASSES and legendre(u, p) == 1
    assert nonres % 2 == 0 and (W == 0 or nonres >= 2)
    for d in divs:
        if d % m == m - 1:
            assert legendre(d, p) == -1 and legendre(N // d, p) == -1
    mag = magnitude(p, a, u, odd, m, W, R)
    assert nonres >= 2 or mag['c_mag'] == 0
    return {'R': R, 'T4': mag['T4'], 'raw': mag['raw_lower'], 'c_mag': mag['c_mag'], 'nonres': nonres,
            'rho_sum': sum((t['rho'] for t in mag['terms']), Fr(0))}


def prime_summary(p, rng, full_check=False, expect_B=None):
    B, rows = gtest.prime_rows(p, rng, full_check)
    if expect_B is not None:
        assert B == expect_B, (p, B)
    box = [r for r in rows if r['in_budget']]
    for r in box:
        if r['L'] == 8:
            r.update(l8(p, r['a'], r['u'], r))
            r['c_comb'] = max(r['c'], r['c_mag'])
        else:
            r['c_comb'] = r['c']
        assert 0 <= r['c'] <= r['c_comb'] <= r['W'] and math.ceil(r['c_comb']) <= r['W'], (p, r['a'], r['u'])
    s = gtest.summarize(p, B, rows)
    l8rows = [r for r in box if r['L'] == 8]
    assert all(r['W'] > 0 for r in l8rows if r['c_mag'] > 0)
    s.update({
        'n2': n2(p), 'W_B': s['Wtot'], 'Cert_old': s['Cert'], 'Cert_comb': sum((r['c_comb'] for r in box), Fr(0)),
        'W_pos_pairs': [[r['a'], r['u'], r['W'], int(r['in_budget'])] for r in rows if r['W'] > 0],
        'l8_rows': len(l8rows), 'l8_Wpos': sum(1 for r in l8rows if r['W'] > 0),
        'l8_old0': sum(1 for r in l8rows if r['c'] == 0),
        'l8_old0_Wpos': sum(1 for r in l8rows if r['c'] == 0 and r['W'] > 0),
        'l8_rescued': sum(1 for r in l8rows if r['c'] == 0 and r['c_mag'] > 0),
        'l8_cmag_pos': sum(1 for r in l8rows if r['c_mag'] > 0),
        'l8_both': sum(1 for r in l8rows if r['c'] > 0 and r['c_mag'] > 0),
        'l8_loses': sum(1 for r in l8rows if r['c'] > 0 and r['c_mag'] == 0),
        'l8_detail': [{'a': r['a'], 'u': r['u'], 'W': r['W'], 'c_old': r['c'], 'c_mag': r['c_mag'], 'T4': r['T4'],
                       'R': r['R'], 'nonres': r['nonres']} for r in l8rows if r['W'] > 0],
        'full_checked': full_check})
    assert s['Cert_old'] <= s['Cert_comb'] <= s['W_B']
    return s, rows


def brief(s):
    keys = ('p', 'B', 'n2', 'W_B', 'S', 'E_log', 'E_all', 'Cert_old', 'Cert_comb', 'Cert_E', 'Cert_0notE', 'Cert_plus',
            'quad', 'group', 'expo', 'l8_rescued', 'full_checked', 'W_pos_pairs', 'l8_detail')
    return {k: s[k] for k in keys}


# ---------------------------------------------------------------- judgement

def quantiles(xs, qs=(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 1)):
    xs = sorted(xs)
    return {str(q): xs[min(len(xs) - 1, int(q * (len(xs) - 1)))] for q in qs}


def spearman(x, y):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    mu = (len(x) + 1) / 2
    num = sum((a - mu) * (b - mu) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mu) ** 2 for a in rx) * sum((b - mu) ** 2 for b in ry))
    return num / den


def judge(sums):
    n = len(sums)
    helog = [s for s in sums if s['E_log'] == 0]
    he = [s for s in sums if s['E_all'] == 0]
    wz = [s for s in sums if s['W_B'] == 0]
    old0 = [s for s in sums if s['Cert_old'] == 0]
    boxres = [s for s in old0 if s['Cert_comb'] > 0]
    miss = [s for s in sums if s['Cert_comb'] == 0 and s['W_B'] > 0]
    denomA = sum(s['l8_old0_Wpos'] for s in sums)
    denomB = sum(s['l8_old0'] for s in sums)
    resc = sum(s['l8_rescued'] for s in sums)
    rateA = resc / denomA if denomA else None
    el_min = min(s['E_log'] for s in sums)
    ea_min = min(s['E_all'] for s in sums)
    res = {
        'N1_H_Elog_exceptions': {'count': len(helog), 'forecast': '0 or 1', 'primes': [brief(s) for s in helog]},
        'N2_H_E_exceptions': {'count': len(he), 'forecast': 0, 'primes': [brief(s) for s in he]},
        'N3_W_B_zero': {'count': len(wz), 'forecast': 0, 'primes': [brief(s) for s in wz]},
        'N4_minima': {'E_log_min': el_min, 'E_all_min': ea_min, 'forecast': 'E_log_min <= 2 and E_all_min <= 2',
                      'pass': el_min <= 2 and ea_min <= 2,
                      'primes_E_log_le1': [brief(s) for s in sums if s['E_log'] <= 1],
                      'primes_E_all_le2': [brief(s) for s in sums if s['E_all'] <= 2]},
        'N5_cmag': {'a_rescued_rows': resc, 'a_denominator_rows_Wpos_old0': denomA, 'a_rate': rateA, 'a_band': list(BAND_A),
                    'a_pass': rateA is not None and BAND_A[0] <= rateA <= BAND_A[1],
                    'a_rate_all_old0_rows': resc / denomB if denomB else None, 'a_denominator_all_old0_rows': denomB,
                    'b_box_rescued_primes': len(boxres), 'b_forecast': '0 or 1', 'b_primes': [brief(s) for s in boxres],
                    'primes_Cert_old_zero': len(old0), 'primes_Cert_old_zero_list': [brief(s) for s in old0],
                    'c_soundness_rows_asserted': sum(s['l8_rows'] for s in sums),
                    'd_comb_zero_W_pos_primes': len(miss), 'd_forecast': 0, 'd_primes': [brief(s) for s in miss]},
        'N6_order_0_le_old_le_comb_le_W_all_box_rows': True,
    }
    res['branches'] = {'R1_H_Elog_exception': len(helog) > 0, 'R2_H_E_exception': len(he) > 0,
                       'R3_W_B_zero': len(wz) > 0, 'R4_box_rescued_by_cmag': len(boxres) > 0,
                       'R5_bound_misses_existing_witness': len(miss) > 0,
                       'R6_rescue_rate_out_of_band': not res['N5_cmag']['a_pass'],
                       'R7_no_exception': not (helog or he or wz or boxres or miss)}
    Wt = sum(s['W_B'] for s in sums)
    Co = sum((s['Cert_old'] for s in sums), Fr(0))
    Cc = sum((s['Cert_comb'] for s in sums), Fr(0))
    order = sorted(sums, key=lambda s: (s['E_log'], s['E_all'], s['p']))
    pos = [s for s in sums if s['W_B'] > 0]
    fails = sum(s['fail'] for s in sums)
    second = sum(s['group'] + s['expo'] for s in sums)
    l8pair = {}
    rescued_rows = []
    for s in sums:
        for d in s['l8_detail']:
            key = f"{d['a']},{d['u']}"
            e = l8pair.setdefault(key, {'rows_Wpos': 0, 'old0_Wpos': 0, 'rescued': 0, 'both': 0})
            e['rows_Wpos'] += 1
            if d['c_old'] == 0:
                e['old0_Wpos'] += 1
                if d['c_mag'] > 0:
                    e['rescued'] += 1
                    rescued_rows.append(d)
            elif d['c_mag'] > 0:
                e['both'] += 1
    for e in l8pair.values():
        e['rate'] = e['rescued'] / e['old0_Wpos'] if e['old0_Wpos'] else None
    n2s = [s['n2'] for s in sums]
    res['measure'] = {
        'primes': n, 'W_B_sum': Wt, 'Cert_old_sum': Co, 'Cert_comb_sum': Cc,
        'Cert_old_over_W_B_pooled': float(Co / Wt) if Wt else None, 'Cert_comb_over_W_B_pooled': float(Cc / Wt) if Wt else None,
        'Cert_old_over_W_B_per_prime': quantiles([float(s['Cert_old'] / s['W_B']) for s in pos]),
        'Cert_comb_over_W_B_per_prime': quantiles([float(s['Cert_comb'] / s['W_B']) for s in pos]),
        'Cert_E_sum': sum(s['Cert_E'] for s in sums), 'Cert_0notE_sum': sum(s['Cert_0notE'] for s in sums),
        'Cert_plus_sum': sum((s['Cert_plus'] for s in sums), Fr(0)),
        'Cert_E_over_Cert_old_pooled': float(sum(s['Cert_E'] for s in sums) / Co) if Co else None,
        'W_B_quantiles': quantiles([s['W_B'] for s in sums]), 'S_quantiles': quantiles([s['S'] for s in sums]),
        'E_log_quantiles': quantiles([s['E_log'] for s in sums]), 'E_all_quantiles': quantiles([s['E_all'] for s in sums]),
        'E_log_hist_small': {str(v): sum(1 for s in sums if s['E_log'] == v) for v in range(0, 6)},
        'E_all_hist_small': {str(v): sum(1 for s in sums if s['E_all'] == v) for v in range(0, 6)},
        'lowest_20_E_log': [brief(s) for s in order[:20]],
        'primes_Cert_comb_gt_Cert_old': sum(1 for s in sums if s['Cert_comb'] > s['Cert_old']),
        'l8_rows': sum(s['l8_rows'] for s in sums), 'l8_Wpos_rows': sum(s['l8_Wpos'] for s in sums),
        'l8_cmag_pos_rows': sum(s['l8_cmag_pos'] for s in sums), 'l8_both_rows': sum(s['l8_both'] for s in sums),
        'l8_loses_rows': sum(s['l8_loses'] for s in sums), 'l8_by_pair': l8pair,
        'rescued_T4_hist': dict(Counter(str(d['T4']) for d in rescued_rows)),
        'rescued_R_hist': dict(Counter(str(d['R']) for d in rescued_rows)),
        'rescued_c_mag_max': max((d['c_mag'] for d in rescued_rows), default=Fr(0)),
        'c_mag_max': max((d['c_mag'] for s in sums for d in s['l8_detail']), default=Fr(0)),
        'failure_rows': fails, 'quad_share_of_failures': sum(s['quad'] for s in sums) / fails if fails else None,
        'group_share_of_second_stage': sum(s['group'] for s in sums) / second if second else None,
        'second_uniform_share': sum(s['second_uniform'] for s in sums) / second if second else None,
        'second_oddH_all_exponent': all(s['second_oddH'] == s['second_oddH_expo'] for s in sums),
        'W_pos_c0_rows_old': sum(len(s['W_pos_c0_rows']) for s in sums),
        'W_pos_c0_by_pair_old': dict(Counter(f"{a},{u}" for s in sums for a, u in s['W_pos_c0_rows']).most_common(15)),
        'n2_hist': dict(sorted(Counter(str(v) for v in n2s).items(), key=lambda t: int(t[0]))),
        'mean_S_by_n2': {str(v): sum(s['S'] for s in sums if s['n2'] == v) / c for v, c in sorted(Counter(n2s).items())},
        'spearman_n2_S': spearman(n2s, [s['S'] for s in sums]),
        'spearman_n2_W_B': spearman(n2s, [s['W_B'] for s in sums]),
        'full_checked_primes': sum(1 for s in sums if s['full_checked']),
    }
    return res


# ---------------------------------------------------------------- modes

def git(*args):
    return subprocess.run(['git', '-C', ROOT, *args], capture_output=True, text=True)


def cohort():
    primes = [p for p in sieve_range(LO, HI) if p % 840 in CLASSES]
    assert len(primes) == NPRIMES, len(primes)
    assert all(isprime(p) for p in primes) and all((p - 1).bit_length() == BUDGET for p in primes)
    h = hashlib.sha256(''.join(f'{p}\n' for p in primes).encode()).hexdigest()
    assert h == PRIME_SHA, h
    head = git('rev-parse', 'HEAD').stdout.strip()
    rel = os.path.relpath(HERE, ROOT).replace(os.sep, '/')
    tracked = git('ls-files', rel).stdout.split()
    for fn in ('n8scan.py', 'contract.md'):
        assert f'{rel}/{fn}' in tracked, fn
    assert git('diff', '--quiet', 'HEAD', '--', rel).returncode == 0, 'working tree differs from HEAD'
    print(f'G0 ok: primes {len(primes)} sha256 {h}, HEAD {head}, {time.time() - T0:.0f}s', flush=True)
    rng = random.Random(SEED)
    sums = []
    fp = gzip.open(os.path.join(HERE, 'n8_rows.csv.gz'), 'wt', encoding='utf-8', newline='\n')
    fp.write('p,a,u,in_budget,inE,L,W,W_G2,G2size,m1K,cls,uniform,eta,c_old,nonconst,c_mag,T4,raw,nonres,c_comb\n')
    for i, p in enumerate(primes):
        s, rows = prime_summary(p, rng, expect_B=BUDGET)
        if s['E_log'] <= 1 or s['E_all'] <= 1 or s['W_B'] <= 2:
            # independent route on every row of this prime: trial division up to sqrt N for all divisors,
            # sympy cyclotomic reduction for the characters (primality of every factor is asserted on every row anyway)
            s2, rows2 = prime_summary(p, random.Random(0), full_check=True, expect_B=BUDGET)
            assert [(r['a'], r['u'], r['W']) for r in rows2] == [(r['a'], r['u'], r['W']) for r in rows]
            assert (s2['Cert_old'], s2['Cert_comb'], s2['E_log'], s2['E_all']) == (s['Cert_old'], s['Cert_comb'], s['E_log'], s['E_all'])
            s['full_checked'] = True
        sums.append(s)
        for r in rows:
            if r['in_budget']:
                tail = f",{r['c_mag']},{r['T4']},{r['raw']},{r['nonres']}" if r['L'] == 8 else ',,,,'
                fp.write(f"{p},{r['a']},{r['u']},1,{int(r['inE'])},{r['L']},{r['W']},{r['W_G2']},{r['G2size']},{int(r['m1K'])},"
                         f"{r['cls']},{int(r['uniform'])},{r['eta']},{r['c']},{r['nonconst']}{tail},{r['c_comb']}\n")
            else:
                fp.write(f"{p},{r['a']},{r['u']},0,1,,{r['W']},,,,,,,,,,,,,\n")
        if (i + 1) % 1000 == 0:
            print(f'{i + 1}/{len(primes)} primes, {time.time() - T0:.0f}s', flush=True)
        assert time.time() - T0 < LIMIT
    fp.close()
    with gzip.open(os.path.join(HERE, 'n8_primes.jsonl.gz'), 'wt', encoding='utf-8', newline='\n') as g:
        for s in sums:
            g.write(json.dumps(ser(s), ensure_ascii=False) + '\n')
    res = judge(sums)
    res['G0'] = {'primes': len(primes), 'prime_list_sha256': h, 'HEAD': head, 'B': BUDGET,
                 'box_rows': sum(s['pairs'] for s in sums), 'extra_E_rows': sum(s['extraE'] for s in sums),
                 'E_rows_in_box_eta_checked': sum(s['E_pairs_box'] for s in sums),
                 'sample_rows_trial_division_and_sympy': sum(s['sampleT'] for s in sums),
                 'full_checked_primes': sum(1 for s in sums if s['full_checked'])}
    res['seconds'] = round(time.time() - T0, 1)
    json.dump(ser(res), open(os.path.join(HERE, 'n8_results.json'), 'w', encoding='utf-8', newline='\n'), ensure_ascii=False, indent=1)
    sys.stdout.reconfigure(encoding='utf-8')
    short = {k: (v if k not in ('measure',) else {kk: vv for kk, vv in v.items() if kk != 'lowest_20_E_log'}) for k, v in res.items()}
    print(json.dumps(ser(short), ensure_ascii=False, indent=1))


def selftest():
    """Outside the range only. Does not touch the cohort."""
    ref = json.load(gzip.open(os.path.join(ROOT, 'Research', 'L8Magnitude69Check20261018', 'check69_rows.json.gz'), 'rt', encoding='utf-8'))
    out = {}
    # 1. the 11 fixed inputs: every L = 8 box row vs check69_rows.json.gz (fixed42), own 2026-10-18 computation
    mine = {}
    boxes = {}
    for p in FIXED11:
        s, rows = prime_summary(p, random.Random(SEED))
        boxes[p] = s
        for r in rows:
            if r['in_budget'] and r['L'] == 8:
                mine[(p, r['a'], r['u'])] = r
    fields = {'W': 'W', 'old_c': 'c', 'R': 'R', 'T4': 'T4', 'raw_lower': 'raw', 'c_mag': 'c_mag', 'nonresidue_multiplicity': 'nonres'}
    mism = []
    for t in ref['fixed42']:
        key = (t['p'], t['a'], t['u'])
        for theirs, ours in fields.items():
            a, b = t[theirs], mine[key][ours]
            if isinstance(a, str):
                a = Fr(a)
            if a != b:
                mism.append([key, theirs, str(a), str(b)])
    out['fixed42'] = {'rows': len(ref['fixed42']), 'rows_mine': len(mine), 'mismatches': mism}
    assert len(ref['fixed42']) == 42 == len(mine) and not mism
    # 2. the 38296729 box: W, c_old, eta per row vs rows73; W of the 81 E rows vs E81; summary values of >>66 / >>69
    s = boxes[38296729]
    s2, rows = prime_summary(38296729, random.Random(SEED))
    byk = {(r['a'], r['u']): r for r in rows}
    m73 = [[t['a'], t['u']] for t in ref['rows73']
           if (byk[(t['a'], t['u'])]['W'], byk[(t['a'], t['u'])]['c'], byk[(t['a'], t['u'])]['eta']) != (t['W'], Fr(t['old_c']), Fr(t['eta']))]
    m81 = [[t['a'], t['u']] for t in ref['E81'] if byk[(t['a'], t['u'])]['W'] != t['W']]
    vals = {'B': s['B'], 'pairs': s['pairs'], 'W_B': s['W_B'], 'S': s['S'], 'E_log': s['E_log'], 'E_all': s['E_all'],
            'Cert_old': s['Cert_old'], 'Cert_comb': s['Cert_comb'], 'ceil_Cert_comb': math.ceil(s['Cert_comb']),
            'Cert_E': s['Cert_E'], 'Cert_0notE': s['Cert_0notE'], 'Cert_plus': s['Cert_plus'], 'n2': s['n2'],
            'l8_detail': s['l8_detail'], 'W_pos_pairs': s['W_pos_pairs']}
    out['box_38296729'] = {'mismatches_rows73': m73, 'mismatches_E81': m81, 'values': vals}
    assert not m73 and not m81
    assert (vals['B'], vals['pairs'], vals['W_B'], vals['S'], vals['E_log'], vals['E_all']) == (26, 73, 6, 4, 0, 1)
    assert vals['Cert_old'] == 1 and vals['Cert_comb'] == Fr(10451, 8192) and vals['ceil_Cert_comb'] == 2
    # 3. 300 development rows (G range, first 300 of dev4709 in file order) vs own 2026-10-18 values
    mism = []
    for t in ref['dev4709'][:300]:
        p, a, u = t['p'], t['a'], t['u']
        r = gtest.row(p, a, u, (p - 1).bit_length(), random.Random(0))
        r.update(l8(p, a, u, r))
        for theirs, ours in fields.items():
            x = Fr(t[theirs]) if isinstance(t[theirs], str) else t[theirs]
            if x != r[ours]:
                mism.append([[p, a, u], theirs, str(x), str(r[ours])])
    out['dev300'] = {'rows': 300, 'mismatches': mism}
    assert not mism
    # 4. timing on the first 50 six-class primes >= 1.2e8 (outside the range; B = 27 as in the range), 5 of them with full_check
    tp = [p for p in sieve_range(HI, HI + 200000) if p % 840 in CLASSES][:50]
    assert len(tp) == 50 and all(isprime(p) for p in tp)
    t1 = time.time()
    ts = [prime_summary(p, random.Random(SEED), expect_B=BUDGET)[0] for p in tp]
    t2 = time.time()
    tf = [prime_summary(p, random.Random(0), full_check=True, expect_B=BUDGET)[0] for p in tp[:5]]
    t3 = time.time()
    for a, b in zip(ts[:5], tf):
        assert (a['W_B'], a['Cert_old'], a['Cert_comb'], a['E_log'], a['E_all']) == (b['W_B'], b['Cert_old'], b['Cert_comb'], b['E_log'], b['E_all'])
    out['timing'] = {'primes': tp, 'seconds_per_prime': (t2 - t1) / 50, 'seconds_per_prime_full_check': (t3 - t2) / 5,
                     'estimate_seconds_33663': (t2 - t1) / 50 * NPRIMES,
                     'summary_50': {'W_B_min': min(s['W_B'] for s in ts), 'E_log_min': min(s['E_log'] for s in ts),
                                    'E_all_min': min(s['E_all'] for s in ts), 'Cert_old_zero': sum(1 for s in ts if s['Cert_old'] == 0),
                                    'l8_old0_Wpos': sum(s['l8_old0_Wpos'] for s in ts), 'l8_rescued': sum(s['l8_rescued'] for s in ts)}}
    out['seconds'] = round(time.time() - T0, 1)
    json.dump(ser(out), open(os.path.join(HERE, 'selftest_results.json'), 'w', encoding='utf-8', newline='\n'), ensure_ascii=False, indent=1)
    sys.stdout.reconfigure(encoding='utf-8')
    print(json.dumps(ser({k: v for k, v in out.items() if k != 'box_38296729'} | {'box_38296729_values': vals}), ensure_ascii=False, indent=1))


if __name__ == '__main__':
    {'selftest': selftest, 'cohort': cohort}[sys.argv[1]]()
