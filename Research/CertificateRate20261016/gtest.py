"""G test (Outcasts 80c02e64, contract.md in this directory): factor-only certificate c_i = (W_i - eta_i)_+, the
guaranteed-exact family E, H_Elog / H_E, and the three failure types of >>62, on the 5,625 primes of the K test.
This is a prespecified reanalysis / implementation audit of an already-used cohort, not a new sample.

Usage (from anywhere):
  python -I gtest.py selftest     # the 10 discovery inputs (not the cohort); compares with ExactFamily65Check / Slack60G2Check
  python -I gtest.py cohort       # the 5,625 primes; refuses to run unless the prime list hash and the git state match
Outputs (cohort): gtest_results.json, gtest_primes.jsonl.gz, gtest_rows.csv.gz next to this file.
"""
import sys, os, json, gzip, math, random, subprocess, time, hashlib
from collections import Counter
from fractions import Fraction as Fr
from sympy import factorint, isprime, cyclotomic_poly, Poly, symbols

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
from chars20261015 import characters, reduce_mod_cyclo, F_poly, add_polys, divisors_direct, divisors_from_factors, ser  # noqa: E402

E = frozenset({1, 2, 3, 4, 5, 6, 7, 9, 10, 12, 14, 15, 18, 20, 21, 30, 42, 60, 63, 126})
SMALL_L = {1, 2, 3, 4, 6}
SEED = 20261016
PRIME_SHA = 'e3cb7f64159ea509e2c856c554d80d383a33ff2fc60be4a1c6b8f506c3b469fc'
DISCOVERY = [66529, 345601, 670849, 1740481, 3678481, 5843041, 8615161, 18490201, 53722321, 84525841]
X = symbols('x')

# ---------------------------------------------------------------- algebra helpers

_cyc = {}


def cyc(L):
    if L not in _cyc:
        _cyc[L] = [int(c) for c in Poly(cyclotomic_poly(L, X), X).all_coeffs()[::-1]]  # low -> high, monic
    return _cyc[L]


def reduce_fast(poly, L):
    """Remainder of sum_k c_k x^k (k taken mod L) modulo Phi_L, coefficients low -> high, trailing zeros stripped."""
    arr = [0] * L
    for k, c in poly.items():
        arr[k % L] += c
    ph = cyc(L)
    d = len(ph) - 1
    for i in range(L - 1, d - 1, -1):
        c = arr[i]
        if c:
            for j in range(d + 1):
                arr[i - d + j] -= c * ph[j]
    out = arr[:d] if d > 0 else [0]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out if out else [0]


def closure(m, gens, start):
    S = set(start)
    frontier = list(S)
    while frontier:
        new = []
        for g in frontier:
            for q in gens:
                h = g * q % m
                if h not in S:
                    S.add(h)
                    new.append(h)
        frontier = new
    return S


_grp = {}


def group(m):
    """Units, G^2, -G^2, coset index of each unit in G/G^2, lambda(m) (from the character module)."""
    if m not in _grp:
        units = [r for r in range(1, m) if math.gcd(r, m) == 1]
        G2 = frozenset(r * r % m for r in units)
        negG2 = frozenset((-s) % m for s in G2)
        idx = {}
        reps = []
        for r in units:
            if r not in idx:
                k = len(reps)
                reps.append(r)
                for s in G2:
                    idx[r * s % m] = k
        assert len(reps) * len(G2) == len(units)
        _grp[m] = {'units': units, 'G2': G2, 'negG2': negG2, 'coset': idx, 'ncosets': len(reps)}
    return _grp[m]


_chr = {}


def charinfo(m):
    if m not in _chr:
        L, units, chars, bykey = characters(m)
        principal = tuple(0 for _ in chars[0]['param'])
        reals = [(c['phase'], c['sign'], c['param'] == principal) for c in chars if c['real']]
        seen, pairs = set(), []
        for c in chars:
            if c['real'] or c['param'] in seen:
                continue
            conj = bykey[tuple((-c['phase'][r]) % L for r in units)]
            seen.update((c['param'], conj['param']))
            pairs.append((c['phase'], conj['phase'], c['sign']))
        assert len(reals) + 2 * len(pairs) == len(units)
        _chr[m] = (L, reals, pairs)
    return _chr[m]


def lam(m):
    out = 1
    for q, e in factorint(m).items():
        out = math.lcm(out, (1 if e == 1 else 2 if e == 2 else 2 ** (e - 2)) if q == 2 else q ** (e - 1) * (q - 1))
    return out


# ---------------------------------------------------------------- one row

def divisor_data(N, m):
    f = factorint(N)
    assert math.prod(q ** e for q, e in f.items()) == N and all(isprime(q) for q in f), N
    odd = sorted((q, e) for q, e in f.items() if q != 2)
    divs = [1]
    for q, e in odd:
        divs = [d * q ** j for d in divs for j in range(e + 1)]
    tau = math.prod(e + 1 for _, e in odd)
    assert len(divs) == tau == len(set(divs)) and all(N % d == 0 for d in divs)
    W1 = sum(1 for d in divs if d % m == m - 1)
    dp = {1 % m: 1}
    for q, e in odd:
        nd = {}
        for r, c in dp.items():
            x = r
            for j in range(e + 1):
                nd[x] = nd.get(x, 0) + c
                x = x * q % m
        dp = nd
    W2 = dp.get(m - 1, 0)
    assert W1 == W2, (N, m, W1, W2)
    return f, odd, divs, tau, W1


def char_eta(odd, m, tau, W, WG2, G2size, cross_check=False):
    """eta, number of non-constant conjugate pairs; asserts orthogonality and the G^2 identity (exact)."""
    L, reals, pairs = charinfo(m)
    phi = len(group(m)['units'])
    real_sum = 0
    for phase, sgn, is_principal in reals:
        F = F_poly(phase, odd, m, L)
        v = reduce_fast(F, L)
        assert len(v) == 1
        if is_principal:
            assert v[0] == tau
        real_sum += sgn * v[0]
        if cross_check:
            assert v == reduce_mod_cyclo(F, L)
    ell_sum, nonconst, Ctot = 0, 0, {}
    for ph1, ph2, sgn in pairs:
        s = add_polys(F_poly(ph1, odd, m, L), F_poly(ph2, odd, m, L))
        red = [sgn * k for k in reduce_fast(s, L)]
        if cross_check:
            assert red == [sgn * k for k in reduce_mod_cyclo(s, L)]
        if len(red) == 1:
            ell_sum += red[0]
        else:
            ell_sum += -sum(abs(k) for k in red)
            nonconst += 1
        Ctot = add_polys(Ctot, s, sgn)
    Cred = reduce_fast(Ctot, L) if Ctot else [0]
    assert len(Cred) == 1
    C_num = Cred[0]
    assert real_sum + C_num == phi * W                       # orthogonality
    assert phi * WG2 == real_sum * G2size                   # R_i = W_{G^2}/|G^2| = (1/phi) sum_{chi^2=1} chi(-1) F_chi
    eta = Fr(C_num - ell_sum, phi)
    assert eta >= 0
    return eta, nonconst, L


def row(p, a, u, B, rng, full_check=False):
    k = a * u
    m = 4 * k
    N = a * p + u
    in_budget = k <= B
    f, odd, divs, tau, W = divisor_data(N, m)
    out = {'p': p, 'a': a, 'u': u, 'k': k, 'm': m, 'in_budget': in_budget, 'inE': k in E, 'W': W}
    if not in_budget:
        assert k in E
        return out
    g = group(m)
    WG2 = sum(1 for d in divs if d % m in g['negG2'])
    counts = [0] * g['ncosets']
    for d in divs:
        counts[g['coset'][d % m]] += 1
    qs = [q % m for q, _ in odd]
    K = closure(m, qs, [1])
    KH = closure(m, qs, g['G2'])
    assert ((m - 1) in KH) == (WG2 > 0)
    m1K = (m - 1) in K
    if W > 0:
        assert m1K and WG2 > 0
    cls = '' if W > 0 else ('quad' if WG2 == 0 else ('group' if not m1K else 'expo'))
    G2size = len(g['G2'])
    sampleT = rng.random() < 0.01 or full_check
    if sampleT:
        assert divisors_direct(N) == divisors_from_factors(f)
    L = lam(m)
    if k in E:
        assert L in SMALL_L
    eta, nonconst, Lc = char_eta(odd, m, tau, W, WG2, G2size, cross_check=sampleT)
    assert Lc == L
    if k in E:
        assert eta == 0 and nonconst == 0                    # the paper claim, checked on every E row of the box
    if nonconst == 0:
        assert eta == 0
    c = max(Fr(0), W - eta)
    out.update({'L': L, 'tau': tau, 'W_G2': WG2, 'G2size': G2size, 'm1K': m1K, 'cls': cls,
                'uniform': len(set(counts)) == 1, 'eta': eta, 'c': c, 'nonconst': nonconst, 'sampleT': sampleT})
    return out


def prime_rows(p, rng, full_check=False):
    B = (p - 1).bit_length()
    pairs = [(a, u) for a in range(1, B + 1) for u in range(1, B + 1) if a * u <= B and math.gcd(a, u) == 1]
    extra = sorted((a, k // a) for k in E if k > B for a in range(1, k + 1) if k % a == 0 and math.gcd(a, k // a) == 1)
    rows = [row(p, a, u, B, rng, full_check) for a, u in pairs]
    rows += [row(p, a, u, B, rng, full_check) for a, u in extra]
    return B, rows


def summarize(p, B, rows):
    box = [r for r in rows if r['in_budget']]
    s = {'p': p, 'B': B, 'pairs': len(box), 'extraE': len(rows) - len(box),
         'Wtot': sum(r['W'] for r in box), 'S': sum(1 for r in box if r['W'] > 0),
         'Cert': sum((r['c'] for r in box), Fr(0)),
         'Cert_E': sum(r['W'] for r in box if r['inE']),
         'Cert_0notE': sum(r['W'] for r in box if not r['inE'] and r['eta'] == 0),
         'Cert_plus': sum((r['c'] for r in box if r['eta'] > 0), Fr(0)),
         'E_all': sum(r['W'] for r in rows if r['inE']),
         'E_pairs_box': sum(1 for r in box if r['inE']),
         'outsideE_eta0_rows': sum(1 for r in box if not r['inE'] and r['eta'] == 0),
         'eta_pos_rows': sum(1 for r in box if r['eta'] > 0),
         'eta_sum': sum((r['eta'] for r in box), Fr(0)),
         'pos_rows': sum(1 for r in box if r['c'] > 0),
         'pos_eta0_rows': sum(1 for r in box if r['c'] > 0 and r['eta'] == 0),
         'W_pos_c0_rows': [[r['a'], r['u']] for r in box if r['W'] > 0 and r['c'] == 0],
         'W_gt_eta_pos': [{'a': r['a'], 'u': r['u'], 'L': r['L'], 'W': r['W'], 'eta': r['eta'], 'nonconst': r['nonconst']}
                          for r in box if r['eta'] > 0 and r['W'] > r['eta']],
         'fail': sum(1 for r in box if r['W'] == 0),
         'quad': sum(1 for r in box if r['cls'] == 'quad'),
         'group': sum(1 for r in box if r['cls'] == 'group'),
         'expo': sum(1 for r in box if r['cls'] == 'expo'),
         'second_uniform': sum(1 for r in box if r['cls'] in ('group', 'expo') and r['uniform']),
         'second_oddH': sum(1 for r in box if r['cls'] in ('group', 'expo') and r['G2size'] % 2 == 1),
         'second_oddH_expo': sum(1 for r in box if r['cls'] == 'expo' and r['G2size'] % 2 == 1),
         'second_G2size': dict(Counter(str(r['G2size']) for r in box if r['cls'] in ('group', 'expo'))),
         'sampleT': sum(1 for r in box if r['sampleT'])}
    s['Cert0'] = s['Cert_E'] + s['Cert_0notE']
    s['E_log'] = s['Cert_E']
    assert s['Cert'] == s['Cert0'] + s['Cert_plus']
    assert s['fail'] == s['quad'] + s['group'] + s['expo']
    assert s['second_oddH'] == s['second_oddH_expo']           # >>62: |G^2| odd => no group obstruction
    return s


# ---------------------------------------------------------------- modes

def selftest():
    rng = random.Random(SEED)
    ref65 = json.load(open(os.path.join(ROOT, 'Research', 'ExactFamily65Check20261016', 'check65_results.json'), encoding='utf-8'))
    ref62 = json.load(open(os.path.join(ROOT, 'Research', 'Slack60G2Check20261015', 'check62_results.json'), encoding='utf-8'))
    boxes = {b['p']: b for b in ref65['part2_log_rows']['boxes']}
    perE = {b['p']: b for b in ref65['part3_E_rows']['per_prime']}
    t0 = time.time()
    sums, bad = [], []
    for full in (False, True):
        for p in DISCOVERY:
            t = time.time()
            B, rows = prime_rows(p, rng, full_check=full)
            s = summarize(p, B, rows)
            b = boxes[p]
            pairs = [('pairs', s['pairs'], b['pairs']), ('Wtot', s['Wtot'], b['Wtot']), ('Cert', s['Cert'], Fr(b['Cert'])),
                     ('Cert0', s['Cert0'], Fr(b['Cert0'])), ('Certplus', s['Cert_plus'], Fr(b['Certplus'])),
                     ('E_Wtot', s['Cert_E'], b['E_Wtot']), ('E_pairs', s['E_pairs_box'], b['E_pairs']),
                     ('outside_E_zero_eta_rows', s['outsideE_eta0_rows'], b['outside_E_zero_eta_rows']),
                     ('outside_E_zero_eta_Wtot', s['Cert_0notE'], b['outside_E_zero_eta_Wtot']),
                     ('positive_rows', s['pos_rows'], b['positive_rows']), ('positive_eta_zero', s['pos_eta0_rows'], b['positive_eta_zero']),
                     ('eta_positive_rows', s['eta_pos_rows'], b['eta_positive_rows']), ('eta_sum', s['eta_sum'], Fr(b['eta_sum'])),
                     ('E_all', s['E_all'], perE[p]['W_E_total'])]
            bad += [(full, p, n, str(x), str(y)) for n, x, y in pairs if x != y]
            if not full:
                sums.append(s)
            print(f'selftest full_check={full} p={p} B={B} rows={len(rows)} {time.time() - t:.2f}s', flush=True)
    tot = [sum(s[k] for s in sums) for k in ('fail', 'quad', 'group', 'expo')]
    uni = sum(s['second_uniform'] for s in sums)
    res = {'mode': 'selftest', 'inputs': DISCOVERY, 'box_mismatches': bad,
           'failure_split': tot, 'failure_split_ref62': [ref62['failure_split'][k] for k in ('total', 'quadratic_W_H_zero', 'W_H_pos_minus1_not_in_K', 'W_H_pos_minus1_in_K_W_zero')],
           'second_stage_uniform': uni, 'second_stage_uniform_ref62': ref62['uniform_coset_lemma']['uniform_failures'],
           'seconds': round(time.time() - t0, 1)}
    res['all_match'] = not bad and res['failure_split'] == res['failure_split_ref62'] and uni == res['second_stage_uniform_ref62']
    sys.stdout.reconfigure(encoding='utf-8')
    print(json.dumps(res, ensure_ascii=False, indent=1))
    json.dump(res, open(os.path.join(HERE, 'selftest_results.json'), 'w', encoding='utf-8', newline='\n'), ensure_ascii=False, indent=1)
    assert res['all_match']


def git(*args):
    return subprocess.run(['git', '-C', ROOT, *args], capture_output=True, text=True)


def cohort():
    t0 = time.time()
    # G0: inputs and code state before any computation
    primes = [p for p in range(16_000_001, 19_000_000) if p % 840 in (1, 121, 169, 289, 361, 529) and isprime(p)]
    h = hashlib.sha256(''.join(f'{p}\n' for p in primes).encode()).hexdigest()
    assert h == PRIME_SHA and len(primes) == 5625, h
    head = git('rev-parse', 'HEAD').stdout.strip()
    rel = os.path.relpath(HERE, ROOT).replace(os.sep, '/')
    tracked = git('ls-files', rel).stdout.split()
    for fn in ('gtest.py', 'chars20261015.py', 'contract.md'):
        assert f'{rel}/{fn}' in tracked, fn
    assert git('diff', '--quiet', 'HEAD', '--', rel).returncode == 0, 'working tree differs from HEAD'
    kref = {r['p']: r for r in json.load(gzip.open(os.path.join(ROOT, 'Research', 'MainTermRatio20261012', 'mainterm_test_rows.json.gz'), 'rt', encoding='utf-8'))}
    assert sorted(kref) == primes
    print(f'G0 ok: primes {len(primes)} sha256 {h}, HEAD {head}', flush=True)
    rng = random.Random(SEED)
    sums = []
    fp = gzip.open(os.path.join(HERE, 'gtest_rows.csv.gz'), 'wt', encoding='utf-8', newline='\n')
    fp.write('p,a,u,in_budget,inE,L,W,W_G2,G2size,m1K,cls,uniform,eta,c,nonconst\n')
    for i, p in enumerate(primes):
        B, rows = prime_rows(p, rng)
        s = summarize(p, B, rows)
        assert s['Wtot'] == kref[p]['W'] and s['S'] == kref[p]['S'] and s['pairs'] == kref[p]['pairs'], (p, s['Wtot'], s['S'])
        sums.append(s)
        for r in rows:
            if r['in_budget']:
                fp.write(f"{p},{r['a']},{r['u']},1,{int(r['inE'])},{r['L']},{r['W']},{r['W_G2']},{r['G2size']},{int(r['m1K'])},"
                         f"{r['cls']},{int(r['uniform'])},{r['eta']},{r['c']},{r['nonconst']}\n")
            else:
                fp.write(f"{p},{r['a']},{r['u']},0,1,,{r['W']},,,,,,,,\n")
        if (i + 1) % 250 == 0:
            print(f'{i + 1}/{len(primes)} primes, {time.time() - t0:.0f}s', flush=True)
        assert time.time() - t0 < 3 * 3600
    fp.close()
    with gzip.open(os.path.join(HERE, 'gtest_primes.jsonl.gz'), 'wt', encoding='utf-8', newline='\n') as g:
        for s in sums:
            g.write(json.dumps(ser(s), ensure_ascii=False) + '\n')
    res = judge(sums)
    res['G0'] = {'primes': len(primes), 'prime_list_sha256': h, 'HEAD': head, 'K_test_Wtot_S_pairs_match': True,
                 'box_rows': sum(s['pairs'] for s in sums), 'extra_E_rows': sum(s['extraE'] for s in sums),
                 'E_rows_in_box_eta_checked': sum(s['E_pairs_box'] for s in sums),
                 'sample_rows_trial_division_and_sympy': sum(s['sampleT'] for s in sums)}
    res['seconds'] = round(time.time() - t0, 1)
    json.dump(ser(res), open(os.path.join(HERE, 'gtest_results.json'), 'w', encoding='utf-8', newline='\n'), ensure_ascii=False, indent=1)
    sys.stdout.reconfigure(encoding='utf-8')
    print(json.dumps(ser(res), ensure_ascii=False, indent=1))


def quantiles(xs, qs=(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 1)):
    ys = sorted(xs)
    return {str(q): ys[min(len(ys) - 1, int(q * (len(ys) - 1) + 0.5))] for q in qs}


def judge(sums):
    n = len(sums)
    cert0 = [s for s in sums if s['Cert'] == 0]
    g2 = [s for s in sums if s['Cert0'] == 0 and s['Cert_plus'] > 0]
    helog = [s for s in sums if s['E_log'] == 0]
    he = [s for s in sums if s['E_all'] == 0]
    wzero = [s for s in sums if s['Wtot'] == 0]
    fail = sum(s['fail'] for s in sums)
    quad = sum(s['quad'] for s in sums)
    grp = sum(s['group'] for s in sums)
    expo = sum(s['expo'] for s in sums)
    second = grp + expo
    brief = lambda s: {k: s[k] for k in ('p', 'Wtot', 'S', 'Cert', 'Cert_E', 'Cert_0notE', 'Cert_plus', 'E_log', 'E_all', 'quad', 'group', 'expo')}
    G1 = 1 - len(cert0) / n
    G2 = len(g2) / n
    qshare = quad / fail
    gshare = grp / second
    res = {
        'G1': {'fraction_Cert_pos': G1, 'pass_ge_0.99': G1 >= 0.99, 'Cert_zero_primes': [brief(s) for s in cert0]},
        'G2': {'fraction_Cert0_zero_Certplus_pos': G2, 'pass_le_0.01': G2 <= 0.01, 'primes': [brief(s) for s in g2]},
        'GE': {'H_Elog_exceptions': len(helog), 'H_E_exceptions': len(he),
               'H_Elog_exception_primes': [brief(s) for s in helog], 'H_E_exception_primes': [brief(s) for s in he]},
        'G3': {'failures': fail, 'quadratic': quad, 'group': grp, 'exponent': expo,
               'quadratic_share': qshare, 'pass_quadratic_in_[0.75,0.87]': 0.75 <= qshare <= 0.87,
               'group_share_of_second_stage': gshare, 'pass_group_le_0.10': gshare <= 0.10,
               'second_stage_uniform_share': sum(s['second_uniform'] for s in sums) / second,
               'second_stage_oddH': sum(s['second_oddH'] for s in sums),
               'second_stage_oddH_all_exponent': all(s['second_oddH'] == s['second_oddH_expo'] for s in sums)},
        'branches': {'D1_Cert_zero': len(cert0) > 0, 'D2_Wtot_zero': len(wzero) > 0, 'D3_Cert0_zero_Certplus_pos': len(g2) > 0,
                     'D4_H_Elog_exception': len(helog) > 0, 'D5_H_E_exception': len(he) > 0,
                     'D6_no_exception': not (cert0 or wzero or g2 or helog or he),
                     'D7_G3_out_of_band': not (0.75 <= qshare <= 0.87 and gshare <= 0.10)},
    }
    Wt = sum(s['Wtot'] for s in sums)
    Ct = sum((s['Cert'] for s in sums), Fr(0))
    order = sorted(sums, key=lambda s: (s['E_log'], s['E_all'], s['p']))
    res['G4'] = {
        'Wtot_sum': Wt, 'Cert_sum': Ct, 'Cert_over_Wtot_pooled': float(Ct / Wt),
        'Cert_over_Wtot_per_prime': quantiles([float(s['Cert'] / s['Wtot']) for s in sums]),
        'Cert_E_sum': sum(s['Cert_E'] for s in sums), 'Cert_0notE_sum': sum(s['Cert_0notE'] for s in sums),
        'Cert_plus_sum': sum((s['Cert_plus'] for s in sums), Fr(0)),
        'Cert_E_over_Cert_pooled': float(sum(s['Cert_E'] for s in sums) / Ct),
        'E_log_quantiles': quantiles([s['E_log'] for s in sums]), 'E_all_quantiles': quantiles([s['E_all'] for s in sums]),
        'E_log_hist_small': {str(v): sum(1 for s in sums if s['E_log'] == v) for v in range(0, 6)},
        'lowest_20_E_log': [brief(s) for s in order[:20]],
        'W_gt_eta_pos_rows': [dict(p=s['p'], **x) for s in sums for x in s['W_gt_eta_pos']],
        'W_pos_c0_rows': sum(len(s['W_pos_c0_rows']) for s in sums),
        'W_pos_c0_by_pair': dict(Counter(f"{a},{u}" for s in sums for a, u in s['W_pos_c0_rows']).most_common(15)),
        'eta_pos_rows': sum(s['eta_pos_rows'] for s in sums),
        'second_uniform_per_box': dict(Counter(str(s['second_uniform']) for s in sums)),
        'second_G2size': dict(sum((Counter(s['second_G2size']) for s in sums), Counter())),
    }
    return res


if __name__ == '__main__':
    {'selftest': selftest, 'cohort': cohort}[sys.argv[1]]()
