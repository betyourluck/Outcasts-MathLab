"""Child-level check on the three root-failing inputs (>>97 section 3, >>98) -- see contract.md.

Own implementation. The square-source arithmetic is the one of Bfs12Inputs20261024/bfs12.py (not imported):
vertex q, sigma(q) = 1 (q = 3 mod 4) / 3, A_q = sigma q, odd tau with A_q tau^2 < 3p, c = A_q tau^2, x = (p + c)/4,
success <=> some d | x^2 has d = -x (Type II) or d = -p x (Type I) mod c. No peer code is imported or executed.

Usage:
  python -I childobs.py selftest   # 1201, 100687129 (outside the three) and the known >>88 stages of the 2521 root
  python -I childobs.py cohort     # 2521, 8615161, 84525841; refuses to run unless contract.md and childobs.py
                                   # are committed, unchanged and contained in the pushed fork branch
stdout is the log (redirect it to selftest.log / run.log).
"""
import sys, os, json, math, time, subprocess
from fractions import Fraction
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
BFS12 = os.path.join(ROOT, 'Research', 'Bfs12Inputs20261024')
INPUTS = [2521, 8615161, 84525841]
JUDGED = [8615161, 84525841]
SELFTEST = [1201, 100687129]
CLASSES = (1, 121, 169, 289, 361, 529)
FORK_REF = 'refs/remotes/fork/luna/p001-failure-atlas'
CAP_INPUT = 15 * 60
CAP_TOTAL = 30 * 60
T0 = time.monotonic()
LOG = []
# numbers fixed in contract.md (counted from Legendre symbols only, before any success was looked at)
KNOWN = {2521: {'ell': 11, 'n_root': 13, 'r_B': 28, 'W': 3, 'W_children': 2},
         8615161: {'ell': 17, 'n_root': 356, 'r_B': 713, 'W': 59, 'W_children': 35},
         84525841: {'ell': 17, 'n_root': 1115, 'r_B': 2231, 'W': 180, 'W_children': 87}}
# >>88 (reproduced in Reply8688Check20261023): stages of the 13 root sources of 2521. >>88 fixes the exponent of the
# nonresidue prime r of x (Omega_- = 1) per branch: Type I uses d with v_r(d) = 1, Type II d with v_r(d) in {0, 2}.
ROOT2521 = {'base_both': [1, 5, 9, 15, 19], 'square_both': [3, 7, 11, 13, 17, 23, 25],
            'mixed': [(21, 'square', 'base')]}  # (tau, Type I, Type II)
# the stages of contract.md (any d | x^2 at each modulus), found by the first self-test (selftest-first.log):
# at tau = 7, 17 Type I hits mod 11 only with v_r(d) = 1 and mod tau^2 only with v_r(d) = 0, so the stage is joint
ROOT2521_DIRECT = {'base_both': [1, 5, 9, 15, 19], 'square_both': [3, 11, 13, 23, 25],
                   'mixed': [(7, 'joint', 'square'), (17, 'joint', 'square'), (21, 'square', 'base')]}
STAGES = ('success', 'base', 'square', 'joint')
OBST = ('real', 'group', 'exponent', 'composite')


def log(*a):
    s = ' '.join(str(x) for x in a)
    LOG.append(s)
    print(s, flush=True)


class OutOfBudget(Exception):
    pass


def check_budget(t_input):
    now = time.monotonic()
    if now - t_input > CAP_INPUT or now - T0 > CAP_TOTAL:
        raise OutOfBudget


# ---------------------------------------------------------------- arithmetic

def trial_factor(n):
    f, q = {}, 2
    while q * q <= n:
        while n % q == 0:
            f[q] = f.get(q, 0) + 1
            n //= q
        q += 1 if q == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def factor(n):
    """G0-4: trial division and sympy.factorint must agree."""
    a = trial_factor(n)
    b = {int(k): int(v) for k, v in sympy.factorint(n).items()}
    assert a == b, (n, a, b)
    assert math.prod(r ** e for r, e in a.items()) == n
    return dict(sorted(a.items()))


def leg(a, r):
    """Legendre symbol (a|r) by Euler's criterion, r odd prime, r does not divide a."""
    k = pow(a % r, (r - 1) // 2, r)
    assert k in (1, r - 1), (a, r)
    return 1 if k == 1 else -1


def chi(a, p):
    return leg(a, p)


def divisors(f):
    ds = [1]
    for r, e in f.items():
        ds = [d * r ** j for d in ds for j in range(e + 1)]
    return sorted(ds)


def val(n, r):
    n, e = abs(n), 0
    while n % r == 0:
        n //= r
        e += 1
    return e


def least_nonresidue(p):
    for q in sympy.primerange(2, p):
        if chi(q, p) == -1:
            return q
    raise AssertionError(p)


def sigma_of(q):
    return 1 if q % 4 == 3 else 3


def allowed_taus(p, q):
    A = sigma_of(q) * q
    taus, t = [], 1
    while A * t * t < 3 * p:
        taus.append(t)
        t += 2
    return taus


def sqrt_mod(a, r):
    """Tonelli-Shanks, r odd prime, a a nonzero square mod r."""
    a %= r
    assert a and leg(a, r) == 1
    q, s = r - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while leg(z, r) == 1:
        z += 1
    m, c, t, x = s, pow(z, q, r), pow(a, q, r), pow(a, (q + 1) // 2, r)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % r
            i += 1
        b = pow(c, 1 << (m - i - 1), r)
        m, c, t, x = i, b * b % r, t * b * b % r, x * b % r
    assert x * x % r == a
    return x


def beta(p, A, r):
    """The odd root of p + A tau^2 = 0 mod r (1 <= beta <= r - 2), or None if there is none."""
    a = (-p * pow(A, -1, r)) % r
    if a == 0 or leg(a, r) == -1:
        return None
    b = sqrt_mod(a, r)
    bb = b if b % 2 == 1 else r - b
    assert bb % 2 == 1 and 1 <= bb <= r - 2 and (p + A * bb * bb) % r == 0
    return bb


_PM1 = {}


def order(a, m):
    """Multiplicative order of a mod the prime m."""
    if m not in _PM1:
        _PM1[m] = list(factor(m - 1))
    a %= m
    assert a
    o = m - 1
    for l in _PM1[m]:
        while o % l == 0 and pow(a, o // l, m) == 1:
            o //= l
    return o


def in_span(vectors, v):
    """GF(2): is the bit vector v in the span of vectors?"""
    basis = {}
    for w in vectors:
        while w:
            h = w.bit_length() - 1
            if h not in basis:
                basis[h] = w
                break
            w ^= basis[h]
    while v:
        h = v.bit_length() - 1
        if h not in basis:
            return False
        v ^= basis[h]
    return True


def kind_at(t, gens, m, mprimes):
    """Classify the target t at the modulus m (odd, > 1) regardless of hits:
    'real'  a real character psi mod m (product of Legendre symbols over a subset of the primes of m)
            has psi(g) = 1 for all generators and psi(t) = -1;
    'group' m prime, t not in <gens>;  'in_group' m prime, t in <gens>;  'composite' m composite, not 'real'."""
    vec = lambda a: sum(1 << i for i, r in enumerate(mprimes) if leg(a, r) == -1)
    if not in_span([vec(g) for g in gens], vec(t)):
        return 'real'
    if mprimes == [m]:
        L = 1
        for g in gens:
            L = math.lcm(L, order(g, m))
        return 'group' if L % order(t, m) else 'in_group'
    return 'composite'


# ---------------------------------------------------------------- one square source

def analyze(p, q, tau, full=True):
    """Success of the source (q, tau), two routes; with full=True also the stage of each branch
    (moduli m1 = sigma q, m2 = tau^2, m3 = c) and the obstruction type at the failing modulus."""
    sig = sigma_of(q)
    A = sig * q
    c = A * tau * tau
    assert tau % 2 == 1 and c < 3 * p and (p + c) % 4 == 0
    x = (p + c) // 4
    assert 0 < x < p and math.gcd(c, x) == 1  # G0-4
    f = factor(x)
    gens = list(f)
    # route 1: the residues of the divisors of x^2 mod c, prime power by prime power
    R = {1 % c}
    for g, e in f.items():
        pw = [pow(g, j, c) for j in range(2 * e + 1)]
        R = {v * w % c for v in R for w in pw}
    # route 2: enumerate the divisors of x^2
    ds = divisors({g: 2 * e for g, e in f.items()})
    assert {d % c for d in ds} == R
    tgt = {'II': (-x) % c, 'I': (-p * x) % c}
    hit = {br: tgt[br] in R for br in tgt}
    cert = None
    for br in ('II', 'I'):  # the order of bfs12.py: Type II first, smallest d
        if hit[br]:
            d = next(d for d in ds if d % c == tgt[br])
            if br == 'II':
                y, z = p * (d + x) // c, p * (x * x // d + x) // c
                assert p * (d + x) % c == 0 and p * (x * x // d + x) % c == 0
            else:
                y, z = (d + p * x) // c, (p * p * (x * x // d) + p * x) // c
                assert (d + p * x) % c == 0 and (p * p * (x * x // d) + p * x) % c == 0
            assert y > 0 and z > 0 and 4 * x * y * z == p * (x * y + x * z + y * z)
            cert = {'branch': br, 'd': d, 'c': c, 'x': x}
            break
    out = {'q': q, 'tau': tau, 'c': c, 'x': x, 'f': f, 'hit': hit, 'cert': cert,
           'nonres': [r for r in f if chi(r, p) == -1]}
    if not full:
        return out
    m1, m2 = A, tau * tau
    tf = list(factor(tau)) if tau > 1 else []
    # tau = 1: m2 = 1 (always a hit) and m3 = m1, so the joint stage cannot occur
    primes = {m1: [3, q] if sig == 3 else [q], m2: tf, c: sorted(set(([3] if sig == 3 else []) + [q] + tf))}
    R1, R2 = {v % m1 for v in R}, {v % m2 for v in R}
    st = {}
    for br, t in tgt.items():
        h = {m1: t % m1 in R1, m2: t % m2 in R2, c: t in R}
        # G0-6: a hit modulo c implies hits modulo m1 and m2; real / group obstructions mean no hit
        if h[c]:
            assert h[m1] and h[m2]
        kinds = {}
        for m in (m1, m2, c):
            if m == 1:
                assert h[m]
                continue
            k = kinds[m] = kind_at(t, gens, m, primes[m])
            if k in ('real', 'group'):
                assert not h[m], (p, q, tau, br, m, k)
        # G0-5: sigma = 1, Type II has no real-character obstruction modulo q
        if sig == 1 and br == 'II':
            assert leg(-x, q) == 1 and kinds[m1] != 'real', (p, q, tau)
        if h[c]:
            st[br] = ('success', None)
        else:
            m = m1 if not h[m1] else (m2 if not h[m2] else c)
            stage = 'base' if m == m1 else ('square' if m == m2 else 'joint')
            k = kinds[m]
            st[br] = (stage, 'exponent' if k == 'in_group' else k)
    out['stage'] = st
    return out


# ---------------------------------------------------------------- one input

def stage_counts(rows):
    cnt = {s: 0 for s in STAGES}
    obs = {o: 0 for o in OBST}
    cross = {s: {o: 0 for o in OBST} for s in STAGES[1:]}
    for s in rows:
        for br in ('II', 'I'):
            stage, k = s['stage'][br]
            cnt[stage] += 1
            if k:
                obs[k] += 1
                cross[stage][k] += 1
    return {'stages': cnt, 'obstructions': obs, 'stage_x_obstruction': cross, 'units': 2 * len(rows)}


def code(sk):
    stage, k = sk
    return stage if k is None else f'{stage}:{k}'


def run_input(p, known, tsv, t_input):
    """known: the bfs12 JSON entry of p (regression, G0-2)."""
    assert sympy.isprime(p) and factor(p) == {p: 1}
    assert p % 840 in CLASSES and p % 8 == 1 and chi(2, p) == 1
    ell = least_nonresidue(p)
    sig = sigma_of(ell)
    A = sig * ell
    taus = allowed_taus(p, ell)
    assert (known['ell'], known['n_root_tau']) == (ell, len(taus))
    res = {'p': p, 'class': p % 840, 'ell': ell, 'sigma': sig, 'A': A, 'n_root_tau': len(taus)}

    # the root, all sources
    root = []
    for t in taus:
        check_budget(t_input)
        root.append(analyze(p, ell, t))
    S0 = [s['tau'] for s in root if s['cert']]
    assert len(S0) == known['sweep']['S0']
    children = sorted({r for s in root for r in s['nonres']})
    from_failing = sorted({r for s in root if not s['cert'] for r in s['nonres']})
    for s in root:
        tsv.append([p, 'root', ell, s['tau'], code(s['stage']['II']), code(s['stage']['I'])])
    res['root'] = {'S0': S0, 'rows': [[s['tau'], s['x'], code(s['stage']['II']), code(s['stage']['I'])] for s in root],
                   'summary': stage_counts([s for s in root if not s['cert']])}
    if p == 2521:
        res['root']['regression_88'] = regression_2521(root)

    # children: beta, reach, all sources
    cs = known['sweep']['child_status']
    assert sorted(int(k) for k in cs) == from_failing  # G0-2
    rB = max(r for r in range(2, 10 * math.isqrt(3 * p) + 10) if A * (r - 2) ** 2 < 3 * p)
    rows, fail_src = [], []
    for r in children:
        check_budget(t_input)
        reach = [(s['tau'], s['f'][r]) for s in root if r in s['f']]
        b = beta(p, A, r)
        # G0-3
        assert chi(r, p) == -1 and A % r and leg(-A, r) == -1 and b is not None
        assert b == min(t for t, _ in reach) and A * b * b < 3 * p
        src = []
        for t in allowed_taus(p, r):
            check_budget(t_input)
            src.append(analyze(p, r, t))
        succ = [s for s in src if s['cert']]
        row = {'r': r, 'sigma': sigma_of(r), 'chi_p': chi(r, p), 'sym_minusA': leg(-A, r), 'beta': b,
               'reach': reach, 'in_B': A * (r - 2) ** 2 < 3 * p, 'n_tau': len(src),
               'status': 'success' if succ else 'fail', 'n_success_sources': len(succ)}
        if succ:
            s = succ[0]
            row['first'] = {'tau': s['tau'], 'branch': s['cert']['branch'], 'd': s['cert']['d'], 'c': s['c'], 'x': s['x']}
        else:
            row['summary'] = stage_counts(src)
            fail_src += src
            for s in src:
                tsv.append([p, 'child', r, s['tau'], code(s['stage']['II']), code(s['stage']['I'])])
        if str(r) in cs:  # G0-2
            k = cs[str(r)]
            assert k['status'] == row['status'], (p, r)
            if succ:
                assert (k['tau'], k['branch']) == (row['first']['tau'], row['first']['branch']), (p, r)
            else:
                assert k['n_tau'] == row['n_tau'], (p, r)
        rows.append(row)
    res['children'] = rows
    res['children_summary'] = {
        'n': len(rows), 'n_success': sum(1 for x in rows if x['status'] == 'success'),
        'success': [x['r'] for x in rows if x['status'] == 'success'],
        'success_in_B': [x['r'] for x in rows if x['status'] == 'success' and x['in_B']],
        'n_in_B': sum(1 for x in rows if x['in_B']), 'max_child': max(children) if children else None}
    res['failing_children'] = stage_counts(fail_src)
    res['failing_children']['n_children'] = sum(1 for x in rows if x['status'] == 'fail')
    res['failing_children']['n_sources'] = len(fail_src)

    # the window W(p): child <=> (-A/r) = -1 (G0-7); success scanned up to the first success
    status = {x['r']: x for x in rows}
    W = [r for r in sympy.primerange(3, rB + 1) if chi(r, p) == -1 and r != ell and A % r]
    wrows = []
    for r in W:
        check_budget(t_input)
        is_child = leg(-A, r) == -1
        assert is_child == (r in status), (p, r)
        if is_child:
            x = status[r]
            wrows.append({'r': r, 'child': True, 'status': x['status'], 'scanned': x['n_tau']})
            continue
        st, n = 'fail', 0
        for t in allowed_taus(p, r):
            check_budget(t_input)
            n += 1
            s = analyze(p, r, t, full=False)
            if s['cert']:
                st = 'success'
                break
        wrows.append({'r': r, 'child': False, 'status': st, 'scanned': n, 'n_tau': len(allowed_taus(p, r))})
    tab = {'child_success': 0, 'child_fail': 0, 'nonchild_success': 0, 'nonchild_fail': 0}
    for w in wrows:
        tab[('child_' if w['child'] else 'nonchild_') + w['status']] += 1
    res['W'] = {'r_B': rB, 'size': len(W), 'n_children': sum(1 for w in wrows if w['child']), 'table': tab,
                'fisher_two_sided': float(fisher_two_sided(tab)), 'rows': wrows}
    if p in KNOWN:
        k = KNOWN[p]
        assert (ell, len(taus), rB, len(W), res['W']['n_children']) == (k['ell'], k['n_root'], k['r_B'], k['W'], k['W_children'])
        assert not S0 and children == from_failing  # G0-1
        if p == 2521:
            assert W == [17, 19, 23]

    # 2521 (and the small self-test input): the reverse direction of G0-3 over all primes below p
    if p < 10 ** 4:
        rev = [r for r in sympy.primerange(3, p) if chi(r, p) == -1 and A % r and leg(-A, r) == -1
               and A * beta(p, A, r) ** 2 < 3 * p]
        assert rev == children, (p, rev, children)
        res['reverse_G0_3'] = {'primes_checked': len(list(sympy.primerange(3, p))), 'equal': True}
    return res


def stage_restricted(p, s):
    """The >>88 form: Omega_- = 1, x = r m; Type I uses only d | x^2 with v_r(d) = 1, Type II only v_r(d) in {0, 2}."""
    q, tau, c, x, f = s['q'], s['tau'], s['c'], s['x'], s['f']
    assert len(s['nonres']) == 1 and f[s['nonres'][0]] == 1
    r = s['nonres'][0]
    ds = divisors({g: 2 * e for g, e in f.items()})
    m1, m2 = sigma_of(q) * q, tau * tau
    out = {}
    for br, t, J in (('II', -x, (0, 2)), ('I', -p * x, (1,))):
        h = {m: any((d - t) % m == 0 for d in ds if val(d, r) in J) for m in (m1, m2, c)}
        out[br] = 'success' if h[c] else ('base' if not h[m1] else ('square' if not h[m2] else 'joint'))
    return out


def table_2521(st):
    """st: tau -> {'I': stage, 'II': stage}"""
    return {'base_both': [t for t in st if st[t]['I'] == 'base' and st[t]['II'] == 'base'],
            'square_both': [t for t in st if st[t]['I'] == 'square' and st[t]['II'] == 'square'],
            'mixed': [(t, st[t]['I'], st[t]['II']) for t in st if st[t]['I'] != st[t]['II']]}


def regression_2521(root):
    """The restricted form must give the >>88 table exactly; the stages of contract.md may differ from it only
    by restricted 'square' -> direct 'joint' (a hit mod tau^2 that uses another exponent of r)."""
    direct = {s['tau']: {br: s['stage'][br][0] for br in ('I', 'II')} for s in root}
    restr = {s['tau']: stage_restricted(2521, s) for s in root}
    got_d, got_r = table_2521(direct), table_2521(restr)
    assert got_r == ROOT2521, got_r
    assert got_d == ROOT2521_DIRECT, got_d
    diff = [(t, br, restr[t][br], direct[t][br]) for t in direct for br in ('I', 'II') if direct[t][br] != restr[t][br]]
    assert all((a, b) == ('square', 'joint') for _, _, a, b in diff), diff
    return {'ok': True, 'restricted_equals_88': got_r, 'direct': got_d, 'differences': diff}


# ---------------------------------------------------------------- statistics and judgments

def fisher_two_sided(tab):
    """Two-sided Fisher exact test on [[child_success, child_fail], [nonchild_success, nonchild_fail]]:
    the sum of the hypergeometric probabilities not larger than the observed one (exact rationals)."""
    a, b, c, d = tab['child_success'], tab['child_fail'], tab['nonchild_success'], tab['nonchild_fail']
    row, col, n = a + b, a + c, a + b + c + d
    if n == 0:
        return Fraction(1)
    w = lambda k: math.comb(col, k) * math.comb(n - col, row - k)
    obs = w(a)
    return Fraction(sum(w(k) for k in range(max(0, row + col - n), min(row, col) + 1) if w(k) <= obs), math.comb(n, row))


def judge(res):
    by = {r['p']: r for r in res}
    if any(p not in by or by[p].get('undetermined') for p in JUDGED):
        return {'status': 'undetermined (an input exceeded the time limit)'}
    tab = {k: sum(by[p]['W']['table'][k] for p in JUDGED) for k in by[JUDGED[0]]['W']['table']}
    pf = fisher_two_sided(tab)
    rc = Fraction(tab['child_success'], tab['child_success'] + tab['child_fail'])
    rn = Fraction(tab['nonchild_success'], tab['nonchild_success'] + tab['nonchild_fail'])
    j1 = 'J1-a' if pf >= Fraction(5, 100) else ('J1-b' if rc > rn else 'J1-c')
    st = {s: sum(by[p]['failing_children']['stages'][s] for p in JUDGED) for s in STAGES}
    ob = {o: sum(by[p]['failing_children']['obstructions'][o] for p in JUDGED) for o in OBST}
    units = sum(by[p]['failing_children']['units'] for p in JUDGED)
    assert st['success'] == 0 and sum(st.values()) == units == sum(ob.values())
    base, joint, real = Fraction(st['base'], units), Fraction(st['joint'], units), Fraction(ob['real'], units)
    j2 = 'J2-a' if base >= Fraction(60, 100) and joint <= Fraction(5, 100) else 'J2-b'
    j3 = 'J3-a' if real >= Fraction(50, 100) else 'J3-b'
    return {'J1': {'pooled_table': tab, 'fisher_two_sided_p': float(pf), 'success_rate_children': float(rc),
                   'success_rate_nonchildren': float(rn), 'branch': j1},
            'J2': {'units': units, 'stages': st, 'share_base': float(base), 'share_square': float(Fraction(st['square'], units)),
                   'share_joint': float(joint), 'branch': j2},
            'J3': {'units': units, 'obstructions': ob, 'share_real': float(real), 'branch': j3}}


def unit_checks():
    """Small known cases for the statistics and the obstruction classifier."""
    # Fisher's tea tasting [[3,1],[1,3]]: 34/70; [[1,9],[11,3]]: 5412/1961256 (R: 0.002759)
    assert fisher_two_sided({'child_success': 3, 'child_fail': 1, 'nonchild_success': 1, 'nonchild_fail': 3}) == Fraction(34, 70)
    assert fisher_two_sided({'child_success': 1, 'child_fail': 9, 'nonchild_success': 11, 'nonchild_fail': 3}) == Fraction(5412, 1961256)
    assert kind_at(3, [2], 7, [7]) == 'real'           # <2> = {1,2,4} = squares mod 7, 3 is not a square
    assert kind_at(4, [3], 13, [13]) == 'group'        # <3> = {1,3,9}, 4 is a square of order 6
    assert kind_at(9, [3], 13, [13]) == 'in_group'
    assert kind_at(2, [4], 15, [3, 5]) == 'real'       # (4|3) = (4|5) = 1, (2|5) = -1
    assert kind_at(8, [2], 15, [3, 5]) == 'composite'  # 8 and 2 have the same symbols (-1, -1) mod 3 and mod 5
    assert in_span([0b011, 0b110], 0b101) and not in_span([0b011, 0b110], 0b100)
    for r in (7, 13, 101, 7963009):
        for a in range(1, 30):
            if a % r and leg(a, r) == 1:
                s = sqrt_mod(a, r)
                assert s * s % r == a % r
    return {'ok': True}


# ---------------------------------------------------------------- main

def git(*args):
    return subprocess.run(['git', '-C', ROOT, *args], capture_output=True, text=True)


def run_all(ps, known_file):
    known = {k['p']: k for k in json.load(open(os.path.join(BFS12, known_file), encoding='utf-8'))['inputs']}
    res, tsv = [], []
    for p in ps:
        t_input = time.monotonic()
        try:
            r = run_input(p, known[p], tsv, t_input)
        except OutOfBudget:
            r = {'p': p, 'undetermined': True}
            log(f'p={p}: over the time limit, undetermined')
        r['seconds'] = round(time.monotonic() - t_input, 2)
        if not r.get('undetermined'):
            cs, fc, w = r['children_summary'], r['failing_children'], r['W']
            log(f"p={p} ell={r['ell']} A={r['A']}: root {r['n_root_tau']} sources, S0={len(r['root']['S0'])}, "
                f"root stages {r['root']['summary']['stages']}, obstructions {r['root']['summary']['obstructions']}")
            log(f"p={p} children {cs['n']} (success {cs['n_success']}: {cs['success']}, in B {cs['success_in_B']}), "
                f"failing children {fc['n_children']} with {fc['n_sources']} sources")
            log(f"p={p} failing children: stages {fc['stages']}, obstructions {fc['obstructions']}")
            log(f"p={p} W: r_B={w['r_B']}, |W|={w['size']}, children {w['n_children']}, table {w['table']}, "
                f"Fisher two-sided {w['fisher_two_sided']:.4g}; {r['seconds']}s")
        res.append(r)
    return res, tsv


def write_tsv(fn, rows):
    with open(os.path.join(HERE, fn), 'w', encoding='utf-8', newline='\n') as fp:
        fp.write('p\trole\tq\ttau\ttype_II\ttype_I\n')
        for row in rows:
            fp.write('\t'.join(str(v) for v in row) + '\n')


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else ''
    out = {'contract': 'contract.md', 'mode': mode, 'caps_seconds': {'input': CAP_INPUT, 'total': CAP_TOTAL},
           'sympy': sympy.__version__, 'python': sys.version.split()[0]}
    out['unit_checks'] = unit_checks()
    log('unit checks ok')
    if mode == 'selftest':
        # the known >>88 stages of the 2521 root (13 sources; no child, no window, nothing used by J1-J3)
        root = [analyze(2521, 11, t) for t in allowed_taus(2521, 11)]
        out['regression_2521_root'] = regression_2521(root)
        log(f"regression 2521 root stages vs >>88: OK {out['regression_2521_root']}")
        res, tsv = run_all(SELFTEST, 'selftest_results.json')
        out['inputs'] = res
        write_tsv('selftest_source_stages.tsv', tsv)
        fn = 'selftest_results.json'
    elif mode == 'cohort':
        head = git('rev-parse', 'HEAD').stdout.strip()
        rel = os.path.relpath(HERE, ROOT).replace(os.sep, '/')
        tracked = git('ls-files', rel).stdout.split()
        for f in ('childobs.py', 'contract.md'):
            assert f'{rel}/{f}' in tracked, f
        assert git('diff', '--quiet', 'HEAD', '--', f'{rel}/childobs.py', f'{rel}/contract.md').returncode == 0, 'working tree differs from HEAD'
        assert git('merge-base', '--is-ancestor', 'HEAD', FORK_REF).returncode == 0, 'HEAD is not pushed to the fork'
        log(f'G0 ok: HEAD {head} (pushed), inputs {INPUTS}')
        out['head'] = head
        res, tsv = run_all(INPUTS, 'bfs12_results.json')
        out['inputs'] = res
        out['judgments'] = judge(res)
        log('judgments:', json.dumps(out['judgments']))
        write_tsv('source_stages.tsv', tsv)
        fn = 'childobs_results.json'
    else:
        raise SystemExit(__doc__)
    out['seconds'] = round(time.monotonic() - T0, 2)
    out['log'] = LOG
    log(f'done: {out["seconds"]}s')
    with open(os.path.join(HERE, fn), 'w', encoding='utf-8', newline='\n') as fp:
        fp.write(json.dumps(out, ensure_ascii=False, indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
