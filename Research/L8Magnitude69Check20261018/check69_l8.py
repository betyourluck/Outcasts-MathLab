"""Independent check of >>69 (Outcasts 80c02e64, 2026-10-17): the phase-discarding L = 8 lower bound c_mag,
the 42 pilot rows, the 4,709 development rows of the G test, and the Qmag counterexample 670849 (1,8).
Contract: contract.md (fixed before computation). Own implementation on top of
Research/CertificateRate20261016/chars20261015.py (same author, 2026-10-15). The participant's scripts are
neither imported nor executed; their JSON / CSV at 2eb1317 are read as data from the directory given as argv[1].

usage: python check69_l8.py <dir with the 2eb1317 files>
"""
import sys, json, csv, gzip, math, cmath, time, hashlib, importlib.util
from fractions import Fraction as Fr
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
THEIRS = Path(sys.argv[1])
T0 = time.monotonic()
LIMIT = 30 * 60

spec = importlib.util.spec_from_file_location('chars20261015', ROOT / 'Research/CertificateRate20261016/chars20261015.py')
ch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ch)
from sympy import totient, isprime, factorint

SIX = {1, 121, 169, 289, 361, 529}
INPUTS = [66529, 345601, 670849, 1740481, 3678481, 5843041, 8615161, 18490201, 53722321, 84525841, 38296729]
P = 38296729
SQ2 = math.sqrt(2)


def checktime():
    assert time.monotonic() - T0 < LIMIT, 'contract time limit'


def legendre(n, p):
    v = pow(n % p, (p - 1) // 2, p)
    return -1 if v == p - 1 else v


def carmichael(n):
    out = 1
    for q, e in factorint(n).items():
        if q == 2:
            lam = 1 if e == 1 else (2 if e == 2 else 2 ** (e - 2))
        else:
            lam = (q - 1) * q ** (e - 1)
        out = math.lcm(out, lam)
    return out


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
            d += w            # 2 w cos(pi/4) = w sqrt2
        else:
            d -= w            # t in (3, 5)
    return g, d


def mul_sqrt2(x, y):
    return x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def magnitude(p, a, u, odd, m, W, R):
    """c_mag from the factorization and the modulus only (plus R, which comes from the real characters)."""
    L, units, chars, bykey = ch.characters(m)
    assert L == 8
    phi = int(totient(m))
    zeta = cmath.exp(1j * math.pi / 4)
    seen = set()
    T4 = 0
    terms = []
    sum_re8 = 0          # sum over order-8 representatives of sign * 2 c0  (integer part of sign*(F + Fbar))
    sum_irr8 = 0         # sum of sign * (c1 - c3)  (sqrt2 coefficient; must cancel to 0 in total)
    for c in chars:
        if c['real'] or c['param'] in seen:
            continue
        conj = bykey[tuple((-c['phase'][r]) % L for r in units)]
        assert conj['param'] != c['param'] and conj['sign'] == c['sign']
        seen.update((c['param'], conj['param']))
        d = L // math.gcd(L, *c['phase'].values())
        F = ch.F_poly(c['phase'], odd, m, L)
        red = pad4(ch.reduce_mod_cyclo(F, L))
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
            Fbar = ch.F_poly(conj['phase'], odd, m, L)
            prod = pad4(ch.reduce_mod_cyclo(mul_dict(F, Fbar, L), L))
            G, D = prod[0], prod[1]
            assert prod == [G, D, 0, -D], (p, a, u, prod)
            # second route: the local cosine formula multiplied in Z[sqrt2]
            gd = (1, 0)
            for q, e in odd:
                gd = mul_sqrt2(gd, local_norm_cos(c['phase'][q % m], e))
            assert gd == (G, D), (p, a, u, gd, (G, D))
            # third route: floating point
            assert abs(abs(Ff) ** 2 - (G + D * SQ2)) < 1e-6
            U = Fr(G) + (Fr(3, 2) if D >= 0 else Fr(7, 5)) * D
            assert U >= G + D * SQ2 - 1e-9
            rho = ceil_sqrt_1024(U)
            assert rho * rho >= U
            assert float(rho) >= abs(Ff) - 1e-9
            # pointwise bound used in the derivation: sign*(F + Fbar) = sign*2 Re F >= -2|F| >= -2 rho
            assert c['sign'] * 2 * Ff.real >= -2 * float(rho) - 1e-9
            sum_re8 += c['sign'] * 2 * red[0]
            sum_irr8 += c['sign'] * (red[1] - red[3])
            terms.append({'gamma': G, 'delta': D, 'U': U, 'rho': rho})
    assert len(seen) == sum(1 for c in chars if not c['real'])
    # exact orthogonality identity: phi W = phi R + T4 + sum_{order 8 reps} sign (F + Fbar)
    assert sum_irr8 == 0, (p, a, u, sum_irr8)
    assert phi * W == phi * R + T4 + sum_re8, (p, a, u, phi * W, phi * R, T4, sum_re8)
    raw = R + Fr(T4, phi) - Fr(2, phi) * sum((t['rho'] for t in terms), Fr(0))
    c_mag = max(Fr(0), raw)
    assert 0 <= c_mag <= W
    return {'T4': T4, 'terms': terms, 'raw_lower': raw, 'c_mag': c_mag}


def l8row(p, a, u, direct=True):
    checktime()
    r = ch.row(p, a, u, with_orbits=False, with_g2=False)
    N, m, phi = r['N'], r['m'], r['phi']
    f = {int(q): e for q, e in r['factors'].items()}
    odd = sorted((q, e) for q, e in f.items() if q != 2)
    L, units, chars, bykey = ch.characters(m)
    # divisors: generated from the factors (and, when direct, by trial division up to sqrt N)
    ds = ch.divisors_from_factors(f)
    if direct:
        assert ds == ch.divisors_direct(N)
    odd_divs = [d for d in ds if d % 2 == 1]
    # residue DP of the odd divisors mod m
    dp = {1: 1}
    for q, e in odd:
        new = {}
        for res, cnt in dp.items():
            for j in range(e + 1):
                k = (res * pow(q, j, m)) % m
                new[k] = new.get(k, 0) + cnt
        dp = new
    assert dp.get(m - 1, 0) == r['W'] == len(r['witnesses'])
    assert all((d + 1) % m == 0 for d in r['witnesses'])
    # R from the real characters only, and its identity with the square-coset count
    R_num = 0
    for c in chars:
        if c['real']:
            F = ch.F_poly(c['phase'], odd, m, L)
            val = sum(v * (1 if k == 0 else -1) for k, v in F.items() if v)
            assert all(k in (0, L // 2) for k, v in F.items() if v)
            R_num += c['sign'] * val
    R = Fr(R_num, phi)
    assert R == r['R'] if 'R' in r else True
    G2 = sorted({(x * x) % m for x in units})
    negG2 = {(-s) % m for s in G2}
    W_H = sum(cnt for res, cnt in dp.items() if res in negG2)
    assert R == Fr(W_H, len(G2)), (p, a, u, R, W_H, len(G2))
    assert R + r['C'] == r['W']
    c_old = max(Fr(0), r['W'] - r['eta'])
    assert r['b'] == min(r['V'], c_old)
    # nonresidue factors (actual multiplicity) and the generated group K
    nr = [(q, e) for q, e in sorted(f.items()) if legendre(q, p) == -1]
    nonres = sum(e for _, e in nr)
    K = {1}
    changed = True
    while changed:
        changed = False
        for q, _ in odd:
            new = {(k * q) % m for k in K}
            if not new <= K:
                K |= new
                changed = True
    out = {'p': p, 'a': a, 'u': u, 'N': N, 'm': m, 'L': L, 'factors': [[q, e] for q, e in sorted(f.items())],
           'W': r['W'], 'witnesses': r['witnesses'], 'R': R, 'C': r['C'], 'ell_over_phi': r['ell_sum'],
           'eta': r['eta'], 'old_c': c_old, 'b': r['b'], 'nonresidue_factors': [[q, e] for q, e in nr],
           'nonresidue_multiplicity': nonres, 'target_in_K': (m - 1) in K, 'W_H': W_H, 'H_card': len(G2),
           'u_legendre': legendre(u, p)}
    if L == 8:
        assert p % 840 in SIX and p > 120 and u % p != 0
        assert legendre(u, p) == 1
        assert nonres % 2 == 0 and (r['W'] == 0 or nonres >= 2)
        for s in r['witnesses']:
            assert legendre(s, p) == -1 and legendre(N // s, p) == -1
        mag = magnitude(p, a, u, odd, m, r['W'], R)
        assert nonres >= 2 or mag['c_mag'] == 0
        out.update(mag)
        out['rescues_old_zero'] = c_old == 0 and mag['c_mag'] > 0
        out['loses_old_positive'] = c_old > 0 and mag['c_mag'] == 0
        out['Qmag_counterexample'] = r['W'] > 0 and nonres >= 2 and (m - 1) in K and mag['c_mag'] == 0
    return out


def S(x):
    return ch.ser(x)


def same(mine, theirs, fields, tag, mism):
    for k in fields:
        a, b = mine[k], theirs['order8_magnitude_terms' if k == 'terms' else k]
        if k in ('R', 'C', 'ell_over_phi', 'eta', 'old_c', 'b', 'raw_lower', 'c_mag'):
            ok = Fr(str(b)) == a
        elif k == 'factors' or k == 'nonresidue_factors':
            ok = [list(x) for x in a] == [list(x) for x in b]
        elif k == 'terms':
            ok = sorted((t['gamma'], t['delta'], str(t['U']), str(t['rho'])) for t in a) == \
                 sorted((t['gamma'], t['delta'], str(Fr(t['norm_squared_upper'])), str(Fr(t['rho']))) for t in b)
        else:
            ok = a == b
        if not ok:
            mism.append((tag, k, S(a), b))


def load(name):
    return json.loads((THEIRS / name).read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


report = {'contract': 'contract.md', 'source_commit': '2eb1317', 'mismatches': []}
mism = report['mismatches']
log = []


def say(s):
    log.append(s)
    print(s, flush=True)


# ---------------------------------------------------------------- 0. lambda(4k) = 8 and the families
lam8 = sorted(k for k in range(1, 5001) if carmichael(4 * k) == 8)
assert lam8 == [8, 24, 40, 120], lam8
L8pairs = sorted((a, k // a) for k in lam8 for a in ch.divisors_direct(k) if math.gcd(a, k // a) == 1)
Ek = sorted(set(ch.divisors_direct(60)) | set(ch.divisors_direct(126)))
Epairs = sorted((a, k // a) for k in Ek for a in ch.divisors_direct(k) if math.gcd(a, k // a) == 1)
assert len(L8pairs) == 18 and len(Epairs) == 81 and not set(L8pairs) & set(Epairs)
assert len(set(L8pairs) | set(Epairs)) == 99
report['families'] = {'lambda8_products': lam8, 'L8_pairs': 18, 'E_pairs': 81, 'union': 99}
say('0. lambda(4k)=8 <=> k in {8,24,40,120} (k<=5000); 18 L8 pairs, 81 E pairs, disjoint, union 99')

# ---------------------------------------------------------------- 2a. 38,296,729: 73 box pairs and 81 E pairs
assert isprime(P) and P % 840 == 289
B, pairs = ch.box_pairs(P)
assert B == 26 and len(pairs) == 73
their73 = load('counterexample73-rows.json')
assert len(their73) == 73 and [(r['a'], r['u']) for r in their73] == pairs
rows73 = []
for (a, u), th in zip(pairs, their73):
    mine = l8row(P, a, u, direct=True)
    fields = ['p', 'a', 'u', 'N', 'm', 'L', 'factors', 'W', 'witnesses', 'R', 'C', 'ell_over_phi', 'eta', 'old_c', 'b',
              'nonresidue_factors', 'nonresidue_multiplicity', 'target_in_K', 'W_H', 'H_card', 'u_legendre']
    if mine['L'] == 8:
        fields += ['T4', 'terms', 'raw_lower', 'c_mag', 'rescues_old_zero', 'loses_old_positive', 'Qmag_counterexample']
    same(mine, th, fields, f'73:{a},{u}', mism)
    rows73.append(mine)
pos = [(r['a'], r['u']) for r in rows73 if r['W'] > 0]
assert pos == [(1, 25), (3, 8), (8, 1), (11, 1)]
Cert0 = sum(r['old_c'] for r in rows73 if r['eta'] == 0)
Certp = sum(r['old_c'] for r in rows73 if r['eta'] > 0)
Elog = sum(r['W'] for r in rows73 if (r['a'], r['u']) in set(Epairs))
assert Cert0 == 0 and Certp == 1 and Elog == 0
assert sum(r['W'] for r in rows73) == 6
report['p38296729_box'] = {'pairs': 73, 'W_positive': pos, 'Wtot': 6, 'Cert0': S(Cert0), 'Certplus': S(Certp), 'E_log': Elog,
                           'mismatch_rows_so_far': len(mism)}
say(f'2a. 38296729 box: 73 pairs, W>0 at {pos}, Cert0={Cert0} Cert+={Certp} E_log={Elog}; mismatches {len(mism)}')

their81 = load('counterexample81-E-rows.json')
assert len(their81) == 81
theirE = {(r['a'], r['u']): r for r in their81}
assert set(theirE) == set(Epairs)
layers = {}
Erows = []
for (a, u) in Epairs:
    checktime()
    N, m = a * P + u, 4 * a * u
    f = factorint(N)
    ds = ch.divisors_from_factors(f)
    assert ds == ch.divisors_direct(N)
    wit = [d for d in ds if (d + 1) % m == 0]
    lay = carmichael(m) // 2
    layers[lay] = layers.get(lay, 0) + 1
    mine = {'p': P, 'a': a, 'u': u, 'N': N, 'm': m, 'factors': [[int(q), e] for q, e in sorted(f.items())], 'W': len(wit), 'witnesses': wit}
    same(mine, theirE[(a, u)], ['p', 'a', 'u', 'N', 'm', 'factors', 'W', 'witnesses'], f'E:{a},{u}', mism)
    Erows.append(mine)
assert layers == {1: 9, 2: 36, 3: 36}
assert [(r['a'], r['u'], r['witnesses']) for r in Erows if r['W']] == [(14, 9, [544319])]
report['p38296729_E'] = {'pairs': 81, 'layers': layers, 'E_all': 1, 'success': [14, 9, 544319], 'mismatch_rows_so_far': len(mism)}
say(f'2a. 38296729 E: 81 pairs, layers {layers}, E_all=1 at (14,9) s=544319; mismatches {len(mism)}')

# ---------------------------------------------------------------- 2b. pilot: 42 L = 8 rows of the 11 fixed inputs
theirfixed = load('l8-fixed-rows.json')
fixed = []
for p in INPUTS:
    assert isprime(p) and p % 840 in SIX
    Bp, prs = ch.box_pairs(p)
    for (a, u) in prs:
        if carmichael(4 * a * u) == 8:
            assert (a, u) in set(L8pairs)
            fixed.append(l8row(p, a, u, direct=True))
assert len(fixed) == 42 == len(theirfixed)
assert [(r['p'], r['a'], r['u']) for r in fixed] == [(r['p'], r['a'], r['u']) for r in theirfixed]
FIELDS8 = ['p', 'a', 'u', 'N', 'm', 'L', 'factors', 'W', 'witnesses', 'R', 'C', 'ell_over_phi', 'eta', 'old_c', 'b',
           'nonresidue_factors', 'nonresidue_multiplicity', 'target_in_K', 'W_H', 'H_card', 'u_legendre',
           'T4', 'terms', 'raw_lower', 'c_mag', 'rescues_old_zero', 'loses_old_positive', 'Qmag_counterexample']
for mine, th in zip(fixed, theirfixed):
    same(mine, th, FIELDS8, f"fixed:{mine['p']},{mine['a']},{mine['u']}", mism)


def summary(rows):
    return {'rows': len(rows), 'W_positive': sum(r['W'] > 0 for r in rows), 'old_positive': sum(r['old_c'] > 0 for r in rows),
            'new_positive': sum(r['c_mag'] > 0 for r in rows), 'rescued_old_zero': sum(r['rescues_old_zero'] for r in rows),
            'lost_old_positive': sum(r['loses_old_positive'] for r in rows),
            'fewer_than_two_nonresidue': sum(r['nonresidue_multiplicity'] < 2 for r in rows),
            'Qmag_counterexamples': sum(r['Qmag_counterexample'] for r in rows)}


sfixed = summary(fixed)
assert sfixed == {'rows': 42, 'W_positive': 4, 'old_positive': 1, 'new_positive': 1, 'rescued_old_zero': 0, 'lost_old_positive': 0,
                  'fewer_than_two_nonresidue': 27, 'Qmag_counterexamples': 3}, sfixed
byk = {(r['p'], r['a'], r['u']): r for r in fixed}
r81 = byk[(P, 8, 1)]
assert r81['W'] == 2 and r81['old_c'] == 1 and r81['c_mag'] == Fr(10451, 8192) and r81['c_mag'] > 1
assert r81['nonresidue_factors'] == [[199, 1], [229, 1]] and r81['factors'] == [[3, 4], [83, 1], [199, 1], [229, 1]]
assert all(legendre(s, P) == -1 for s in r81['witnesses']) and r81['witnesses'] == [1791, 19007]
assert math.ceil(r81['c_mag']) == 2 <= r81['W']
r38 = byk[(P, 3, 8)]
assert r38['W'] == 2 and r38['old_c'] == 0 and r38['c_mag'] == 0
r670 = byk[(670849, 1, 8)]
assert r670['factors'] == [[3, 1], [11, 1], [29, 1], [701, 1]] and r670['W'] == 2 and r670['nonresidue_multiplicity'] == 2
assert r670['target_in_K'] and r670['c_mag'] == 0 and r670['raw_lower'] == Fr(-8937, 8192) and r670['Qmag_counterexample']
assert r670['R'] == 1 and r670['T4'] == 0
r665 = byk[(66529, 8, 1)]
assert r665['factors'] == [[3, 2], [13, 1], [4549, 1]] and r665['W'] == 0 and r665['nonresidue_multiplicity'] == 2 and r665['target_in_K']
qmag_rows = [(r['p'], r['a'], r['u']) for r in fixed if r['Qmag_counterexample']]
report['pilot'] = {'summary': sfixed, 'Qmag_counterexample_rows': qmag_rows,
                   'p38296729_8_1': S({k: r81[k] for k in ['W', 'old_c', 'c_mag', 'raw_lower', 'T4', 'R']}),
                   'p38296729_3_8': S({k: r38[k] for k in ['W', 'old_c', 'c_mag', 'raw_lower', 'T4', 'R']}),
                   'p670849_1_8': S({k: r670[k] for k in ['W', 'old_c', 'c_mag', 'raw_lower', 'T4', 'R', 'nonresidue_multiplicity', 'target_in_K']}),
                   'p66529_8_1': S({k: r665[k] for k in ['W', 'nonresidue_multiplicity', 'target_in_K', 'c_mag']}),
                   'mismatch_rows_so_far': len(mism)}
say(f'2b. pilot 42 rows: {sfixed}; Qmag counterexamples {qmag_rows}; mismatches {len(mism)}')

# ---------------------------------------------------------------- 2c. development rows from our own gtest_rows.csv.gz
GZ = ROOT / 'Research/CertificateRate20261016/gtest_rows.csv.gz'
gzsha = sha(GZ)
assert gzsha == '2851fb4611be7461ddc92d63d94ba9939b071c2fba236d9f353c01b44a43b8b1'
theirdev = list(csv.DictReader((THEIRS / 'l8-development-results.csv').open(encoding='utf-8', newline='')))
cohort = missed = 0
dev = []
with gzip.open(GZ, 'rt', newline='') as fh:
    for rec in csv.DictReader(fh):
        cohort += 1
        if int(rec['in_budget']) and int(rec['W']) > 0 and Fr(rec['c']) == 0:
            missed += 1
            if int(rec['L']) != 8:
                continue
            p, a, u = int(rec['p']), int(rec['a']), int(rec['u'])
            mine = l8row(p, a, u, direct=True)
            assert mine['W'] == int(rec['W']) and mine['eta'] == Fr(rec['eta']) and mine['old_c'] == 0
            assert mine['R'] * int(rec['G2size']) == int(rec['W_G2']) and mine['H_card'] == int(rec['G2size'])
            dev.append(mine)
            if len(dev) % 500 == 0:
                say(f'   dev rows {len(dev)} ({time.monotonic() - T0:.0f}s)')
assert cohort == 587689 and missed == 16809
assert len(dev) == 4709 == len(theirdev)
for mine, th in zip(dev, theirdev):
    assert (mine['p'], mine['a'], mine['u']) == (int(th['p']), int(th['a']), int(th['u']))
    th2 = dict(th)
    th2['rescues_old_zero'] = th['rescues_old_zero'] == 'True'
    th2['W'] = int(th['W']); th2['T4'] = int(th['T4']); th2['nonresidue_multiplicity'] = int(th['nonresidue_multiplicity'])
    same(mine, th2, ['W', 'R', 'eta', 'old_c', 'T4', 'raw_lower', 'c_mag', 'nonresidue_multiplicity', 'rescues_old_zero'],
         f"dev:{mine['p']},{mine['a']},{mine['u']}", mism)
sdev = summary(dev)
assert sdev == {'rows': 4709, 'W_positive': 4709, 'old_positive': 0, 'new_positive': 355, 'rescued_old_zero': 355,
                'lost_old_positive': 0, 'fewer_than_two_nonresidue': 0, 'Qmag_counterexamples': 4354}, sdev
first = next(r for r in dev if r['c_mag'] > 0)
assert (first['p'], first['a'], first['u']) == (16000609, 3, 8) and first['c_mag'] == Fr(363, 2048) and first['T4'] == 16
assert first['W'] == 1 and first['witnesses'] == [84959] and first['R'] == 1 and first['eta'] == 1
norms = sorted((t['gamma'], t['delta']) for t in first['terms'])
assert norms == [(0, 0)] * 4 + [(32, -16)] * 2 + [(32, 16)] * 2, norms
# pre-registered observation (contract section 3): c_mag > 0 implies T4 > 0 or R >= 2
obs = [(r['p'], r['a'], r['u'], S(r['R']), r['T4']) for r in fixed + dev if r['c_mag'] > 0 and not (r['T4'] > 0 or r['R'] >= 2)]
pos_rows = [r for r in dev if r['c_mag'] > 0]
report['development'] = {'gzip_sha256': gzsha, 'cohort_rows': cohort, 'old_misses_all_L': missed, 'summary': sdev,
                         'first_rescue': S({k: first[k] for k in ['p', 'a', 'u', 'N', 'factors', 'W', 'witnesses', 'R', 'eta', 'old_c', 'T4', 'raw_lower', 'c_mag']}),
                         'first_rescue_order8_norms': norms,
                         'positive_rows_T4_distribution': dict(sorted(((str(k), v) for k, v in
                             __import__('collections').Counter(r['T4'] for r in pos_rows).items()), key=lambda kv: int(kv[0]))),
                         'positive_rows_R_distribution': dict(sorted(((S(k), v) for k, v in
                             __import__('collections').Counter(r['R'] for r in pos_rows).items()), key=lambda kv: Fr(kv[0]))),
                         'positive_rows_by_pair': dict(sorted(((f"{a},{u}", v) for (a, u), v in
                             __import__('collections').Counter((r['a'], r['u']) for r in pos_rows).items()))),
                         'all_rows_by_pair': dict(sorted(((f"{a},{u}", v) for (a, u), v in
                             __import__('collections').Counter((r['a'], r['u']) for r in dev).items()))),
                         'c_mag_max': S(max(r['c_mag'] for r in dev)),
                         'observation_T4pos_or_R_ge_2_violations': obs,
                         'mismatch_rows_so_far': len(mism)}
say(f'2c. development: {sdev}; first rescue 16000609 (3,8) c_mag=363/2048; observation violations {len(obs)}; mismatches {len(mism)}')

# ---------------------------------------------------------------- 2d. soundness on every L = 8 row (asserted inside l8row); counts
all8 = fixed + dev
assert all(0 <= r['c_mag'] <= r['W'] for r in all8)
report['soundness'] = {'L8_rows_checked': len(all8), 'c_mag_le_W': True, 'u_legendre_all_1': all(r['u_legendre'] == 1 for r in all8),
                       'witness_symbols_minus1': True, 'even_nonresidue_multiplicity': True}

# ---------------------------------------------------------------- 2e. manifest and replay record
sums = {}
for line in (THEIRS / 'SHA256SUMS').read_text().splitlines():
    h, name = line.split()
    sums[name] = h
assert len(sums) == 15
manifest = {}
for name, h in sums.items():
    raw = (THEIRS / name).read_bytes()
    direct = hashlib.sha256(raw).hexdigest() == h
    crlf = (b'\r' not in raw) and hashlib.sha256(raw.replace(b'\n', b'\r\n')).hexdigest() == h
    manifest[name] = 'blob' if direct else ('blob-as-CRLF' if crlf else 'MISMATCH')
replay = load('portable-replay.json')
replay_ok = replay['final_source_sha256'] == sha(THEIRS / 'check-l8-norm-20261017.py') == sums['check-l8-norm-20261017.py']
theirres = load('l8-norm-results.json')
assert theirres['fixed_L8_pilot'] == sfixed and theirres['development_L8_old_misses'] == sdev
assert theirres['download_sha256'] == gzsha and theirres['development_old_misses_all_L'] == 16809
nblob = sum(v == 'blob' for v in manifest.values())
ncrlf = sum(v == 'blob-as-CRLF' for v in manifest.values())
assert nblob + ncrlf == 15, manifest
report['manifest'] = {'SHA256SUMS_matched_from_blob': nblob, 'matched_after_LF_to_CRLF': ncrlf, 'of': 15, 'detail': manifest,
                      'replay_source_hash_ok': replay_ok, 'their_results_json_summaries_match': True}
say(f"2e. SHA256SUMS {nblob}/15 from blobs, {ncrlf} after LF->CRLF ({[k for k, v in manifest.items() if v != 'blob']}), replay source hash ok={replay_ok}")

report['mismatch_count'] = len(mism)
report['seconds'] = round(time.monotonic() - T0, 1)
(HERE / 'check69_results.json').write_text(json.dumps(S(report), indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
(HERE / 'check69_rows.json').write_text(json.dumps(S({'rows73': rows73, 'E81': Erows, 'fixed42': fixed, 'dev4709': dev}), ensure_ascii=False) + '\n', encoding='utf-8')
say(f"done: mismatches {len(mism)}, {report['seconds']} s")
