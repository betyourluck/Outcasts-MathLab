"""Independent check of >>76 (Outcasts 80c02e64, 2026-10-18): four-layer classification of the 3,391 square-source
candidates with the generated group K and the square-class group K G^2, the Type I / II quadratic equivalence,
the S+ / S- split with Delta(S+), the unique fixed point for sigma = 1, t = 1, and the edge constraint.
Contract: contract.md (fixed before computation). Own implementation; the reproduction code in >>76 is not executed.
Rows (candidates, factors, S, both selectors) come from our own Research/Reply73Check20261018/check73_rows.json.gz.
"""
import sys, json, gzip, math, time
from pathlib import Path
from collections import defaultdict
from sympy import factorint, primitive_root, totient

sys.stdout.reconfigure(encoding='utf-8')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
T0 = time.monotonic()
rows = json.load(gzip.open(ROOT / 'Research/Reply73Check20261018/check73_rows.json.gz', 'rt', encoding='utf-8'))
assert len(rows) == 3391


def checktime():
    assert time.monotonic() - T0 < 20 * 60, 'contract time limit'


def leg(a, p):
    v = pow(a % p, (p - 1) // 2, p)
    return -1 if v == p - 1 else v


def kron_minus_sigma_q(sq, d):
    """psi(d) = (-sq | d) (Kronecker) for odd d > 0 coprime to sq: Jacobi symbol via factorization of d."""
    out = 1
    for r, e in factorint(d).items():
        if e % 2:
            out *= leg(-sq, r) if r != 2 else (1 if (-sq) % 8 in (1, 7) else -1)
    return out


_phi_cache = {}


def phi_fac(c):
    if c not in _phi_cache:
        ph = int(totient(c))
        _phi_cache[c] = (ph, factorint(ph), sorted(factorint(c)))
    return _phi_cache[c]


def in_K_sylow(gens, targets, c):
    """membership in <gens> inside (Z/c)^x by projection to each Sylow subgroup of G (G abelian)."""
    ph, phf, _ = phi_fac(c)
    res = [True] * len(targets)
    for ell, k in phf.items():
        cof = ph // ell ** k
        H = {1}
        for g in gens:
            pg = pow(g, cof, c)
            if pg == 1:
                continue
            # closure of H under multiplication by pg (pg has order a power of ell, at most ell^k)
            new = set(H)
            acc = pg
            while acc not in H:
                new |= {(h * acc) % c for h in H}
                acc = acc * pg % c
            H = new
        for i, tgt in enumerate(targets):
            if res[i] and pow(tgt, cof, c) not in H:
                res[i] = False
    return res


def closure(gens, c):
    K = {1}
    todo = [1]
    for v in todo:
        for g in gens:
            w = v * g % c
            if w not in K:
                K.add(w)
                todo.append(w)
    return K


def in_KG2(gens, targets, c):
    primes_c = phi_fac(c)[2]
    def mask(u):
        return sum(1 << j for j, ell in enumerate(primes_c) if leg(u, ell) == -1)
    span = {0}
    for g in gens:
        mg = mask(g)
        span |= {s ^ mg for s in span}
    return [mask(t) in span for t in targets]


layers = defaultdict(lambda: defaultdict(int))
ident = {'quad_equiv': 0, 'II_group_iff_minus1': 0, 'I_group_iff_minus4': 0, 'sigma3_t1': 0, 'split': 0,
         'reduction': 0, 'edge': 0, 'omega1_form': 0, 'direct_closure_rows': 0, 'fixed_point_rows': 0}
equality_rows = []
sigma3_quadratic_rows = []
sigma3_typeI_residual = []
omega = defaultdict(int)
log = []


def say(s):
    log.append(s)
    print(s, flush=True)


for r in rows:
    checktime()
    p, q, t, sg, c, x = r['p'], r['q'], r['t'], r['sigma'], r['c'], r['x']
    fs = [(int(a), int(e)) for a, e in r['factors']]
    assert c % 2 == 1 and 4 * x == p + c and math.gcd(x, c) == 1
    gens = [a % c for a, _ in fs]
    tI, tII = (-p * x) % c, (-x) % c
    kK = in_K_sylow(gens, [tII, tI, (c - 1) % c, (-4) % c], c)
    kQ = in_KG2(gens, [tII, tI, (c - 1) % c], c)
    ph = phi_fac(c)[0]
    if ph <= 30000:
        K = closure(gens, c)
        assert [tII in K, tI in K, (c - 1) in K, (-4) % c in K] == kK, (p, q, t)
        G2 = {u * u % c for u in range(1, c) if math.gcd(u, c) == 1}
        primes_c = phi_fac(c)[2]
        reps = {}
        for k in K:                                   # one representative of K per class of G / G^2
            reps.setdefault(tuple(leg(k, ell) for ell in primes_c), k)
        KG2 = {a * b % c for a in reps.values() for b in G2}
        assert [tII in KG2, tI in KG2, (c - 1) in KG2] == kQ, (p, q, t)
        ident['direct_closure_rows'] += 1
    else:
        K = None
    # (i) quadratic equivalence; (ii) group conditions
    assert kQ[0] == kQ[1] == kQ[2]
    ident['quad_equiv'] += 1
    assert kK[0] == kK[2]
    ident['II_group_iff_minus1'] += 1
    assert kK[1] == kK[3]
    ident['I_group_iff_minus4'] += 1
    hitII, hitI = r['hit_II'], r['hit_I']
    assert (not hitII) or kK[0]
    assert (not hitI) or kK[1]
    # (iii) sigma = 3, t = 1
    if sg == 3 and t == 1:
        assert x % 3 == 1 and leg(x, q) == -1
        has2 = any(a % 3 == 2 for a, _ in fs)
        assert kK[0] == kQ[0] == has2
        ident['sigma3_t1'] += 1
        if not kQ[0]:
            sigma3_quadratic_rows.append((p, q))
            assert all(a % 3 == 1 for a, _ in fs)
        elif not kK[1]:
            sigma3_typeI_residual.append((p, q))
    # S+ / S- split
    divs = [(1, 1)]
    for a, e in fs:
        ca = leg(a, p)
        divs = [(d * a ** k, cd * ca ** k) for d, cd in divs for k in range(e + 1)]
    assert sorted({d % c for d, _ in divs}) == r['S']
    Sp = {d % c for d, cd in divs if cd == 1}
    Sm = {d % c for d, cd in divs if cd == -1}
    assert not (Sp & Sm) and len(Sp) == len(Sm) and len(Sp) + len(Sm) == len(r['S'])   # residues mod c; divisors may collide mod c
    assert Sm == {x * pow(s, -1, c) % c for s in Sp}
    sq = sg * q
    # psi(d) = chi(d) on divisors; psi(-1) = -1 since -sigma q < 0 and -sigma q = 1 mod 4; psi(-p) = +1
    for d, cd in divs[:12]:
        assert kron_minus_sigma_q(sq, d) == cd
    assert (-sq) % 4 == 1
    assert leg(-sq % p, p) == -1          # psi(p) = (-sq | p) = -1  (reciprocity, -sq = 1 mod 4)
    ident['split'] += 1
    # reduction of the two non-crossing conditions
    negSm = {(-s) % c for s in Sm}
    SpSp = {a * b % c for a in Sp for b in Sp}
    DeltaSp = {a * pow(b, -1, c) % c for a in Sp for b in Sp}
    negpSp = {(-p * s) % c for s in Sp}
    negpSm = {(-p * s) % c for s in Sm}
    assert hitII == bool(Sp & negSm) == (tII in SpSp)
    assert hitI == bool(Sp & negpSp) == (((-p) % c) in DeltaSp)
    assert not (Sp & {(-s) % c for s in Sp}) and not (Sm & negSm)
    assert not (Sp & negpSm) and not (Sm & negpSp)
    assert len(Sp & negpSp) == len(Sm & negpSm)
    ident['reduction'] += 1
    # Omega_-, edge constraint
    nr = [(a, e) for a, e in fs if leg(a, p) == -1]
    om = sum(e for _, e in nr)
    assert om % 2 == 1
    omega[om] += 1
    for a, e in nr:
        Spr = {d % c for d, cd in divs if cd == 1 and (x // d) % a == 0}
        assert {a * s % c for s in Spr} <= Sm
        if not hitII:
            assert ((-a) % c) not in {u * pow(v, -1, c) % c for u in Sp for v in Spr}
    ident['edge'] += 1
    if om == 1:
        a0 = nr[0][0]
        m_plus = x // a0
        assert Sp == {d % c for d in (dd for dd, _ in divs) if m_plus % d == 0} and Sm == {a0 * s % c for s in Sp}
        assert hitII == (((-a0) % c) in DeltaSp) and hitI == (((-p) % c) in DeltaSp)
        ident['omega1_form'] += 1
    # sigma = 1, t = 1: unique fixed point and the cardinality bound
    if sg == 1 and t == 1:
        assert K is not None
        Kp = {u for u in K if leg(u, q) == 1}       # ker psi inside K (psi = (-q|.) = (.|q) for q = 3 mod 4)
        assert len(Kp) * 2 == len(K) and len(Kp) % 2 == 1
        assert Sp <= Kp and Sm <= K - Kp
        fixed = [u for u in Kp if u * u % q == (-x) % q]
        assert len(fixed) == 1
        u0 = fixed[0]
        if not hitII:
            assert u0 not in Sp and u0 not in negSm and not (Sp & negSm)
            assert 2 * len(Sp) <= len(Kp) - 1
            if 2 * len(Sp) == len(Kp) - 1:
                equality_rows.append((p, q, len(Kp), len(Sp)))
        ident['fixed_point_rows'] += 1
    # layer table
    L = layers[(sg, 't=1' if t == 1 else 't>1')]
    L['n'] += 1
    L['II_ok'] += hitII
    L['I_ok'] += hitI
    if not hitII:
        L['II_quad' if not kQ[0] else ('II_group' if not kK[0] else 'II_exp')] += 1
    if not hitI:
        L['I_quad' if not kQ[1] else ('I_group' if not kK[1] else 'I_exp')] += 1
    if not hitII and not hitI:
        L['both'] += 1
        if kK[0] and kK[1]:
            L['both_pure_exp'] += 1

cols = ['n', 'II_ok', 'I_ok', 'II_quad', 'I_quad', 'II_group', 'II_exp', 'I_group', 'I_exp', 'both', 'both_pure_exp']
table = {f'{sg},{tt}': {k: layers[(sg, tt)][k] for k in cols} for sg, tt in [(1, 't=1'), (3, 't=1'), (1, 't>1'), (3, 't>1')]}
total = {k: sum(v[k] for v in table.values()) for k in cols}
expected = {'1,t=1': [230, 51, 61, 0, 0, 0, 179, 6, 163, 156, 150], '3,t=1': [217, 15, 22, 34, 34, 0, 168, 2, 159, 187, 151],
            '1,t>1': [1522, 6, 22, 330, 330, 19, 1167, 95, 1075, 1497, 1060], '3,t>1': [1422, 4, 8, 553, 553, 21, 844, 39, 822, 1411, 805]}
mism = []
for key, exp in expected.items():
    mine = [table[key][k] for k in cols]
    if mine != exp:
        mism.append((key, cols, mine, exp))
exp_total = [3391, 76, 113, 917, 917, 40, 2358, 142, 2219, 3251, 2166]
if [total[k] for k in cols] != exp_total:
    mism.append(('total', cols, [total[k] for k in cols], exp_total))
if dict(omega) != {1: 3314, 3: 77}:
    mism.append(('omega', dict(omega), {1: 3314, 3: 77}))
if sorted(equality_rows) != [(2521, 11, 5, 2), (18490201, 211, 3, 1)]:
    mism.append(('equality_rows', equality_rows))
if len(sigma3_quadratic_rows) != 34:
    mism.append(('sigma3_quadratic', len(sigma3_quadratic_rows)))
if sorted(sigma3_typeI_residual) != [(66529, 241), (66529, 337)]:
    mism.append(('sigma3_typeI_residual', sigma3_typeI_residual))
# pre-registered consistency: both_pure_exp = both - (rows where some branch is quad or group obstructed)
for key in table:
    L = layers[(int(key[0]), key[2:])]
    assert L['both_pure_exp'] <= L['both']
report = {'contract': 'contract.md', 'source': '>>76 (post 359c9fa2), no PR materials', 'rows': len(rows),
          'table': table, 'total': total, 'omega_minus': dict(omega), 'equality_rows': equality_rows,
          'sigma3_t1_quadratic_rows': len(sigma3_quadratic_rows), 'sigma3_t1_typeI_residual_rows': sigma3_typeI_residual,
          'identities_checked_rows': ident, 'mismatches': mism, 'mismatch_count': len(mism),
          'seconds': round(time.monotonic() - T0, 1)}
(HERE / 'check76_results.json').write_text(json.dumps(report, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
for key in table:
    say(f"{key}: " + ' '.join(f"{k}={table[key][k]}" for k in cols))
say('total: ' + ' '.join(f"{k}={total[k]}" for k in cols))
say(f"omega_- {dict(omega)}; equality rows {equality_rows}; sigma3 quad rows {len(sigma3_quadratic_rows)}; sigma3 Type I residual {sigma3_typeI_residual}")
say(f"identities: {ident}")
say(f"done: mismatches {len(mism)} {mism}, {report['seconds']} s")
