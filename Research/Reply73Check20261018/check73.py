"""Independent check of >>73 (Outcasts 80c02e64, 2026-10-18): reciprocity on the square-source candidates of 12 fixed
primes, the three equivalent selectors for both branches, the sigma = 1, t = 1 group / quadratic classification
(230 / 179 / 6), the Type I group condition <=> 2 in K, and the finite parts of the closed-set corollary of >>70.
Contract: contract.md (fixed before computation). Own implementation; the participant's script is neither imported
nor executed; its summary.json at 4977cf0 is read as data from argv[1].

usage: python check73.py <dir with the 4977cf0 files>
"""
import sys, json, math, time, hashlib, itertools
from fractions import Fraction as Fr
from pathlib import Path
from sympy import factorint, isprime, primerange

sys.stdout.reconfigure(encoding='utf-8')
HERE = Path(__file__).resolve().parent
THEIRS = Path(sys.argv[1])
T0 = time.monotonic()
PS = [2521, 66529, 345601, 670849, 1740481, 3678481, 5843041, 8615161, 18490201, 38296729, 53722321, 84525841]
SIX = {1, 121, 169, 289, 361, 529}


def checktime():
    assert time.monotonic() - T0 < 600, 'contract time limit'


def leg(a, p):
    v = pow(a % p, (p - 1) // 2, p)
    return -1 if v == p - 1 else v


def kron(a, r):
    if r == 2:
        return 0 if a % 2 == 0 else (1 if a % 8 in (1, 7) else -1)
    return leg(a, r)


def trial_factor(n):
    out = []
    r = 2
    while r * r <= n:
        e = 0
        while n % r == 0:
            n //= r
            e += 1
        if e:
            out.append((r, e))
        r += 1 if r == 2 else 2
    if n > 1:
        out.append((n, 1))
    return out


def divisors(fs, scale=1):
    out = [1]
    for r, e in fs:
        out = [d * r ** k for d in out for k in range(scale * e + 1)]
    return sorted(out)


def closure(gens, q):
    K = {1}
    todo = [1]
    for v in todo:
        for g in gens:
            w = v * g % q
            if w not in K:
                K.add(w)
                todo.append(w)
    return K


def order(a, q):
    k, v = 1, a % q
    while v != 1:
        v = v * a % q
        k += 1
    return k


log = []


def say(s):
    log.append(s)
    print(s, flush=True)


report = {'contract': 'contract.md', 'source_commit': '4977cf0', 'mismatches': []}
mism = report['mismatches']

# ------------------------------------------------------------------ 1. candidates, reciprocity, selectors
rows = []
canon = []
n_checks = n_two = 0
for p in PS:
    assert isprime(p) and p % 840 in SIX
    for q in primerange(11, 400):
        if leg(q, p) != -1:
            continue
        sigma = 1 if q % 4 == 3 else 3
        for t in range(1, 16, 2):
            checktime()
            c = sigma * q * t * t
            if c >= 3 * p:
                continue
            assert (p + c) % 4 == 0
            x = (p + c) // 4
            assert 0 < x < p and leg(x, p) == -1 and math.gcd(x, p * c) == 1
            fs = sorted((int(r), e) for r, e in factorint(x).items())
            assert fs == trial_factor(x) and math.prod(r ** e for r, e in fs) == x
            for r, e in fs:
                # reciprocity chi(r) = (-sigma q | r), with the premises of the derivation
                if r == 2:
                    assert p % 8 == 1 and (t * t) % 8 == 1 and (-sigma * q) % 8 == 1
                    n_two += 1
                else:
                    assert (t * sigma * q) % r != 0 and (p + sigma * q * t * t) % r == 0
                assert leg(r, p) == kron(-sigma * q, r), (p, q, t, r)
                n_checks += 1
            ds = divisors(fs)
            ds2 = divisors(fs, 2)
            # d | x^2  <=>  d = a b with a, b | x
            assert set(ds2) == {a * b for a in ds for b in ds}
            S = {d % c for d in ds}
            D2 = {d % c for d in ds2}
            hitII = (-x) % c in D2
            hitI = (-p * x) % c in D2
            # bounded-exponent residue DP
            dp = {1}
            for r, e in fs:
                nxt = set()
                for v in dp:
                    w = v
                    for _ in range(2 * e + 1):
                        nxt.add(w)
                        w = w * r % c
                dp = nxt
            assert dp == D2
            # non-crossing form
            negS = {(-s) % c for s in S}
            negpS = {(-p * s) % c for s in S}
            assert hitII == bool(S & negS) and hitI == bool(S & negpS)
            rec = {'p': p, 'q': q, 't': t, 'sigma': sigma, 'c': c, 'x': x, 'factors': fs, 'S': sorted(S), 'hit_I': hitI, 'hit_II': hitII}
            if sigma == 1 and t == 1:
                gens = [r % q for r, _ in fs]
                K = closure(gens, q)
                assert len(K) == math.lcm(*(order(g, q) for g in gens))
                assert leg(x, q) == -1 and ((q - 1) // 2) % 2 == 1
                ox = order(x, q)
                assert ox % 2 == 0 and (ox // 2) % 2 == 1 and pow(x, ox // 2, q) == q - 1
                assert (q - 1) in K and (-x) % q in K
                assert p % q == 4 * x % q
                gI = (-p * x) % q in K
                assert gI == (4 in K) == (2 in K)
                assert ((q - 1) // len(K)) % 2 == 1
                G2 = {a * a % q for a in range(1, q)}
                KG2 = {a * b % q for a in K for b in G2}
                assert len(KG2) == q - 1 and (-x) % q in KG2 and (-p * x) % q in KG2
                assert len(S) % 2 == 0 and {x * pow(s, -1, q) % q for s in S} == S
                assert all(leg(s, q) == -leg(x * pow(s, -1, q) % q, q) for s in S)
                assert all(s in K for s in S)
                if not hitII:
                    assert not (S & negS) and 2 * len(S) <= len(K) - 2
                rec.update({'K': sorted(K), 'K_card': len(K), 'group_I': gI, 'group_II': True, 'quadratic_I': True, 'quadratic_II': True,
                            'two_in_K': 2 in K, 'x_even': x % 2 == 0})
                canon.append(rec)
            rows.append(rec)
assert len(rows) == 3391 and n_checks == 9441 and n_two == 1787
assert len(canon) == 230
failII = [r for r in canon if not r['hit_II']]
obsI = [r for r in canon if not r['group_I']]
assert len(failII) == 179 and all(r['group_II'] for r in failII)
assert len(obsI) == 6
assert [(r['p'], r['q'], r['x']) for r in obsI] == [(2521, 211, 683), (2521, 331, 713), (2521, 379, 725), (3678481, 331, 919703),
                                                   (18490201, 211, 4622603), (53722321, 211, 13430633)]
assert all(not r['hit_I'] for r in obsI)
assert all(not r['two_in_K'] and not r['x_even'] for r in obsI)
assert all(r['two_in_K'] for r in canon if r['x_even'])
their = json.loads((THEIRS / 'summary.json').read_text(encoding='utf-8'))
ts = their['summary']
mine_summary = {'inputs': PS, 'sources': len(rows), 'factorwise_checks': n_checks, 'factor_two_checks': n_two, 'reciprocity_violations': 0,
                'selector_intersections_and_DP_all_agree': True, 'canonical_sigma1_sources': len(canon),
                'canonical_TypeII_failure': len(failII), 'canonical_TypeII_group_failure': 0, 'canonical_TypeI_group_failure': len(obsI),
                'canonical_TypeI_group_pass_iff_2_in_K': True, 'canonical_quadratic_failure_both_branches': 0,
                'group_BFS_and_order_LCM_agree': True, 'failure_cardinality_bound_all_agree': True}
for k, v in mine_summary.items():
    if ts.get(k) != v:
        mism.append(('summary', k, v, ts.get(k)))
for mine, th in zip(obsI, their['group_I_obstruction_rows']):
    for k in ['p', 'q', 't', 'sigma', 'c', 'x', 'factors', 'S', 'hit_I', 'hit_II', 'K', 'group_I', 'group_II', 'quadratic_I', 'quadratic_II']:
        a = [list(f) for f in mine[k]] if k == 'factors' else mine[k]
        if a != th[k]:
            mism.append((f"obstruction:{mine['p']},{mine['q']}", k, a, th[k]))
script_sha = hashlib.sha256((THEIRS / 'check-reciprocity-20261018.py').read_bytes()).hexdigest()
report['script_sha256_matches_summary'] = (script_sha == their['script_sha256'])
extra = {'canonical_hit_II': sum(r['hit_II'] for r in canon), 'canonical_hit_I': sum(r['hit_I'] for r in canon),
         'canonical_both_fail': sum((not r['hit_I']) and (not r['hit_II']) for r in canon),
         'canonical_x_even': sum(r['x_even'] for r in canon), 'canonical_two_in_K': sum(r['two_in_K'] for r in canon),
         'canonical_TypeI_fail_with_group_pass': sum((not r['hit_I']) and r['group_I'] for r in canon),
         'all_hit_II': sum(r['hit_II'] for r in rows), 'all_hit_I': sum(r['hit_I'] for r in rows),
         'all_both_fail': sum((not r['hit_I']) and (not r['hit_II']) for r in rows),
         'S_card_max_canonical': max(len(r['S']) for r in canon),
         'TypeII_fail_slack_min': min(r['K_card'] - 2 - 2 * len(r['S']) for r in failII),
         'per_prime_candidates': {str(p): sum(1 for r in rows if r['p'] == p) for p in PS}}
report['part1'] = {'summary': mine_summary, 'obstruction_rows': [(r['p'], r['q'], r['x'], r['K_card'], len(r['S'])) for r in obsI],
                   'measurements': extra, 'mismatches_so_far': len(mism)}
say(f"1. {len(rows)} candidates, {n_checks} factor checks ({n_two} with r=2), violations 0; canonical {len(canon)}: II-fail {len(failII)} "
    f"(all group-pass), I-group-obstruction {len(obsI)}, (group I <=> 2 in K) on all; mismatches {len(mism)}; script sha ok={report['script_sha256_matches_summary']}")

# ------------------------------------------------------------------ 2. finite parts of the closed-set corollary (>>70 section 3)
paper = {}
for p in PS:
    checktime()
    ell = next(r for r in primerange(2, p) if leg(r, p) == -1)
    M = math.isqrt(p) + (0 if math.isqrt(p) ** 2 == p else 1)
    assert ell >= 11 and ell <= M
    sig = 1 if ell % 4 == 3 else 3
    F = lambda j: (p + sig * ell) // 4 + sig * ell * j * (j + 1)
    assert (p + sig * ell) % 4 == 0
    assert all(F(j) % ell != 0 for j in range(ell))
    others = [r for r in primerange(2, 400) if leg(r, p) == -1 and r != ell]
    rootinfo = {}
    for r in others:
        roots = [j for j in range(r) if F(j) % r == 0]
        assert len(roots) in (0, 2) and (-sig * ell * p) % r != 0
        assert F((r - 1) // 2) % r != 0
        if roots:
            a, b = roots
            assert (a + b + 1) % r == 0 and a < b
            rootinfo[r] = (a, b)
            # prefix error of one pair, all N up to 3r
            for N in range(0, 3 * r + 1):
                cnt = sum(1 for j in range(N) if j % r in (a, b))
                assert abs(Fr(cnt) - Fr(2 * N, r)) < 1
    # CRT pairs for two primes with roots: 4 common roots, 2 pairs, each with prefix error < 1
    crt_checked = 0
    for r1, r2 in list(itertools.combinations(sorted(rootinfo)[:6], 2)):
        d = r1 * r2
        common = [j for j in range(d) if F(j) % r1 == 0 and F(j) % r2 == 0]
        assert len(common) == 4 and all((d - 1 - j) in common for j in common)
        pairs = {tuple(sorted((j, d - 1 - j))) for j in common}
        assert len(pairs) == 2
        for a, b in pairs:
            for N in range(0, 2 * d + 1, max(1, d // 50)):
                cnt = sum(1 for j in range(N) if j % d in (a, b))
                assert abs(Fr(cnt) - Fr(2 * N, d)) < 1
        crt_checked += 1
    # inclusion-exclusion bound for C = {ell} + next n nonresidue primes, N = 4^n
    iebound = []
    for n in range(1, 6):
        C = others[:n]
        N = 4 ** n
        prodC = math.prod(C)
        Z = sum(1 for j in range(N) if math.gcd(F(j), prodC) == 1)
        active = [r for r in C if r in rootinfo]
        h = len(active)
        assert all(r >= ell + 2 * (i + 1) for i, r in enumerate(sorted(active)))
        dens = math.prod(Fr(r - 2, r) for r in active)
        assert dens >= Fr(11, 11 + 2 * n)
        if h == 0:
            assert Z == N
        else:
            assert abs(Fr(Z) - N * dens) < Fr(3 ** h - 1, 2)
        assert Z > 0
        iebound.append({'n': n, 'C': C, 'N': N, 'Z': Z, 'h': h, 'main': str(N * dens), 'error_bound': str(Fr(3 ** h - 1, 2))})
    paper[str(p)] = {'ell': ell, 'ceil_sqrt_p': M, 'sigma': sig, 'primes_with_roots_below_400': len(rootinfo), 'primes_checked': len(others),
                     'crt_pairs_checked': crt_checked, 'inclusion_exclusion': iebound}
# arithmetic lemma
for n in range(0, 201):
    assert (1 + Fr(2 * n, 11)) * 3 ** n <= 4 ** n
    assert Fr(3 * (2 * n + 13), 2 * n + 11) <= 4
report['part2_paper_finite'] = paper
say('2. corollary finite parts: ell >= 11 and ell <= ceil(sqrt p) for 12 primes; roots 0/2 paired, no fixed point; pair prefix error < 1; '
    'inclusion-exclusion bound and Z > 0 for n = 1..5; (1+2n/11)3^n <= 4^n for n <= 200')

report['mismatch_count'] = len(mism)
report['seconds'] = round(time.monotonic() - T0, 1)
(HERE / 'check73_results.json').write_text(json.dumps(report, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
(HERE / 'check73_rows.json').write_text(json.dumps(rows, ensure_ascii=False) + '\n', encoding='utf-8')
say(f"done: mismatches {len(mism)}, {report['seconds']} s")
