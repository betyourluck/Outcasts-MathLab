"""Independent check of >>57 (Outcasts 80c02e64, 2026-10-13): factor-only certificate K_i, complex-pair numerators,
box sums, and the Lean-free parts of PR #4 b84bbb7 Research/ResidualCharacters20261013/.
Own implementation from the contract (contract.md); the participant's scripts are not imported or executed.
Usage: python -I check57_cert.py <dir with the participant's residual-audit-results.json and residual-audit-rows.json>
"""
import sys, json, time, math, itertools, cmath
from fractions import Fraction as Fr
from pathlib import Path
from sympy import factorint, totient, primitive_root, jacobi_symbol, cyclotomic_poly, Poly, symbols, isprime

x = symbols('x')
T0 = time.time()


def divisors_direct(n):
    out = set()
    for d in range(1, math.isqrt(n) + 1):
        if n % d == 0:
            out.update((d, n // d))
    return sorted(out)


def divisors_from_factors(f):
    out = [1]
    for q, e in f.items():
        out = [d * q ** k for d in out for k in range(e + 1)]
    return sorted(out)


_char_cache = {}


def characters(m):
    """All Dirichlet characters mod m: (L, units, [dict(param, phase: unit->Z/L, sign=chi(-1), real)]).
    CRT into cyclic blocks; unit r gets coordinates (j_1..j_n) with r = prod g_i^{j_i};
    a character is a tuple (k_1..k_n) with chi(r) = zeta_L^{sum k_i j_i L/o_i}."""
    if m in _char_cache:
        return _char_cache[m]
    blocks = []  # (modulus, table r mod n -> coords, orders)
    for q, e in sorted(factorint(m).items()):
        n = q ** e
        if q == 2:
            if e == 1:
                blocks.append((n, {1: ()}, ()))
            elif e == 2:
                blocks.append((n, {1: (0,), 3: (1,)}, (2,)))
            else:
                k = 2 ** (e - 2)
                tab = {}
                for b in range(2):
                    for j in range(k):
                        tab[((-1) ** b * pow(5, j, n)) % n] = (b, j)
                assert len(tab) == n // 2
                blocks.append((n, tab, (2, k)))
        else:
            g = int(primitive_root(n))
            ph = int(totient(n))
            tab = {pow(g, j, n): (j,) for j in range(ph)}
            assert len(tab) == ph
            blocks.append((n, tab, (ph,)))
    orders = tuple(o for _, _, os_ in blocks for o in os_)
    units = [r for r in range(1, m) if math.gcd(r, m) == 1]
    coords = {r: tuple(c for n, tab, _ in blocks for c in tab[r % n]) for r in units}
    phi = int(totient(m))
    assert len(units) == phi == (math.prod(orders) if orders else 1)
    L = math.lcm(*orders) if orders else 1
    chars = []
    for param in itertools.product(*(range(o) for o in orders)):
        phase = {r: sum(k * j * (L // o) for k, j, o in zip(param, coords[r], orders)) % L for r in units}
        for r in units:  # multiplicativity
            for s in units:
                assert phase[r * s % m] == (phase[r] + phase[s]) % L
        sgn = phase[m - 1]
        assert sgn in (0, L // 2) or L == 1
        chars.append({'param': param, 'phase': phase, 'sign': 1 if sgn == 0 else -1,
                      'real': all(v in (0, L // 2) for v in phase.values())})
    assert len({tuple(c['phase'][r] for r in units) for c in chars}) == phi
    _char_cache[m] = (L, units, chars)
    return _char_cache[m]


def reduce_mod_cyclo(coeffs_by_power, L):
    """coeffs_by_power: dict power->int. Returns coefficients (low->high) of the remainder modulo Phi_L."""
    if L == 1:
        return [sum(coeffs_by_power.values())]
    p = Poly(sum(c * x ** k for k, c in coeffs_by_power.items()), x, domain='ZZ')
    r = p.rem(Poly(cyclotomic_poly(L, x), x, domain='ZZ'))
    cs = [int(c) for c in r.all_coeffs()[::-1]]
    while len(cs) > 1 and cs[-1] == 0:
        cs.pop()
    return cs if cs else [0]


def F_poly(chi, odd_factors, m, L):
    """F_chi = prod_(q^e) (1 + chi(q) + ... + chi(q)^e) as dict power->coeff (exponents mod L)."""
    poly = {0: 1}
    for q, e in odd_factors:
        ph = chi['phase'][q % m]
        new = {}
        for k, c in poly.items():
            for j in range(e + 1):
                kk = (k + j * ph) % L
                new[kk] = new.get(kk, 0) + c
        poly = new
    return poly


def row(p, a, u):
    N = a * p + u
    m = 4 * a * u
    phi = int(totient(m))
    f = factorint(N)
    assert math.gcd(N, a * u) == 1
    odd = sorted((q, e) for q, e in f.items() if q != 2)
    tau = math.prod(e + 1 for _, e in odd)
    L, units, chars = characters(m)
    # distinguished odd real character chi_{-au} = Jacobi(-au | r) on units
    jac = {r: int(jacobi_symbol((-a * u) % r, r)) for r in units}
    dist = [c for c in chars if c['real'] and all((1 if c['phase'][r] == 0 else -1) == jac[r] for r in units)]
    assert len(dist) == 1
    chi0 = dist[0]
    assert chi0['sign'] == -1
    for q, _ in odd:  # chi0(q) = (p|q) on odd prime factors
        assert jacobi_symbol((-a * u) % q, q) == jacobi_symbol(p % q, q)
    polys = {c['param']: reduce_mod_cyclo(F_poly(c, odd, m, L), L) for c in chars}
    principal = tuple(0 for _ in chi0['param'])
    assert polys[principal] == [tau]
    for c in chars:
        if c['real']:
            assert len(polys[c['param']]) == 1 and polys[c['param']][0] >= 0
    tnum = polys[chi0['param']][0]
    assert 0 <= tnum <= tau
    A_num = sum(polys[c['param']][0] for c in chars if c['real'] and c['sign'] == 1 and c['param'] != principal)
    B_num = sum(polys[c['param']][0] for c in chars if c['real'] and c['sign'] == -1 and c['param'] != chi0['param'])
    # complex conjugate pairs
    seen = set()
    pairs = []
    Ctot = {}
    ell_sum = 0
    bykey = {tuple(c['phase'][r] for r in units): c for c in chars}
    for c in chars:
        if c['real'] or c['param'] in seen:
            continue
        conj = bykey[tuple((-c['phase'][r]) % L for r in units)]
        seen.update((c['param'], conj['param']))
        raw = F_poly(c, odd, m, L)
        rawc = F_poly(conj, odd, m, L)
        s = {}
        for k, v in list(raw.items()) + list(rawc.items()):
            s[k] = s.get(k, 0) + v
        red = [c['sign'] * k for k in reduce_mod_cyclo(s, L)]
        ell = red[0] if len(red) == 1 else -sum(abs(k) for k in red)
        ell_sum += ell
        for k, v in s.items():
            Ctot[k] = Ctot.get(k, 0) + c['sign'] * v
        pairs.append({'param': list(c['param']), 'conj': list(conj['param']), 'sign': c['sign'],
                      'signed_pair_poly': red, 'ell': ell})
    Cred = reduce_mod_cyclo(Ctot, L) if Ctot else [0]
    assert len(Cred) == 1, (p, a, u, Cred)
    C_num = Cred[0]
    ds = divisors_direct(N)
    assert ds == divisors_from_factors(f)
    W = sum(1 for d in ds if (d + 1) % m == 0)
    total = tau - tnum + A_num - B_num + C_num
    assert total % phi == 0 and total // phi == W, (p, a, u, total, W)
    # floating cross-check of orthogonality
    zeta = cmath.exp(2j * cmath.pi / L)
    Wf = 0
    for c in chars:
        F = 1
        for q, e in odd:
            z = zeta ** c['phase'][q % m]
            F *= sum(z ** k for k in range(e + 1))
        Wf += c['sign'] * F
    assert abs(Wf / phi - W) < 1e-6
    U = Fr(tau, phi)
    t = Fr(tnum, phi)
    V = U - t
    A = Fr(A_num, phi)
    B = Fr(B_num, phi)
    C = Fr(C_num, phi)
    H = A - B + Fr(ell_sum, phi)
    K = min(V, max(Fr(0), -H))
    lower = max(Fr(0), V + H)
    r = max(Fr(0), V - W)
    delta = t - min(t, max(Fr(0), U - W))
    assert 0 <= r <= K <= V
    assert W == V - r + max(Fr(0), W - U) + delta
    const_nums = sorted(pr['signed_pair_poly'][0] for pr in pairs if len(pr['signed_pair_poly']) == 1)
    return {'p': p, 'a': a, 'u': u, 'N': N, 'm': m, 'phi': phi, 'L': L,
            'factors': {str(q): e for q, e in sorted(f.items())}, 'tau_odd': tau,
            'U': U, 't': t, 'V': V, 'A': A, 'B': B, 'C': C, 'ell_sum': Fr(ell_sum, phi), 'H': H, 'K': K,
            'lower': lower, 'W': W, 'r': r, 'delta': delta,
            'n_real_nonprincipal_zero': sum(1 for c in chars if c['real'] and c['param'] != principal and polys[c['param']][0] == 0),
            'complex_pair_numerators': const_nums,
            'complex_pairs_nonconstant': sum(1 for pr in pairs if len(pr['signed_pair_poly']) > 1),
            'pairs': pairs}


def ser(o):
    if isinstance(o, Fr):
        return str(o)
    if isinstance(o, dict):
        return {k: ser(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [ser(v) for v in o]
    return o


KEYS = ['U', 't', 'V', 'W', 'rest', 'delta', 'factor_only_K', 'factor_only_W_lower',
        'other_even_real', 'other_odd_real', 'signed_complex']
MINE = {'U': 'U', 't': 't', 'V': 'V', 'W': 'W', 'rest': 'r', 'delta': 'delta', 'factor_only_K': 'K',
        'factor_only_W_lower': 'lower', 'other_even_real': 'A', 'other_odd_real': 'B', 'signed_complex': 'C'}


def main():
    their_dir = Path(sys.argv[1])
    theirs = json.loads(their_dir.joinpath('residual-audit-results.json').read_text())
    their_rows = json.loads(their_dir.joinpath('residual-audit-rows.json').read_text())
    their_by = {(r['p'], r['a'], r['u']): r for r in their_rows}
    out = {'contract': 'contract.md', 'boxes': [], 'row_matches': {}, 'mismatches': []}
    all_rows = []
    for p in (66529, 18490201):
        assert isprime(p)
        B = math.ceil(math.log2(p))
        t0 = time.time()
        pairs = [(a, u) for a in range(1, B + 1) for u in range(1, B + 1) if a * u <= B and math.gcd(a, u) == 1]
        rows = [row(p, a, u) for a, u in pairs]
        all_rows += rows
        assert time.time() - t0 < 120
        sums = {k: sum((Fr(rw[MINE[k]]) for rw in rows), Fr(0)) for k in KEYS}
        th = next(s for s in theirs['summaries'] if s['p'] == p)
        cmp = {k: {'mine': str(sums[k]), 'theirs': th['sums'][k], 'match': Fr(th['sums'][k]) == sums[k]} for k in KEYS}
        pos_pairs = [[rw['a'], rw['u']] for rw in rows if rw['lower'] > 0]
        box = {'p': p, 'B': B, 'pairs': len(pairs), 'pairs_theirs': th['pairs'], 'sums': cmp,
               'certificate_lower_sumV_minus_sumK': str(sums['V'] - sums['factor_only_K']),
               'sum_other_real_signed': str(sums['other_even_real'] - sums['other_odd_real']),
               'positive_certificate_pairs': pos_pairs,
               'positive_certificate_pairs_theirs': th['positive_certificate_pairs'],
               'positive_certificate_pairs_match': pos_pairs == th['positive_certificate_pairs'],
               'all_sums_match': all(v['match'] for v in cmp.values()) and th['pairs'] == len(pairs)}
        out['boxes'].append(box)
        for rw in rows:
            tr = their_by[(p, rw['a'], rw['u'])]
            if not all(Fr(tr[k]) == Fr(rw[MINE[k]]) for k in KEYS):
                out['mismatches'].append({'p': p, 'a': rw['a'], 'u': rw['u'],
                                          'mine': ser({k: rw[MINE[k]] for k in KEYS}),
                                          'theirs': {k: tr[k] for k in KEYS}})
        out['row_matches'][str(p)] = {'rows': len(rows),
                                      'all_fields_match': not any(mm['p'] == p for mm in out['mismatches'])}
    r114 = next(rw for rw in all_rows if (rw['p'], rw['a'], rw['u']) == (18490201, 1, 14))
    tr114 = their_by[(18490201, 1, 14)]
    their_nums = sorted(pr['signed_pair_poly'][0] for pr in tr114['complex_pairs'] if len(pr['signed_pair_poly']) == 1)
    claim = sorted([6, -2, 0, 0, 0, 0, -18, -18])
    out['row_18490201_1_14'] = {
        'tau': r114['tau_odd'], 'phi': r114['phi'], 'U': str(r114['U']), 't': str(r114['t']), 'W': r114['W'],
        'r': str(r114['r']), 'real_nonprincipal_zero': r114['n_real_nonprincipal_zero'],
        'complex_pair_numerators_sorted': r114['complex_pair_numerators'],
        'complex_pairs_nonconstant': r114['complex_pairs_nonconstant'], 'C': str(r114['C']),
        'thread_claim_multiset': claim,
        'match_thread': r114['complex_pair_numerators'] == claim and r114['C'] == Fr(-4, 3),
        'match_their_rows': r114['complex_pair_numerators'] == their_nums}
    s = row(1201, 3, 1)
    ts = theirs['sharpness_pair']
    out['row_1201_3_1'] = {k: str(s[MINE[k]]) for k in KEYS}
    out['row_1201_3_1']['match'] = all(Fr(ts[k]) == Fr(s[MINE[k]]) for k in KEYS)
    b1 = next(b for b in out['boxes'] if b['p'] == 66529)
    b2 = next(b for b in out['boxes'] if b['p'] == 18490201)
    out['thread_table_match'] = {
        '66529': {'sumV': b1['sums']['V']['mine'] == '55/8', 'sum_r': b1['sums']['rest']['mine'] == '39/8',
                  'sumK': b1['sums']['factor_only_K']['mine'] == '39/8',
                  'lower': b1['certificate_lower_sumV_minus_sumK'] == '2', 'Wtot': b1['sums']['W']['mine'] == '4',
                  'A_minus_B': b1['sum_other_real_signed'] == str(Fr(241, 48) - Fr(319, 48)),
                  'complex_signed': b1['sums']['signed_complex']['mine'] == '-5/4'},
        '18490201': {'sumV': b2['sums']['V']['mine'] == '12571/1320', 'sum_r': b2['sums']['rest']['mine'] == '3831/440',
                     'sumK': b2['sums']['factor_only_K']['mine'] == '4007/440',
                     'lower': b2['certificate_lower_sumV_minus_sumK'] == '5/12', 'Wtot': b2['sums']['W']['mine'] == '3',
                     'A_minus_B': b2['sum_other_real_signed'] == '-259/120',
                     'complex_signed': b2['sums']['signed_complex']['mine'] == '-2881/660'}}
    out['seconds'] = round(time.time() - T0, 2)
    Path('check57_results.json').write_text(json.dumps(ser(out), indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
    Path('check57_rows.json').write_text(
        json.dumps(ser([{k: v for k, v in rw.items() if k != 'pairs'} for rw in all_rows + [s]]), indent=0,
                   ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(ser({k: v for k, v in out.items() if k != 'mismatches'}), indent=1, ensure_ascii=False))
    print('mismatched rows:', len(out['mismatches']))
    for mm in out['mismatches'][:10]:
        print(mm)


if __name__ == '__main__':
    main()
