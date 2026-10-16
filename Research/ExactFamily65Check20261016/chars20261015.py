"""Shared machinery for the independent check of >>60 / >>61 (Outcasts 80c02e64, 2026-10-15).
Own implementation from contract.md; descends from Research/Residual57Check20261013/check57_cert.py (same author, 2026-10-13).
The participants' scripts are neither imported nor executed.
"""
import math, itertools, cmath
from fractions import Fraction as Fr
from sympy import factorint, totient, primitive_root, jacobi_symbol, cyclotomic_poly, Poly, symbols, isprime, mobius

x = symbols('x')


def ceil_log2(p):
    return (p - 1).bit_length()


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
    """All Dirichlet characters mod m: (L, units, chars). chars[i] = dict(param, phase: unit -> Z/L, sign = chi(-1), real).
    CRT into cyclic blocks; unit r has coordinates (j_1..j_n) with r = prod g_i^{j_i};
    character (k_1..k_n) has chi(r) = zeta_L^{sum k_i j_i L/o_i}. Multiplicativity is asserted for every pair of units."""
    if m in _char_cache:
        return _char_cache[m]
    blocks = []
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
        for r in units:
            for s in units:
                assert phase[r * s % m] == (phase[r] + phase[s]) % L
        sgn = phase[m - 1]
        assert sgn in (0, L // 2) or L == 1
        chars.append({'param': param, 'phase': phase, 'sign': 1 if sgn == 0 else -1,
                      'real': all(v in (0, L // 2) for v in phase.values())})
    keys = {tuple(c['phase'][r] for r in units) for c in chars}
    assert len(keys) == phi
    bykey = {tuple(c['phase'][r] for r in units): c for c in chars}
    _char_cache[m] = (L, units, chars, bykey)
    return _char_cache[m]


def reduce_mod_cyclo(coeffs_by_power, L):
    """coeffs_by_power: dict power -> int (powers of zeta_L). Returns the remainder modulo Phi_L (coefficients low -> high)."""
    if L == 1:
        return [sum(coeffs_by_power.values())]
    p = Poly(sum(c * x ** (k % L) for k, c in coeffs_by_power.items()), x, domain='ZZ')
    r = p.rem(Poly(cyclotomic_poly(L, x), x, domain='ZZ'))
    cs = [int(c) for c in r.all_coeffs()[::-1]]
    while len(cs) > 1 and cs[-1] == 0:
        cs.pop()
    return cs if cs else [0]


def F_poly(phase, odd_factors, m, L):
    """F_chi = prod_(q^e || N_odd) (1 + chi(q) + ... + chi(q)^e) as dict power(mod L) -> coeff."""
    poly = {0: 1}
    for q, e in odd_factors:
        ph = phase[q % m]
        new = {}
        for k, c in poly.items():
            for j in range(e + 1):
                kk = (k + j * ph) % L
                new[kk] = new.get(kk, 0) + c
        poly = new
    return poly


def add_polys(a, b, sign=1):
    s = dict(a)
    for k, v in b.items():
        s[k] = s.get(k, 0) + sign * v
    return s


def ramanujan_mobius(d, j):
    return sum(h * int(mobius(d // h)) for h in range(1, d + 1) if d % h == 0 and math.gcd(d, j) % h == 0)


def ramanujan_direct(d, j):
    z = sum(cmath.exp(2j * cmath.pi * k * j / d) for k in range(1, d + 1) if math.gcd(k, d) == 1)
    assert abs(z.imag) < 1e-6 and abs(z.real - round(z.real)) < 1e-6
    return int(round(z.real))


def orbits_of(chars, bykey, units, L, odd, m):
    """>>60: complete power-conjugacy orbits of the complex characters. For each orbit: order d, size, sign,
    D(j) (phase DP in Z[X]/(X^d-1) with the actual exponents), T_O = sign * sum_j D(j) c_d(j),
    and the direct constant sum_{k} chi^k(-1) F_{chi^k} reduced mod Phi_L (must equal T_O)."""
    seen = set()
    orbits = []
    for c in chars:
        if c['real'] or c['param'] in seen:
            continue
        g = math.gcd(L, *c['phase'].values())
        d = L // g
        assert d > 2
        members = []
        for k in range(1, d + 1):
            if math.gcd(k, d) == 1:
                mem = bykey[tuple((k * c['phase'][r]) % L for r in units)]
                members.append(mem)
        assert len({mm['param'] for mm in members}) == len(members) == int(totient(d))
        assert all(mm['sign'] == c['sign'] for mm in members)
        seen.update(mm['param'] for mm in members)
        D = [1] + [0] * (d - 1)
        for q, e in odd:
            b = c['phase'][q % m] // g
            assert c['phase'][q % m] % g == 0
            local = [0] * d
            for j in range(e + 1):
                local[(j * b) % d] += 1
            new = [0] * d
            for i, ci in enumerate(D):
                if ci:
                    for j, lj in enumerate(local):
                        if lj:
                            new[(i + j) % d] += ci * lj
            D = new
        assert sum(D) == math.prod(e + 1 for _, e in odd)
        cd = [ramanujan_mobius(d, j) for j in range(d)]
        assert cd == [ramanujan_direct(d, j) for j in range(d)]
        T = c['sign'] * sum(Dj * cj for Dj, cj in zip(D, cd))
        direct = {}
        for mm in members:
            direct = add_polys(direct, F_poly(mm['phase'], odd, m, L), mm['sign'])
        red = reduce_mod_cyclo(direct, L)
        assert red == [T], (d, D, T, red)
        orbits.append({'order': d, 'size': len(members), 'sign': c['sign'], 'D': D, 'trace': T,
                       'representative_param': list(c['param'])})
    assert len(seen) == sum(1 for c in chars if not c['real'])
    return orbits


def row(p, a, u, with_orbits=True, with_g2=True):
    N = a * p + u
    m = 4 * a * u
    phi = int(totient(m))
    f = factorint(N)
    assert math.gcd(N, a * u) == 1
    odd = sorted((q, e) for q, e in f.items() if q != 2)
    tau = math.prod(e + 1 for _, e in odd)
    L, units, chars, bykey = characters(m)
    jac = {r: int(jacobi_symbol((-a * u) % r, r)) for r in units}
    dist = [c for c in chars if c['real'] and all((1 if c['phase'][r] == 0 else -1) == jac[r] for r in units)]
    assert len(dist) == 1
    chi0 = dist[0]
    assert chi0['sign'] == -1
    for q, _ in odd:
        assert jacobi_symbol((-a * u) % q, q) == jacobi_symbol(p % q, q)
    raw = {c['param']: F_poly(c['phase'], odd, m, L) for c in chars}
    polys = {k: reduce_mod_cyclo(v, L) for k, v in raw.items()}
    principal = tuple(0 for _ in chi0['param'])
    assert polys[principal] == [tau]
    for c in chars:
        if c['real']:
            assert len(polys[c['param']]) == 1 and polys[c['param']][0] >= 0
    tnum = polys[chi0['param']][0]
    assert 0 <= tnum <= tau
    A_num = sum(polys[c['param']][0] for c in chars if c['real'] and c['sign'] == 1 and c['param'] != principal)
    B_num = sum(polys[c['param']][0] for c in chars if c['real'] and c['sign'] == -1 and c['param'] != chi0['param'])
    # complex conjugate pairs (>>57 convention)
    seen = set()
    pairs = []
    Ctot = {}
    ell_sum = 0
    for c in chars:
        if c['real'] or c['param'] in seen:
            continue
        conj = bykey[tuple((-c['phase'][r]) % L for r in units)]
        seen.update((c['param'], conj['param']))
        s = add_polys(raw[c['param']], raw[conj['param']])
        red = [c['sign'] * k for k in reduce_mod_cyclo(s, L)]
        ell = red[0] if len(red) == 1 else -sum(abs(k) for k in red)
        ell_sum += ell
        Ctot = add_polys(Ctot, s, c['sign'])
        pairs.append({'sign': c['sign'], 'signed_pair_poly': red, 'ell': ell})
    Cred = reduce_mod_cyclo(Ctot, L) if Ctot else [0]
    assert len(Cred) == 1, (p, a, u, Cred)
    C_num = Cred[0]
    ds = divisors_direct(N)
    assert ds == divisors_from_factors(f)
    witnesses = [d for d in ds if (d + 1) % m == 0]
    W = len(witnesses)
    total = tau - tnum + A_num - B_num + C_num
    assert total % phi == 0 and total // phi == W, (p, a, u, total, W)
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
    ell = Fr(ell_sum, phi)
    H = A - B + ell
    K = min(V, max(Fr(0), -H))
    lower = max(Fr(0), V + H)
    r_ = max(Fr(0), V - W)
    delta = t - min(t, max(Fr(0), U - W))
    eta = C - ell
    b = V - K
    assert 0 <= r_ <= K <= V
    assert eta >= 0
    # >>60 identities
    assert b == min(V, max(Fr(0), W - eta))
    assert (b > 0) == (V > 0 and W > eta)
    assert phi * V >= 2 * W
    odd_divs = [d for d in ds if d % 2 == 1]
    assert len(odd_divs) == tau and all(math.gcd(d, m) == 1 for d in odd_divs)
    n_minus = sum(1 for d in odd_divs if jac[d % m] == -1)
    assert tau - tnum == 2 * n_minus
    assert V == Fr(2 * n_minus, phi)
    out = {'p': p, 'a': a, 'u': u, 'N': N, 'm': m, 'phi': phi, 'L': L,
           'factors': {str(q): e for q, e in sorted(f.items())}, 'tau_odd': tau,
           'U': U, 't': t, 'V': V, 'A': A, 'B': B, 'C': C, 'ell_sum': ell, 'H': H, 'K': K,
           'lower': lower, 'W': W, 'witnesses': witnesses, 'r': r_, 'delta': delta, 'eta': eta, 'b': b,
           'n_minus': n_minus,
           'n_real': sum(1 for c in chars if c['real']),
           'n_real_nonprincipal_zero': sum(1 for c in chars if c['real'] and c['param'] != principal and polys[c['param']][0] == 0),
           'complex_pair_numerators': sorted(pr['signed_pair_poly'][0] for pr in pairs if len(pr['signed_pair_poly']) == 1),
           'complex_pairs_nonconstant': sum(1 for pr in pairs if len(pr['signed_pair_poly']) > 1)}
    if with_orbits:
        orbs = orbits_of(chars, bykey, units, L, odd, m)
        assert sum(o['trace'] for o in orbs) == C_num
        out['orbits'] = orbs
        out['orbit_C'] = Fr(sum(o['trace'] for o in orbs), phi)
        out['orbit_b'] = V - min(V, max(Fr(0), -(A - B + out['orbit_C'])))
        assert out['orbit_b'] == min(V, W)
    if with_g2:
        G2 = sorted({(r * r) % m for r in units})
        negG2 = {(-s) % m for s in G2}
        W_G2 = sum(1 for d in odd_divs if d % m in negG2)
        R = V + A - B
        assert R == Fr(sum(c['sign'] * polys[c['param']][0] for c in chars if c['real']), phi)
        # coset counts of odd divisors in G/G^2
        cosets = {}
        for r in units:
            key = frozenset((r * s) % m for s in G2)
            cosets.setdefault(key, 0)
        for d in odd_divs:
            key = frozenset((d * s) % m for s in G2)
            cosets[key] += 1
        assert len(cosets) == phi // len(G2) and sum(cosets.values()) == tau
        span = set(G2)
        changed = True
        while changed:
            changed = False
            for q, _ in odd:
                new = {(g * q) % m for g in span}
                if not new <= span:
                    span |= new
                    changed = True
        # bundled identity (>>61 §3): for every complex psi,
        #   sum_{chi real} (psi chi)(-1) F_{psi chi}  ==  [G:G^2] * sum_{s odd | N, s in -G^2} psi(-s)   (as elements of Z[zeta_L])
        index = phi // len(G2)
        bundled_ok = True
        for c in chars:
            if c['real']:
                continue
            lhs = {}
            for ch in chars:
                if ch['real']:
                    prod = bykey[tuple((c['phase'][r] + ch['phase'][r]) % L for r in units)]
                    lhs = add_polys(lhs, raw[prod['param']], prod['sign'])
            rhs = {}
            for d in odd_divs:
                if d % m in negG2:
                    k = c['phase'][(-d) % m]
                    rhs[k] = rhs.get(k, 0) + index
            if reduce_mod_cyclo(lhs, L) != reduce_mod_cyclo(rhs, L):
                bundled_ok = False
        out.update({'G2_size': len(G2), 'G2': G2, 'negG2': sorted(negG2), 'W_G2': W_G2, 'R': R,
                    'R_identity': R == Fr(W_G2, len(G2)), 'C_identity': C == W - R,
                    'coset_counts': sorted(cosets.values()),
                    'all_cosets_equal': len(set(cosets.values())) == 1,
                    'negG2_counts': {str(s): sum(1 for d in odd_divs if d % m == s) for s in sorted(negG2)},
                    'minus1_in_span': (m - 1) in span,
                    'span_identity': ((m - 1) in span) == (W_G2 > 0),
                    'bundled_identity': bundled_ok})
    return out


def ser(o):
    if isinstance(o, Fr):
        return str(o)
    if isinstance(o, dict):
        return {k: ser(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [ser(v) for v in o]
    return o


def box_pairs(p):
    B = ceil_log2(p)
    return B, [(a, u) for a in range(1, B + 1) for u in range(1, B + 1) if a * u <= B and math.gcd(a, u) == 1]
