"""Independent check of >>66 (Outcasts 80c02e64, ジェミー, 2026-10-16): scan of the 179,468 primes p < 10^8 in the six
classes for H_Elog / H_E over the family E (81 pairs), the counterexample p = 38,296,729, and the layer table of section 1.
Own implementation from contract.md (committed before computation). The code in >>66 is read but not executed.
Reuses Research/CertificateRate20261016/gtest.py and chars20261015.py (same author) for the character-level items.
Usage: python -I scan66.py   (writes scan66_results.json and scan66_primes.csv.gz next to this file)
"""
import sys, os, json, gzip, math, random, time, hashlib
from collections import Counter
from fractions import Fraction as Fr
from sympy import factorint, isprime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
GT = os.path.join(ROOT, 'Research', 'CertificateRate20261016')
sys.path.insert(0, GT)
import gtest  # noqa: E402
from chars20261015 import row as fullrow, divisors_direct, divisors_from_factors, ser  # noqa: E402

CLASSES = (1, 121, 169, 289, 361, 529)
E = sorted(gtest.E)
PAIRS = sorted(((a, k // a) for k in E for a in range(1, k + 1) if k % a == 0 and math.gcd(a, k // a) == 1),
               key=lambda t: (t[0] * t[1], t[0]))
SEED = 20261016
LIMIT = 10 ** 8
T0 = time.time()


def sieve(n):
    s = bytearray([1]) * n
    s[0] = s[1] = 0
    for i in range(2, math.isqrt(n - 1) + 1):
        if s[i]:
            s[i * i::i] = bytes(len(range(i * i, n, i)))
    return s


def wcount(N, m):
    f = factorint(N)
    assert math.prod(q ** e for q, e in f.items()) == N, N
    odd = [(q, e) for q, e in f.items() if q != 2]
    divs = [1]
    for q, e in odd:
        divs = [d * q ** j for d in divs for j in range(e + 1)]
    wit = sorted(d for d in divs if d % m == m - 1)
    dp = {1: 1}
    for q, e in odd:
        nd = {}
        for r, c in dp.items():
            x = r
            for j in range(e + 1):
                nd[x] = nd.get(x, 0) + c
                x = x * q % m
        dp = nd
    assert len(wit) == dp.get(m - 1, 0), (N, m)
    if wit:
        s = wit[0]
        assert N % s == 0 and (s + 1) % m == 0            # certificate, independent of the factorization
    return f, wit


def full_check(N, m, f, wit):
    assert all(isprime(q) for q in f)
    alld = divisors_direct(N)
    assert alld == divisors_from_factors(f)
    assert [d for d in alld if d % m == m - 1] == wit


def scan():
    sv = sieve(LIMIT)
    primes = [p for p in range(LIMIT) if sv[p] and p % 840 in CLASSES]
    del sv
    assert all(isprime(p) for p in primes)
    h = hashlib.sha256(''.join(f'{p}\n' for p in primes).encode()).hexdigest()
    print(f'primes {len(primes)} sha256 {h} ({time.time() - T0:.0f}s)', flush=True)
    rng = random.Random(SEED)
    recs = []
    nfull = nsample = 0
    out = gzip.open(os.path.join(HERE, 'scan66_primes.csv.gz'), 'wt', encoding='utf-8', newline='\n')
    out.write('p,B,E_log,E_all,witness_pairs\n')
    for i, p in enumerate(primes):
        B = (p - 1).bit_length()
        rows = []
        for a, u in PAIRS:
            N, m = a * p + u, 4 * a * u
            f, wit = wcount(N, m)
            samp = rng.random() < 0.01
            rows.append((a, u, N, m, f, wit, samp))
        E_log = sum(len(r[5]) for r in rows if r[0] * r[1] <= B)
        E_all = sum(len(r[5]) for r in rows)
        for a, u, N, m, f, wit, samp in rows:
            if samp or E_log <= 2:
                full_check(N, m, f, wit)
                nsample += samp
                nfull += 1
        wp = [(a, u, len(wit), wit[0]) for a, u, N, m, f, wit, samp in rows if wit]
        recs.append({'p': p, 'B': B, 'E_log': E_log, 'E_all': E_all, 'wp': wp})
        out.write(f"{p},{B},{E_log},{E_all},{';'.join(f'{a}:{u}:{w}:{s}' for a, u, w, s in wp)}\n")
        if (i + 1) % 10000 == 0:
            print(f'{i + 1}/{len(primes)} primes, {time.time() - T0:.0f}s', flush=True)
        assert time.time() - T0 < 2 * 3600
    out.close()
    return primes, h, recs, nfull, nsample


def item4():
    p = 38296729
    B, rows = gtest.prime_rows(p, random.Random(SEED), full_check=True)
    s = gtest.summarize(p, B, rows)
    exp = {(1, 25): ({'199': 1, '96223': 1}, [199], '1/10', '9/10', '-9/5', '27/10', '0', '0'),
           (3, 8): ({'5': 1, '7': 1, '43': 1, '97': 1, '787': 1}, [3935, 381695], '1', '1', '-3', '4', '0', '0'),
           (8, 1): ({'3': 4, '83': 1, '199': 1, '229': 1}, [1791, 19007], '5/2', '-1/2', '-3/2', '1', '1', '1'),
           (11, 1): ({'3': 1, '5': 1, '73': 1, '96179': 1}, [219], '4/5', '1/5', '-2', '11/5', '0', '0')}
    got = {}
    for r in rows:
        if r['in_budget'] and r['W'] > 0:
            fr = fullrow(p, r['a'], r['u'], with_orbits=False, with_g2=True)
            oddf = {q: e for q, e in fr['factors'].items() if q != '2'}
            got[(r['a'], r['u'])] = (oddf, fr['witnesses'], str(fr['R']), str(fr['C']), str(fr['ell_sum']), str(fr['eta']),
                                    str(max(Fr(0), fr['W'] - fr['eta'])), str(fr['b']))
    beyond = [(r['a'], r['u'], r['W']) for r in rows if not r['in_budget'] and r['W'] > 0]
    N149 = 14 * p + 9
    return {'p': p, 'mod840': p % 840, 'B': B, 'pairs': s['pairs'], 'Wtot': s['Wtot'], 'S': s['S'],
            'all_W_pos_outside_E_eta_pos': all(not r['inE'] and r['eta'] > 0 for r in rows if r['in_budget'] and r['W'] > 0),
            'E_in_budget_rows': s['E_pairs_box'], 'E_log': s['E_log'], 'Cert0': s['Cert0'], 'Cert_plus': s['Cert_plus'], 'Cert': s['Cert'],
            'rows_match': got == exp, 'rows': {f'{k[0]},{k[1]}': v for k, v in got.items()},
            'beyond_E_witness_rows': beyond, 'E_all': s['E_all'],
            '14_9': {'N': N149, 'factors': {str(q): e for q, e in factorint(N149).items()}, 's': 157 * 3467,
                     'divides': N149 % (157 * 3467) == 0, 's_mod_504': (157 * 3467) % 504}}


def item5():
    layer = {}
    for k in E:
        L = gtest.lam(4 * k)
        m = 4 * k
        g = gtest.group(m)
        G2 = g['G2']
        expo = 1
        for x in G2:
            o, y = 1, x
            while y != 1:
                y = y * x % m
                o += 1
            expo = math.lcm(expo, o)
        layer[k] = {'lambda': L, 'exp_G2': expo, 'G2size': len(G2), 'pairs': 2 ** len(factorint(k)) if k > 1 else 1}
    lay = lambda k: {2: 1, 4: 2, 6: 3}[layer[k]['lambda']]
    per = {}
    for k in E:
        d = per.setdefault(lay(k), {'k': [], 'pairs': 0, 'G2sizes': set(), 'exp_G2': set()})
        d['k'].append(k)
        d['pairs'] += layer[k]['pairs']
        d['G2sizes'].add(layer[k]['G2size'])
        d['exp_G2'].add(layer[k]['exp_G2'])
    exp_half = all(layer[k]['exp_G2'] * 2 == layer[k]['lambda'] for k in E)
    tab = {}
    grp = []
    for q in gtest.DISCOVERY:
        B, rr = gtest.prime_rows(q, random.Random(SEED))
        for r in rr:
            if not r['in_budget']:
                continue
            if r['cls'] == 'group':
                grp.append({'p': q, 'a': r['a'], 'u': r['u'], 'k': r['k'], 'layer': lay(r['k']) if r['inE'] else None})
            if r['inE']:
                t = tab.setdefault(lay(r['k']), Counter())
                t['rows'] += 1
                t['Wpos'] += r['W'] > 0
                if r['cls']:
                    t[r['cls']] += 1
    return {'per_layer': {str(k): {'k': v['k'], 'pairs': v['pairs'], 'G2sizes': sorted(v['G2sizes']), 'exp_G2': sorted(v['exp_G2'])}
                          for k, v in sorted(per.items())},
            'exp_G2_equals_lambda_half_on_E': exp_half,
            'claimed_pairs_per_layer': [8, 32, 41],
            'discovery_rows_per_layer': {str(k): dict(v) for k, v in sorted(tab.items())},
            'claimed_rows': {'1': [90, 3, 87, 0, 0], '2': [192, 13, 153, 6, 20], '3': [144, 15, 99, 0, 30]},
            'group_rows': grp}


def item6():
    p, a, u = 345601, 4, 1
    N, m = a * p + u, 16
    f = factorint(N)
    odd = [q for q in f if q != 2]
    K = gtest.closure(m, [q % m for q in odd], [1])
    G2 = gtest.group(m)['G2']
    W = sum(1 for d in divisors_from_factors(f) if d % m == m - 1)
    o5 = next(k for k in range(1, 17) if pow(5, k, 16) == 1)
    return {'N': N, 'factors': {str(q): e for q, e in f.items()}, 'minus1_in_K': (m - 1) in K, 'W': W, 'order_of_5_mod_16': o5,
            '5sq_mod16': 25 % 16, 'in_G2_minus_1': (25 % 16) in G2 and 25 % 16 != 1, 'exponent_of_5': f.get(5, 0), 'G2': sorted(G2)}


def item7(recs):
    by = {r['p']: r for r in recs}
    ref = [json.loads(l) for l in gzip.open(os.path.join(GT, 'gtest_primes.jsonl.gz'), 'rt', encoding='utf-8')]
    bad = [(s['p'], s['E_log'], s['E_all'], by[s['p']]['E_log'], by[s['p']]['E_all']) for s in ref
           if (s['E_log'], s['E_all']) != (by[s['p']]['E_log'], by[s['p']]['E_all'])]
    return {'primes': len(ref), 'mismatches': bad}


def main():
    primes, h, recs, nfull, nsample = scan()
    t_scan = time.time() - T0
    low = [r for r in recs if r['p'] < 2 * 10 ** 7]
    exc_log = [r for r in recs if r['E_log'] == 0]
    exc_all = [r for r in recs if r['E_all'] == 0]
    res = {'item1_counts': {'primes': len(primes), 'below_2e7': len(low), 'sha256': h, 'claims': [179468, 39391],
                            'match': [len(primes), len(low)] == [179468, 39391]},
           'item2_H_Elog_exceptions': [ser(r) for r in exc_log], 'item2_match': [r['p'] for r in exc_log] == [38296729],
           'item2_below_2e7': [r['p'] for r in exc_log if r['p'] < 2 * 10 ** 7],
           'item3_H_E_exceptions': [ser(r) for r in exc_all], 'item3_match': not exc_all,
           'full_checks': {'rows_with_trial_division_and_primality': nfull, 'of_which_1pct_sample': nsample}}
    res['item4_38296729'] = ser(item4())
    res['item5_layers'] = ser(item5())
    res['item6_345601_4_1'] = ser(item6())
    res['item7_vs_gtest'] = item7(recs)
    meas = {}
    for lo, hi in ((0, 10 ** 6), (10 ** 6, 10 ** 7), (10 ** 7, 10 ** 8)):
        rr = [r for r in recs if lo <= r['p'] < hi]
        ea = sorted(r['E_all'] for r in rr)
        el = sorted(r['E_log'] for r in rr)
        meas[f'[{lo},{hi})'] = {'primes': len(rr), 'E_log_min': el[0], 'E_log_median': el[len(el) // 2],
                                 'E_all_min': ea[0], 'E_all_median': ea[len(ea) // 2],
                                 'E_log_le2': sum(1 for x in el if x <= 2), 'E_all_le2': sum(1 for x in ea if x <= 2)}
    res['measure_by_range'] = meas
    res['measure_E_all_le2'] = [ser(r) for r in recs if r['E_all'] <= 2]
    res['measure_E_log_le1'] = [ser(r) for r in recs if r['E_log'] <= 1]
    res['seconds'] = {'scan': round(t_scan, 1), 'total': round(time.time() - T0, 1)}
    json.dump(res, open(os.path.join(HERE, 'scan66_results.json'), 'w', encoding='utf-8', newline='\n'), ensure_ascii=False, indent=1)
    sys.stdout.reconfigure(encoding='utf-8')
    short = {k: v for k, v in res.items() if k not in ('measure_E_all_le2', 'measure_E_log_le1')}
    print(json.dumps(short, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
