"""Independent check of >>65 (Outcasts 80c02e64, 2026-10-16): the guaranteed-exact family E and the 610 / 810 rows.
Own implementation from contract.md (committed before computation). The participant's scripts are read but not executed.
Usage: python check65_exact.py <path to their ExactFamily20261016 directory (extracted from aa0c3d1)> <path to check60_rows.json>
Writes check65_results.json and check65_rows.json next to this file.
"""
import json, math, os, sys, time
from fractions import Fraction as Fr
from sympy import factorint, isprime, primefactors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from chars20261015 import row, ser, box_pairs, ceil_log2  # noqa: E402

PRIMES = [66529, 345601, 670849, 1740481, 3678481, 5843041, 8615161, 18490201, 53722321, 84525841]
SMALL_L = {1, 2, 3, 4, 6}
T0 = time.time()


def lam_formula(m):
    out = 1
    for q, e in factorint(m).items():
        if q == 2:
            v = 1 if e == 1 else (2 if e == 2 else 2 ** (e - 2))
        else:
            v = q ** (e - 1) * (q - 1)
        out = math.lcm(out, v)
    return out


def lam_brute(m):
    units = [r for r in range(1, m) if math.gcd(r, m) == 1]
    ph = len(units)
    for e in sorted(d for d in range(1, ph + 1) if ph % d == 0):
        if all(pow(r, e, m) == 1 for r in units):
            return e
    raise AssertionError(m)


def part1(theirs_dir):
    out = {}
    KMAX = 5000
    lam = {}
    for k in range(1, KMAX + 1):
        f = lam_formula(4 * k)
        b = lam_brute(4 * k)
        assert f == b, (k, f, b)
        lam[k] = f
        assert f % 2 == 0
    E_scan = [k for k in range(1, KMAX + 1) if lam[k] in SMALL_L]
    div60 = [k for k in range(1, 61) if 60 % k == 0]
    div126 = [k for k in range(1, 127) if 126 % k == 0]
    E = sorted(set(div60) | set(div126))
    out['lambda_two_ways_agree_up_to'] = KMAX
    out['E_scan'] = E_scan
    out['E_divisors_60_or_126'] = E
    out['E_scan_equals_divisors'] = E_scan == E
    out['E_size'] = len(E)
    # lambda | 4 <=> 4k | 240 ; lambda | 6 <=> 4k | 504, checked on the scan range
    out['lambda_div4_iff_4k_div240'] = all((4 % lam[k] == 0) == (240 % (4 * k) == 0) for k in lam)
    out['lambda_div6_iff_4k_div504'] = all((6 % lam[k] == 0) == (504 % (4 * k) == 0) for k in lam)
    # structure claims of >>65
    def v(n, q):
        c = 0
        while n % q == 0:
            n //= q
            c += 1
        return c
    out['div60_structure'] = all(set(primefactors(k)) <= {2, 3, 5} and v(k, 3) <= 1 and v(k, 5) <= 1 and v(4 * k, 2) <= 4 for k in div60) \
        and max(v(4 * k, 2) for k in div60) == 4
    out['div126_structure'] = all(set(primefactors(k)) <= {2, 3, 7} and v(k, 3) <= 2 and v(k, 7) <= 1 and v(4 * k, 2) <= 3 for k in div126) \
        and max(v(4 * k, 2) for k in div126) == 3 and max(v(k, 3) for k in div126) == 2
    # pair count
    pairs = sorted((a, k // a) for k in E for a in range(1, k + 1) if k % a == 0 and math.gcd(a, k // a) == 1)
    out['E_pairs'] = len(pairs)
    out['sum_2_pow_omega'] = sum(2 ** len(primefactors(k)) for k in E)
    # compare their lambda checks (k = 1..126)
    th = json.load(open(os.path.join(theirs_dir, 'exact-family-lambda-checks.json'), encoding='utf-8'))
    bad = [r for r in th if not (r['lambda'] == lam[r['k']] and r['in_E'] == (r['k'] in E))]
    out['their_lambda_checks'] = {'rows': len(th), 'k_range': [min(r['k'] for r in th), max(r['k'] for r in th)], 'mismatches': bad}
    return E, pairs, out


def dec(s):
    return Fr(s) if isinstance(s, str) else Fr(s)


def part2(E, theirs_dir, saved_path):
    rows = []
    for p in PRIMES:
        B, ps = box_pairs(p)
        t = time.time()
        for a, u in ps:
            r = row(p, a, u, with_orbits=False, with_g2=False)
            r['B'] = B
            r['k'] = a * u
            r['forced_E'] = a * u in E
            r['all_pairs_constant'] = r['complex_pairs_nonconstant'] == 0
            r['c'] = max(Fr(0), r['W'] - r['eta'])
            # 2e: E => all pairs constant => eta = 0 => c = W
            if r['forced_E']:
                assert r['L'] in SMALL_L and r['all_pairs_constant'] and r['eta'] == 0 and r['c'] == r['W'], (p, a, u)
            if r['all_pairs_constant']:
                assert r['eta'] == 0
            rows.append(r)
        assert time.time() - t < 120
    out = {'rows': len(rows)}
    # 2a: against my saved 10-15 rows
    saved = json.load(open(saved_path, encoding='utf-8'))
    sv = {(s['p'], s['a'], s['u']): s for s in saved if s['p'] in PRIMES}
    assert len(sv) == len(rows)
    flds = ['W', 'eta', 'C', 'ell_sum', 'V', 'L', 'lower', 'K', 'b']
    bad = []
    for r in rows:
        s = sv[(r['p'], r['a'], r['u'])]
        for f in flds:
            if dec(s[f]) != Fr(r[f]):
                bad.append((r['p'], r['a'], r['u'], f, s[f], str(r[f])))
    out['vs_saved_20261015'] = {'rows': len(sv), 'fields': flds, 'mismatches': bad}
    # 2b: against their log rows and baseline rows
    key = {(r['p'], r['a'], r['u']): r for r in rows}
    th = json.load(open(os.path.join(theirs_dir, 'exact-family-log-rows.json'), encoding='utf-8'))
    assert len(th) == len(rows) and {(t['p'], t['a'], t['u']) for t in th} == set(key)
    bad = []
    for t in th:
        r = key[(t['p'], t['a'], t['u'])]
        mine = {'m': r['m'], 'L': r['L'], 'W': r['W'], 'eta': r['eta'], 'c': r['c'],
                'forced_E': r['forced_E'], 'all_pairs_constant': r['all_pairs_constant']}
        for f, v in mine.items():
            tv = dec(t[f]) if f in ('eta', 'c') else t[f]
            if tv != v:
                bad.append((t['p'], t['a'], t['u'], f, t[f], ser(v)))
    out['vs_their_log_rows'] = {'rows': len(th), 'fields': 7, 'mismatches': bad}
    bl = json.load(open(os.path.join(theirs_dir, 'baseline-log-rows.json'), encoding='utf-8'))
    bmap = {'W': 'W', 'V': 'V', 'signed_complex': 'C', 'complex_triangle_lower': 'ell_sum', 'factor_only_W_lower': 'lower'}
    bad = []
    for t in bl:
        r = key[(t['p'], t['a'], t['u'])]
        for f, g in bmap.items():
            if dec(t[f]) != Fr(r[g]):
                bad.append((t['p'], t['a'], t['u'], f, t[f], str(r[g])))
    out['vs_their_baseline_rows'] = {'rows': len(bl), 'fields': list(bmap), 'mismatches': bad}
    # 2c / 2d
    pos = [r for r in rows if r['c'] > 0]
    out['positive_rows'] = len(pos)
    out['positive_eta_zero'] = sum(1 for r in pos if r['eta'] == 0)
    out['positive_forced_E'] = sum(1 for r in pos if r['forced_E'])
    sel = ('p', 'a', 'u', 'm', 'L', 'W', 'eta', 'c', 'forced_E', 'all_pairs_constant')
    out['positive_outside_E_eta_zero'] = [ser({k: r[k] for k in sel}) for r in pos if not r['forced_E'] and r['eta'] == 0]
    out['positive_eta_nonzero'] = [ser({k: r[k] for k in sel}) for r in pos if r['eta'] > 0]
    out['E_rows'] = sum(1 for r in rows if r['forced_E'])
    out['E_success_rows'] = sum(1 for r in rows if r['forced_E'] and r['W'] > 0)
    out['outside_E_rows'] = sum(1 for r in rows if not r['forced_E'])
    out['outside_E_eta_zero_rows'] = sum(1 for r in rows if not r['forced_E'] and r['eta'] == 0)
    out['successful_rows'] = sum(1 for r in rows if r['W'] > 0)
    out['successful_rows_eta_positive'] = sum(1 for r in rows if r['W'] > 0 and r['eta'] > 0)
    boxes = []
    for p in PRIMES:
        rs = [r for r in rows if r['p'] == p]
        boxes.append({'p': p, 'B': rs[0]['B'], 'pairs': len(rs), 'Wtot': sum(r['W'] for r in rs),
                      'positive_rows': sum(1 for r in rs if r['c'] > 0),
                      'positive_eta_zero': sum(1 for r in rs if r['c'] > 0 and r['eta'] == 0),
                      'E_pairs': sum(1 for r in rs if r['forced_E']),
                      'E_Wtot': sum(r['W'] for r in rs if r['forced_E']),
                      'outside_E_zero_eta_rows': sum(1 for r in rs if not r['forced_E'] and r['eta'] == 0),
                      'outside_E_zero_eta_Wtot': sum(r['W'] for r in rs if not r['forced_E'] and r['eta'] == 0),
                      'Cert': sum((r['c'] for r in rs), Fr(0)),
                      'Cert0': sum((r['c'] for r in rs if r['eta'] == 0), Fr(0)),
                      'Certplus': sum((r['c'] for r in rs if r['eta'] > 0), Fr(0)),
                      'eta_positive_rows': sum(1 for r in rs if r['eta'] > 0),
                      'eta_sum': sum((r['eta'] for r in rs), Fr(0))})
    out['boxes'] = boxes
    out['every_box_has_E_witness_in_budget'] = all(b['E_Wtot'] > 0 for b in boxes)
    res = json.load(open(os.path.join(theirs_dir, 'exact-family-results.json'), encoding='utf-8'))
    bad = []
    n = 0
    for tb, mb in zip(res['box_summary'], boxes):
        assert tb['p'] == mb['p']
        for f, tv in tb.items():
            n += 1
            mv = mb[f]
            if (dec(tv) if isinstance(tv, str) else tv) != mv:
                bad.append((tb['p'], f, tv, ser(mv)))
    out['vs_their_box_summary'] = {'items': n, 'mismatches': bad}
    tt = res['totals']
    strd = lambda lst: [{k: str(v) for k, v in d.items()} for d in lst]
    out['vs_their_totals'] = {
        'positive_rows': tt['positive_rows'] == out['positive_rows'],
        'positive_eta_zero': tt['positive_eta_zero'] == out['positive_eta_zero'],
        'E_rows': tt['E_rows'] == out['E_rows'],
        'E_success_rows': tt['E_success_rows'] == out['E_success_rows'],
        'positive_eta_nonzero': strd(tt['positive_eta_nonzero']) == strd(out['positive_eta_nonzero']),
        'outside_E_eta_zero_successes': strd(tt['outside_E_eta_zero_successes']) == strd(out['positive_outside_E_eta_zero']),
    }
    return rows, out


def part3(E, pairs, theirs_dir):
    rows = []
    for p in PRIMES:
        B = ceil_log2(p)
        for a, u in pairs:
            r = row(p, a, u, with_orbits=False, with_g2=False)
            N = a * p + u
            f = factorint(N)
            assert math.prod(q ** e for q, e in f.items()) == N and all(isprime(q) for q in f)
            assert r['L'] in SMALL_L
            assert r['complex_pairs_nonconstant'] == 0 and r['eta'] == 0, (p, a, u)  # 3b
            rows.append({'p': p, 'a': a, 'u': u, 'k': a * u, 'N': N, 'm': 4 * a * u, 'L': r['L'],
                         'factors': [[q, e] for q, e in sorted(f.items())], 'W': r['W'], 'witnesses': r['witnesses'],
                         'within_log_budget': a * u <= B, 'eta': r['eta'], 'C': r['C'],
                         'complex_pairs': (r['phi'] - r['n_real']) // 2})
    out = {'rows': len(rows), 'all_eta_zero_all_pairs_constant': True}
    th = json.load(open(os.path.join(theirs_dir, 'exact-family-extended-rows.json'), encoding='utf-8'))
    key = {(r['p'], r['a'], r['u']): r for r in rows}
    assert len(th) == len(rows) and {(t['p'], t['a'], t['u']) for t in th} == set(key)
    flds = ['k', 'N', 'm', 'L', 'factors', 'W', 'witnesses', 'within_log_budget']
    bad = [(t['p'], t['a'], t['u'], f, t[f], key[(t['p'], t['a'], t['u'])][f]) for t in th for f in flds if t[f] != key[(t['p'], t['a'], t['u'])][f]]
    out['vs_their_extended_rows'] = {'rows': len(th), 'fields': flds, 'mismatches': bad}
    out['complex_pairs_checked'] = sum(r['complex_pairs'] for r in rows)
    summ = []
    for p in PRIMES:
        rs = [r for r in rows if r['p'] == p]
        succ = [(r['a'], r['u']) for r in rs if r['W'] > 0]
        succ_log = [(r['a'], r['u']) for r in rs if r['W'] > 0 and r['within_log_budget']]
        summ.append({'p': p, 'pairs': len(rs), 'success_pairs': succ, 'success_pairs_in_log_budget': succ_log,
                     'W_E_total': sum(r['W'] for r in rs), 'W_E_log': sum(r['W'] for r in rs if r['within_log_budget'])})
    out['per_prime'] = summ
    out['H_Elog_all_pass'] = all(s['success_pairs_in_log_budget'] for s in summ)
    out['H_E_all_pass'] = all(s['success_pairs'] for s in summ)
    res = json.load(open(os.path.join(theirs_dir, 'exact-family-results.json'), encoding='utf-8'))
    bad = []
    for ts, ms in zip(res['full_E_summary'], summ):
        assert ts['p'] == ms['p']
        tset = sorted((x['a'], x['u']) for x in ts['successes'])
        if ts['pairs'] != ms['pairs'] or tset != sorted(ms['success_pairs']):
            bad.append((ts['p'], tset, ms['success_pairs']))
        for x in ts['successes']:
            if x != {kk: (key[(x['p'], x['a'], x['u'])][kk]) for kk in x}:
                bad.append(('success_row', x['p'], x['a'], x['u']))
    out['vs_their_full_E_summary'] = {'primes': len(res['full_E_summary']), 'mismatches': bad}
    out['their_H_Elog'] = res['H_Elog']
    out['their_H_E'] = res['H_E']
    return rows, out


def main():
    theirs, saved = sys.argv[1], sys.argv[2]
    E, pairs, o1 = part1(theirs)
    t1 = time.time() - T0
    rows2, o2 = part2(E, theirs, saved)
    t2 = time.time() - T0
    rows3, o3 = part3(E, pairs, theirs)
    t3 = time.time() - T0
    assert t3 < 600
    res = {'part1_E': o1, 'part2_log_rows': o2, 'part3_E_rows': o3,
           'seconds': {'part1': round(t1, 1), 'part1_2': round(t2, 1), 'total': round(t3, 1)}}
    json.dump(ser(res), open(os.path.join(HERE, 'check65_results.json'), 'w', encoding='utf-8', newline='\n'), ensure_ascii=False, indent=1)
    keep = ('p', 'a', 'u', 'k', 'm', 'L', 'W', 'eta', 'c', 'C', 'ell_sum', 'V', 'forced_E', 'all_pairs_constant', 'complex_pairs_nonconstant')
    json.dump({'log_rows': [ser({k: r[k] for k in keep}) for r in rows2], 'E_rows': [ser(r) for r in rows3]},
              open(os.path.join(HERE, 'check65_rows.json'), 'w', encoding='utf-8', newline='\n'), ensure_ascii=False)
    sys.stdout.reconfigure(encoding='utf-8')
    print(json.dumps(ser(res), ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
