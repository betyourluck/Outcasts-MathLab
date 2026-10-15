"""Check of >>62 (境界測量士, Outcasts 80c02e64, 2026-10-15; PR #4 ddf9066 Research/SquareCoset20261015/) on the rows
already computed under contract.md (check60_rows.json, 610 rows of the 10 boxes). Not pre-registered: >>62 arrived after the
contract; this script only re-counts the saved rows with K = <odd prime factors of N> (participant's K), which contract.md
did not define (it used K*H = <G^2, {q}>).
Usage: python check62_k.py   (from this directory)
"""
import json, math
from collections import Counter
from pathlib import Path

rows = [r for r in json.load(open('check60_rows.json', encoding='utf-8')) if r['p'] != 1201]
assert len(rows) == 610


def units(m):
    return [x for x in range(1, m) if math.gcd(x, m) == 1]


def generated(m, gens, start):
    S = set(start)
    changed = True
    while changed:
        changed = False
        for q in gens:
            new = {g * q % m for g in S}
            if not new <= S:
                S |= new
                changed = True
    return S


for r in rows:
    m = r['m']
    qs = [int(q) % m for q in r['factors'] if int(q) != 2]
    U = units(m)
    H = {x * x % m for x in U}
    K = generated(m, qs, [1])
    KH = generated(m, qs, H)
    r['K'] = sorted(K)
    r['minus1_in_K'] = (m - 1) in K
    r['KH_eq_G'] = len(KH) == len(U)
    r['H_odd'] = len(H) % 2 == 1
    assert ((m - 1) in KH) == (r['W_G2'] > 0)          # first stage <=> -1 in K*H
    assert r['G2_size'] == len(H)

fail = [r for r in rows if r['W'] == 0]
quad = [r for r in fail if r['W_G2'] == 0]
group = [r for r in fail if r['W_G2'] > 0 and not r['minus1_in_K']]
expo = [r for r in fail if r['W_G2'] > 0 and r['minus1_in_K']]
per = {}
for r in fail:
    if r['W_G2'] > 0:
        per.setdefault(r['p'], [0, 0])
        per[r['p']][0 if not r['minus1_in_K'] else 1] += 1
claim_per = {66529: [1, 6], 345601: [0, 9], 670849: [0, 10], 1740481: [0, 11], 3678481: [2, 8], 5843041: [1, 11],
             8615161: [1, 13], 18490201: [2, 11], 53722321: [0, 9], 84525841: [0, 13]}
uniform = [r for r in fail if r['all_cosets_equal']]
odd_rows = [r for r in rows if r['H_odd']]
r10 = next(r for r in rows if (r['p'], r['a'], r['u']) == (66529, 10, 1))
r9 = next(r for r in rows if (r['p'], r['a'], r['u']) == (66529, 9, 1))
box = [r for r in rows if r['p'] == 5843041]
oddH = [r for r in box if r['H_odd']]
out = {
    'basis': 'check60_rows.json (610 rows computed under contract.md); K = <odd prime factors of N> added here',
    'failure_split': {'total': len(fail), 'quadratic_W_H_zero': len(quad), 'W_H_pos_minus1_not_in_K': len(group),
                      'W_H_pos_minus1_in_K_W_zero': len(expo), 'claim': [566, 458, 7, 101],
                      'match': [len(fail), len(quad), len(group), len(expo)] == [566, 458, 7, 101]},
    'group_obstruction_rows': [{'p': r['p'], 'a': r['a'], 'u': r['u'], 'm': r['m'], 'H_size': r['G2_size'],
                                'W_G2': r['W_G2'], 'K': r['K']} for r in group],
    'per_prime_second_stage': {str(p): {'mine': v, 'claim': claim_per[p], 'match': v == claim_per[p]} for p, v in sorted(per.items())},
    'per_prime_all_match': all(per[p] == claim_per[p] for p in claim_per),
    'odd_H_lemma': {'odd_H_rows': len(odd_rows),
                    'violations_W_H_pos_but_minus1_not_in_K': sum(1 for r in odd_rows if r['W_G2'] > 0 and not r['minus1_in_K']),
                    'odd_H_second_stage_failures': sum(1 for r in fail if r['H_odd'] and r['W_G2'] > 0),
                    'all_of_them_have_minus1_in_K': all(r['minus1_in_K'] for r in fail if r['H_odd'] and r['W_G2'] > 0),
                    'claim_second_stage_failures': 53},
    'uniform_coset_lemma': {'uniform_failures': len(uniform), 'all_KH_eq_G': all(r['KH_eq_G'] for r in uniform),
                            'all_minus1_in_K': all(r['minus1_in_K'] for r in uniform), 'claim': 81},
    'KH_eq_G_rows': {'count': sum(r['KH_eq_G'] for r in rows), 'claim': 125,
                     'implies_minus1_in_K_on_all_rows': all(r['minus1_in_K'] for r in rows if r['KH_eq_G'])},
    'group_obstruction_rows_all_even_H': all(not r['H_odd'] for r in group),
    'example_66529_10_1': {'factors': r10['factors'], 'm': r10['m'], 'H_size': r10['G2_size'], 'W_G2': r10['W_G2'], 'K': r10['K'],
                           'minus1_in_K': r10['minus1_in_K'], 'R': r10['R'], 'C': r10['C'], 'W': r10['W'],
                           'claim_match': r10['K'] == [1, 11, 21, 31] and r10['W_G2'] == 2 and r10['R'] == '1' and r10['C'] == '-1' and r10['W'] == 0},
    'example_66529_9_1': {'factors': r9['factors'], 'm': r9['m'], 'H_size': r9['G2_size'], 'W_G2': r9['W_G2'], 'W': r9['W'],
                          'minus1_in_K': r9['minus1_in_K'], 'pow_83_3_mod_36': pow(83, 3, 36), '83_cubed_divides_N': r9['N'] % 83 ** 3 == 0,
                          'claim_match': r9['G2_size'] == 3 and pow(83, 3, 36) == 35 and r9['N'] % 83 ** 3 != 0 and r9['W'] == 0 and r9['minus1_in_K']},
    'refutation_5843041': {'pairs': len(box), 'odd_H_pairs': len(oddH), 'odd_H_pairs_all_W_zero': all(r['W'] == 0 for r in oddH),
                           'successes': [{'a': r['a'], 'u': r['u'], 'H_size': r['G2_size'], 'witnesses': r['witnesses']} for r in box if r['W'] > 0],
                           'claim_match': len(box) == 63 and len(oddH) == 35 and all(r['W'] == 0 for r in oddH)
                                          and [(r['a'], r['u'], r['G2_size']) for r in box if r['W'] > 0] == [(1, 5, 2), (1, 17, 8), (5, 1, 2)]
                                          and [r['witnesses'] for r in box if r['W'] > 0] == [[7359], [127023], [279, 11439]]},
    'control_m20_N11': {'minus_H': sorted((-h) % 20 for h in {x * x % 20 for x in units(20)}), 'K': sorted(generated(20, [11], [1]))},
}
Path('check62_results.json').write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
print(json.dumps(out, indent=1, ensure_ascii=False))
