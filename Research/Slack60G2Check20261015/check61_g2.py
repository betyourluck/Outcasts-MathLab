"""Independent check of >>61 (ジェミー, Outcasts 80c02e64, 2026-10-14): square-subgroup decomposition
R_i = V_i + A_i - B_i = W_{G^2,i}/|G^2|, the 10-prime / 610-pair table, the span criterion, the (1,14) coset picture,
the |G^2| distribution for p = 84525841 and the bundled identity of §3.
Own implementation (chars20261015.py); >>61's embedded code is read, not executed.
Usage: python -I check61_g2.py   (reads check60_rows.json? no: recomputes, so that this file is self-contained)
"""
import json, time, math
from collections import Counter
from fractions import Fraction as Fr
from pathlib import Path
from chars20261015 import row, box_pairs, ser

T0 = time.time()
PRIMES = [66529, 345601, 670849, 1740481, 3678481, 5843041, 8615161, 18490201, 53722321, 84525841]
# >>61 table: rows, Wtot, fail, stage1 (W_G2=0), stage2 (W_G2>0, W=0), all-real-zero among stage2, sum(A-B), sum(C)
TABLE = {
    66529: (43, 4, 40, 33, 7, 4, Fr(-13, 8), Fr(-5, 4)),
    345601: (49, 5, 45, 36, 9, 9, Fr(-7, 4), Fr(-857, 180)),
    670849: (53, 4, 50, 40, 10, 10, Fr(1, 4), Fr(-151, 45)),
    1740481: (57, 6, 53, 42, 11, 9, Fr(-85, 12), Fr(-1013, 180)),
    3678481: (61, 7, 56, 46, 10, 7, Fr(47, 60), Fr(-4)),
    5843041: (63, 4, 60, 48, 12, 9, Fr(-689, 120), Fr(-581, 60)),
    8615161: (67, 14, 57, 43, 14, 8, Fr(-17, 72), Fr(-6239, 990)),
    18490201: (69, 3, 66, 53, 13, 7, Fr(-259, 120), Fr(-2881, 660)),
    53722321: (73, 8, 67, 58, 9, 8, Fr(-127, 120), Fr(-13, 5)),
    84525841: (75, 3, 72, 59, 13, 10, Fr(-2309, 720), Fr(-14437, 3960)),
}
COLS = ['rows', 'Wtot', 'fail', 'stage1_fail', 'stage2_fail', 'stage2_all_real_zero', 'sum_A_minus_B', 'sum_C']


def main():
    out = {'contract': 'contract.md', 'table': [], 'identities': {}, 'seconds': None}
    allrows = []
    for p in PRIMES:
        t0 = time.time()
        B, pairs = box_pairs(p)
        rows = [row(p, a, u, with_orbits=False) for a, u in pairs]
        assert time.time() - t0 < 120
        allrows += rows
        stage1 = [r for r in rows if r['W'] == 0 and r['W_G2'] == 0]
        stage2 = [r for r in rows if r['W'] == 0 and r['W_G2'] > 0]
        mine = {'rows': len(rows), 'Wtot': sum(r['W'] for r in rows), 'fail': sum(1 for r in rows if r['W'] == 0),
                'stage1_fail': len(stage1), 'stage2_fail': len(stage2),
                'stage2_all_real_zero': sum(1 for r in stage2 if r['n_real_nonprincipal_zero'] == r['n_real'] - 1),
                'sum_A_minus_B': sum((r['A'] - r['B'] for r in rows), Fr(0)),
                'sum_C': sum((r['C'] for r in rows), Fr(0))}
        # >>61 counts the (1,14) type by equal coset counts; check the two criteria agree on every row
        agree = all((r['n_real_nonprincipal_zero'] == r['n_real'] - 1) == r['all_cosets_equal'] for r in rows)
        claim = dict(zip(COLS, TABLE[p]))
        match = {k: mine[k] == claim[k] for k in COLS}
        out['table'].append({'p': p, 'B': B, 'mine': mine, 'claim': claim, 'match': match, 'all_match': all(match.values()),
                             'stage1_all_C_zero': all(r['C'] == 0 for r in stage1),
                             'stage2_all_C_negative': all(r['C'] < 0 for r in stage2),
                             'R_identity_all_rows': all(r['R_identity'] for r in rows),
                             'C_identity_all_rows': all(r['C_identity'] for r in rows),
                             'span_identity_all_rows': all(r['span_identity'] for r in rows),
                             'bundled_identity_all_rows': all(r['bundled_identity'] for r in rows),
                             'type_criteria_agree': agree,
                             'G2_size_distribution': dict(sorted(Counter(r['G2_size'] for r in rows).items())),
                             'seconds': round(time.time() - t0, 2)})
    tot = {k: sum(t['mine'][k] for t in out['table']) for k in ('rows', 'fail', 'stage1_fail', 'stage2_fail', 'stage2_all_real_zero')}
    out['totals'] = dict(tot, stage1_share=str(Fr(tot['stage1_fail'], tot['fail'])),
                         stage1_pct=round(100 * tot['stage1_fail'] / tot['fail'], 1),
                         stage2_pct=round(100 * tot['stage2_fail'] / tot['fail'], 1),
                         claim={'rows': 610, 'fail': 566, 'stage1_fail': 458, 'stage2_fail': 108, 'stage2_all_real_zero': 81,
                                'stage1_pct': 80.9, 'stage2_pct': 19.1})
    out['totals']['match'] = all(out['totals']['claim'][k] == (tot[k] if k in tot else out['totals'][k])
                                 for k in out['totals']['claim'])
    out['identities'] = {
        'rows': len(allrows),
        'R_eq_WG2_over_G2size': all(r['R_identity'] for r in allrows),
        'C_eq_W_minus_R': all(r['C_identity'] for r in allrows),
        'R_nonnegative': all(r['R'] >= 0 for r in allrows),
        'B_le_V_plus_A': all(r['B'] <= r['V'] + r['A'] for r in allrows),
        'span_criterion': all(r['span_identity'] for r in allrows),
        'bundled_identity': all(r['bundled_identity'] for r in allrows),
        'stage1_implies_C_zero': all(r['C'] == 0 for r in allrows if r['W'] == 0 and r['W_G2'] == 0),
        'failure_implies_C_le_0': all(r['C'] <= 0 for r in allrows if r['W'] == 0),
        'C_negative_iff_stage2': all((r['C'] < 0) == (r['W'] == 0 and r['W_G2'] > 0) for r in allrows if r['W'] == 0),
        'type_criteria_agree': all((r['n_real_nonprincipal_zero'] == r['n_real'] - 1) == r['all_cosets_equal'] for r in allrows)}
    # (1,14) picture
    r = next(rw for rw in allrows if (rw['p'], rw['a'], rw['u']) == (18490201, 1, 14))
    out['row_18490201_1_14'] = {'phi': r['phi'], 'G2': r['G2'], 'negG2': r['negG2'], 'G2_size': r['G2_size'],
                                'tau_odd': r['tau_odd'], 'coset_counts': r['coset_counts'], 'negG2_counts': r['negG2_counts'],
                                'W_G2': r['W_G2'], 'R': r['R'], 'C': r['C'],
                                'claim_match': (r['phi'] == 24 and r['G2'] == [1, 9, 25] and r['negG2'] == [31, 47, 55]
                                                and r['tau_odd'] == 32 and r['coset_counts'] == [4] * 8
                                                and r['negG2_counts'] == {'31': 1, '47': 3, '55': 0}
                                                and r['R'] == Fr(4, 3) and r['C'] == Fr(-4, 3))}
    dist = next(t for t in out['table'] if t['p'] == 84525841)['G2_size_distribution']
    le4 = sum(v for k, v in dist.items() if k <= 4)
    out['G2_distribution_84525841'] = {'distribution': dist, 'le4': le4, 'pct_le4': round(100 * le4 / 75),
                                       'claim_match': dist.get(1) == 9 and dist.get(2) == 16 and dist.get(3) == 16 and dist.get(4) == 10 and le4 == 51}
    out['seconds'] = round(time.time() - T0, 2)
    Path('check61_results.json').write_text(json.dumps(ser(out), indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(ser(out), indent=1, ensure_ascii=False))


if __name__ == '__main__':
    main()
