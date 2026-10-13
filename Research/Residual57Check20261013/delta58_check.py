"""Check of >>58 (ParallaxSol, 2026-10-13) on the rows produced by check57_cert.py:
delta_i = C_i - sum_j ell_ij >= 0 and V_i - K_i = min{V_i, (W_i - delta_i)_+} on every row; box sums of delta."""
import json
from fractions import Fraction as Fr
rows = json.load(open('check57_rows.json', encoding='utf-8'))
out = {}
for p in (66529, 18490201, 1201):
    rs = [r for r in rows if r['p'] == p]
    tot = Fr(0); npos = 0; bad = 0; nonconst = 0; max_d = (Fr(0), None)
    for r in rs:
        V, K, W, C, ell = (Fr(r[k]) for k in ('V', 'K', 'W', 'C', 'ell_sum'))
        d = C - ell
        assert d >= 0
        if d > 0: npos += 1
        if d > max_d[0]: max_d = (d, (r['a'], r['u']))
        nonconst += r['complex_pairs_nonconstant']
        if V - K != min(V, max(Fr(0), W - d)): bad += 1
        tot += d
    out[str(p)] = {'rows': len(rs), 'sum_delta58': str(tot), 'rows_with_positive_delta58': npos,
                   'nonconstant_complex_pairs_total': nonconst, 'max_delta58': [str(max_d[0]), max_d[1]],
                   'identity_V_minus_K_eq_min_V_Wminusdelta_fails': bad,
                   'rows_with_W_gt_delta58': sum(1 for r in rs if Fr(r['W']) > Fr(r['C']) - Fr(r['ell_sum']))}
json.dump(out, open('delta58_results.json', 'w', encoding='utf-8'), indent=1)
print(json.dumps(out, indent=1))
