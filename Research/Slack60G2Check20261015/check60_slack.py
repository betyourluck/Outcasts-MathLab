"""Independent check of >>60 (Outcasts 80c02e64, 2026-10-14; PR #4 4e514ac Research/CertificateSlack20261014/):
8 fixed boxes (498 pairs), per-row 16 fields, slack identities, 113-row summaries, complete power-conjugacy orbits.
Own implementation (chars20261015.py); the participant's JSON is read only after the computation, for comparison.
Usage: python -I check60_slack.py <dir with slack-certificate-results.json, slack-certificate-rows-summary.json, selected-full-row.json>
"""
import sys, json, time, math
from collections import Counter
from fractions import Fraction as Fr
from pathlib import Path
from sympy import isprime
from chars20261015 import row, box_pairs, ser

T0 = time.time()
INPUTS = [84525841, 8615161, 53722321, 345601, 670849, 1740481, 3678481, 5843041]
OLD = [66529, 18490201]
FIELDS = {'tau_odd': 'tau_odd', 'phi': 'phi', 'root_order': 'L', 'U': 'U', 't': 't', 'V': 'V',
          'other_even_real': 'A', 'other_odd_real': 'B', 'signed_complex': 'C', 'complex_triangle_lower': 'ell_sum',
          'factor_only_K': 'K', 'factor_only_W_lower': 'lower', 'W': 'W', 'witnesses': 'witnesses',
          'approximation_slack': 'eta', 'capped_certificate_gap': 'b'}


def same(theirs, mine):
    if isinstance(mine, Fr):
        return Fr(str(theirs)) == mine
    return theirs == mine


def box_summary(p, rows):
    return {'p': p, 'pairs': len(rows),
            'sum_V': sum((r['V'] for r in rows), Fr(0)), 'sum_K': sum((r['K'] for r in rows), Fr(0)),
            'gap': sum((r['b'] for r in rows), Fr(0)), 'Wtot': sum(r['W'] for r in rows),
            'successful_pairs': [[r['a'], r['u']] for r in rows if r['W'] > 0],
            'sum_slack': sum((r['eta'] for r in rows), Fr(0)),
            'positive_slack_rows': sum(1 for r in rows if r['eta'] > 0),
            'nonconstant_pairs': sum(r['complex_pairs_nonconstant'] for r in rows),
            'W_above_slack_pairs': [[r['a'], r['u']] for r in rows if r['W'] > r['eta']],
            'positive_certificate_pairs': [[r['a'], r['u']] for r in rows if r['b'] > 0],
            'sum_A_minus_B': sum((r['A'] - r['B'] for r in rows), Fr(0)),
            'sum_C': sum((r['C'] for r in rows), Fr(0))}


def orbit_key(o):
    return (o['order'], o['size'], o['sign'], o['trace'], tuple(sorted(o['D'])))


def main():
    their_dir = Path(sys.argv[1])
    rows_all = {}
    boxes = []
    for p in INPUTS + OLD:
        assert isprime(p) and p % 840 in {1, 121, 169, 289, 361, 529}
        t0 = time.time()
        B, pairs = box_pairs(p)
        rows = [row(p, a, u) for a, u in pairs]
        assert time.time() - t0 < 120, (p, time.time() - t0)
        rows_all[p] = rows
        s = box_summary(p, rows)
        s['B'] = B
        s['seconds'] = round(time.time() - t0, 2)
        boxes.append(s)
    r1201 = row(1201, 3, 1)
    rows_all[1201] = [r1201]
    boxes.append(dict(box_summary(1201, [r1201]), B=None))
    rows_new = [r for p in INPUTS for r in rows_all[p]]
    assert len(rows_new) == 498
    rows_113 = [r for p in OLD for r in rows_all[p]] + [r1201]
    assert len(rows_113) == 113
    # ---- comparison with the participant (read only now) ----
    theirs = json.loads(their_dir.joinpath('slack-certificate-results.json').read_text(encoding='utf-8'))
    their_rows = json.loads(their_dir.joinpath('slack-certificate-rows-summary.json').read_text(encoding='utf-8'))
    selected = json.loads(their_dir.joinpath('selected-full-row.json').read_text(encoding='utf-8'))
    their_by = {(r['p'], r['a'], r['u']): r for r in their_rows}
    assert len(their_rows) == 498 and len(their_by) == 498
    out = {'contract': 'contract.md', 'fixed_inputs': INPUTS, 'boxes': [], 'row_mismatches': [], 'orbit_examples': {}}
    thread_table = {84525841: (75, Fr(2), 3), 8615161: (67, Fr(85, 12), 14), 53722321: (73, Fr(7, 3), 8),
                    345601: (49, Fr(17, 6), 5), 670849: (53, Fr(1, 4), 4), 1740481: (57, Fr(11, 3), 6),
                    3678481: (61, Fr(23, 4), 7), 5843041: (63, Fr(3), 4)}
    for s in boxes:
        p = s['p']
        rec = dict(s)
        if p in INPUTS:
            th = next(b for b in theirs['completed_boxes'] if b['p'] == p)
            cmp = {k: same(th[k], s[k]) for k in ['pairs', 'sum_V', 'sum_K', 'gap', 'Wtot', 'successful_pairs',
                                                     'sum_slack', 'positive_certificate_pairs']}
            cmp['B'] = th['B'] == s['B']
            tt = thread_table[p]
            rec['thread_table_match'] = (s['pairs'], s['gap'], s['Wtot']) == tt
            rec['json_match'] = cmp
            rec['all_match'] = all(cmp.values()) and rec['thread_table_match']
        else:
            th = next(b for b in theirs['prior113_slack_identity_checks'] if b['p'] == p)
            cmp = {k: same(th[k], s[k]) for k in ['sum_slack', 'positive_slack_rows', 'nonconstant_pairs', 'W_above_slack_pairs']}
            cmp['rows'] = th['rows'] == s['pairs']
            rec['json_match'] = cmp
            rec['all_match'] = all(cmp.values())
        out['boxes'].append(rec)
    for r in rows_new:
        tr = their_by[(r['p'], r['a'], r['u'])]
        bad = {k: (tr[k], ser(r[v])) for k, v in FIELDS.items() if not same(tr[k], r[v])}
        if bad:
            out['row_mismatches'].append({'p': r['p'], 'a': r['a'], 'u': r['u'], 'fields': ser(bad)})
    out['rows_compared'] = len(rows_new)
    out['fields_per_row'] = len(FIELDS)
    out['all_rows_match'] = not out['row_mismatches']
    # identities over all 611 rows are asserted inside row(); record counts
    allrows = rows_new + rows_113
    out['identity_checks'] = {
        'rows': len(allrows),
        'slack_identity_b_eq_min_V_Wminus_eta_plus': True,
        'positivity_criterion': True,
        'eta_nonnegative': True,
        'phi_V_ge_2W': True,
        'V_eq_2nminus_over_phi': True,
        'orbit_total_eq_C': all(r['orbit_C'] == r['C'] for r in allrows),
        'orbit_gap_eq_min_V_W': all(r['orbit_b'] == min(r['V'], r['W']) for r in allrows),
        'orbits_total': sum(len(r['orbits']) for r in allrows),
        'rows_where_orbit_grouping_turns_zero_certificate_positive':
            sum(1 for r in allrows if r['b'] == 0 and r['orbit_b'] > 0),
        'rows_with_W_positive': sum(1 for r in allrows if r['W'] > 0),
        'rows_with_b_positive': sum(1 for r in allrows if r['b'] > 0)}
    # orbit examples
    for (p, a, u), claim in {(84525841, 1, 19): {'traces': [0, 6, 0, 4, 6, 0, 4, 0], 'sum': 20, 'C': Fr(5, 9),
                                                  'old_eta': Fr(14, 9), 'old_b': Fr(0), 'orbit_b': Fr(2, 9), 'V': Fr(2, 9), 'W': 1},
                             (18490201, 1, 14): {'traces': [6, -2, 0, 0, 0, 0, -18, -18], 'sum': -32, 'C': Fr(-4, 3),
                                                 'old_eta': Fr(0), 'old_b': Fr(0), 'orbit_b': Fr(0), 'V': Fr(4, 3), 'W': 0}}.items():
        r = next(rw for rw in rows_all[p] if (rw['a'], rw['u']) == (a, u))
        mine = {'m': r['m'], 'orbits': len(r['orbits']), 'traces_sorted': sorted(o['trace'] for o in r['orbits']),
                'sum': sum(o['trace'] for o in r['orbits']), 'C': r['C'], 'old_eta': r['eta'], 'old_b': r['b'],
                'orbit_b': r['orbit_b'], 'V': r['V'], 'W': r['W'],
                'orbit_table': [{k: o[k] for k in ('order', 'size', 'sign', 'trace', 'D')} for o in r['orbits']]}
        th = next(o for o in theirs['orbit_checks'] if (o['p'], o['a'], o['u']) == (p, a, u))
        their_keys = Counter((o['order'], o['size'], o['sign'], o['trace_numerator'], tuple(sorted(o['phase_coefficients'])))
                             for o in th['orbits'])
        my_keys = Counter(orbit_key(o) for o in r['orbits'])
        mine['thread_match'] = (mine['traces_sorted'] == sorted(claim['traces']) and mine['sum'] == claim['sum']
                                and all(mine[k] == claim[k] for k in ('C', 'old_eta', 'old_b', 'orbit_b', 'V', 'W')))
        mine['json_match'] = (my_keys == their_keys and Fr(th['signed_complex']) == r['C'] and Fr(th['old_slack']) == r['eta']
                              and Fr(th['old_gap']) == r['b'] and Fr(th['orbit_gap']) == r['orbit_b'] and th['W'] == r['W']
                              and Fr(th['V']) == r['V'] and Fr(th['orbit_slack']) == 0)
        out['orbit_examples']['%d_%d_%d' % (p, a, u)] = mine
    # selected full row (84525841, 1, 19) extra fields
    sel = selected[0]
    r = next(rw for rw in rows_all[84525841] if (rw['a'], rw['u']) == (1, 19))
    out['selected_full_row_match'] = {k: same(sel[k], r[v]) for k, v in FIELDS.items()}
    out['selected_full_row_match']['complex_pairs_count'] = len(sel['complex_pairs']) == (r['phi'] - r['n_real']) // 2
    out['selected_full_row_match']['real_terms_count'] = len(sel['real_terms']) == r['n_real']
    out['participant_flags'] = {k: theirs[k] for k in ('stopped_on_first_zero_gap', 'participant_script_executed',
                                                        'new_general_existence_proved', 'Hlog_or_ESC_refuted', 'seconds')}
    out['seconds'] = round(time.time() - T0, 2)
    Path('check60_results.json').write_text(json.dumps(ser(out), indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
    slim = [{k: v for k, v in rw.items() if k not in ('orbits', 'G2', 'negG2', 'coset_counts', 'negG2_counts')}
            for rw in allrows]
    Path('check60_rows.json').write_text(json.dumps(ser(slim), indent=0, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(ser({k: v for k, v in out.items() if k not in ('row_mismatches',)}), indent=1, ensure_ascii=False))
    print('row mismatches:', len(out['row_mismatches']))
    for mm in out['row_mismatches'][:10]:
        print(mm)


if __name__ == '__main__':
    main()
