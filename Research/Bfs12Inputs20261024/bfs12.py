"""Square-source BFS of >>86 section 4 on 12 fixed inputs (>>89 section 4, >>90) -- see contract.md.

Own implementation. The arithmetic is the one of our check8688.py, which reproduced the >>86 BFS
order line by line (53 sources, 39 vertices for p = 2521). No peer code is imported or executed.

Usage:
  python -I bfs12.py selftest   # inputs outside the 12, plus the known BFS order of 2521; writes selftest_results.json
  python -I bfs12.py cohort     # the 12 inputs; refuses to run unless contract.md and bfs12.py are committed and unchanged
stdout is the log (redirect it to selftest.log / run.log).
"""
import sys, os, json, math, time, subprocess
from collections import deque
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
INPUTS = [2521, 66529, 345601, 670849, 1740481, 3678481, 5843041, 8615161, 18490201, 38296729, 53722321, 84525841]
SELFTEST = [1201, 100000609, 100687129]
CLASSES = (1, 121, 169, 289, 361, 529)
CAP_SOURCES = 500_000   # per input and per run (BFS / sweep); beyond it the run stops with type 3
CAP_SECONDS = 1800      # per input and per run
CAP_TOTAL = 4 * 3600    # whole process; inputs not reached are type 3
T0 = time.monotonic()
LOG = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    LOG.append(s)
    print(s, flush=True)


# ---------------------------------------------------------------- arithmetic

def trial_factor(n):
    f, q = {}, 2
    while q * q <= n:
        while n % q == 0:
            f[q] = f.get(q, 0) + 1
            n //= q
        q += 1 if q == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def factor(n):
    a = trial_factor(n)
    b = {int(k): int(v) for k, v in sympy.factorint(n).items()}
    assert a == b, (n, a, b)
    assert math.prod(r ** e for r, e in a.items()) == n
    return dict(sorted(a.items()))


def chi(a, p):
    """Legendre symbol (a|p) by Euler's criterion, p odd prime."""
    k = pow(a % p, (p - 1) // 2, p)
    assert k in (0, 1, p - 1)
    return 0 if k == 0 else (1 if k == 1 else -1)


def divisors(f):
    ds = [1]
    for r, e in f.items():
        ds = [d * r ** j for d in ds for j in range(e + 1)]
    return sorted(ds)


def val(n, r):
    assert n != 0
    n, e = abs(n), 0
    while n % r == 0:
        n //= r
        e += 1
    return e


def least_nonresidue(p):
    for q in sympy.primerange(2, p):
        if chi(q, p) == -1:
            return q
    raise AssertionError(p)


def sigma_of(q):
    return 1 if q % 4 == 3 else 3


def allowed_taus(p, q):
    A = sigma_of(q) * q
    taus, t = [], 1
    while A * t * t < 3 * p:
        taus.append(t)
        t += 2
    return taus


# ---------------------------------------------------------------- one square source

def test_source(p, q, tau):
    """x = (p + sigma q tau^2)/4; Type II d = -x, Type I d = -p x (mod c) over all d | x^2, two routes."""
    A = sigma_of(q) * q
    c = A * tau * tau
    assert tau % 2 == 1 and c < 3 * p and (p + c) % 4 == 0
    x = (p + c) // 4
    assert 0 < x < p and math.gcd(c, x) == 1
    f = factor(x)
    f2 = {r: 2 * e for r, e in f.items()}
    tgt = {'II': (-x) % c, 'I': (-p * x) % c}
    # route 1: enumerate the divisors of x^2
    ds = divisors(f2)
    hits = {br: [d for d in ds if d % c == t] for br, t in tgt.items()}
    # route 2: the set of residues of the divisors of x^2 mod c, built prime power by prime power
    R = {1 % c}
    for r, e in f2.items():
        pw = [pow(r, j, c) for j in range(e + 1)]
        R = {v * w % c for v in R for w in pw}
    for br in tgt:
        assert (tgt[br] in R) == bool(hits[br]), (p, q, tau, br)
    cert = None
    for br in ('II', 'I'):  # the order of the >>86 implementation (Type II first, smallest d)
        if hits[br]:
            d = hits[br][0]
            if br == 'II':
                assert p * (d + x) % c == 0 and p * (x * x // d + x) % c == 0
                y, z = p * (d + x) // c, p * (x * x // d + x) // c
            else:
                assert (d + p * x) % c == 0 and (p * p * (x * x // d) + p * x) % c == 0
                y, z = (d + p * x) // c, (p * p * (x * x // d) + p * x) // c
            assert y > 0 and z > 0 and 4 * x * y * z == p * (x * y + x * z + y * z)
            cert = {'branch': br, 'c': c, 'x': x, 'd': d, 'y': y, 'z': z}
            break
    nonres = [r for r in f if chi(r, p) == -1]
    assert all(r < p for r in nonres)
    return {'q': q, 'tau': tau, 'c': c, 'x': x, 'f': f, 'ok': {br: bool(hits[br]) for br in tgt},
            'n_hits': {br: len(hits[br]) for br in tgt}, 'cert': cert, 'nonres': nonres}


class Cache:
    def __init__(self, p):
        self.p, self.d, self.byq = p, {}, {}

    def get(self, q, tau):
        k = (q, tau)
        if k not in self.d:
            self.d[k] = test_source(self.p, q, tau)
            self.byq.setdefault(q, []).append(self.d[k])
        return self.d[k]

    def tested_at(self, q):
        return self.byq.get(q, [])


def out_of_budget(n_done, t_start):
    return n_done >= CAP_SOURCES or time.monotonic() - t_start > CAP_SECONDS or time.monotonic() - T0 > CAP_TOTAL


# ---------------------------------------------------------------- the procedure of >>86 section 4

def bfs(p, cache, keep_order=10_000):
    t0 = time.monotonic()
    ell = least_nonresidue(p)
    depth, parent = {ell: 0}, {}
    queue, order, processed = deque([ell]), [], []
    while queue:
        q = queue.popleft()
        processed.append(q)
        for tau in allowed_taus(p, q):
            if out_of_budget(len(order), t0):
                return {'status': 3, 'ell': ell, 'checked': len(order), 'vertices': len(depth), 'processed': len(processed),
                        'depth_reached': depth[q], 'seconds': round(time.monotonic() - t0, 2)}
            s = cache.get(q, tau)
            order.append((q, tau))
            if s['cert']:
                path, r = [], q
                while r != ell:
                    pq, pt, v = parent[r]
                    path.append({'parent': pq, 'tau': pt, 'child': r, 'v_child': v})
                    r = pq
                path.reverse()
                assert len(path) == depth[q]
                return {'status': 1, 'ell': ell, 'D': depth[q], 'checked': len(order), 'vertices': len(depth),
                        'processed': len(processed), 'path': path, 'hit': {'q': q, 'tau': tau, **s['cert']},
                        'hit_ok_both': s['ok'], 'order': order if len(order) <= keep_order else None,
                        'vertex_list': sorted(depth) if len(depth) <= keep_order else None,
                        'seconds': round(time.monotonic() - t0, 2)}
            for r in sorted(s['nonres']):
                if r not in depth:
                    depth[r] = depth[q] + 1
                    parent[r] = (q, tau, s['f'][r])
                    queue.append(r)
    # the queue is exhausted: every vertex was processed with every allowed tau and failed
    C = sorted(depth)
    assert sorted(processed) == C
    for q in C:
        rows = cache.tested_at(q)
        assert sorted(s['tau'] for s in rows) == allowed_taus(p, q)
        assert all(s['cert'] is None for s in rows)
        assert all(r in depth for s in rows for r in s['nonres'])
    return {'status': 2, 'ell': ell, 'C': C, 'size': len(C), 'checked': len(order), 'max_depth': max(depth.values()),
            'seconds': round(time.monotonic() - t0, 2)}


# ---------------------------------------------------------------- measurements (not used for the judgments)

def vertex_status(p, cache, q, t0, counter):
    """Scan the allowed tau of q in increasing order up to the first success."""
    for tau in allowed_taus(p, q):
        if out_of_budget(counter[0], t0):
            return {'status': 'undetermined'}
        s = cache.get(q, tau)
        counter[0] += 1
        if s['cert']:
            return {'status': 'success', 'tau': tau, 'branch': s['cert']['branch']}
    return {'status': 'fail', 'n_tau': len(allowed_taus(p, q))}


def sweep(p, cache):
    t0 = time.monotonic()
    counter = [0]
    ell = least_nonresidue(p)
    root = [cache.get(ell, t) for t in allowed_taus(p, ell)]
    counter[0] += len(root)
    succ = [{'tau': s['tau'], 'II': s['ok']['II'], 'I': s['ok']['I']} for s in root if s['cert']]
    failing = [s for s in root if not s['cert']]
    children = sorted({r for s in failing for r in s['nonres']})
    status = {}
    for r in children:
        status[r] = vertex_status(p, cache, r, t0, counter)
    undetermined = [r for r in children if status[r]['status'] == 'undetermined']
    e1p = [{'tau': s['tau'], 'children': s['nonres'], 'successful_children': [r for r in s['nonres'] if status[r]['status'] == 'success']}
           for s in failing]
    any_child = any(status[r]['status'] == 'success' for r in children)
    if succ or any_child:
        e1 = True
    elif undetermined:
        e1 = None
    else:
        e1 = False
    return {'ell': ell, 'n_root': len(root), 'S0': len(succ), 'root_successes': succ,
            'root_II': sum(1 for s in root if s['ok']['II']), 'root_I': sum(1 for s in root if s['ok']['I']),
            'n_failing_root': len(failing), 'n_children': len(children),
            'n_children_success': sum(1 for r in children if status[r]['status'] == 'success'),
            'n_children_fail': sum(1 for r in children if status[r]['status'] == 'fail'),
            'n_children_undetermined': len(undetermined),
            'E1prime_holds': sum(1 for e in e1p if e['successful_children']),
            'E1prime_fails': [e['tau'] for e in e1p if not e['successful_children']
                              and all(status[r]['status'] == 'fail' for r in e['children'])],
            'E1_by_sweep': e1,
            'child_status': {str(r): status[r] for r in children} if len(children) <= 5000 else None,
            'sources': counter[0], 'seconds': round(time.monotonic() - t0, 2)}


def hensel(p, cache):
    """>>90: A r^(2k) <= 3p  =>  some allowed odd tau has r^k | x(tau). Refined: A (r^k - 2)^2 < 3p suffices
    (the odd one of b, r^k - b is at most r^k - 2). Brute force over all allowed tau of every tested vertex,
    for the odd primes r that divide some tested x at that vertex and have k_ref >= 2."""
    out = {'pairs': 0, 'k90_ge2': 0, 'kref_gt_k90': 0, 'beyond_kref': 0, 'examples_beyond': [], 'rows': []}
    for q in sorted(cache.byq):
        A = sigma_of(q) * q
        taus = allowed_taus(p, q)
        rs = sorted({r for s in cache.tested_at(q) for r in s['f'] if r != 2 and A * (r * r - 2) ** 2 < 3 * p})
        for r in rs:
            k90 = 0
            while A * r ** (2 * (k90 + 1)) <= 3 * p:
                k90 += 1
            kref = 0
            while A * (r ** (kref + 1) - 2) ** 2 < 3 * p:
                kref += 1
            vmax = max(val((p + A * t * t) // 4, r) for t in taus)
            assert kref >= k90 and kref >= 2, (q, r, k90, kref)
            assert vmax >= kref, ('Hensel claim fails', p, q, r, k90, kref, vmax)
            out['pairs'] += 1
            out['k90_ge2'] += k90 >= 2
            out['kref_gt_k90'] += kref > k90
            if vmax > kref:
                out['beyond_kref'] += 1
                if len(out['examples_beyond']) < 20:
                    tau = next(t for t in taus if val((p + A * t * t) // 4, r) == vmax)
                    out['examples_beyond'].append({'q': q, 'r': r, 'k90': k90, 'kref': kref, 'vmax': vmax, 'tau': tau})
            out['rows'].append((q, r, k90, kref, vmax))
    if len(out['rows']) > 5000:
        out['rows'] = None
    return out


# ---------------------------------------------------------------- one input

def run_input(p, with_sweep=True):
    assert sympy.isprime(p) and factor(p) == {p: 1}
    assert p % 840 in CLASSES and p % 8 == 1 and chi(2, p) == 1
    ell = least_nonresidue(p)
    cache = Cache(p)
    res = {'p': p, 'class': p % 840, 'B': (p - 1).bit_length(), 'ell': ell, 'sigma': sigma_of(ell),
           'n_root_tau': len(allowed_taus(p, ell))}
    res['bfs'] = bfs(p, cache)
    b = res['bfs']
    if b['status'] == 1:
        log(f"p={p} ell={ell} BFS: type 1, D={b['D']}, checked {b['checked']}, vertices {b['vertices']}, "
            f"hit ({b['hit']['q']},{b['hit']['tau']}) {b['hit']['branch']} d={b['hit']['d']}, path {b['path']}, {b['seconds']}s")
    elif b['status'] == 2:
        log(f"p={p} ell={ell} BFS: type 2 (closed failure set), |C|={b['size']}, checked {b['checked']}, {b['seconds']}s")
    else:
        log(f"p={p} ell={ell} BFS: type 3 (budget), checked {b['checked']}, vertices {b['vertices']}, {b['seconds']}s")
    if with_sweep:
        sw = res['sweep'] = sweep(p, cache)
        log(f"p={p} sweep: S0={sw['S0']}/{sw['n_root']} (II {sw['root_II']}, I {sw['root_I']}), failing root {sw['n_failing_root']}, "
            f"children {sw['n_children']} (success {sw['n_children_success']}, fail {sw['n_children_fail']}, undetermined {sw['n_children_undetermined']}), "
            f"E1' holds for {sw['E1prime_holds']}/{sw['n_failing_root']}, E1 by sweep {sw['E1_by_sweep']}, {sw['seconds']}s")
        # two traversals, one judgment: the BFS depth and the sweep must agree (asserted only when both are determined)
        if b['status'] in (1, 2) and sw['E1_by_sweep'] is not None:
            if sw['S0'] > 0:
                assert b['status'] == 1 and b['D'] == 0
            else:
                assert (b['status'] == 1 and b['D'] == 1) == sw['E1_by_sweep']
        h = res['hensel'] = hensel(p, cache)
        log(f"p={p} Hensel: pairs {h['pairs']} (k90>=2: {h['k90_ge2']}, kref>k90: {h['kref_gt_k90']}), exceptions 0, "
            f"lifts beyond kref {h['beyond_kref']}")
    res['sources_tested'] = len(cache.d)
    return res


def judge(res):
    """J1-J3 and the branches R1-R5 of contract.md."""
    st = {r['p']: r['bfs']['status'] for r in res}
    D = {r['p']: r['bfs'].get('D') for r in res}
    new = [r['p'] for r in res if r['p'] != 2521]
    det = [p for p in st if st[p] != 3]
    j2 = [p for p in st if st[p] == 2 or (st[p] == 1 and D[p] >= 2)]
    d_new = {k: sum(1 for p in new if st[p] == 1 and D[p] == k) for k in sorted({D[p] for p in new if st[p] == 1})}
    br = []
    if any(st[p] == 2 for p in st):
        br.append('R1')
    if any(st[p] == 1 and D[p] >= 2 for p in st):
        br.append('R2')
    if all(st[p] == 1 and D[p] <= 1 for p in det):
        if any(st[p] == 1 and D[p] == 1 for p in new):
            br.append('R3')
        elif all(st[p] == 1 and D[p] == 0 for p in new if st[p] != 3):
            br.append('R4')
    if any(st[p] == 3 for p in st):
        br.append('R5')
    return {'J1_status': st, 'J1_all_type1': all(v == 1 for v in st.values()),
            'J2_inputs_with_D_ge2_or_type2': j2, 'J2_E1_rejected': bool(j2),
            'J3_D_counts_new11': d_new, 'J3_prediction_met': d_new.get(0, 0) >= 2 and d_new.get(1, 0) >= 2,
            'D': D, 'branches': br}


def git(*args):
    return subprocess.run(['git', '-C', ROOT, *args], capture_output=True, text=True)


def main():
    mode = sys.argv[1]
    out = {'contract': 'contract.md', 'mode': mode, 'caps': {'sources': CAP_SOURCES, 'seconds': CAP_SECONDS, 'total': CAP_TOTAL},
           'sympy': sympy.__version__, 'python': sys.version.split()[0]}
    if mode == 'selftest':
        res = [run_input(p) for p in SELFTEST]
        # regression: the BFS order of 2521 reproduced in Reply8688Check20261023 (= the >>86 JSON)
        known = json.load(open(os.path.join(ROOT, 'Research', 'Reply8688Check20261023', 'check8688_results.json'), encoding='utf-8'))['bfs']['p2521']
        c2521 = Cache(2521)
        b = bfs(2521, c2521)
        ok = ([list(o) for o in b['order']] == [list(o) for o in known['order']] and b['vertex_list'] == known['vertices']
              and (b['hit']['d'], b['hit']['x'], b['hit']['c']) == (known['hit']['d'], known['hit']['x'], known['hit']['c']) and b['D'] == 1)
        log(f"regression 2521 BFS (order {len(b['order'])}, vertices {len(b['vertex_list'])}, D={b['D']}): {'OK' if ok else 'DIFF'}")
        assert ok
        out.update({'inputs': res, 'regression_2521': {'ok': ok, 'checked': b['checked'], 'vertices': b['vertices'], 'D': b['D']}})
        fn = 'selftest_results.json'
    elif mode == 'cohort':
        head = git('rev-parse', 'HEAD').stdout.strip()
        rel = os.path.relpath(HERE, ROOT).replace(os.sep, '/')
        tracked = git('ls-files', rel).stdout.split()
        for f in ('bfs12.py', 'contract.md'):
            assert f'{rel}/{f}' in tracked, f
        assert git('diff', '--quiet', 'HEAD', '--', f'{rel}/bfs12.py', f'{rel}/contract.md').returncode == 0, 'working tree differs from HEAD'
        log(f'G0 ok: HEAD {head}, inputs {INPUTS}')
        out['head'] = head
        res = []
        for p in INPUTS:
            res.append(run_input(p))
        out['inputs'] = res
        out['judgments'] = judge(res)
        log('judgments:', json.dumps(out['judgments'], default=str))
        fn = 'bfs12_results.json'
    else:
        raise SystemExit(__doc__)
    out['seconds'] = round(time.monotonic() - T0, 2)
    out['log'] = LOG
    log(f'done: {out["seconds"]}s')
    with open(os.path.join(HERE, fn), 'w', encoding='utf-8', newline='\n') as fp:
        fp.write(json.dumps(out, ensure_ascii=False, indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
