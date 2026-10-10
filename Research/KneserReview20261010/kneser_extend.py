"""Measurement (no judgement): how often HKlog fails among hard primes p < 10^8, for K and for A,
and whether Hlog (some actual witness in the coprime box) still holds there.  The author scanned
only up to the first failure (17,816 inputs); this counts all failures in the same population."""
import json, sys, time
import kneser_check as k

def main():
    t0 = time.time(); N = 10**8
    hard = [q for q in k.sieve(N) if q % 840 in k.SQ840]
    K_fail, A_fail, both_fail, hlog_fail, false_pos = [], [], [], [], []
    for p in hard:
        rows = [k.conditions(p, a, u) for a, u in k.coprime_pairs(k.clog2(p))]
        Kp = any(r["K"] for r in rows); Ap = any(r["A"] for r in rows); W = any(r["witnesses"] for r in rows)
        false_pos += [(p, r["a"], r["u"]) for r in rows if (r["K"] or r["A"]) and not r["witnesses"]]
        if not Kp: K_fail.append(p)
        if not Ap: A_fail.append(p)
        if not Kp and not Ap: both_fail.append(p)
        if not W: hlog_fail.append(p)
    out = {"population": "primes p < 10^8 with p mod 840 in {1,121,169,289,361,529}", "count": len(hard),
           "K_fail": K_fail, "A_fail": A_fail, "both_fail": both_fail, "hlog_fail": hlog_fail,
           "false_positives": false_pos, "seconds": round(time.time() - t0, 1)}
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
