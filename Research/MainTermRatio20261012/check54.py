"""Independent check of the table in Outcasts >>54 (M, R, M - R, W_tot for 66529, 8615161, 53722321) with exact rationals,
plus 1201 and 2521 (the other two inputs >>54 says it completed).  Own implementation: trial division, divisor enumeration."""
import json, sys
from fractions import Fraction
from math import gcd, isqrt

def factor(n, small):
    f = {}
    for q in small:
        if q * q > n: break
        while n % q == 0: f[q] = f.get(q, 0) + 1; n //= q
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def phi(m):
    r = m; n = m; q = 2
    while q * q <= n:
        if n % q == 0:
            r -= r // q
            while n % q == 0: n //= q
        q += 1
    if n > 1: r -= r // n
    return r

def sieve(n):
    b = bytearray([1]) * (n + 1); b[0] = b[1] = 0
    for i in range(2, isqrt(n) + 1):
        if b[i]: b[i * i::i] = bytearray(len(b[i * i::i]))
    return [i for i in range(n + 1) if b[i]]

def main():
    small = sieve(40000); out = {}
    for p in (1201, 2521, 66529, 8615161, 53722321):
        B = (p - 1).bit_length(); M = Fraction(0); Ep = Fraction(0); Em = Fraction(0); W = 0; pairs = 0; Tchi = Fraction(0)
        for a in range(1, B + 1):
            for u in range(1, B // a + 1):
                if gcd(a, u) != 1: continue
                m = 4 * a * u; N = a * p + u; f = factor(N, small); pairs += 1
                assert all(q * q <= N or q == N or True for q in f) and small[-1] ** 2 > N or max(f) < small[-1] ** 2
                odd = {q: e for q, e in f.items() if q % 2 == 1}
                tau = 1; chi = 1
                for q, e in odd.items():
                    tau *= e + 1
                    leg = pow(p % q, (q - 1) // 2, q) == 1
                    chi *= (e + 1) if leg else (1 if e % 2 == 0 else 0)
                divs = [1]
                for q, e in odd.items(): divs = [d * q ** k for d in divs for k in range(e + 1)]
                w = sum(1 for s in divs if (s + 1) % m == 0)
                U = Fraction(tau, phi(m)); E = w - U
                M += U; W += w; Tchi += Fraction(chi, phi(m))
                if E > 0: Ep += E
                else: Em += -E
        R = Ep + Em
        out[str(p)] = {"B": B, "pairs": pairs, "M": str(M), "R": str(R), "M_minus_R": str(M - R), "W_tot": W, "E_plus": str(Ep), "E_minus": str(Em),
                       "T_chi": str(Tchi), "M_minus_Tchi": str(M - Tchi), "identity_M_minus_R_eq_W_minus_2Eplus": M - R == W - 2 * Ep}
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
