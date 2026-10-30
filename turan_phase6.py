"""
turan_phase6.py  --  the smooth-phase construction + rigorous check.

What we learned from the free-search tuples (n = 50..300):
  * every power sum has the SAME size:  |S_k| = T  for all k = 1..n
  * the angle of S_k follows one smooth curve  phi(k/n), the same curve for every n
  * our old family keeps the angle FIXED on the first third; that is what costs ~0.003

The construction here:
  S_k = T * exp(i * phi(k/n))        for k = 1 .. n-L-1      (phi = 24 straight pieces, finer near 0)
  S_m for the last L+1 indices       chosen so that b_n = 0  (exact; this needs L < n/2)
  then the n-tuple is {1} + roots of  Z^(n-1) + b_1 Z^(n-2) + ... + b_(n-1)

Rigorous check (python-flint balls, works for ANY profile, no special formula):
  b_l from exp(-sum S_k x^k / k) / (1 - x), take exact midpoints  ->  this DEFINES the tuple
  power sums of that tuple from log((1 - x) * B(x))   ->   R_n <= max_k |power sum k|
  (the tuple contains 1, so scaling it to max modulus 1 can only make the sums smaller)

Run:   python turan_phase6.py           n = 300 .. 30000     (a few minutes)
       python turan_phase6.py big       also 100000 and 300000 (long, a few GB RAM)
       python turan_phase6.py certs [n ...]   no optimization: rebuild the tuples of smooth_table.json
                                              from their params and only write the certificates
       python turan_phase6.py run n1 n2 ...   like big, for these n only (warm start from the nearest n in
                                              smooth_table.json); adds new rows, never changes existing ones
Needs: numpy scipy python-flint  and  smooth_start.json in the same folder
Saves: smooth_table.json   (all parameters + rigorous bounds)
       results/certs/n_<n>.json   (certificate: the exact b_0..b_(n-1); check it with verify_tuple.py)
"""
import json, math, os, sys, time
from decimal import Decimal, getcontext
from fractions import Fraction
import numpy as np
from scipy.optimize import minimize
from flint import acb, acb_poly, acb_series, fmpq, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
M = 24                                    # number of straight pieces in phi


# ---------------------------------------------------------------- fast power series (FFT)
def _mul(a, b, N):
    m = 1 << int(np.ceil(np.log2(len(a) + len(b))))
    return np.fft.ifft(np.fft.fft(a, m) * np.fft.fft(b, m))[:N]


def _inv(f, N):
    g = np.array([1 / f[0]], complex); m = 1
    while m < N:
        m = min(2 * m, N)
        g = (2 * np.pad(g, (0, m - len(g))) - _mul(_mul(f[:m], g, m), g, m))[:m]
    return g


def _log(f, N):
    q = _mul(np.arange(1, N) * f[1:N], _inv(f, N), N - 1)
    out = np.zeros(N, complex); out[1:] = q / np.arange(1, N); return out


def _exp(h, N):
    g = np.ones(1, complex); m = 1
    while m < N:
        m = min(2 * m, N)
        g = np.pad(g, (0, m - len(g)))
        hm = h[:m] if len(h) >= m else np.pad(h, (0, m - len(h)))
        g = _mul(g, hm - _log(g, m) + np.eye(1, m, 0)[0], m)
    return g[:N]


# ---------------------------------------------------------------- the construction (floats)
def knots(top):
    return top * (np.arange(M + 1) / M) ** 2


def choose_L(n):
    return (n - 1) // 2 - 1                          # last block = last ~half, and 2(n-L) > n


def profile(params, n):
    """Returns S_1..S_n and the float bound."""
    L = choose_L(n)
    K = n - L - 1
    T, ph = params[0], params[1:]
    phi = np.interp(np.arange(1, K + 1) / n, knots(K / n), ph)
    Sf = T * np.exp(1j * phi)
    h = np.zeros(n + 1, complex); h[1:K + 1] = -Sf / np.arange(1, K + 1)
    b = np.cumsum(_exp(h, n + 1))                    # coefficients of exp(-sum S_k x^k/k)/(1-x)
    m = np.arange(n - L, n + 1)
    w = np.abs(b[n - m]) / m
    xF = abs(b[n]) / w.sum()
    SF = (b[n] / w.sum()) * np.conj(b[n - m]) / np.abs(b[n - m])   # makes b_n exactly 0
    return np.concatenate([Sf, SF]), max(T, xF), T, xF


def optimize(n, p0):
    cons = {"type": "ineq", "fun": lambda p: p[0] - profile(p, n)[3]}
    res = minimize(lambda p: p[0], p0, method="SLSQP", constraints=[cons],
                   options=dict(maxiter=300, ftol=1e-14))
    return res.x


# ---------------------------------------------------------------- rigorous check
def _ex(x):
    f = Fraction(float(x)); return fmpq(f.numerator, f.denominator)


def _to_frac(x):
    man, e = x.man_exp(); return Fraction(int(man)) * Fraction(2) ** int(e)


def _up(fr, d=12):
    return Fraction(math.ceil(fr * 10 ** d), 10 ** d)


def _dec(fr, d=12):
    getcontext().prec = 30
    return f"{Decimal(fr.numerator) / Decimal(fr.denominator):.{d}f}"


def _save_cert(path, b, U, prec):
    """Certificate: b_0..b_(n-1) exactly. Each real/imag part is [numerator, k] = numerator / 2^k."""
    def q(x):
        f = _to_frac(x)                                    # x is an exact midpoint
        assert f.denominator & (f.denominator - 1) == 0
        return [str(f.numerator), f.denominator.bit_length() - 1]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(dict(n=len(b), script="turan_phase6.py", prec_bits=prec, rigorous_upper_bound=_dec(U),
                   tuple="{1} + roots of Z^(n-1) + b_1 Z^(n-2) + ... + b_(n-1)",
                   format="b[l] = [re, im]; each part [numerator, k] means numerator / 2^k exactly",
                   b=[[q(c.real), q(c.imag)] for c in b]), open(path, "w"))


def rigorous(S, prec=320, cert=None):
    n = len(S)
    ctx.prec = prec; ctx.cap = n + 2
    Sb = [acb(_ex(c.real), _ex(c.imag)) for c in S]
    G = acb_series([acb(0)] + [Sb[k - 1] / k for k in range(1, n + 1)], prec=n + 1)
    e = (-G).exp().coeffs()
    b, run = [], acb(0)
    for l in range(n):
        run += e[l]; b.append(run.mid())               # exact midpoints define the tuple
    lg = acb_series((acb_poly(b) * acb_poly([1, -1])).coeffs(), prec=n + 1).log().coeffs()
    best = Fraction(0)
    for k in range(1, n + 1):
        best = max(best, _to_frac((-k * lg[k]).abs_upper().upper()))
    if cert:
        _save_cert(cert, b, _up(best), prec)
    return _up(best)


# ---------------------------------------------------------------- main
def cert_path(n):
    return os.path.join(HERE, "results", "certs", f"n_{n}.json")


def certs(ns):
    """Rebuild the tuples of smooth_table.json from their saved params (no optimization) and save certificates.
    Does NOT overwrite smooth_table.json. The bound must come out identical to the stored one."""
    tab = json.load(open(os.path.join(HERE, "smooth_table.json")))
    ns = ns or sorted(map(int, tab))
    print(f"{'n':>7} {'RIGOROUS R_n <=':>17} {'smooth_table.json':>18} {'same?':>6}")
    for n in ns:
        t0 = time.time()
        p = np.array([float(x) for x in tab[str(n)]["params"]])
        S, fb, T, xF = profile(p, n)
        U = rigorous(S, cert=cert_path(n))
        old = tab[str(n)]["rigorous_upper_bound"]
        print(f"{n:>7} {_dec(U):>17} {old:>18} {'yes' if _dec(U) == old else 'NO':>6}   ({time.time() - t0:.0f}s)", flush=True)
    print("\nsaved results/certs/n_<n>.json")


def run(ns):
    """Like 'big' but for the given n: warm start from the n in smooth_table.json nearest in ratio,
    optimize, rigorous check, certificate, and ADD the row to smooth_table.json (existing rows are never changed)."""
    path = os.path.join(HERE, "smooth_table.json")
    tab = json.load(open(path))
    print(f"{'n':>7} {'start from':>10} {'smooth (float)':>15} {'RIGOROUS R_n <=':>17}")
    for n in ns:
        if str(n) in tab:
            print(f"{n:>7}  already in smooth_table.json, skipped (use 'certs' to rebuild its certificate)", flush=True)
            continue
        t0 = time.time()
        near = min(tab, key=lambda m: (abs(math.log(int(m) / n)), int(m)))
        p = optimize(n, np.array([float(x) for x in tab[near]["params"]]))
        S, fb, T, xF = profile(p, n)
        U = rigorous(S, cert=cert_path(n))
        tab[str(n)] = dict(M=M, L=choose_L(n), params=[repr(float(x)) for x in p],
                           float_bound=float(fb), rigorous_upper_bound=_dec(U))
        tab = dict(sorted(tab.items(), key=lambda kv: int(kv[0])))
        json.dump(tab, open(path, "w"), indent=1)
        print(f"{n:>7} {near:>10} {fb:>15.12f} {_dec(U):>17}   ({time.time() - t0:.0f}s)", flush=True)
    print("\ndone: new rows (if any) in smooth_table.json, certificates in results/certs/")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "certs":
        return certs([int(x) for x in sys.argv[2:]])
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        return run([int(x) for x in sys.argv[2:]])
    ns = [300, 1000, 3000, 10_000, 30_000]
    if len(sys.argv) > 1 and sys.argv[1] == "big":
        ns += [100_000, 300_000]
    p = np.array([float(x) for x in json.load(open(os.path.join(HERE, "smooth_start.json")))["params"]])
    old = {300: 0.694769708, 1000: 0.691609453, 3000: 0.690246947, 10000: 0.689520872,
           30000: 0.689203546, 100000: 0.689035050, 300000: 0.688961575}
    print("Smooth-phase construction.  'old family' = Stage 2 of turan_phase5.  Record (asymptotic): 0.6873358063")
    print(f"{'n':>7} {'old family':>11} {'smooth (float)':>15} {'RIGOROUS R_n <=':>17}")
    out = {}
    for n in ns:
        t0 = time.time()
        p = optimize(n, p)                                 # warm start from the previous n
        S, fb, T, xF = profile(p, n)
        U = rigorous(S, cert=cert_path(n))
        out[str(n)] = dict(M=M, L=choose_L(n), params=[repr(float(x)) for x in p],
                           float_bound=float(fb), rigorous_upper_bound=_dec(U))
        json.dump(out, open(os.path.join(HERE, "smooth_table.json"), "w"), indent=1)
        print(f"{n:>7} {old.get(n, float('nan')):>11.9f} {fb:>15.12f} {_dec(U):>17}   ({time.time() - t0:.0f}s)", flush=True)
    print("\nsaved smooth_table.json")


if __name__ == "__main__":
    main()
