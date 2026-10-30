r"""
exact_check.py  --  EXACT check: Python integers and fractions only.
No floating point arithmetic, no ball/interval arithmetic, no flint/numpy/mpmath, no code from other scripts.

Mode 1, certificates (results/certs/n_<n>.json, written by turan_phase6.py):
    The certificate holds b_0..b_(n-1), each real/imag part = numerator / 2^k.  The n-tuple is
    z_1 = 1 and the roots of Z^(n-1) + b_1 Z^(n-2) + ... + b_(n-1), i.e.
        W(x) = (1 - x)(b_0 + b_1 x + ... + b_(n-1) x^(n-1)) = (1 - z_1 x)...(1 - z_n x).
    All numbers are put over one denominator 2^D:  w_l = W_l / 2^D,  W_l Gaussian integers.
    Checks b_0 = 1 and W(1) = 0 exactly. Since w_0 = 1, Newton's identities need no division:
        s_k = -(k w_k + w_1 s_(k-1) + ... + w_(k-1) s_1),     s_k = N_k / 2^(kD),  N_k Gaussian integer.
    For EVERY k = 1..n:   q^2 (Re(N_k)^2 + Im(N_k)^2) <= p^2 4^(kD),   where c = p/q is the bound to check
    (the "rigorous_upper_bound" of the certificate, or --bound).
    Cost ~ n^3 (the numbers N_k have k*D bits), memory ~ n^2 D bits.

Mode 2, small-n tuples (small_n_tuples.json, written by turan_phase5.py):
    The tuple is the exact binary value of each stored float: Fraction(float(s)) (the file stores repr strings).
    Checks z_1 = 1 exactly and n numbers, then s_k = sum_j z_j^k exactly for k = 1..n.
    Saves results/small_n_exact.json.

Run:   python exact_check.py results\certs\n_300.json [more files] [--bound 0.69101890823]
       python exact_check.py small [nmax]          (all tuples, or only n <= nmax)
Needs: nothing outside the Python standard library
"""
import hashlib, json, math, os, sys, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))


def sqrt_decimals(A, sh, d):
    """sqrt(A / 2^sh) with d decimals: (rounded down, rounded up), exact integers"""
    x = (A * 10 ** (2 * d)) >> sh                          # floor(A 10^(2d) / 2^sh)
    lo = math.isqrt(x)
    exact = lo * lo == x and (A * 10 ** (2 * d)) == x << sh
    hi = lo if exact else lo + 1
    return f"{lo // 10 ** d}.{str(lo % 10 ** d).zfill(d)}", f"{hi // 10 ** d}.{str(hi % 10 ** d).zfill(d)}"


# ---------------------------------------------------------------- mode 1
def certificate(path, bound):
    t0 = time.time()
    raw = open(path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    cert = json.loads(raw)
    n = cert["n"]
    c = Fraction(bound if bound is not None else cert["rigorous_upper_bound"])
    p, q = c.numerator, c.denominator
    pairs = cert["b"]
    if len(pairs) != n:
        sys.exit(f"{path}: {len(pairs)} coefficients, expected n = {n}")
    D = max(int(k) for re, im in pairs for _, k in (re, im))
    Br = [int(re[0]) << (D - int(re[1])) for re, im in pairs]       # b_l = (Br[l] + i Bi[l]) / 2^D
    Bi = [int(im[0]) << (D - int(im[1])) for re, im in pairs]
    if Br[0] != 1 << D or Bi[0] != 0:
        sys.exit(f"{path}: b_0 is not exactly 1")
    Wr = [Br[0]] + [Br[l] - Br[l - 1] for l in range(1, n)] + [-Br[n - 1]]
    Wi = [Bi[0]] + [Bi[l] - Bi[l - 1] for l in range(1, n)] + [-Bi[n - 1]]
    if sum(Wr) != 0 or sum(Wi) != 0:
        sys.exit(f"{path}: W(1) is not exactly 0")

    Nr, Ni = [0] * (n + 1), [0] * (n + 1)
    p2, q2 = p * p, q * q
    best, kbest, failed = -1, 0, []
    for k in range(1, n + 1):
        # Horner in 2^D:  acc = sum_{i=1}^{k-1} W_i N_(k-i) 2^((i-1)D),  numbers stay ~ (k-i)D bits long
        ar = ai = 0
        for i in range(k - 1, 0, -1):
            a, b, x, y = Wr[i], Wi[i], Nr[k - i], Ni[k - i]
            ar = (ar << D) + a * x - b * y
            ai = (ai << D) + a * y + b * x
        sh = (k - 1) * D
        Nr[k] = -((k * Wr[k]) << sh) - ar
        Ni[k] = -((k * Wi[k]) << sh) - ai
        A = Nr[k] * Nr[k] + Ni[k] * Ni[k]                  # |s_k|^2 = A / 4^(kD)
        if q2 * A > p2 << (2 * k * D):
            failed.append(k)
        if kbest == 0 or (A << (2 * kbest * D)) > (best << (2 * k * D)):
            best, kbest = A, k
    lo, hi = sqrt_decimals(best, 2 * kbest * D, 20)
    return dict(path=path, sha=sha, n=n, D=D, k=kbest, lo=lo, hi=hi, cstr=bound if bound is not None else cert["rigorous_upper_bound"], failed=failed, t=time.time() - t0)


def mode1(paths, bound):
    print("exact_check.py, certificates: integers only (no floats, no balls); every k = 1..n checked")
    print(f"{'n':>7} {'D':>5} {'k_max':>6} {'max|s_k| >= (20 dec)':>23} {'max|s_k| <= (20 dec)':>23} "
          f"{'bound c':>16} {'result':>8} {'time':>8}  sha256")
    for path in paths:
        r = certificate(path, bound)
        res = "PASSED" if not r["failed"] else f"FAILED at k = {r['failed'][:5]}{' ...' if len(r['failed']) > 5 else ''} ({len(r['failed'])} k)"
        print(f"{r['n']:>7} {r['D']:>5} {r['k']:>6} {r['lo']:>23} {r['hi']:>23} {r['cstr']:>16} {res:>8} {r['t']:>7.1f}s  {r['sha'][:16]}", flush=True)


# ---------------------------------------------------------------- mode 2
def dyadic(s):
    """exact value of the binary float written as the string s: (integer m, e) with value m / 2^e"""
    f = Fraction(float(s))
    e = f.denominator.bit_length() - 1
    assert f.denominator == 1 << e
    return f.numerator, e


def small_tuple(n, rec):
    t0 = time.time()
    z = rec["z"]
    if len(z) != n:
        sys.exit(f"n = {n}: {len(z)} numbers in the tuple")
    parts = [(dyadic(re), dyadic(im)) for re, im in z]
    E = max(e for pr in parts for _, e in pr)
    Zr = [m << (E - e) for (m, e), _ in parts]                # z_j = (Zr[j] + i Zi[j]) / 2^E
    Zi = [m << (E - e) for _, (m, e) in parts]
    if Zr[0] != 1 << E or Zi[0] != 0:
        sys.exit(f"n = {n}: z_1 is not exactly 1")
    Pr, Pi = [1] * n, [0] * n                               # Z_j^k
    best, kbest = -1, 0
    for k in range(1, n + 1):
        for j in range(n):
            Pr[j], Pi[j] = Pr[j] * Zr[j] - Pi[j] * Zi[j], Pr[j] * Zi[j] + Pi[j] * Zr[j]
        sr, si = sum(Pr), sum(Pi)                           # s_k = (sr + i si) / 2^(kE)
        A = sr * sr + si * si
        if kbest == 0 or (A << (2 * kbest * E)) > (best << (2 * k * E)):
            best, kbest = A, k
    sh = 2 * kbest * E
    lo, hi = sqrt_decimals(best, sh, 20)
    _, up12 = sqrt_decimals(best, sh, 12)
    st = Fraction(rec["bound"])
    stored_ok = st.denominator ** 2 * best <= st.numerator ** 2 << sh
    return dict(n=n, E=E, k_max=kbest, max_floor20=lo, max_ceil20=hi, rigorous_upper_bound=up12,
                stored_bound=rec["bound"], stored_bound_valid=stored_ok, seconds=round(time.time() - t0, 2))


def mode2(only=None):
    data = json.load(open(os.path.join(HERE, "small_n_tuples.json")))
    print("exact_check.py, small_n_tuples.json: exact s_k = sum z_j^k (z_j = exact binary floats), every k = 1..n")
    print(f"{'n':>4} {'k_max':>6} {'max|s_k| >= (20 dec)':>23} {'max|s_k| <= (20 dec)':>23} {'R_n <= (12 dec)':>16} "
          f"{'stored':>12} {'stored ok?':>10} {'time':>7}")
    out = {}
    for key in sorted(data, key=int):
        n = int(key)
        if only and n > only:
            continue
        r = small_tuple(n, data[key])
        out[key] = r
        print(f"{n:>4} {r['k_max']:>6} {r['max_floor20']:>23} {r['max_ceil20']:>23} {r['rigorous_upper_bound']:>16} "
              f"{r['stored_bound']:>12} {'yes' if r['stored_bound_valid'] else 'NO':>10} {r['seconds']:>6.1f}s", flush=True)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, "results", "small_n_exact.json"), "w"), indent=1)
    print("\nsaved results/small_n_exact.json")


def main():
    args = sys.argv[1:]
    if args and args[0] == "small":
        return mode2(int(args[1]) if len(args) > 1 else None)
    bound = None
    if "--bound" in args:
        i = args.index("--bound")
        bound = args[i + 1]
        args = args[:i] + args[i + 2:]
    mode1(args, bound)


if __name__ == "__main__":
    main()
