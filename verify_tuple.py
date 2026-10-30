r"""
verify_tuple.py  --  independent check of a certificate  results/certs/n_<n>.json

Shares NO code with turan_phase*.py: it only reads the json file.

The certificate holds exact numbers b_0 .. b_(n-1) (each real/imag part = integer / 2^k). The n-tuple is
    z_1 = 1   and   z_2 .. z_n = the roots of  Z^(n-1) + b_1 Z^(n-2) + ... + b_(n-1),
in other words  W(x) = (1 - x) B(x) = (1 - z_1 x)(1 - z_2 x)...(1 - z_n x)  with  B(x) = b_0 + b_1 x + ... + b_(n-1) x^(n-1).
The roots are never computed: the power sums s_k = z_1^k + ... + z_n^k follow from W alone.

What it checks
  1. the file has n numbers and b_0 = 1 exactly
  2. W is built exactly (Python integers); W(0) = 1 and W(1) = 0 exactly, so z_1 = 1 is in the tuple
  3. s_k for EVERY k = 1..n from Newton's identities
         s_k = -k w_k - (w_1 s_(k-1) + w_2 s_(k-2) + ... + w_(k-1) s_1)        (w_l = coefficient of x^l in W)
     in python-flint balls, solved in one of two ways:
       newton  : this recursion, term by term                                  O(n^2), for n up to ~2000
       inverse : s_1 + s_2 x + s_3 x^2 + ... = -W'(x) / W(x), with 1/W from our own Newton iteration
                 that uses only polynomial products (no series log / exp / division from the library)
     The precision is raised until every ball is tight.
  4. prints a rigorous interval for max_k |s_k| and compares it with the bound written in the certificate

Run:   python verify_tuple.py results\certs\n_1000.json           (newton if n <= 2000, else inverse)
       python verify_tuple.py results\certs\n_1000.json both      (both methods, and checks that they agree)
       python verify_tuple.py results\certs\n_300.json results\certs\n_1000.json      (several files)
Needs: python-flint
"""
import json, sys, time
from fractions import Fraction
from flint import acb, acb_poly, fmpq, ctx

TIGHT = Fraction(1, 2 ** 80)          # every |s_k| ball must be narrower than this
NEWTON_MAX = 10_000                   # newton is O(n^2): never used above this n


# ---------------------------------------------------------------- exact part (Python integers only)
def load(path):
    """W as exact Gaussian integers over a common 2^K:  w_l = (wr[l] + i wi[l]) / 2^K,  l = 0..n."""
    d = json.load(open(path))
    n, b = d["n"], d["b"]
    assert len(b) == n, f"expected {n} coefficients, found {len(b)}"
    parts = [(int(num), int(k)) for c in b for num, k in c]
    K = max(k for _, k in parts)
    ints = [num << (K - k) for num, k in parts]
    br, bi = ints[0::2], ints[1::2]                         # b_l = (br[l] + i bi[l]) / 2^K
    assert br[0] == 1 << K and bi[0] == 0, "b_0 is not exactly 1"
    wr = [br[0]] + [br[l] - br[l - 1] for l in range(1, n)] + [-br[n - 1]]
    wi = [bi[0]] + [bi[l] - bi[l - 1] for l in range(1, n)] + [-bi[n - 1]]
    assert sum(wr) == 0 and sum(wi) == 0, "W(1) is not exactly 0"
    return n, d.get("rigorous_upper_bound"), (wr, wi, K)


def exact(x):
    """an exact arb (e.g. from .upper()) as a Fraction"""
    m, e = (int(t) for t in x.man_exp())
    return Fraction(m << e) if e >= 0 else Fraction(m, 1 << -e)


def dec(fr, d, up):
    q = fr * 10 ** d
    q = -((-q.numerator) // q.denominator) if up else q.numerator // q.denominator
    return f"{q // 10 ** d}.{q % 10 ** d:0{d}d}"


# ---------------------------------------------------------------- balls
def balls(re, im, K):
    den = 1 << K
    return [acb(fmpq(a, den), fmpq(c, den)) for a, c in zip(re, im)]


def coeffs(p, lo, hi):
    """coefficients lo..hi-1 of p (acb_poly drops trailing zeros, so pad with 0)"""
    c = p.coeffs()
    return [c[i] if i < len(c) else acb(0) for i in range(lo, hi)]


def sums_newton(w, n):
    s = [acb(0)] * (n + 1)
    for k in range(1, n + 1):
        acc = -k * w[k]
        for i in range(1, k):
            acc -= w[i] * s[k - i]
        s[k] = acc
    return s


def inverse(w, n):
    """1/W mod x^n.  If W g = 1 + x^m E then g - x^m g E is the inverse mod x^(2m)."""
    g = [acb(1)]                                            # w_0 = 1 exactly
    while len(g) < n:
        m = len(g)
        m2 = min(2 * m, n)
        E = coeffs(acb_poly(w[:m2]) * acb_poly(g), m, m2)
        g += [-c for c in coeffs(acb_poly(g) * acb_poly(E), 0, m2 - m)]
    return g


def sums_inverse(w, dw, n):
    """s_1 + s_2 x + ... + s_n x^(n-1) = -W'(x) / W(x)  mod x^n"""
    P = coeffs(acb_poly(dw) * acb_poly(inverse(w, n)), 0, n)
    return [acb(0)] + [-c for c in P]


def run(method, n, W, prec):
    wr, wi, K = W
    ctx.prec = prec
    w = balls(wr, wi, K)
    if method == "newton":
        s = sums_newton(w, n)
    else:
        dw = balls([l * wr[l] for l in range(1, n + 1)], [l * wi[l] for l in range(1, n + 1)], K)
        s = sums_inverse(w, dw, n)
    assert len(s) == n + 1
    iv = []                                                 # [lower, upper] for |s_k|, k = 1..n
    for k in range(1, n + 1):
        a = abs(s[k])
        iv.append((max(Fraction(0), exact(a.lower())), exact(a.upper())))
    return iv, all(x.is_exact() for x in w)


def solve(method, n, W, prec):
    """raise the precision until every |s_k| ball is narrower than TIGHT"""
    t0 = time.time()                                        # time includes the retries
    while True:
        try:
            iv, wexact = run(method, n, W, prec)
            width = max(u - l for l, u in iv)
        except ValueError:                                  # a ball became infinite
            iv, width = None, None
        if width is not None and width < TIGHT:
            return iv, prec, wexact, time.time() - t0
        if prec >= 1 << 16:
            sys.exit(f"{method}: balls still too wide at {prec} bits")
        lost = prec if width is None else prec + width.numerator.bit_length() - width.denominator.bit_length()
        print(f"  {method}: {prec} bits not enough, retrying", flush=True)
        prec = max(2 * prec, lost + 128)


# ---------------------------------------------------------------- main
def main():
    paths = [a for a in sys.argv[1:] if a.endswith(".json")]
    choice = [a for a in sys.argv[1:] if not a.endswith(".json")]
    for path in paths:
        check(path, choice[0] if choice else None)


def check(path, method):
    t_all = time.time()
    n, claimed, W = load(path)
    method = method or ("newton" if n <= 2000 else "inverse")
    methods = ["newton", "inverse"] if method == "both" else [method]
    print(f"verify_tuple.py  {path}")
    if n > NEWTON_MAX and "newton" in methods:              # O(n^2): days at n = 300000
        print(f"  n > {NEWTON_MAX}: newton skipped, inverse only")
        methods = ["inverse"]
    print(f"  n = {n}: b_0 = 1 exactly, W(1) = 0 exactly  ->  z_1 = 1 is in the tuple (n numbers in all)")
    res = {}
    for m in methods:
        iv, prec, wexact, dt = solve(m, n, W, 128 if m == "newton" else 320)
        lo, up = max(l for l, _ in iv), max(u for _, u in iv)
        kmax = 1 + max(range(n), key=lambda i: iv[i][1])
        ok = claimed is not None and up <= Fraction(claimed)
        print(f"  {m:<7}: all k = 1..{n} done | {prec} bits | W exact in balls: {'yes' if wexact else 'no'} | {dt:.1f}s")
        print(f"           max_k |s_k| in [{dec(lo, 15, False)}, {dec(up, 15, True)}]   (largest at k = {kmax})")
        print(f"           RIGOROUS max_k |s_k| <= {dec(up, 12, True)}   certificate: {claimed}   -> "
              f"{'CONFIRMED' if ok else 'NOT confirmed'}")
        res[m] = iv
    if len(res) == 2:
        bad = [k + 1 for k, ((a, b), (c, d)) in enumerate(zip(res["newton"], res["inverse"])) if max(a, c) > min(b, d)]
        print(f"  newton vs inverse: " + ("intervals overlap for every k" if not bad else f"DISAGREE at k = {bad[:10]}"))
    print(f"  total time for n = {n}: {time.time() - t_all:.1f}s", flush=True)


if __name__ == "__main__":
    main()
