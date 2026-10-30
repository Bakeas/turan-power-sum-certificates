# Certified upper bounds in Turán's pure power sum problem (Erdős #519)

Explicit $n$-tuples of complex numbers, stored as exact certificates, that give **rigorous upper bounds** for

$$R_n \;=\; \min_{z_1 = 1,\ z_2, \dots, z_n \in \mathbb{C}}\ \max_{1 \le k \le n}\ \bigl|z_1^k + z_2^k + \dots + z_n^k\bigr|$$

for 13 values of $n$ between 300 and 1 000 000, together with two independent checkers.

**Author:** Vasileios Bakeas, Department of Mechanical Engineering, University of Thessaly, Volos, Greece (bakeasb@gmail.com)
**Status:** not peer reviewed. Made with substantial help from AI models; see [AI use](#ai-use).

## Background

This is Turán's pure power sum problem: [Erdős problem #519](https://www.erdosproblems.com/519), and constant C42 ($=\limsup R_n$) in the crowdsourced [Optimization Constants repository](https://github.com/teorth/optimizationproblems).

- Lower bound: $R_n > q$ for every $n$, for an absolute constant $q > 1/2$ (Biró 2000, improving Atkinson's 1/6 and Biró's 1/2).
- Upper bound, published: $\limsup R_n < 0.69368$ (computation of Harcos, reported in Biró 2000; short proof in Harcos's note).
- Cheer and Goldston (1996) computed near-optimal configurations for $n \le 55$.
- Unrefereed 2026 certificates in the Optimization Constants repository claim smaller values of $\limsup R_n$, the smallest being $0.6873358063$ (Koshta, PR #263). These are asymptotic statements whose analytic reduction, by the authors' own account, still needs human review.

**Normalization.** Some authors require $1 = z_1 \ge |z_2| \ge \dots \ge |z_n|$. This gives the same $R_n$. One direction is trivial. Conversely, take any tuple with $z_1 = 1$ and $\lambda = \max_j |z_j| \ge 1$, divide it by $\lambda$, and rotate so that a point of largest modulus becomes 1. Each $|s_k|$ is then multiplied by $\lambda^{-k} \le 1$. So every bound below holds in both normalizations.

## Results

### Table 1: smooth-phase construction

Each bound is the exact value of $\max_{1\le k\le n} |s_k|$ for one explicit tuple, over **every** $k = 1, \dots, n$, rounded **up** to 12 decimals.

| $n$ | $R_n \le$ | checks |
|---:|:---|:---:|
| 300 | 0.691018908239 | P V E |
| 500 | 0.689542310941 | P V E |
| 1 000 | 0.688245154929 | P V E |
| 2 000 | 0.687468053119 | P V E |
| 3 000 | 0.687169029139 | P V E |
| 5 000 | 0.686902345661 | P V E |
| 10 000 | 0.686673342516 | P V |
| 20 000 | 0.686540042277 | P V |
| 30 000 | 0.686489957913 | P V |
| 50 000 | 0.686446081102 | P V |
| 100 000 | 0.686409189526 | P V |
| 300 000 | 0.686380380702 | P V |
| 1 000 000 | 0.686368028005 | P V |

The checks are:
- **P**: the generator's own check (`turan_phase6.py`).
- **V**: the independent ball-arithmetic checker (`verify_tuple.py`).
- **E**: the exact integer checker (`exact_check.py`).

They are described [below](#how-the-bounds-are-certified). The maximum is attained at $k = n$ for every row except $n = 10\,000$ (there at $k = 5554$; the $|s_k|$ are almost all equal). [`results/final_table.md`](results/final_table.md) has the full table with enclosures.

### Table 2: small $n$ (free numerical search, checked exactly)

| $n$ | $R_n \le$ (this repo) | Cheer–Goldston 1996 |
|---:|:---|:---|
| 10 | 0.742655049145 | 0.74265504914476718905 |
| 15 | 0.728108891701 | 0.72810889170059241298 |
| 20 | 0.720169916107 | 0.72016991610628663886 |
| 25 | 0.715082407017 | 0.71508240701664073135 |
| 30 | 0.711505714456 | 0.71150571445562842882 |
| 40 | 0.706751927398 | 0.70675192739708693051 |
| 50 | 0.703694627087 | 0.70369462708689897930 |
| 55 | 0.702532663224 | 0.70253266322361406713 |
| 70 | 0.699927441737 | – |
| 100 | 0.696818594535 | – |
| 150 | 0.694131485985 | – |
| 200 | 0.692648924528 | – |
| 300 | 0.691018903643 | – |

The tuples are in `small_n_tuples.json` as floating-point numbers. The exact binary values of those floats are the certified tuples (`exact_check.py small`). For $n \le 55$ they reproduce the values of Cheer and Goldston to 12 decimals and are **not new**. At $n = 300$ this free search and the smooth-phase construction differ by only $4.6\cdot 10^{-9}$.

## What this shows, and what it does not

- These are rigorous upper bounds for $R_n$ at the listed $n$. As far as we could find, they are the first rigorous upper bounds for specific $n > 55$.
- For every $n$ from 3 000 to 1 000 000 in Table 1, the certified value is below $0.6873358063$, the smallest claimed value of $\limsup R_n$. This finite-$n$ statement does not depend on the unreviewed 2026 asymptotic arguments.
- These bounds say nothing directly about $\limsup R_n$ (constant C42), because $R_n$ is not known to be monotone in $n$. They do not change the C42 entry.
- The values in Table 1 decrease steadily. A floating-point fit $L + C n^{-a}$ to them gives $L \approx 0.68636$. **This is a heuristic extrapolation, not a proven bound.**

## The construction (`turan_phase6.py`)

The construction extends the idea of Harcos's short proof of the bound 0.693677. Write

$$W(x) = \prod_{j=1}^{n} (1 - z_j x) = \exp\Bigl(-\sum_{k \ge 1} \frac{s_k}{k} x^k\Bigr).$$

Then $z = 1$ is in the tuple exactly when $W(1) = 0$.

1. Fix the first $K = n - \lfloor (n-1)/2 \rfloor$ power sums.
2. Every product of two of the remaining power sums $s_{K+1}, \dots, s_n$ has index larger than $n$. So the condition $W(1) = 0$ is **linear** in them, and they are chosen with equal modulus and aligned phases. This choice gives them the smallest possible common modulus.
3. Harcos takes the first block constant. Here it is
$$s_k = c\, e^{i\varphi(k/n)}, \qquad k = 1, \dots, K,$$
   with a constant modulus $c$ and a slowly turning phase $\varphi$. The phase is piecewise linear with 24 pieces, with knots finer near 0.
4. The 26 parameters ($c$ and 25 knot values) are optimized for each $n$ (SLSQP), so that both blocks end up with the same modulus. The parameters for every $n$ are in `smooth_table.json`.

All of this happens in floating point and only **guides** the search. The certified object is the polynomial stored in the certificate:

1. The coefficients $b_0, \dots, b_{n-1}$ of $B(x) = W(x)/(1-x)$ are computed in ball arithmetic.
2. Each one is replaced by the exact dyadic midpoint of its ball.
3. $b_n$ is dropped.

The tuple is then **defined** by $W(x) = (1-x)\,B(x)$. It consists of $z_1 = 1$ (exactly, because of the factor $1 - x$) together with the $n-1$ roots of $Z^{n-1} + b_1 Z^{n-2} + \dots + b_{n-1}$. The roots are never computed. Every power sum follows from $W$ through Newton's identities,

$$s_k = -k\,w_k - \sum_{i=1}^{k-1} w_i\, s_{k-i},$$

where $w_i$ is the coefficient of $x^i$ in $W$.

**Certificate format** (`results/certs/n_<n>.json`). The fields are:
- `n`;
- `rigorous_upper_bound`;
- `b`: a list of $n$ pairs `[re, im]`. Each part is `[numerator, k]`, meaning `numerator / 2^k` exactly.

## How the bounds are certified

- **P**: `turan_phase6.py` computes $s_k = -k\,[x^k]\log W(x)$ with python-flint (Arb) ball arithmetic at 320 bits.
- **V**: `verify_tuple.py` was written separately and shares no code with the generator; it reads only the certificate.
  - It rebuilds $W$ exactly with Python integers and checks $b_0 = 1$ and $W(1) = 0$ exactly.
  - It then computes all $s_k$ from Newton's identities in ball arithmetic, in one of two ways: term by term (`newton`, $O(n^2)$), or as $-W'/W$ with its own power-series inverse built only from polynomial products (`inverse`). The second avoids the library's series division and logarithm.
  - It raises the precision until every ball is narrower than $2^{-80}$.
- **E**: `exact_check.py` uses only Python integers: no floating point and no external libraries.
  - It puts all coefficients over a common denominator $2^D$. Since $w_0 = 1$, Newton's identities need no division.
  - For every $k$ it checks $q^2 |N_k|^2 \le p^2\, 4^{kD}$, where $s_k = N_k / 2^{kD}$ and $p/q$ is the claimed bound.
  - The cost grows like $n^3$. It was run for $n \le 5000$; $n = 5000$ took 3.7 hours.
- For every $n \le 5000$, the exact 20-digit interval from E lies inside the ball from V.

**Different machine and operating system.** The certificates and all logs in `results/` were produced on Windows 10 (Python 3.12, python-flint 0.9.0). The checks were then repeated on Linux (Python 3.13, python-flint 0.9.0):
- V was re-run for $n = 300$, 1 000 (both methods) and 5 000.
- E was re-run for $n = 300$ and 1 000, and in `small` mode.
- The output was identical to the Windows output in every digit.

**The generator is not bit-reproducible.** Re-running `turan_phase6.py certs 300 1000 5000` on Linux gave the same bounds to all 12 decimals. However, the coefficients differed from the Windows certificates by about $3\cdot 10^{-16}$, because of floating-point FFTs. The certificates are the certified objects. The generator is included to document how they were found.

## How to check

You need Python 3 and, for V, `pip install python-flint==0.9.0`. E needs only the standard library.

```bash
# 1. the certificates are the published ones (Windows: certutil -hashfile results\certs\n_300.json SHA256)
cd results/certs && sha256sum -c --ignore-missing SHA256SUMS.txt && cd ../..

# 2. checker V (seconds)
python verify_tuple.py results/certs/n_300.json results/certs/n_1000.json both
python verify_tuple.py results/certs/n_3000.json results/certs/n_20000.json

# 3. checker E (exact): n = 300 takes seconds, n = 1000 about 1.5 minutes, n = 3000 about 40 minutes
python exact_check.py results/certs/n_300.json results/certs/n_500.json

# 4. small n (exact, seconds; rewrites results/small_n_exact.json)
python exact_check.py small
```

V prints `CONFIRMED` and E prints `PASSED` when a certificate satisfies its bound.

**The larger certificates** are for $n$ = 30 000, 50 000, 100 000, 300 000 and 1 000 000 (6.5 MB to 218 MB). They are attached as a zip file to the [GitHub release](https://github.com/Bakeas/turan-power-sum-certificates/releases). Extract them into `results/certs/`. Their SHA-256 sums are already in `results/certs/SHA256SUMS.txt`. V takes about 2 minutes at $n = 10^6$.

**Re-running the generator** (optional; needs numpy, scipy and python-flint):
- `python turan_phase6.py certs 300 1000` rebuilds certificates from the parameters in `smooth_table.json`. The bounds come out the same, but the bytes may differ (see above).
- `python turan_phase6.py run <n>` optimizes a new $n$ and adds it to the table.
- Running it with no argument, or with `big`, restarts the optimization from scratch and **overwrites** `smooth_table.json`.

## Files

| File | Contents |
|---|---|
| `results/certs/n_<n>.json` | certificates for $n \le 20\,000$ (the larger ones are in the release) |
| `results/certs/SHA256SUMS.txt` | SHA-256 of all 13 certificates |
| `results/final_table.md` | full results table, generated from the logs |
| `results/*.txt` | raw outputs: `run_table`, `certs_big` (generator), `verify_*` (V), `exact_*` (E) |
| `results/small_n_exact.json` | output of `exact_check.py small` |
| `verify_tuple.py`, `exact_check.py` | checkers V and E |
| `turan_phase6.py` | generator (construction and check P) |
| `smooth_table.json`, `smooth_start.json` | parameters and bounds for each $n$; starting point of the optimization |
| `small_n_tuples.json` | small-$n$ tuples. These come from an earlier search script, which is not included; the tuples are certified as stored |

## AI use

This project was carried out with substantial help from AI models: **Claude Opus 5.5** and **Claude Sonnet 5.5** (Anthropic). They were used for:
- the literature search;
- designing the construction;
- writing all code;
- drafting this text.

The author planned and directed the work, ran all computations, and checked the results.

Every number above can be re-checked with the short checkers in this repository, or with any other software, directly from the certificates. The mathematics has not yet been reviewed by a professional mathematician. Corrections and comments are very welcome, through a GitHub issue or by email.

## How to cite

See [`CITATION.cff`](CITATION.cff), or use:

```bibtex
@misc{bakeas2026turan,
  author       = {Bakeas, Vasileios},
  title        = {Certified upper bounds in {T}ur\'an's pure power sum problem ({E}rd\H{o}s \#519)},
  year         = {2026},
  howpublished = {\url{https://github.com/Bakeas/turan-power-sum-certificates}},
  note         = {Version 1.0.0}
}
```

## References

- F. V. Atkinson, *On sums of powers of complex numbers*, Acta Math. Acad. Sci. Hungar. 12 (1961), 185–188.
- A. Biró, *On a problem of Turán concerning sums of powers of complex numbers*, Acta Math. Hungar. 65 (1994), 209–216.
- A. Biró, *An improved estimate in a power sum problem of Turán*, Indag. Math. (N.S.) 11 (2000), 343–358.
- A. Biró, *An upper estimate in Turán's pure power sum problem*, Indag. Math. (N.S.) 11 (2000), 499–508.
- G. Harcos, *On a result of András Biró*, note, <https://www.renyi.hu/~gharcos/biro.pdf>.
- A. Y. Cheer and D. A. Goldston, *Turán's pure power sum problem*, Math. Comp. 65 (1996), 1349–1358.
- P. Turán, *On a New Method of Analysis and its Applications*, Wiley, 1984.
- Erdős Problems, problem #519, <https://www.erdosproblems.com/519>.
- T. Tao et al., Optimization Constants repository, constant C42, <https://github.com/teorth/optimizationproblems> (PRs #85, #184, #251, #263; unrefereed).
- The FLINT team, FLINT/Arb and python-flint, <https://flintlib.org>.

## License

- Code (`*.py`): MIT, see [`LICENSE`](LICENSE).
- Data (certificates, tables, logs and JSON files): CC BY 4.0, see [`LICENSE-DATA.md`](LICENSE-DATA.md).
