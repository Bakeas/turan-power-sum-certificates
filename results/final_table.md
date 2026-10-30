# Final table: rigorous upper bounds for R_n

R_n = min over complex z_1 = 1, z_2, ..., z_n of max_(1<=k<=n) |z_1^k + ... + z_n^k|.
Every bound below comes from one explicit n-tuple; the bound holds for both normalizations
("z_1 = 1" and "z_1 = 1, all |z_j| <= 1"), see README.md. Values with 12 decimals are rounded UP.
Generated from the log files in results/ and smooth_table.json (no numbers typed by hand).

## Table 1. Smooth-phase construction (turan_phase6.py), n = 300 .. 1000000

Checks: P = turan_phase6.py (python-flint balls, log series); V = verify_tuple.py (python-flint balls, Newton identities,
other algorithm, no shared code); E = exact_check.py (exact integer arithmetic, no floating point, no python-flint).
Enclosure = interval containing the exact max_k |s_k| of the certified tuple (exact: 20 decimals; balls: 15 decimals).
Certificates: results/certs/n_<n>.json (SHA-256 in results/certs/SHA256SUMS.txt).

| n | R_n <= | k of max | max_k \|s_k\| in | checks |
|---|---|---|---|---|
| 300 | 0.691018908239 | 300 | [0.69101890823890394806, 0.69101890823890394807] (exact) | P V E |
| 500 | 0.689542310941 | 500 | [0.68954231094014107070, 0.68954231094014107071] (exact) | P V E |
| 1000 | 0.688245154929 | 1000 | [0.68824515492870161684, 0.68824515492870161685] (exact) | P V E |
| 2000 | 0.687468053119 | 2000 | [0.68746805311861023480, 0.68746805311861023481] (exact) | P V E |
| 3000 | 0.687169029139 | 3000 | [0.68716902913847906781, 0.68716902913847906782] (exact) | P V E |
| 5000 | 0.686902345661 | 5000 | [0.68690234566085135502, 0.68690234566085135503] (exact) | P V E |
| 10000 | 0.686673342516 | 5554 | [0.686673342515948, 0.686673342515949] (balls) | P V |
| 20000 | 0.686540042277 | 20000 | [0.686540042276604, 0.686540042276605] (balls) | P V |
| 30000 | 0.686489957913 | 30000 | [0.686489957912097, 0.686489957912098] (balls) | P V |
| 50000 | 0.686446081102 | 50000 | [0.686446081101166, 0.686446081101167] (balls) | P V |
| 100000 | 0.686409189526 | 100000 | [0.686409189525917, 0.686409189525918] (balls) | P V |
| 300000 | 0.686380380702 | 300000 | [0.686380380701450, 0.686380380701451] (balls) | P V |
| 1000000 | 0.686368028005 | 1000000 | [0.686368028004727, 0.686368028004728] (balls) | P V |

For comparison (not rigorous finite-n results): best claimed asymptotic constant limsup R_n < 0.6873358063 (Koshta 2026, unrefereed);
Harcos/Biro 2000: limsup R_n < 0.69368.

## Table 2. Small n: free numerical search (turan_phase5.py Stage 1), checked exactly (exact_check.py small)

The tuple is the exact binary value of each stored float (small_n_tuples.json). "exact max" = max_k |s_k| rounded down to 20 decimals;
"R_n <=" = the same, rounded up to 12 decimals. Cheer-Goldston 1996 (Math. Comp. 65, Table 1, stated as a Theorem) is given where
available; n <= 55 is NOT new. All 8 C-G values agree with ours to 12 decimals; ours is larger by 7e-16 .. 7e-15.

| n | k of max | exact max_k \|s_k\| (20 dec, down) | R_n <= | Cheer-Goldston 1996 |
|---|---|---|---|---|
| 10 | 7 | 0.74265504914476792846 | 0.742655049145 | 0.74265504914476718905 |
| 15 | 6 | 0.72810889170059334140 | 0.728108891701 | 0.72810889170059241298 |
| 20 | 5 | 0.72016991610628733583 | 0.720169916107 | 0.72016991610628663886 |
| 25 | 18 | 0.71508240701664253140 | 0.715082407017 | 0.71508240701664073135 |
| 30 | 19 | 0.71150571445563006928 | 0.711505714456 | 0.71150571445562842882 |
| 40 | 21 | 0.70675192739709053152 | 0.706751927398 | 0.70675192739708693051 |
| 50 | 23 | 0.70369462708690375231 | 0.703694627087 | 0.70369462708689897930 |
| 55 | 35 | 0.70253266322362112834 | 0.702532663224 | 0.70253266322361406713 |
| 70 | 22 | 0.69992744173619015413 | 0.699927441737 | - |
| 100 | 27 | 0.69681859453459515817 | 0.696818594535 | - |
| 150 | 30 | 0.69413148598456627924 | 0.694131485985 | - |
| 200 | 45 | 0.69264892452795809241 | 0.692648924528 | - |
| 300 | 100 | 0.69101890364279679326 | 0.691018903643 | - |

At n = 300 the free-search tuple (0.691018903643) beats the smooth-phase tuple (0.691018908239) by only 4.6e-9.
