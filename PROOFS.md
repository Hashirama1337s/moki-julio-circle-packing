# Proof: optimal packings of equal circles in thin fixed rectangles (the zig-zag range) — Moki&Julio, 2026-09-23

## Theorem
Let N ≥ 2 equal circles of radius r lie in the rectangle of width 1 and height h. Put k = (N − 1)² and

  r_zz(N) = ( (1 + k h) − √( (1 + k h)² − 1 − k h² ) ) / 2,

the root of (N − 1)·√(4 h r − h²) + 2 r = 1. **If h(2 − √3) ≤ r_zz(N) ≤ h/2, the largest possible radius is exactly r_zz(N)**, and
it is attained by the zig-zag packing (circles alternately touching the top and the bottom wall).

## Proof
*Upper bound.* Take any packing with r > h/4 (if r ≤ h/4 there is nothing to prove, since r_zz > h/4 in the stated range). The
centres lie in a strip of height h − 2r < 2r, so two centres differ in y by at most h − 2r, and their distance is at least 2r;
hence they differ in x by at least s(r) = √(4r² − (h − 2r)²) = √(4 h r − h²) > 0. Sorting the centres by x, the N − 1 gaps are each
at least s(r) and together at most (1 − 2r), the width available to centres. So (N − 1)·s(r) + 2r ≤ 1. The left side strictly
increases with r, so r ≤ r_zz(N).
*Attainment.* Put x_i = −1/2 + r + i·s and y_i = ±(h/2 − r) alternately (i = 0…N−1) with r = r_zz, s = s(r_zz). Neighbours are at
distance √(s² + (h − 2r)²) = 2r; circles two apart are 2s ≥ 2r apart exactly when s ≥ r, i.e. r ≥ h(2 − √3); the last centre sits
at 1/2 − r by the definition of r_zz; r ≤ h/2 keeps every circle inside. ∎

**The bound is Füredi's, not ours.** Füredi (Discrete Comput. Geom. 6 (1991) 95–106, pp. 96–97) calls the case w ≤ √3/2 trivial
and states it for a closed w × x rectangle R_{w,x}: points at mutual distance ≥ 1 have |x_i − x_j| ≥ √(1 − w²), so R_{w,x} holds
at most 1 + x/√(1 − w²) of them. Scaled to centres in the (1 − 2r) × (h − 2r) box at distance 2r, that is exactly the upper bound
above. (Correction in v1.5: versions 1.0–1.4 said we had found no statement of it for fixed rectangles. That was wrong.)
What is ours is the translation to Packomania's fixed-aspect containers, the window in which the
zig-zag attains the bound, and the table entries it settles. Packomania marks none of the zig-zag entries below as proven (only
single-row entries are bold). Known special cases we do not claim: two circles in any rectangle (classical), and the 2 × 1
rectangle (= crc_500) for N = 3–6, described in the literature as proved by zig-zag arguments. (Ruda 1970, sometimes cited for
n ≤ 8, concerns minimum area with the aspect ratio left FREE — a different problem.)

**Independent check (Grok, xAI):** re-derived the root and every edge case (r ≤ h/4; pairs two or more apart; N = 2), fetched
Packomania's crc_100 / crc_200 pages and matched five published radii to r_zz to ≥ 16 digits (for crc_100 N = 11 the root is
exactly 1/(11 + √119)); verdict CONFIRMED. His fingerprint check: every entry in the range has 2N + 1 contacts, all circles on the
boundary, no rattler; the pattern breaks exactly where r_zz drops below h(2 − √3) (crc_100 N = 37, crc_200 N = 18).

## What it proves in Packomania's tables (check it yourself: `python proofs/proof_zigzag_check.py`, standard library only;
## data with Packomania's published radii and the status of every entry: `proofs/PROOF_ZIGZAG.csv`)
| table | height h | proven optimal here | already known (not claimed) |
|---|---|---|---|
| crc_100 | 0.1 | N = 11–36 (26) | — |
| crc_200 | 0.2 | N = 6–17 (12) | N = 5 (single row, bold on Packomania) |
| crc_300 | 0.3 | N = 4–11 (8) | — |
| crc_400 | 0.4 | N = 3–8 (6) | — |
| crc_500 | 0.5 | — | N = 2 (single row), N = 3–6 (2 × 1 rectangle, prior zig-zag proofs) |
| crc_600 | 0.6 | N = 3–5 (3) | N = 2 (two circles) |
| crc_700 | 0.7 | N = 3–4 (2) | N = 2 (two circles) |
| crc_800 | 0.8 | N = 3 (1) | N = 2 (two circles) |

**58 table entries newly proven optimal** (67 in the zig-zag range in total). At all 67, Packomania's published radius equals
r_zz(N) to 28 digits: the table values are right, and this theorem is why. Beyond the range (e.g. crc_100 from N = 37) the zig-zag
no longer fits and the problem is open.
