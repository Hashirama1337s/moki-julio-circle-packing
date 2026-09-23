# Moki&Julio — 1,945 new best-known packings of equal circles

**1,945 packings of equal circles that beat the best-known records** listed on
[Packomania](https://www.packomania.com/) (E. Specht's record tables), in ten containers whose tables had not changed since
2010–2013. Every packing is supplied with an exact certificate and is verified by **two independently written exact checkers**
(pure rational arithmetic, no floating point in any decision).

| Packomania table | container | new records | largest gain in radius |
|---|---|---|---|
| `crt` | isosceles right triangle, legs 1 | 143 | +0.055 % (N = 79) |
| `ccq` | circular quadrant, radius 1 | 433 | +0.098 % (N = 500) |
| `crc_100` | rectangle 1 × 0.1 | 8 | +0.011 % (N = 87) |
| `crc_200` | rectangle 1 × 0.2 | 132 | +0.080 % (N = 222) |
| `crc_300` | rectangle 1 × 0.3 | 135 | +0.113 % (N = 20) |
| `crc_400` | rectangle 1 × 0.4 | 158 | +0.175 % (N = 272) |
| `crc_500` | rectangle 1 × 0.5 | 156 | +0.090 % (N = 143) |
| `crc_600` | rectangle 1 × 0.6 | 174 | +0.135 % (N = 122) |
| `crc_700` | rectangle 1 × 0.7 | 315 | +0.590 % (N = 190) |
| `crc_800` | rectangle 1 × 0.8 | 291 | +0.371 % (N = 10) |
| **total** | | **1,945** | |

Per-N radii (30 digits, old and new) are in [`RESULTS_TABLE.md`](RESULTS_TABLE.md); one row per record in
[`MANIFEST.csv`](MANIFEST.csv) (published radius, new radius, relative gain, kind, precision).

![before and after, rectangle 1 x 0.8, N = 10](figures/ba_crc_800_10.png)

## What's new in version 1.1 (2026-09-23)

- **1,945 records** (v1.0: 1,933). **12 sizes where we had no record before**; **576 of the v1.0 records improved
  further** (gain over our own v1.0 radius; largest `crc_700` N = 190 +0.590 %, `crc_700` N = 188 +0.558 %, `crc_700` N = 187 +0.511 %, `crc_700` N = 329 +0.269 %, `crc_800` N = 205 +0.254 %). Every file is re-verified from scratch by both
  exact checkers (1945 / 1945), and Packomania's printed radius for every claimed size was re-derived from its own coordinate file in the
  same frame and in mirrored frames (0 mismatches).
- **How:** (1) every certificate driven to its exact local peak in 80-digit arithmetic (a mixed-precision sequential-LP step, then
  Newton on the identified contacts); (2) **neighbour transplants** — seed size N from our packing at N−1 (one circle into the
  largest hole) or N+1 (remove the circle with fewest contacts), then polish; this jumped over weak basins, e.g. `crc_700` N = 187
  from +0.003 % to +0.514 % over the published radius (picture below); (3) a flex walk along load-bearing flexes found by the certificate below.
- **Local optimality certificates** (`local_optimality/<table>.zip`, one JSON per packing: kept constraints and contact forces).
  For **1,051** of the 1,945 packings, two independently written checkers (`checkers/lopt.py`, float-rigorous;
  `checkers/verify_lopt.py`, exact rational in ℚ(√2)) prove: every feasible packing whose load-bearing circles and radius lie
  within an explicit distance ρ of ours (median ρ = 6.5e-09 r) has radius at most ours + Δ (Δ ≤ 6e-26 r in every case,
  median 4e-44 r), and a true local maximum lies within t₀ of ours. In plain words: no small nudge beats these packings.
  Of the other 894: 878 carry a first-order flex (a sliding or buckling motion the first-order theorem cannot
  exclude) and 16 could not be certified by this method; none of them is claimed locally optimal. Verdict, ρ
  and Δ for every packing: [`LOCAL_OPTIMALITY.csv`](LOCAL_OPTIMALITY.csv). Theorem and proof: docstring of `checkers/lopt.py`.

![before and after, rectangle 1 x 0.7, N = 187](figures/ba_crc_700_187.png)

## Verify it yourself

```
python verify.py
```
Standard library only. It unzips `certificates/`, and for every record runs both checkers against the published radius in
`MANIFEST.csv`. A record counts only if **both** say `IMPROVES`. Expected output: `TOTAL: 1945 / 1945 verified by both checkers`.

- **Certificate format** (`certificates/<table>.zip`, one file per N): first line `r <radius>`, then N lines `x y` (circle centres,
  decimal). Same coordinate frame as Packomania's own files: `crt` has its right angle at the origin and legs along the axes;
  `ccq` is the unit quarter disc at the origin; `crc_k` is width 1 × height 0.k centred at the origin.
- **Checker A** (`checkers/certify_circ.py`): every circle inside the container, every pair at distance ≥ 2r, all in exact
  fractions (the only irrational terms, √2 on the triangle's hypotenuse and the quadrant's arc, are handled by a sign check and
  squaring). Claim floor: IMPROVES only if the radius exceeds the published one by more than one part in 10¹⁰ — published values can
  be low by ~10⁻²⁴ from last-digit rounding, which is a tie, not a record.
- **Checker B** (`checkers/verify_exact_crt.py`, `checkers/verify_exact_rect_quad.py`): written independently, blind to checker A.
- `packomania/<table>.zip` holds the same packings in Packomania's submission format (`.pck`: radius, name, coordinates).

## How they were found

Search over the positions of N centres maximising the smallest clearance (circles to walls, circles to each other):
1. **Local solver** — sequential linear programming with a trust region (HiGHS): linearise every near-active distance and wall,
   solve the LP, accept the step only if the true minimum clearance grows (`solver/slp.py`, `solver/slp_circ.py`).
2. **Global search** — monotonic basin hopping (Grosso, Jamali, Locatelli, Schoen): perturb all circles a little or one region a
   lot, re-solve, keep improvements, kick on stagnation; many chains in parallel; seeded with the published record, then with our own
   best (`solver/campaign*.py`, `solver/quickpass.py`, `solver/deeppass.py`).
3. **High precision** — mixed-precision Gauss–Newton on the contact graph to residual < 10⁻⁵⁰; certificates written at 45 digits
   where the refinement converges (`solver/refine*.py`).
4. **Positive controls before any claim** — all published triangle records re-derived from their own coordinates to 10⁻³⁰; the
   search re-found 14 of 14 known small records; both checkers reject planted overlaps, points outside, and inflated radii.

**Kinds of improvement** (our circles matched one-to-one to the published ones; loose "rattler" circles ignored):
1,526 new arrangements (a held circle moved ≥ 5 % of a radius and the gain is ≥ 10⁻⁶), 174 refinements of the published
arrangement, 245 small gains we do not claim as new structures.

## Prior art we checked — and what we therefore do NOT claim

Packomania is not the only record: two papers improved some of these tables without their results reaching the site.
- **Lai, Hao, Yue & Zhou**, *Computers & Operations Research* 181 (2025) 107099 — better values for N = 151–200 in the triangle and
  the quadrant. Our results there beat the website but **not** the paper (compared at the paper's 8 printed decimals, after checking
  the paper rounds rather than truncates), except **triangle N = 199**, which beats both. The other sizes in that window are
  excluded.
- **López & Beasley**, *EJOR* 214 (2011) 512–525 (and López's thesis, Table 4.5) — improvements on the 5 × 1 and 10 × 1
  rectangles (= `crc_200`, `crc_100`) at 15 values of N whose radii were never published; those N are excluded.
91 packings that beat the website were excluded this way; they are listed at the end of `RESULTS_TABLE.md`.

## Credit, citation, licence

Finders: **Moki&Julio**. Search, certification and both checkers were built with AI assistance (Anthropic's Claude and xAI's
Grok). Record tables and published coordinates: E. Specht, packomania.com.
Data (certificates, `.pck` files, tables): **CC BY 4.0** — reuse freely, credit "Moki&Julio". Code: MIT (see `LICENSE`).
Please cite via [`CITATION.cff`](CITATION.cff).
