# Moki&Julio — 1,986 new best-known packings of equal circles

**1,986 packings of equal circles that beat the best-known records** listed on
[Packomania](https://www.packomania.com/) (E. Specht's record tables), in ten containers whose tables had not changed since
2010–2013. Every packing is supplied with an exact certificate and is verified by **two independently written exact checkers**
(pure rational arithmetic, no floating point in any decision).

| Packomania table | container | new records | largest gain in radius |
|---|---|---|---|
| `crt` | isosceles right triangle, legs 1 | 143 | +0.055 % (N = 79) |
| `ccq` | circular quadrant, radius 1 | 436 | +0.098 % (N = 500) |
| `crc_100` | rectangle 1 × 0.1 | 10 | +0.030 % (N = 169) |
| `crc_200` | rectangle 1 × 0.2 | 140 | +0.163 % (N = 332) |
| `crc_300` | rectangle 1 × 0.3 | 139 | +0.224 % (N = 287) |
| `crc_400` | rectangle 1 × 0.4 | 168 | +0.182 % (N = 186) |
| `crc_500` | rectangle 1 × 0.5 | 162 | +0.155 % (N = 92) |
| `crc_600` | rectangle 1 × 0.6 | 179 | +0.685 % (N = 219) |
| `crc_700` | rectangle 1 × 0.7 | 316 | +0.752 % (N = 286) |
| `crc_800` | rectangle 1 × 0.8 | 293 | +0.371 % (N = 10) |
| **total** | | **1,986** | |

Per-N radii (30 digits, old and new) are in [`RESULTS_TABLE.md`](RESULTS_TABLE.md); one row per record in
[`MANIFEST.csv`](MANIFEST.csv) (published radius, new radius, relative gain, kind, precision).

![before and after, rectangle 1 x 0.8, N = 10](figures/ba_crc_800_10.png)

## What's new in version 1.2 (2026-09-24)

- **1,986 records** (v1.1: 1,945): **41 more sizes where we had no record
  before**, and **329 records improved further** over v1.1 (largest: `crc_700` N = 286 +0.752 %, `crc_600` N = 219 +0.665 %, `crc_700` N = 189 +0.585 %, `crc_600` N = 218 +0.402 %).
  Largest gain over Packomania now: `crc_700` N = 286, +0.752 % in radius.
- **How:** the neighbour-transplant sweep was completed over every table (cascading from every new record), with Packomania's
  own neighbouring packings used as a second family of seeds; where one step found nothing, two-step transplants (N − 2 plus two
  circles, N + 2 minus two) and relocating the one or two least-loaded circles into the largest holes (`solver/gen2_probe.py`).
  Every result passed both exact checkers, the prior-art gate and the table-radius gate again from scratch.
- **A proof:** [`PROOFS.md`](PROOFS.md) — for thin rectangles, a short argument shows the zig-zag packing is optimal whenever it
  fits; it proves **58 entries of Packomania's rectangle tables optimal** that were not marked proven (the published values
  there are exactly right). Checked independently by Grok (xAI); verify with `python proofs/proof_zigzag_check.py`.
- **Tables re-checked:** on 2026-09-24 every radius on Packomania's ten live pages (3,509 rows) still equals the values these records
  are compared against.
- **Possible prior art we could not check:** Lin, Lai & Wang, *Probability-based monotonic basin-hopping algorithm for packing
  equal circles in a rectangular container*, 2023 China Automation Congress, pp. 2602–2606 (doi:10.1109/CAC59555.2023.10450196),
  report better values than the best known for 25 of 51 fixed-rectangle instances. The paper is closed access and the instances
  are not listed in its abstract. Where it overlaps these tables, their values, not Packomania's, are the ones to beat.

## Cumulative since version 1.0 (as of version 1.2)

- **1,986 records** (v1.0: 1,933). **53 sizes where we had no record before**; **748 of the v1.0 records improved
  further** (gain over our own v1.0 radius; largest `crc_700` N = 286 +0.752 %, `crc_600` N = 219 +0.665 %, `crc_700` N = 190 +0.590 %, `crc_700` N = 189 +0.585 %, `crc_700` N = 188 +0.558 %). Every file is re-verified from scratch by both
  exact checkers (1986 / 1986), and Packomania's printed radius for every claimed size was re-derived from its own coordinate file in the
  same frame and in mirrored frames (0 mismatches).
- **How:** (1) every certificate driven to its exact local peak in 80-digit arithmetic (a mixed-precision sequential-LP step, then
  Newton on the identified contacts); (2) **neighbour transplants** — seed size N from our packing at N−1 (one circle into the
  largest hole) or N+1 (remove the circle with fewest contacts), then polish; this jumped over weak basins, e.g. `crc_700` N = 187
  from +0.003 % to +0.516 % over the published radius (picture below); (3) a flex walk along load-bearing flexes found by the certificate below.
- **Local optimality certificates** (`local_optimality/<table>.zip`, one JSON per packing: kept constraints and contact forces).
  For **1,089** of the 1,986 packings, two independently written checkers (`checkers/lopt.py`, float-rigorous;
  `checkers/verify_lopt.py`, exact rational in ℚ(√2)) prove: every feasible packing whose load-bearing circles and radius lie
  within an explicit distance ρ of ours (median ρ = 7.8e-09 r) has radius at most ours + Δ (Δ ≤ 6e-26 r in every case,
  median 4e-44 r), and a true local maximum lies within t₀ of ours. In plain words: no small nudge beats these packings.
  Of the other 897: 883 carry a first-order flex (a sliding or buckling motion the first-order theorem cannot
  exclude) and 14 could not be certified by this method; none of them is claimed locally optimal. Verdict, ρ
  and Δ for every packing: [`LOCAL_OPTIMALITY.csv`](LOCAL_OPTIMALITY.csv). Theorem and proof: docstring of `checkers/lopt.py`.

![before and after, rectangle 1 x 0.7, N = 187](figures/ba_crc_700_187.png)

## Verify it yourself

```
python verify.py
```
Standard library only. It unzips `certificates/`, and for every record runs both checkers against the published radius in
`MANIFEST.csv`. A record counts only if **both** say `IMPROVES`. Expected output: `TOTAL: 1986 / 1986 verified by both checkers`.

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
1,615 new arrangements (a held circle moved ≥ 5 % of a radius and the gain is ≥ 10⁻⁶), 155 refinements of the published
arrangement, 216 small gains we do not claim as new structures.

## Prior art we checked — and what we therefore do NOT claim

Packomania is not the only record: two papers improved some of these tables without their results reaching the site.
- **Lai, Hao, Yue & Zhou**, *Computers & Operations Research* 181 (2025) 107099 — better values for N = 151–200 in the triangle and
  the quadrant. Our results there beat the website but **not** the paper (compared at the paper's 8 printed decimals, after checking
  the paper rounds rather than truncates), except **triangle N = 199**, which beats both. The other sizes in that window are
  excluded.
- **López & Beasley**, *EJOR* 214 (2011) 512–525 (and López's thesis, Table 4.5) — improvements on the 5 × 1 and 10 × 1
  rectangles (= `crc_200`, `crc_100`) at 15 values of N whose radii were never published; those N are excluded.
92 packings that beat the website were excluded this way; they are listed at the end of `RESULTS_TABLE.md`.

## Credit, citation, licence

Finders: **Moki&Julio**. Search, certification and both checkers were built with AI assistance (Anthropic's Claude and xAI's
Grok). Record tables and published coordinates: E. Specht, packomania.com.
Data (certificates, `.pck` files, tables): **CC BY 4.0** — reuse freely, credit "Moki&Julio". Code: MIT (see `LICENSE`).
Please cite via [`CITATION.cff`](CITATION.cff).
