# LOCAL — the local stage for N = 16 (design, 2026-09-26 20:20 ADT; results appended below)

Status words: m16 stays "best known" until the global stage, this local stage and two independent checkers all pass.

## 0. What the global stage gives, what the local stage must add

Global (TREE_FORMAT.md §4, cluster run at m_t = m16 (1 − ε)): every m_t-feasible configuration has, after a symmetry
and relabelling, its 15 frame points inside a box of radius ρ_B (barycentric sup-norm) around one of the three listed
configurations (two distinct best-known optima; the third is the mirror image of the second). Tonight: ρ_B = 5e-3,
ε = 1e-5.

Local (this file): inside a ball of radius R around an accurate frame c̃, (i) every configuration whose value is at
least m_lo lies within δ of c̃ (δ tiny), and (ii) its value is at most m_hi. If the global boxes lie inside the local
balls and m_t ≤ m_lo, then m16 ∈ [m_lo, m_hi] and every optimal packing is, up to symmetry, within δ of one of the
three frames (the 16th point, the rattler, is not located by the certificate; see README.md §6). That is Markót's type of result (enclosure of the value and of all
optimisers), with our δ and m_hi − m_lo expected near 1e-25.

## 1. Coordinates and constraints (frame only; the rattler only adds constraints, so it can be ignored)

Skew coordinates of a point: (u, v) = (A/G, B/G) — the first two barycentric coordinates. Triangle: u ≥ 0, v ≥ 0,
1 − u − v ≥ 0. Squared distance: Q(du, dv) = du² + du·dv + dv². Frame x ∈ R^30 (15 points).

Tight set T at the frame: 20 pairs and 13 walls (33 constraints on 30 coordinates; float check `local16.py`):
- pair k = (i, j):  g_k(x; m) = Q(x_i − x_j) − m²  ≥ 0;  gradient w.r.t. x_i: (2du + dv, du + 2dv), w.r.t. x_j: minus that;
  exact expansion g_k(x̃ + d; m) = g_k(x̃; m) + ∇g_k(x̃)·d + Q(d_i − d_j), and 0 ≤ Q(d_i − d_j) ≤ 12‖d‖∞².
- wall k: w_k(x) ∈ {u_i, v_i, 1 − u_i − v_i} ≥ 0, linear (no second-order term).

## 2. The certificate (all exact rational arithmetic; floats only propose numbers)

Data: a rational frame c̃ (Newton-refined to ~1e-30, then rounded; wall coordinates exactly on their walls), and
nonnegative rational vectors found by float LPs and then re-checked exactly:
- a stress λ ≥ 0 on T with Σ_pairs λ_k = Λ > 0, residual r = Jᵀλ (J = 33 × 30 Jacobian at c̃, exact);
- for every coordinate i ∈ {0..29} and sign σ = ±1, a vector y^{iσ} ≥ 0 on T with residual
  ε^{iσ} = Jᵀ y^{iσ} + σ e_i (exact, tiny).

**Lemma L (localisation).** Let x = c̃ + d be any frame inside the triangle with all 20 tight pairs at distance ≥ m.
Put a_k = ∇g_k(c̃)·d. Then a_k ≥ −Q_k(d) − g_k(c̃; m) for every k ∈ T (walls: Q = 0). For each (i, σ):

    σ d_i = −Σ_k y_k a_k + ε·d  ≤  Σ_k y_k (Q_k(d) + g_k(c̃; m)) + ‖ε‖₁‖d‖
          ≤  12 (Σ_pairs y_k) ‖d‖² + Σ_k y_k g_k(c̃; m) + ‖ε‖₁ ‖d‖.

With A = 12 max_{iσ} Σ_pairs y_k, B(m) = max_{iσ} Σ_k y_k g_k(c̃; m), E = max ‖ε^{iσ}‖₁:

    ‖d‖∞ ≤ A‖d‖∞² + B(m) + E‖d‖∞.

So if ‖d‖∞ ≤ R := (1 − E)/(2A), then ‖d‖∞ ≤ δ(m) := 2B(m)/(1 − E) (and if B(m) < 0, no such x exists).

**Lemma V (value).** For the same x with value m(x) (so every g_k(x; m(x)) ≥ 0):
Σ λ_k g_k(x; m(x)) ≥ 0 gives Φ − Λ m(x)² + r·d + Σ λ_k Q_k(d) ≥ 0 with Φ = Σ_pairs λ_k Q(c̃_i − c̃_j) + Σ_walls λ_k w_k(c̃), so

    m(x)² ≤ ( Φ + ‖r‖₁ δ + 12 Λ δ² ) / Λ =: m_hi².

**Lower bound.** m_lo² := the exact minimum of Q over all 120 pairs of the 16-point configuration c̃ + a rational
rattler position (all points checked inside the triangle). This is attained, so m16 ≥ m_lo.

**Theorem (conditional on the global stage).** If every global box lies in the R-ball of its c̃ (‖P_q/G − c̃‖∞ +
ρ_q/G ≤ R on every frame point) and the global target m_t ≤ m_lo, then m16 ∈ [m_lo, max_q m_hi(q)] and every optimal
packing has its frame within δ(m_lo) of some c̃ (after a symmetry), the rattler anywhere feasible.

Why this beats the rectangle proofs' Lemma 1: that lemma uses one fixed 2N-row basis and the smallest stress entry
(ρ0 ≈ η/(24Λ) ≈ 3e-7 here, `local16.py`); the LP picks, for every coordinate separately, the cheapest nonnegative
combination of the 33 constraint gradients, which is expected to give R ~ 1e-4 (to be measured).

## 3. Closing the gap

The global boxes (ρ_B = 5e-3) must shrink below R. The global run's feasible set at m16 (1 − ε) has size ≈ C·ε with a
per-cluster constant C (measured tonight on the two tile sets of the optima). Choose ε so that C·ε plus the centre
offset is below R with margin, rerun the global cluster stage with boxes of radius ρ_B ≈ R/2 (int64 grid step 2.4e-8
is far below 1e-4, so no big integers are needed at that scale), then replay.

## 4. Prototype plan (tonight)

1. `local_cert.py`: Newton refinement in mpmath (60 digits) of the 18 free unknowns (30 coordinates − 13 wall-fixed,
   + m²) on the 20 pair equations (over-determined by 2: consistency shows as a tiny residual); rational rounding;
   exact J; float LPs for λ and the 60 vectors y; exact residuals and the numbers A, B, E, R, δ, m_lo, m_hi per
   cluster; JSON output with every rational so a checker can recompute.
2. Measure C: run the hard tile sets at m16 (1 − ε) for ε = 1e-6, 1e-7 with boxes of radius R/2; record nodes.
3. If the hard tile sets accept, run the full global cluster stage at that ε (2 processes, ~10 min) and replay it.

---------------------------------------------------------------------------------------------------------------------
## RESULTS (appended, dated)

### L1 — 20:16: local certificates (`local_cert.py`, output `out/local_cert16.json`; exact rationals, floats only propose)

| cluster | tight | Newton residual | A | E | R (skew sup-norm) | δ(m_lo) | m_hi − m_lo |
|---|---|---|---|---|---|---|---|
| 0 | 20 pairs + 13 walls | 2.4e-62 | 4676.8 | 2.9e-14 | 1.069e-4 | 1.2e-31 | 3.6e-34 |
| 1 | 20 + 13 | 1.9e-62 | 3393.1 | 1.3e-14 | 1.474e-4 | 8.7e-32 | 3.6e-34 |
| 2 (mirror of 1) | 20 + 13 | 2.4e-62 | 3040.3 | 2.9e-14 | 1.645e-4 | 7.8e-32 | 3.6e-34 |

m_lo = 0.2162272693097818217346349753963486148701 (exact 16-point rational configuration: frame c̃ rounded to 2^-110,
rattler at its float cage optimum rounded to 2^-40; the same value for clusters 0 and 1 — the two distinct optima have
the same value to all printed digits). m_hi = 0.2162272693097818217346349753963489711416 for all three.
(A first pass used the rattler's grid position from the cluster file; in clusters 1 and 2 that position was closer
than m* to a frame point, which lowered m_lo to 0.21622726924 and inflated δ to 1.6e-8; fixed by placing the rattler at
its cage optimum.)

**The Lemma-1 radius of the rectangle proofs would have been ρ0 ≈ 3e-7; the LP localisation gives R ≈ 1.1e-4, 360×
larger.** Gap to the global boxes: 5e-3 vs 1.07e-4 → the global boxes must shrink by ~50×.

### L2 — 20:26–20:29: shrinking the global boxes to ρ = 4000 grid units (9.54e-5 < R_min − offset 1e-8)

`out/clusters16_rho4000.json` (same three clusters, ρ = 4000). The two cluster tile sets at m16 (1 − ε):
- tile set of cluster 0: ε = 1e-6 → accepted, 329 nodes; ε = 3e-6 → capped at 60,000 (feasible set wider than the box,
  as the LP bound 2B ≈ 73ε·… predicts).
- tile set of clusters 1 + 2 (mirror-symmetric): first capped (60,000 nodes at 1e-6; > 9 min at 3e-7) because the
  engine's "never split the rattler" rule was global and the two clusters of this tile set have different rattlers,
  so neither was exempt and the rattler was bisected down to 1e-5. Fix (heuristic only, 20:28): exempt the rattler of
  the cluster the node is closest to. Then ε = 1e-6 → accepted, 15,793 nodes, 22 s; ε = 3e-7 → accepted, 14,863
  nodes, 18 s.

### L3 seal — 20:30 (written before launching): full global cluster stage at ε = 3e-7, ρ = 4000

`run_acc2.sh 0|1`: m_t = m16 (1 − 3e-7) rounded down, all 341,662 representatives, 2 processes, trees
`trees/n16_acc2_part{0,1}.jsonl.gz`. Prediction: every leaf discarded or accepted (PROVED MODULO CLUSTERS) in < 20 min
wall per part; replay PASS. If it passes, the chain global (L3) + local (L1) gives m16 ∈ [m_lo, m_hi]
(width 3.6e-34) with all optimal frames within 1.2e-31 of the three c̃ — to be kept as "best known" until a second
independent checker has passed both stages.

### L3 result — 20:31 (runs), 20:37 (replay)

| part | representatives | status | nodes | accepted leaves | wall s |
|---|---|---|---|---|---|
| 0/2 | 170,831 | PROVED (no cluster tile set here) | 173,295 | 0 | 72.4 |
| 1/2 | 170,831 | PROVED MODULO CLUSTERS | 188,869 | 49 | 102.7 |

m_t = 216227204441601/10^15 (= m16 (1 − 3e-7) rounded down). Replay (`replay_tri.py --join`, streaming, hardened):
**PASS** — 362,164 nodes, 28,363,056 strips, 351,864 discards, 10,251 splits, 49 accepted leaves, all 2,042,975 tile
subsets covered. The tree headers carry exactly the clusters of `out/clusters16_rho4000.json` (checked). Prediction
(< 20 min per part, replay PASS): met.

Exact local re-check (`check_local.py out/local_cert16.json out/clusters16_rho4000.json 216227204441601/10^15`,
standard library only, recomputes every number from the stored rationals): **PASS** for all three clusters —
recomputed = stored; box offset + ρ = 9.5377e-5 ≤ R (1.069e-4 / 1.474e-4 / 1.645e-4); m_t ≤ m_lo; m_hi² − m_lo² =
1.54e-34.

Accept-fate mutations on the N = 4 example: ρ shrunk, a cluster point moved, an unknown cluster cited — all rejected;
"make the exempt point a frame point" was accepted, correctly (that region was inside its box anyway, so the stronger
claim is still true).

### Combined statement (NOT a public claim; keep m16 as "best known")

Given (i) the global cluster stage L3 (tile/orbit coverage + trees, replayed) and (ii) the local certificates L1
(exact re-check):

- **m16 ∈ [m_lo, m_hi]**, m_lo = 0.21622726930978182173463497539634861487…, m_hi − m_lo = 3.6e-34
  (both ends are exact rationals in `out/local_cert16.json`; m_lo is attained by an exactly checked configuration);
- **every optimal packing of 16 points** (equivalently, 16 equal circles of radius r16 = m/(2(1 + √3 m)) in the unit
  triangle) is, up to the triangle's symmetries and relabelling, one of the two known frames within 1.2e-31 in every
  coordinate, with the 16th point (the rattler) anywhere it fits. The two frames have the same value to 34 digits; the
  certificate does not decide whether their exact values coincide (both are enclosed in the same 3.6e-34 interval).

This is the Markót-type result (enclosure of the value and of all optimisers), here with widths ~1e-31 rather than
1e-12..1e-15, because the local stage is an exact-rational LP certificate instead of interval Newton.

### What still stands between this and a claim

1. **A second, independently written checker of the global stage** (from TREE_FORMAT.md alone; the reviewer has it)
   run on `trees/n16_acc2_part{0,1}.jsonl.gz` (and, for the ε-bound, `n16_pos_*`).
2. **An independent check of the local stage** from LOCAL.md §2 alone (inputs: `out/local_cert16.json`,
   `out/clusters16_rho4000.json`); ours is `check_local.py` (same author as the certificate — not independent).
3. Adversarial read of the two lemmas (L and V) and of the global→local link (frame labels = tile order; symmetry
   applied before the box claim; the skew coordinates (u, v) = (A/G, B/G) and Q = du² + du dv + dv² equal the
   barycentric distance ½(dA² + dB² + dC²)/G²).
4. Housekeeping: the engine's `--timeout` is only checked between tile sets (a tile set can overrun it), and the
   split heuristic changed at 20:28 (decisions unchanged; replay is the arbiter).
