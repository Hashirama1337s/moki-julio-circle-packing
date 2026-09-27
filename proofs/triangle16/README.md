# 16 equal circles in an equilateral triangle: a computer-assisted proof of optimality

Moki&Julio. Version 3, 2026-09-26, revised after two independent reviews (the second found no soundness defect; its
fixes are applied). Published 2026-09-27 in repository version 2.7. All file paths are relative to this folder (`proofs/triangle16/`).

## 1. Result

**Problem (point form).** Place 16 points in the closed equilateral triangle T of side 1 so that the minimum pairwise
distance is as large as possible. Let m16 be that maximum. It exists because T^16 is compact.

**Circle form.** Sixteen non-overlapping circles of radius r fit in T if and only if their centres, at mutual distance
≥ 2r, fit in the inner triangle whose sides are parallel to those of T at distance r. That inner triangle has side
1 − 2√3·r. Scaling it to side 1 turns the centre problem into the point problem with m = 2r/(1 − 2√3·r). So the
largest radius r16 of 16 equal circles in the triangle of side 1 is

    r16 = m16 / (2(1 + √3·m16)),      equivalently  m16 = 2·r16 / (1 − 2√3·r16).

The side of the smallest equilateral triangle holding 16 unit circles is 1/r16 = 2/m16 + 2√3. From the Theorem's
bounds (outward rounding, 30 places) it lies in [12.713628774150546014842082430654, 12.713628774150546014842082430655].

**Theorem.** Let m_lo and m_hi be the positive square roots of the rationals m_lo² and m_hi² listed in §3. They are
copied from `out/local_cert16.json`. Then

    m_lo ≤ m16 ≤ m_hi,        m_hi − m_lo < 3.6 × 10⁻³⁴.

Rounded outward (lower end down, upper end up), with each rounding checked by exact integer square roots:

    0.2162272693097818217346349753963486148700 ≤ m16 ≤ 0.2162272693097818217346349753963489711416
    0.078655749492482286188348450524856119     ≤ r16 ≤ 0.078655749492482286188348450524856215

Take any 16 points in T with minimum distance ≥ m_lo; every optimal configuration qualifies. For such a configuration
there are a symmetry of T, a relabelling of the points and a cluster q ∈ {0, 1, 2} (see §4) such that the following
holds. Each of its 15 frame points lies within δ_q of the matching point of the certified rational frame c̃_q, measured
in the sup-norm on the 30 skew coordinates (u_i, v_i) of §3. Here δ_q ≤ 1.21 × 10⁻³¹ (exact values in §3). In the plane
each frame point is then within √3·δ_q of its certified position. The certificate says nothing about the 16th point
beyond the facts stated in §6.

**Not claimed.** The certificate encloses both known frames (clusters 0 and 1) in the same interval. It does not decide
whether their exact values are equal, so it does not say which of them is optimal (or whether both are). It gives no
information about the 16th point beyond §6, and it proves nothing for any N other than 16.

**Comparison with the best-known value.** Graham & Lubachevsky (1995) print m16 = 0.216227269309782 (15 decimals). The
enclosure above, 0.21622726930978182…, rounds to that value. Their value gives r = 0.07865574949248233… in circle form.
Rounded to 16 decimals, that and the r16 interval above agree: 0.0786557494924823.

## 2. Prior work, and what is new here

- **The configurations are not new.** Melissen & Schuur (1995) found the best packings of 16 circles in an equilateral
  triangle and conjectured that they are optimal. They computed the value by solving a high-degree equation derived from
  the contact graph, and noted that, if the conjectures for 13–15 hold, N = 16 is the first N whose optimum is not
  symmetric. Graham & Lubachevsky (1995) give m16 = 0.216227269309782 and two distinct best configurations. They label
  them t16a33.1 and t16a33.2. Each has 33 bonds and one free point, and the two differ in the placement of four points.
- **Our two frames match that description.** Neither paper prints coordinates, so the match rests on these properties.
  Our two frames (clusters 0 and 1) have the same value to the 15 printed digits, 20 point–point plus 13 point–side
  contacts (33) and one free point each. They also differ in exactly four frame points. We do not claim which of our
  frames carries which label.
- **Our contribution** is the enclosure and the localisation stated in §1: an enclosure of m16 of width < 3.6 × 10⁻³⁴,
  and the localisation of every optimal frame to within ~10⁻³¹ of the certified frames. The packings themselves come
  from the literature above.
- **Status before this work.** Optimality is proven for N ≤ 15 (for example N = 13 by Joós, Aequationes Math. 95, 2021), for
  N = 20 (Payan 1997, one less than a triangular number) and for every triangular number (Oler 1961). For N = 16 we found no proof of optimality and no computer-assisted
  attempt. Sources searched:
  - Wikipedia, "Circle packing in an equilateral triangle", read 2026-09-20 and 2026-09-26;
  - E. Friedman's *Circles in Triangles* page, read 2026-09-26;
  - an arXiv / OpenAlex / GitHub / thesis-repository search for 2022–2026;
  - works citing Melissen & Schuur, Graham & Lubachevsky and Joós.

  This is a search result, not a proof that no such work exists.

## 3. Coordinates and exact constants

- **Barycentric coordinates.** A point of T is (a, b, c) with a + b + c = 1 and a, b, c ≥ 0. The squared distance is
  (da² + db² + dc²)/2. The global stage uses the integer scale A = G·a, B = G·b, C = G·c with
  G = 41,943,040 = 5·2²³. Only region bounds are integers. Points are real, so the global stage covers all real
  configurations, not only grid points.
- **Skew coordinates** (local stage): (u, v) = (a, b). The squared distance is Q(du, dv) = du² + du·dv + dv², which
  equals the barycentric formula because dc = −du − dv. A frame is x ∈ ℝ³⁰ (15 points), and ‖d‖∞ is the maximum of the
  30 coordinate differences. Under ‖d‖∞ ≤ δ, the third coordinate c moves by at most 2δ, and the Euclidean displacement
  of each point is at most √3·δ because Q(du, dv) ≤ 3δ².
- **Floating point only proposes numbers:** search order, Newton refinement of the frames and LP multipliers. Every
  inequality the proof uses is decided in exact integer or rational arithmetic by the checkers (§7).

Exact constants. Every value below is copied verbatim from the certificate JSON (`data[q].<field>`) or from the tree
headers:

```
m_t    = 216227204441601/1000000000000000     localisation target (trees n16_acc2_*), Tc = ceil(2 m_t^2 G^2) = 164501730816426
m_R    = 216248892036713/1000000000000000     refutation target   (trees n16_pos_*),  Tc = 164534731528431
rho    = 4000 grid units;  rho/G = 25/262144 = 9.5367431640625e-5 exactly

m_lo^2 = data[q].m_lo_sq (identical for q = 0, 1, 2)
       = 78780725062457142723775750437274854467981958748910328075087958753/1684996666696914987166688442938726917102321526408785780068975640576
         (the denominator is 2^220)
m_hi^2 = max_q data[q].m_hi_sq = data[0].m_hi_sq
       = 189460364811230778998685563791506840018731043812836804000181479039011730203219217207926803273080568531714300184609003215679280030437717778058101121099385396650547396340458956937/4052261297735108269868151188595835382615305236095432576454628717238084918878894019225311571888031359594454393616101391511144895802482728218239293085932720848234956218998357753856

R_0     = 365375409332715031609072699439616544708531407547/3417565876866576205690153644220314335076199465746432
delta_0 = 14243233175053878650161252140381182304882219612487/118571099379008312433120888636726253181249779798780799398704089889065957727404032
R_1     = 60895901555453516224721775082993128220639901357/413249729118746453931010903380737837277313952120832
delta_1 = 27502051362623804858098825784602288254240717203681/316189598344027415674804661474051795459312346306813624052860516445660578182070272
R_2     = 487167212443620386689065017189928713496880269529/2962239869408221204339548273513094990381531105591296
delta_2 = 24661143194409709060499408076643219115965333705185/316189598344022390109926413097008194301675368037067929017427254159434302355406848
```

The certificate's `summary` block holds floating-point previews only; no checker reads it. Its `m_lo` string
(…6148701) is rounded to nearest, so it is slightly above the true m_lo and is not a lower bound. The decimals in §1 are
the directed ones. Its `global_box_fits: false` refers to the superseded boxes of half-width 5·10⁻³
(`out/clusters16.json`), not to the boxes used here.

| q | A_q | E_q | R_q | δ_q = δ_q(m_lo) | m_hi,q² − m_lo² |
|---|---|---|---|---|---|
| 0 | ≈ 4676.79 | ≈ 2.93e-14 | ≈ 1.069110e-4 | ≈ 1.201240e-31 | ≈ 1.5407e-34 |
| 1 | ≈ 3393.08 | ≈ 1.27e-14 | ≈ 1.473586e-4 | ≈ 8.697962e-32 | ≈ 1.5407e-34 |
| 2 | ≈ 3040.27 | ≈ 2.86e-14 | ≈ 1.644591e-4 | ≈ 7.799480e-32 | ≈ 1.5407e-34 |

Comparisons used by the chain, all decided exactly:
- m_t² ≤ m_lo² (m_lo − m_t ≈ 6.49 × 10⁻⁸);
- A_q > 0, E_q < 1 and Λ_q > 0 for every q;
- the box-in-ball test of §5.

## 4. The three clusters

The localisation run lists three clusters. They are in the header field `accept` of both `n16_acc2` trees, identical in
id, tile set, points, frame mask and ρ to `out/clusters16_rho4000.json` (see §8 for the comparison).

Each cluster q has:
- a tile set S_q: 16 indices into the sorted list of 25 tiles, as in TREE_FORMAT.md §3;
- 16 integer points P_q[0..15]: a best-known configuration rounded to the grid;
- a frame mask marking 15 points as frame and 1 as non-frame;
- ρ = 4000.

| q | what it is | tile set S_q | non-frame label | local certificate |
|---|---|---|---|---|
| 0 | **Frame A**: one of the two best-known optima | 0 1 3 4 7 8 10 12 13 15 17 19 20 22 23 24 | 7 | `data[0]` |
| 1 | **Frame B**: the other best-known optimum. Swap coordinates a ↔ b in frame A: 11 of its 15 frame points coincide with frame B's, and 4 do not. | 0 1 3 5 7 8 10 12 14 15 19 20 21 22 23 24 | 7 | `data[1]` |
| 2 | **Mirror image of frame B** under a ↔ c, the reflection that fixes the vertex (0, 1, 0) | same as cluster 1 | 10 | `data[2]`: its own rationals, not derived from `data[1]` |

In each cluster all three vertices are occupied and 10 of the 15 frame points lie on the boundary. Each certificate uses
20 point–point contacts and 13 point–side contacts.

The three boxes cover every symmetry image of the two optima, for the following reasons:
- **Frame A.** Only the identity maps A's tile set to itself. So each of the 6 images of A lies in a different tile set of
  the same orbit, and the permutation that maps that tile set back to the root maps the image back onto A itself.
- **Frame B.** Its tile set is mapped to itself by a ↔ c. So B and its mirror image lie in the same root with the same
  tile labels, and each needs its own box.

This explains why three clusters are needed. The proof does not rely on it: the tree certifies directly that every
configuration at level m_t lands in one of the three boxes.

For each q the local certificate `data[q]` consists of:
- the rational frame c̃_q (15 points; coordinates are dyadic with denominators up to 2^110; points on a side lie exactly
  on it);
- the lists of 20 pairs and 13 walls;
- a stress λ_q ∈ ℚ³³, λ_q ≥ 0;
- 60 vectors y_q^{iσ} ∈ ℚ³³, all ≥ 0;
- a rational position for the 16th point.

The frame labels map to certificate indices through `data[q].frame_indices`, which is the frame mask of cluster q in
label order.

## 5. What each stage proves, and the chain

Configuration space: all X = (x_1, …, x_16) ∈ T^16. The value is m(X) = min over i < j of |x_i − x_j|. The value is
invariant under the 6 symmetries of T (the permutations of (a, b, c)) and under relabelling.

**Stage G, localisation run** (`trees/n16_acc2_part0.jsonl.gz`, `trees/n16_acc2_part1.jsonl.gz`). It proves that every X
with m(X) ≥ m_t = 216227204441601/10^15 lies in one of the three cluster boxes, and it rules out every other X at that
level. Precisely: there are a symmetry π, a cluster q and a labelling of π(X) by the tiles of S_q such that every frame
point i of q satisfies |G·π(X)_i,u − P_q[i][u]| ≤ 4000 for each u ∈ {A, B, C}. It covers the whole of T^16, up to the 6
symmetries and relabelling.

- **Tiles.** The 25 closed sub-triangles of side 1/5 cover T. Each has diameter 1/5 = 0.2 < m_t, so no two points of X
  lie in one tile. Choosing, for each point, any tile that contains it gives a set of 16 distinct tiles and a labelling
  in tile order. Points on a shared edge may choose either tile. No boundary-assignment rule is needed, because the
  tiles are closed and their diameter is strictly below m_t.
- **Symmetry reduction.** There are C(25,16) = 2,042,975 tile subsets, forming 341,662 orbits under the 6 symmetries (by
  Burnside's count). The trees contain exactly 341,662 roots, one tile subset per orbit. Each checker verifies that
  every one of the 2,042,975 subsets is mapped by some symmetry onto a proved root.
- **Trees.** Every node step is an exact integer test: a strip reduction or a discard (max squared vertex distance
  < Tc), an empty region, or a split into two closed halves. In total there are 362,164 nodes, 351,864 discards, 10,251
  splits and 49 accepted leaves, with no witness and no undecided leaf. An accepted leaf is checked vertex by vertex:
  every vertex of each frame region lies in the cube [P_q[i] − ρ, P_q[i] + ρ]³.

**Stage L, local certificate** (`out/local_cert16.json`, one entry per cluster). Notation for a fixed q:
- J is the exact 33 × 30 Jacobian of the 20 pair and 13 wall constraints at c̃.
- g_k(c̃; m) = Q(c̃_i − c̃_j) − m² for a pair, and the wall's value for a wall. Walls are exactly 0 at c̃.
- ε^{iσ} = Jᵀ y^{iσ} + σ e_i and r = Jᵀ λ.
- A = 12 · max_{iσ} Σ_pairs y_k^{iσ}, E = max_{iσ} ‖ε^{iσ}‖₁ and B(m) = max_{iσ} Σ_k y_k^{iσ} g_k(c̃; m).
- Λ = Σ_pairs λ_k and Φ = Σ_pairs λ_k Q(c̃_i − c̃_j).

- **Lemma L (localisation).** Hypotheses: A > 0 and E < 1. Let x = c̃ + d be 15 points in T with |x_i − x_j| ≥ m for
  the 20 listed pairs, and let ‖d‖∞ ≤ R := (1 − E)/(2A). Then ‖d‖∞ ≤ δ(m) := 2·max(B(m), 0)/(1 − E).

  Proof. For each (i, σ) we have σ d_i = ε·d − Σ_k y_k ∇g_k·d. Each pair satisfies
  ∇g_k·d ≥ −g_k(c̃; m) − Q(d_i − d_j), and Q(d_i − d_j) ≤ 12‖d‖∞². Walls are linear, so a wall term is exactly
  ∇g_k·d = g_k(x) − g_k(c̃) ≥ −g_k(c̃): it has no quadratic part and adds nothing to A, and it adds only −g_k(c̃)
  (zero when c̃ lies on that wall) to B(m). Hence t = ‖d‖∞ satisfies
  t ≤ A t² + B(m) + E t. On 0 ≤ t ≤ R we have A t² ≤ (1 − E) t / 2, so t ≤ 2B(m)/(1 − E).

  B(m) decreases in m because y ≥ 0. So for m ≥ m_lo, δ(m) ≤ δ(m_lo) = δ_q.
- **Lemma V (value).** Hypothesis: Λ > 0. For x as above with value m(x) and ‖d‖∞ ≤ δ_q, the sum
  Σ λ_k g_k(x; m(x)) is ≥ 0. Expanding it gives

      m(x)² ≤ (Φ + ‖r‖₁ δ_q + 12 Λ δ_q²) / Λ = m_hi,q².

- **Lower bound.** Take c̃_q together with the stored rational 16th point. All 16 points lie in T, and the minimum of Q
  over the 120 pairs is exactly m_lo². Hence m16 ≥ m_lo.
- **Link to stage G.** Two tests are decided exactly for each q:
  - **Box in ball.** For every frame label i,
    max(|P_q[i]_A/G − ũ_i|, |P_q[i]_B/G − ṽ_i|) + ρ/G ≤ R_q. The left side is at most 9.5378 × 10⁻⁵ for every q,
    against R_q ≥ 1.0691 × 10⁻⁴.
  - **Level.** m_t² ≤ m_lo².

**The chain.**

1. m16 ≥ m_lo, from the explicit configuration above.
2. Let X be any configuration with m(X) ≥ m_lo, for example an optimal one. Since m_lo ≥ m_t, stage G gives π, q and a
   labelling such that each frame point of π(X) lies in its box.
3. That box lies in the R_q ball around c̃_q in (u, v), so d := (frame of π(X)) − c̃_q has ‖d‖∞ ≤ R_q. Leaf ⊂ box is
   checked on all 49 accepted leaves by both tree checkers. Box ⊂ ball is checked once per cluster by both local
   checkers. Together they place every accepted leaf inside its ball.
4. The frame of π(X) lies in T and meets the 20 pair constraints with m = m(X) ≥ m_lo, because they are a subset of all
   120 pairs. Lemma L gives ‖d‖∞ ≤ δ_q.
5. Lemma V gives m(X)² ≤ m_hi,q² ≤ m_hi². Applied to an optimal X, this gives m16 ≤ m_hi.

The 16th point enters only in step 1.

**Stage R, refutation run** (`trees/n16_pos_part0.jsonl.gz`, `trees/n16_pos_part1.jsonl.gz`). This stage is not used in
the chain; it is an independent check. No 16 points in T have all distances ≥ m_R = 216248892036713/10^15. Every leaf is
a discard, over 341,662 roots and 400,088 nodes, again covering all 2,042,975 tile subsets modulo the 6 symmetries. So
m16 < 0.216248892036713. This coarser bound uses no cluster and no local certificate, and it is consistent with
m_hi < m_R.

## 6. The 16th point (the rattler)

In each cluster one label is non-frame: label 7 in clusters 0 and 1, and label 10 in cluster 2.

- **Stage G** places no constraint on that point (TREE_FORMAT.md §4).
- **Lemmas L and V** involve only the 15 frame points, the 20 pair constraints among them and the 13 side constraints.
  Leaving out the pairs that contain the 16th point can only raise a minimum over pairs. So the upper bound holds for
  every position of the 16th point.
- **Consequence.** For an optimal configuration, the certificate says only what the definition already says about the
  16th point: it lies in T at distance ≥ m16 from each of the other 15 points. It does not locate that point, and it
  makes no claim about the set of positions the point can take.
- **Lower bound.** The certificate places the 16th point at an explicit rational position (`data[q].rattlers_uv`, as
  (u, v)):
  - q = 0: (141904580249/2^39, 140741430287/2^39)
  - q = 1: (140741430287/2^39, 141904580249/2^39)
  - q = 2: (33388725419/2^36, 141904580249/2^39)

  At that position its distance to the nearest frame point is ≈ 0.22359, about 0.0074 above m_lo. It is not in contact
  with any frame point.

## 7. Checks

- **Checker A (the engine's own replay).**
  - `replay_tri.py` uses the standard library only and shares no geometry code with the search program `tri_engine.py`.
    It passes both runs: stage R with 400,088 nodes and stage G with 362,164 nodes and 49 accepted leaves. In both runs
    all 2,042,975 subsets are covered.
  - `check_local.py` recomputes every stored rational of the local certificate and decides the link to the boxes of
    `out/clusters16_rho4000.json` and to m_t. It reports PASS for all three clusters.
- **Checker B (written independently from the specification alone).**
  - `verify_tri_tree_blind.py` was written from TREE_FORMAT.md, and `check_local_blind.py` from LOCAL.md §2, without
    reading the engine code.
  - The tree checker passes both runs. Its node, op, discard, split, accept and root counts equal checker A's. Equal
    counts corroborate the result but are not a step of the proof.
  - The local checker passes all three clusters, including the box-in-ball test on every frame point and m_t² ≤ m_lo².
  - Three tampered copies are rejected: an inflated box, an altered multiplier and a shifted frame point.
- **Adjustments, all disclosed.**
  - (a) The tree checker first treated the informational end-of-file status line as a hard requirement. The
    specification says not to trust that line, and the requirement was relaxed to match the specification.
  - (b) LOCAL.md does not fix a JSON layout. `adapt_local_cert.py` therefore renames and re-nests fields for checker B
    without changing any value. Its output is reproducible byte for byte from `out/local_cert16.json` (SHA-256 in §8).
    Checker A reads `out/local_cert16.json` directly.
  - (c) As first built, `out/clusters16_rho4000.json` carried a stale field `m_t` = 216035227139987/10^15 (the
    midpoint of the best-known and runner-up values, used when the file was made), not the binding localisation target.
    Checker B reads m_t from every source and requires them to agree, so its first recorded local PASS used that weaker
    value. The published file now carries `m_t` = 216227204441601/10^15; that single value is the only change (same
    2,797 bytes; SHA-256 in §8). On these bytes checker B passes with no extra argument and with the binding argument,
    and rejects any other value ("m_t values disagree"). Checker A always used the binding value. The old bytes are kept
    as `out/clusters16_rho4000.stale.json` (SHA-256 bfb292e2…ceb2) for the record only; they are not an input.
- **Link between the files.** Neither checker compares the tree headers with the cluster file. That comparison is the
  separate one-line step in §8: the `accept` list in both `n16_acc2` headers equals the clusters of
  `out/clusters16_rho4000.json` in id, tile set, points, frame mask and ρ, and both headers' `m_t` equals the file's `m_t`.
- **Lemma hypotheses.** The conditions A > 0, E < 1, Λ > 0, δ evaluated at m_lo, and m_t ≤ m_lo are each hard checks in
  both local checkers.

## 8. Reproduce and verify

The checkers need only Python 3 and its standard library. Regenerating the certificates (`tri_engine.py`,
`local_cert.py`) is not needed to verify them.

| file | bytes | SHA-256 |
|---|---|---|
| `trees/n16_pos_part0.jsonl.gz` (stage R, part 0/2) | 12,198,941 | 3f1268610089e0a2c9ce91f9761d97df039d0215bc68667392693aa4737793f5 |
| `trees/n16_pos_part1.jsonl.gz` (stage R, part 1/2) | 114,915,793 | 9e5f263e0415c9aafee5227d27cb59c7f8209593f00d730cea31950e19f8a82a |
| `trees/n16_acc2_part0.jsonl.gz` (stage G, part 0/2) | 12,658,913 | f73dad5a9615f5467bf7087fb9f2052eefa16d50d8c005a9f47fba37051db542 |
| `trees/n16_acc2_part1.jsonl.gz` (stage G, part 1/2) | 47,407,927 | ecf47e2bb9cad0687a1ec5faa1fd1ca0bd265e385f19d9be2e5903496483d178 |
| `out/local_cert16.json` (stage L) | 165,024 | a6486475dcfe59f165d786ca75f1de4d8f97c592c2eb977ff61dec609cdcf147 |
| `out/clusters16_rho4000.json` (the three boxes) | 2,797 | 1a093a5434609b754510a339add85f388b5979c2c2cd31590f3476a0d436bef2 |
| `out/local_cert16_for_blind.json` (derived: adapter output) | 179,251 | 5b5f70fbfed42bde0b9c5ea91cda17ba30ec62803ad8872e7e42643610fefed5 |
| `replay_tri.py` (checker A, trees) | 12,479 | 464b752f328361a60845e936c8a9150b2cb653fe23070d4e7cd5fcd3cbcf4625 |
| `check_local.py` (checker A, local) | 5,266 | e38689f84bfed57be8b0fa2f70b29d32755d2461c90ee6a82afd053f860b41bb |
| `verify_tri_tree_blind.py` (checker B, trees) | 38,058 | d332cea1fc134fcd4ec0b31591a3d1f40f9cef305432ff391f2d98d341777b81 |
| `check_local_blind.py` (checker B, local) | 19,117 | 38ec4aab50a11532a2995f75c5e4d05789777eecc4f530c00ef18d2e9763243e |
| `adapt_local_cert.py` (rename-only adapter) | 2,357 | d0dfecc1c259c276bc5e6f3d7bbc903ab1af5240f761516face26bc4bafacc2a |

Specifications: `TREE_FORMAT.md` (global stage) and `LOCAL.md` (local stage). Generators, not needed for verification:
`tri_engine.py` and `local_cert.py` (they import only a process-priority helper from `../../solver/chx/`).
GitHub rejects files over 100 MB, so `trees/n16_pos_part1.jsonl.gz` is stored as two pieces (`.001`, `.002`); step 0a below
rebuilds it byte for byte, and the SHA-256 in the table is that of the rebuilt file. The pieces:

| piece | bytes | SHA-256 |
|---|---|---|
| `trees/n16_pos_part1.jsonl.gz.001` | 57,457,897 | c5b210726d5daaf878f8e88f422a28d3abbf700e2541d22afd311eabf9241e98 |
| `trees/n16_pos_part1.jsonl.gz.002` | 57,457,896 | 32cf4b71084c994ef4297e492631e0c70054f5b3e641591b6da8dd055d448924 |

The two specifications are published as used, except that absolute paths were replaced by relative ones and one informal
phrase in LOCAL.md was replaced by the precise statement of §6.

```
# 0a. rebuild the one split file (byte-for-byte concatenation of its two pieces)
python -c "import shutil;o=open('trees/n16_pos_part1.jsonl.gz','wb');[shutil.copyfileobj(open('trees/n16_pos_part1.jsonl.gz.%03d'%k,'rb'),o) for k in (1,2)];o.close()"

# 0. hashes (compare with the table)
python -c "import hashlib,sys;[print(hashlib.sha256(open(p,'rb').read()).hexdigest(),p) for p in sys.argv[1:]]" trees/n16_pos_part0.jsonl.gz trees/n16_pos_part1.jsonl.gz trees/n16_acc2_part0.jsonl.gz trees/n16_acc2_part1.jsonl.gz out/local_cert16.json out/clusters16_rho4000.json

# 1. checker A (about 11 and 15 minutes for the two tree runs)
python replay_tri.py --join trees/n16_pos_part0.jsonl.gz trees/n16_pos_part1.jsonl.gz
python replay_tri.py --join trees/n16_acc2_part0.jsonl.gz trees/n16_acc2_part1.jsonl.gz
python check_local.py out/local_cert16.json out/clusters16_rho4000.json 216227204441601/1000000000000000

# 2. checker B
python verify_tri_tree_blind.py trees/n16_pos_part0.jsonl.gz trees/n16_pos_part1.jsonl.gz --expect proved
python verify_tri_tree_blind.py trees/n16_acc2_part0.jsonl.gz trees/n16_acc2_part1.jsonl.gz --expect cluster
python adapt_local_cert.py out/local_cert16.json cert_B.json out/clusters16_rho4000.json boxes_B.json
python check_local_blind.py cert_B.json boxes_B.json 216227204441601/1000000000000000

# 3. the tree headers carry exactly the boxes and the m_t the local stage checked (prints True)
python -c "import gzip,json;K=('id','combo','points','frame','rho');f=lambda L:[{k:c[k] for k in K} for c in L];C=json.load(open('out/clusters16_rho4000.json'));H=[json.loads(gzip.open(p,'rt').readline()) for p in ('trees/n16_acc2_part0.jsonl.gz','trees/n16_acc2_part1.jsonl.gz')];print(all(f(h['accept'])==f(C['clusters']) and h['header']['m_t']==C['m_t'] for h in H))"
```

Expected results:

| run | checker | result | nodes | ops | discards | splits | accepts | roots | coverage |
|---|---|---|---|---|---|---|---|---|---|
| stage R | A | PASS | 400,088 | 44,096,735 | 370,875 | 29,213 | 0 | 341,662 | 2,042,975 / 2,042,975 |
| stage R | B | PASS | same counts as checker A | | | | | | |
| stage G | A | PASS | 362,164 | 28,363,056 | 351,864 | 10,251 | 49 | 341,662 | 2,042,975 / 2,042,975 |
| stage G | B | PASS | same counts as checker A | | | | | | |

`check_local.py` prints, for each of the three clusters, recomputed == stored: True; box offset + ρ ≈ 9.5377e-5 ≤ R_q;
m_t ≤ m_lo. It ends with PASS. `check_local_blind.py` prints `VERDICT: PASS`. Step 3 prints `True`.

## 9. Also obtained (not part of this claim)

N = 17: m17 < 13209124868233/62500000000000 ≈ 0.2113459979, a refutation run passed by checker A only; its files are not
in the table above. The best known value is (3 − √3)/6 = 0.211324865405187. N = 18 and 19 are not decided. For N = 18
there is a near-tie 4 × 10⁻⁷ below the best known.

## 10. References

- R. L. Graham, B. D. Lubachevsky, Dense packings of equal disks in an equilateral triangle: from 22 to 34 and beyond,
  Electron. J. Combin. 2 (1995) #A1 (arXiv math/0406252).
- J. B. M. Melissen, P. C. Schuur, Packing 16, 17 or 18 circles in an equilateral triangle, Discrete Math. 145 (1995)
  333–342.
- A. Joós, Packing 13 circles in an equilateral triangle, Aequationes Math. 95 (2021) 35–65,
  doi 10.1007/s00010-020-00753-y.
- J. Oler, A finite packing problem, Canad. Math. Bull. 4 (1961) 153–155.
- C. Payan, Empilement de cercles égaux dans un triangle équilatéral. À propos d'une conjecture d'Erdős–Oler, Discrete Math.
  165/166 (1997) 555–565.
- M. C. Markót, T. Csendes, A new verified optimization technique for the "packing circles in a unit square" problems,
  SIAM J. Optim. 16 (2005) 193–219.
- M. C. Markót, Improved interval methods for solving circle packing problems in the unit square, J. Global Optim. 81
  (2021) 773–803.
- E. Friedman, Packing Center, Circles in Triangles, https://erich-friedman.github.io/packing/cirintri/ (read
  2026-09-26; lists N = 1..15).
- Wikipedia, "Circle packing in an equilateral triangle" (read 2026-09-20 and 2026-09-26).

The global stage follows the interval branch-and-bound approach of Markót & Csendes (2005) and Markót (2021), adapted to
the triangle (tiles, barycentric integer grid, exact tests). The local stage replaces their interval Newton step with an
exact-rational certificate (Lemmas L and V).
