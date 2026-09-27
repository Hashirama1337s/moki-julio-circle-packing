# TREE_FORMAT — specification of the branch-and-bound certificate files (Moki&Julio, 2026-09-26)

Purpose: someone who has NOT seen our code must be able to write, from this file alone, an independent checker of the
global stage for N points in an equilateral triangle (target: N = 16). Everything a checker needs is here: the frame,
the exact claim, why each reduction is sound, every test as an exact integer formula, the file format field by field,
the coverage argument, and small worked examples with their files. Use exact integer / rational arithmetic only.

Files referred to (all relative to this folder, `proofs/triangle16/`):

| file | what | bytes | SHA-256 |
|---|---|---|---|
| `example/ex4_proved.jsonl.gz` | N = 4 refutation, 1 node | 379 | e6381915d46b6905984caa006b44c0edd29b354301198405ba58964edbdd44ea |
| `example/ex4_cluster.jsonl.gz` | N = 4 cluster (acceptance) run, 1 node | 585 | 9cfe7890f6975161f95a2958d728afa5a52f7ed3deb3740f5cb4b9bab501840e |
| `example/ex9_split.jsonl.gz` | N = 9 refutation, 3 nodes, 1 split | 1,790 | f84aabf88ddc4cd8a1603c0aad8fe18d5aa3b90dd65f3d4ab5bcfdd6a4ad8447 |
| `trees/n16_pos_part0.jsonl.gz` | N = 16 refutation, part 0 of 2 | 12,198,941 | 3f1268610089e0a2c9ce91f9761d97df039d0215bc68667392693aa4737793f5 |
| `trees/n16_pos_part1.jsonl.gz` | N = 16 refutation, part 1 of 2 | 114,915,793 | 9e5f263e0415c9aafee5227d27cb59c7f8209593f00d730cea31950e19f8a82a |
| `trees/n16_acc_part0.jsonl.gz` | N = 16 cluster run, part 0 of 2 | 12,719,807 | 0988361e81b5a310bb19c04859b575e669e78aa45eaa1838df85e2146522cd87 |
| `trees/n16_acc_part1.jsonl.gz` | N = 16 cluster run, part 1 of 2 | 186,318,505 | 1f3a79ec7e3dc7ca5b5c337ae793cccabcf1d4f34726b41baab9fd5b43eb2eb8 |
| `trees/n16_acc2_part0.jsonl.gz` | N = 16 cluster run with the SMALL boxes (feeds the local stage), part 0 | 12,658,913 | f73dad5a9615f5467bf7087fb9f2052eefa16d50d8c005a9f47fba37051db542 |
| `trees/n16_acc2_part1.jsonl.gz` | same, part 1 | 47,407,927 | ecf47e2bb9cad0687a1ec5faa1fd1ca0bd265e385f19d9be2e5903496483d178 |

Our own checker is `replay_tri.py`. A blind checker should NOT read it before it is finished.

---------------------------------------------------------------------------------------------------------------------

## 1. Frame and coordinates

- The container is the closed equilateral triangle of side 1. A point is given by barycentric coordinates
  (a, b, c), a + b + c = 1, a, b, c ≥ 0. (Concretely: a point is p = a·V1 + b·V2 + c·V3 for the three vertices; which
  vertex is which never matters, because only the formula below is used.)
- **Squared distance.** For two points with coordinate differences (da, db, dc) (so da + db + dc = 0):

      |p − p'|² = (da² + db² + dc²) / 2.

  (Check: the vertices (1,0,0), (0,1,0) give (1 + 1 + 0)/2 = 1; the height from (0,0,1) to (½,½,0) gives
  (¼ + ¼ + 1)/2 = ¾ = (√3/2)².)
- **Integer grid.** A positive integer G is fixed per file (header field `G`). Integer coordinates are
  A = G·a, B = G·b, C = G·c, so A + B + C = G. All bounds in the files are integers in this scale. Points themselves
  are real; only region bounds are integers.
- **Target.** A rational m_t > 0 (header field `m_t`, a string "p/q"). A configuration is **m_t-feasible** if it is N
  points in the closed triangle with all pairwise distances ≥ m_t. For integer-scaled differences
  (dA, dB, dC) = G·(da, db, dc):

      distance ≥ m_t   ⇔   dA² + dB² + dC² ≥ 2·m_t²·G².

  Define **Tc = ⌈2·m_t²·G²⌉** (an integer; header field `Tc` must equal it — recompute it). For an integer S,
  S < 2·m_t²·G² ⇔ S < Tc. All "too close" tests below are of the form "some integer S < Tc".

## 2. Regions ("tri-boxes")

A region is 6 integers [A0, A1, B0, B1, C0, C1] and denotes the closed real set

    R = { (A, B, C) real : A0 ≤ A ≤ A1, B0 ≤ B ≤ B1, C0 ≤ C ≤ C1, A + B + C = G }.

It is a convex polygon (a hexagon or fewer sides, possibly a segment, a point, or empty). Bounds need not be "tight"
(the set is what it is; only the set matters). Its **vertices** are all integer points and are found by enumerating
the ≤ 12 candidates

    (A, B, G−A−B)   for A ∈ {A0, A1}, B ∈ {B0, B1}   kept if C0 ≤ G−A−B ≤ C1
    (A, G−A−C, C)   for A ∈ {A0, A1}, C ∈ {C0, C1}   kept if B0 ≤ G−A−C ≤ B1
    (G−B−C, B, C)   for B ∈ {B0, B1}, C ∈ {C0, C1}   kept if A0 ≤ G−B−C ≤ A1

R is empty iff no candidate is kept (a nonempty bounded convex polygon in the plane A+B+C = G has a vertex, and each
vertex lies on two of the six lines, which is one of the candidates).

**Maximum squared distance between two regions.** For regions R, R' (both nonempty):

    maxS(R, R') = max over vertices v of R, v' of R' of  (vA−v'A)² + (vB−v'B)² + (vC−v'C)².

This is the exact maximum of the convex function over R × R' (a convex function on a product of polytopes attains its
maximum at a pair of vertices). "R is too close to R' everywhere" means maxS(R, R') < Tc: then every point of R is at
distance < m_t from every point of R'.

**Box containment.** "Every point of R lies in the cube box [c − ρ, c + ρ]³ (per coordinate)" is checked on the
vertices of R: |vA − cA| ≤ ρ, |vB − cB| ≤ ρ, |vC − cC| ≤ ρ for every vertex v (convexity).

## 3. Tiles, tile sets, symmetry

- Header field `k` (integer ≥ 1). Condition to check: **k²·m_t² > 1**, i.e. 1/k < m_t. G must be divisible by k;
  put g = G/k.
- **Tiles**: all integer triples (i, j, l) ≥ 0 with i + j + l ∈ {k − 1, k − 2}, sorted lexicographically. There are
  k² of them (k(k+1)/2 upright with sum k − 1, k(k−1)/2 inverted with sum k − 2). The **tile index** is the position in
  this sorted list; the file's top-level `tiles` list must equal it. Tile (i, j, l) is the region

      [i·g, (i+1)·g, j·g, (j+1)·g, l·g, (l+1)·g]

  (upright: vertices (i+1,j,l)g, (i,j+1,l)g, (i,j,l+1)g; inverted: (i+1,j+1,l)g, (i+1,j,l+1)g, (i,j+1,l+1)g).
- **Tile lemma.** (a) The tiles cover the triangle: for a point with k·a = i + α, k·b = j + β, k·c = l + γ
  (integer parts i, j, l, fractional parts α, β, γ ∈ [0,1)), α + β + γ = k − (i+j+l) is an integer in {0, 1, 2}. If it
  is 1 or 2 the point lies in tile (i, j, l) (sum k − 1 or k − 2); if it is 0 the point is a grid vertex
  (i, j, l)/k with i + j + l = k; lower any positive index by 1 to get an upright tile that contains it. (b) Each tile
  has diameter 1/k (its vertices are pairwise at distance 1/k), so two points in one closed tile are at distance
  ≤ 1/k < m_t. Hence in an m_t-feasible configuration every point lies in at least one tile and no tile contains two
  points: choosing one containing tile per point gives an injective map, i.e. a set S of N distinct tiles.
- **Tile set and labelling.** A **tile set** (a "combo") is a sorted list of N distinct tile indices
  S = [s_0 < s_1 < … < s_{N−1}]. Its **root regions** are R_i = tile s_i. Point i of a configuration "belongs to" S
  if point i lies in tile s_i. By the tile lemma every m_t-feasible configuration can be labelled so that point i lies in
  tile s_i for some tile set S (relabel the points in the order of their chosen tiles).
- **Symmetry.** The six permutations π of the three barycentric coordinates are the symmetries of the triangle; they
  preserve the triangle and the distance formula, so π maps m_t-feasible configurations to m_t-feasible
  configurations. π maps tile (i, j, l) to the tile with the permuted triple (the sum is unchanged), hence tile sets to
  tile sets. A checker must verify: **for every N-subset S of the k² tiles there is a permutation π such that the sorted
  image π(S) is a PROVED root of the certificate** (definition in §6). Roots that are not lexicographically minimal
  in their orbit are allowed (they are just extra work); roots may appear in any order and in any part file.
- **Pigeonhole case.** If k² < N there is no injective map, so no m_t-feasible configuration exists; the file then
  contains a record `{"pigeonhole": [k², N]}` and no roots.

## 4. What a certificate claims

**Refutation run** (files `n16_pos_*`, examples `ex4_proved`, `ex9_split`): *no m_t-feasible configuration exists*,
i.e. every N points in the closed unit triangle have two points at distance < m_t. Consequently the maximin distance
m_N is < m_t.

**Cluster run** (files `n16_acc_*`, example `ex4_cluster`): the header lists clusters q = (id, combo S_q,
points P_q[0..N−1] (integer triples, each summing to G), frame F_q[0..N−1] (booleans), radius ρ_q (integer)).
Claim: *for every m_t-feasible configuration X there are a permutation π and a cluster q such that, labelling the
points of π(X) by the tiles of S_q (point i in tile S_q[i]), every frame point is close to the cluster point:*

    for every i with F_q[i] true and every coordinate u ∈ {A, B, C}:   | G·π(X)_i,u − P_q[i][u] | ≤ ρ_q.

Non-frame points (F_q[i] false; for N = 16 the rattler) are unconstrained by the claim. No uniqueness or optimality is
claimed by this stage; it only says that everything at least m_t-good sits in the listed boxes.

## 5. Node semantics and the exact tests

A **node** holds one region per point, R_0 … R_{N−1}, and stands for the set of labelled configurations with point i
in R_i for every i. Invariant to be maintained by every step (this is what makes the certificate sound):

> (I) every m_t-feasible labelled configuration that lies in the node's regions *before* a step still lies in them
> *after* it.

Steps inside a node, applied in the order they appear in the file:

**(a) Reduction op** `[i, side, t, j]` (integers; 0 ≤ i, j < N, **i ≠ j**, side ∈ {0,…,5}). Coordinate
q = side div 2 (0 = A, 1 = B, 2 = C); side even = low side, odd = high side. The **strip** is

    low  side:  S = R_i with its q-upper bound replaced by t      (the points of R_i with coordinate q ≤ t)
    high side:  S = R_i with its q-lower bound replaced by t      (the points of R_i with coordinate q ≥ t)

Check: current q-lower bound of R_i ≤ t ≤ current q-upper bound of R_i; and either S is empty, or
**maxS(S, R_j) < Tc**. Effect: set R_i's q-lower bound (low side) or q-upper bound (high side) to t.
Soundness: if point i were in S, then (since point j is in R_j) its distance to point j would be < m_t, so no
m_t-feasible configuration has point i in S; removing S keeps (I) (the new region is the closed set R_i ∩ {coord ≥ t}
resp. ≤ t, which contains R_i minus S).

**(b) Discard** fate `["D", i, j]` (i ≠ j): check **maxS(R_i, R_j) < Tc** on the regions after all ops of the node
(or that some region is empty). Soundness: points i and j would be at distance < m_t.

**(c) Empty** fate `["E", i]`: check that region R_i (after the node's ops) is empty. No configuration lies in the node.

**(d) Split** fate `["S", a, q, c, lo_id, hi_id]`: 0 ≤ a < N, q ∈ {0,1,2}, the q-lower bound of R_a ≤ c ≤ the q-upper
bound of R_a (after the node's ops). Children: node `lo_id` = all regions copied, R_a's q-upper bound set to c;
node `hi_id` = all regions copied, R_a's q-lower bound set to c. The two closed halves cover R_a, so (I) passes to the
union of the children. Both child ids must appear later in the same file as nodes, exactly once each.

**(e) Accept** fate `["A", q]` (cluster runs only): cluster q exists in the header, **its combo equals this node's
root combo**, and for every i with F_q[i] true every vertex of R_i (after the node's ops) lies in the box
[P_q[i] − ρ_q, P_q[i] + ρ_q] (all three coordinates). Then every configuration in the node satisfies the cluster
claim of §4 for this q (with the π that led to this root).

**(f) Witness** fate `["F", points]`: N integer triples. Not part of a proof: it shows the target is NOT refutable.
Check: every triple has nonnegative entries summing to G, and every pair has dA² + dB² + dC² ≥ Tc. A root containing an
F leaf is not proved.

**(g) Undecided** fate `["U"]`: the search gave up here. The root is not proved.

## 6. File format (gzip-compressed JSON Lines, UTF-8, one JSON object per line)

All numbers are JSON integers (arbitrary size; use big integers) except inside the header's informational fields.

1. **Line 1 — header object** with keys
   - `header`: object. Fields a checker MUST use: `n` (N), `m_t` (string "p/q"), `k`, `G`, `Tc`. Informational
     (ignore): `m_star`, `eps`, `tiles` (count), `rule`, `rounds`, `cap`, `part` ("i/p"), `priority`.
   - `tiles`: list of k² triples [i, j, l] (must equal the sorted list of §3).
   - `accept` (cluster runs only): list of cluster objects `{"id": int, "combo": [N sorted tile indices],
     "points": [N triples], "frame": [N booleans], "rho": int, …}`; extra keys (`value_float`,
     `rattlers_by_cage_slack`, `slack_max`) are informational.
2. Then either one `{"pigeonhole": [k², N]}` record, or a sequence of blocks, one per root:
   - **root record** `{"id": r, "root": S}` — S a sorted list of N distinct tile indices. The root node's regions are
     the tiles S[0], …, S[N−1] (region i = tile S[i]).
   - **node records** `{"id": x, "ops": [[i, side, t, j], …], "fate": [...]}`. The first node record of a block has
     x = r (the root node itself). Every other node id is a child id announced by an earlier split in the same block.
     Records appear in depth-first order: a node's record comes after its parent's record. Node ids are unique within a
     file; they restart from 0 in each part file.
3. Last line `{"end": status, "nodes": count}` — informational; do not trust it.

A **root is PROVED** iff: its root node and every announced child appear exactly once, every op and fate checks, and
no fate is F or U. (In a cluster run, A fates are allowed; in a refutation run there are no clusters, so an A fate
must be rejected.)

**Multi-part runs.** A run may be split into p part files (header `part` "i/p"). They must agree on n, m_t, k, G, Tc,
the tile list and (cluster runs) the cluster list. The coverage check of §3 is done over the union of their proved
roots.

## 7. Checker algorithm (reference outline)

```
read header(s); recompute Tc = ceil(2 m_t^2 G^2), check k^2 m_t^2 > 1, G % k == 0, tile list
if pigeonhole record: require k^2 < N  → claim holds
proved = set()
for each part file, stream records:
    on root record: push root node (regions = tiles of S), root_status[S] = {ok: true, open: 1}
    on node record: pop its regions (must be pending); open -= 1
        for op in ops: check (a), update the region
        check the fate (b)–(g); on split push two children (open += 2); F or U → ok = false
for each root: PROVED iff ok and open == 0          (all announced children were seen)
for every N-subset S of range(k^2): some permutation of (i,j,l) maps S to a PROVED root  → coverage
report: refutation run → "no m_t-feasible configuration"; cluster run → the cluster claim of §4
```

Memory: stream the files; keep only pending nodes (depth-first order keeps this small). Our checker needs 11 and 15 min
for the two N = 16 pairs in pure Python; the coverage loop is C(25,16) = 2,042,975 subsets × 6 permutations.

## 8. Why the conclusions follow (proof sketch to check against)

Let X be m_t-feasible. By the tile lemma choose a tile set S and labelling. By coverage there is π with π(S) a PROVED
root; π(X), labelled by the tiles of π(S), lies in that root node's regions. By (I) through ops and splits, π(X) lies
in the regions of some leaf of the root's tree (follow the child whose half contains the split coordinate of point a;
if on the boundary either child). That leaf cannot be D or E (they contain no feasible configuration). In a refutation
run every leaf is D or E: contradiction, so X does not exist. In a cluster run the leaf is A with some cluster q whose
combo is π(S), and the box check gives the claim of §4.

## 9. Worked examples (full files)

### 9.1 `example/ex4_proved.jsonl.gz` — N = 4 is refuted at m_t = m4 (1 + 1e-4)

```
{"header":{"n":4,"m_star":"0.57735026918962576451","eps":"1e-4","m_t":"115481600843309/200000000000000","k":2,"G":33554432,"Tc":750750065388662,"tiles":4,"rule":0,"rounds":30,"cap":10000000,"part":"0/1","priority":16384},"tiles":[[0,0,0],[0,0,1],[0,1,0],[1,0,0]]}
{"id":0,"root":[0,1,2,3]}
{"id":0,"ops":[[0,0,4428131,1],[0,2,4428131,1],[0,5,12349085,1],[0,0,9980955,2],[0,3,12349085,2],[0,1,11224392,3],[1,1,1243437,0],[1,3,1201663,0],[1,4,32310995,0],[2,1,1243437,0],[2,2,32310995,0],[2,5,1201663,0],[3,0,32310995,0],[0,2,11147422,1],[0,5,11224392,1],[0,0,11147422,2],[0,3,11182618,2],[0,1,11182618,3],[0,4,11189239,3]],"fate":["D",0,1]}
{"end":"PROVED","nodes":1}
```

Walk-through (you should reproduce every number):
- m_t = 115481600843309/2·10^14 = 0.577406…, k = 2: k²m_t² = 4·0.3334 > 1 ✓. G = 33554432 = 2·2^24, g = 16777216.
  Tc = ⌈2·m_t²·G²⌉ = 750750065388662.
- Tiles (sorted): 0 = (0,0,0) inverted (the central triangle, region [0,g,0,g,0,g]); 1 = (0,0,1) = [0,g,0,g,g,2g];
  2 = (0,1,0) = [0,g,g,2g,0,g]; 3 = (1,0,0) = [g,2g,0,g,0,g]. Only one 4-subset exists: {0,1,2,3}; it is the root.
- First op [0,0,4428131,1]: region 0 (central tile), low side of A, t = 4428131, against region 1 (tile (0,0,1), the
  corner tile at C = G). Strip = R_0 with A-upper bound 4428131 = {A ≤ 4428131} ∩ central tile. Its vertices and the
  vertices of tile 1 give maxS < Tc (every point of the strip is within m_t of every point of the corner tile, whose
  diameter is ½). New R_0 = [4428131, g, 0, g, 0, g].
- The 19 ops push all four regions; after the last op, maxS(R_0, R_1) < Tc → D. The root is PROVED; coverage: the
  single subset is itself → no 4 points have all distances ≥ 0.577406… (the known optimum is 1/√3 = 0.577350…).

### 9.2 `example/ex4_cluster.jsonl.gz` — acceptance at m_t = m4 (1 − 1e-3)

```
{"header":{"n":4,"m_star":"0.57735026918962576451","eps":"-1e-3","m_t":"144193229730109/250000000000000","k":2,"G":33554432,"Tc":749099488619231,"tiles":4,"rule":0,"rounds":30,"cap":10000000,"part":"0/1","priority":16384},"tiles":[[0,0,0],[0,0,1],[0,1,0],[1,0,0]],"accept":[{"id":0,"combo":[0,1,2,3],"points":[[11184810,11184810,11184812],[0,0,33554432],[0,33554432,0],[33554432,0,0]],"frame":[false,true,true,true],"rho":111848}]}
{"id":0,"root":[0,1,2,3]}
{"id":0,"ops":[[0,0,4395894,1],[0,2,4395894,1],[0,5,12381322,1],[0,0,9920729,2],[0,3,12381322,2],[0,1,11252381,3],[1,1,1331652,0],[1,3,1286637,0],[1,4,32222780,0],[2,1,1331652,0],[2,2,32222780,0],[2,5,1286637,0],[3,0,32222780,0],[0,2,11094685,1],[0,5,11252381,1],[0,0,11094685,2],[0,3,11207366,2],[0,1,11207366,3],[0,4,11139851,3],[1,1,67364,0],[1,3,67364,0],[1,4,33487068,0],[2,1,112681,0],[2,2,33441751,0],[2,5,112530,0],[3,0,33441751,0],[3,5,112530,0],[0,0,11139851,1],[0,2,11139851,1],[0,5,11207215,1],[0,0,11140003,2],[0,3,11207215,2],[0,1,11207214,3],[1,1,67211,0],[1,3,67211,0],[1,4,33487221,0],[2,1,67211,0],[2,2,33487221,0],[3,0,33487221,0],[0,2,11140004,1],[0,5,11207214,1],[0,0,11140004,2],[0,3,11207214,2],[1,1,67210,0],[1,3,67210,0],[1,4,33487222,0],[2,1,67210,0],[2,2,33487222,0],[3,0,33487222,0]],"fate":["A",0]}
{"end":"PROVED_MOD_CLUSTERS","nodes":1}
```

The cluster is the known optimum (three vertices + the centroid); point 0 (the centroid, in the inverted tile) is
marked non-frame to exercise the exemption. After the 49 ops, regions 1, 2, 3 lie within ρ = 111848 (= G/300) of the
three vertices → A. Claim proved: every 4 points with all distances ≥ 0.576773… have, after a symmetry, their three
points in the corner tiles within G/300 of the three vertices (the fourth anywhere).

### 9.3 `example/ex9_split.jsonl.gz` — a split (N = 9, m_t = 333366666666667/10^15, k = 3, G = 50331648)

One root (all 9 tiles), node 0: 44 ops, fate `["S", 0, 1, 8391960, 1, 2]` (region 0, coordinate B, split at 8391960);
node 1 (the B ≤ 8391960 half): 97 ops, `["D", 2, 7]`; node 2 (the B ≥ 8391960 half): 121 ops, `["D", 5, 6]`.
Note the header of this older file says `"rule":1` (a search heuristic; irrelevant to checking).

## 10. The N = 16 certificates and the expected results

| run | files | N | m_t | k | G | Tc | expected |
|---|---|---|---|---|---|---|---|
| refutation | `trees/n16_pos_part0/1` | 16 | 216248892036713/1000000000000000 | 5 | 41943040 | 164534731528431 | every one of the C(25,16) = 2,042,975 subsets covered by a PROVED root (341,662 roots, 400,088 nodes) → m16 < 0.216248892036713 |
| cluster | `trees/n16_acc_part0/1` | 16 | 6757034594909/31250000000000 | 5 | 41943040 | 164498539497368 | all subsets covered; 1,788 A leaves in 3 clusters (ρ = 209715 each; one non-frame point per cluster) → the cluster claim of §4 |
| cluster, small boxes (**the one the local stage needs**) | `trees/n16_acc2_part0/1` | 16 | 216227204441601/1000000000000000 | 5 | 41943040 | 164501730816426 | all subsets covered; 49 A leaves in the same 3 clusters with ρ = 4000 each (`out/clusters16_rho4000.json`, SHA-256 bfb292e2f6c21ae1b114c23fc06aadfa5f3554a3a3242d6e4a998eda3464ceb2) → the cluster claim of §4 |

Our checker's output for these (for comparison only after you have your own): refutation — nodes 400,088, ops
44,096,735, discards 370,875, splits 29,213, roots 341,662; cluster — nodes 442,394, ops 61,287,170, discards 390,240,
splits 50,366, accepts 1,788; small-box cluster run — see LOCAL.md L3 for our replay's counts (nodes 362,164,
accepts 49). The cluster list is also in `out/clusters16.json` (SHA-256
9e7fb0f9eb29eb72cc66c0c0427a5f6e25084d49068b0b9a4a2cf1c81d60808c); the headers of both cluster part files contain it.

## 11. Edge cases a checker must handle

- An op whose strip is empty is valid (nothing is removed; the bound may still move — it only intersects R_i with a
  half-plane that already contains it).
- Regions can become degenerate (a segment or a point); vertices then repeat — use a set.
- A D fate may cite a pair (i, j) in either order; i = j must be rejected (and i = j in an op).
- A tile set may appear as a root more than once; it counts as PROVED if any of its copies is PROVED (our runs never
  repeat a root).
- Integers: coordinates ≤ G, dA² + dB² + dC² ≤ 2G² ≈ 3.5·10^15 for these files; any exact integer type works.
- Ops refer to the regions as modified by the earlier ops of the same node (apply them in file order).
- The claim is about closed sets and "≥ m_t"; all tests are strict "< Tc" for too-close, "≥ Tc" for witnesses.
