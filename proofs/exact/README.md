# Exact optimality certificates: replay

**Status:** reviewed by Grok (xAI) in three adversarial rounds (2026-09-25), verdict PUBLISHABLE; published in v1.7. A second replay
checker, `verify_proof_exact.py`, written by Grok without seeing ours, also passes both cases and rejects all 68
corruptions (added in v1.9). Authors: Moki&Julio. Branch-and-bound engine: Grok's engine
(`prove_small.py`, class `Engine`), run through a logging subclass in `exact_prove.py`.

Two entries, both in the centred frame (rectangle [-1/2, 1/2] x [-h/2, h/2]):

| case | N | h | optimal radius r* | files |
|---|---|---|---|---|
| `crc_800_4` | 4 | 4/5 | 7/20 - sqrt2/10 | `crc_800_4_cert.json`, `crc_800_4_tree.json.gz` |
| `crc_600_6` | 6 | 3/5 | 1/5 - sqrt2/30 | `crc_600_6_cert.json`, `crc_600_6_tree.json.gz` |

## Replay (standard library only, a few seconds)

    py -3.11 replay_exact.py                 # both cases
    py -3.11 replay_exact.py crc_600_6       # one case

It prints what it verified and `PASS <case>` or `FAIL <case>: <reason>`. Exit code 0 iff every case passes.
`replay_exact.py` does not import `exact_prove.py` or `prove_small.py`. It uses `fractions.Fraction` and its own exact
arithmetic in Q(sqrt2) (pairs a + b sqrt2 with an exact sign test). It trusts no number in the files.

Second, independent checker (Grok, standard library only, written from this README and PROOFS_EXACT.md without seeing
`replay_exact.py`):

    py -3.11 verify_proof_exact.py crc_800_4_cert.json crc_800_4_tree.json.gz
    py -3.11 verify_proof_exact.py crc_600_6_cert.json crc_600_6_tree.json.gz

One labelled line in it was changed after its first run, at Grok's instruction: the root cell cover may OVERLAP (closed cells;
only a gap breaks covering); its first version demanded an exact partition and so rejected crc_600_6, whose cells overlap by
one grid unit because G = 32768 is not divisible by 3.

## Mutation test (does the checker catch errors?)

    py -3.11 mutate_check.py                 # from any directory; about a minute

This applies 34 corruptions per case (68 in all), one at a time: certificate numbers nudged by as little as 1e-70,
false witnesses, a strip pushed one grid step too far, a wrong permutation, a gap in the cell cover, an undecided leaf,
and others. Each corrupted pair of files is written to a temporary directory; the files here are never modified.
`replay_exact.py` must reject every corruption and PASS the unmodified files. Exit code 0 iff it does. Last run:
68/68 rejected, the unmodified files PASS.

## What the files hold

- `<case>_cert.json`, the Lemma 1 certificate: h, r* = p - q sqrt2 as (p, q), c* exactly in Q(sqrt2), the tight set T
  (contact graph), the stress lambda in Q(sqrt2) with a rational lower bound for each entry, Lambda_hi, the basis B (2N rows
  of T), M and M^-1 in Q(sqrt2), the rational bound normMinv_hi >= ||M^-1||_inf, rho, snap_err_hi, rho', S2_LO and S2_HI,
  r_t, the grid G, and the four snapped mirror images of c* in the engine's point frame.
- `<case>_tree.json.gz`, the branch-and-bound tree. The header gives r_t, G, W, H, d, the chart and the cuts, and it is
  bound to the certificate by `cert_sha256`. Every node has `id`, `parent`, `box` (grid integers
  [X0, X1, Y0, Y1] per point), `ops` (the reduction steps, in order), `tbox` (the box after the reductions), and `fate`.
  - `cells` (root only): the cell cover, the surviving assignments (children), and a witness for every rejected assignment.
  - `split`: `["bisect", i, axis, mid]` or `["hole", i, j, [hx0, hx1, hy0, hy1]]`, with the children ids.
  - `discard`: a witness. `["pair", i, j, fx, fy]` or `["forbid", i, j]` means the range maximum of |p_i - p_j|^2 over the
    two boxes is below d^2. `["empty", i, axis]` means an empty interval. `["xorder", i]` or `["xsum"]` means a cut
    fails on the whole box.
  - `accept`: `mirror` (index into the certificate's images) and `perm`. Point box i lies inside the closed rho'-square
    around image point `perm[i]`.
  - `ops`: `o_hi i`, `o_lo i`, `s1`, `sN`, `y1` are clips by the non-strict cuts x_i <= x_{i+1}, x_0 + x_{N-1} <= 0 and
    y_0 <= 0 (grid: X_0 + X_{N-1} <= G, Y_0 <= G/2). `["sh", i, j, side, v]` moves one side of box i to v. The strip it
    removes has a pair range maximum below d^2 against box j.

## Regenerate (the files above are the output of)

    python exact_prove.py crc_800_4
    python exact_prove.py crc_600_6

Negative controls: add a relative offset such as `1/1000`. They must end `UNDECIDED` and they write nothing. The full
argument is in `../../PROOFS_EXACT.md`.
