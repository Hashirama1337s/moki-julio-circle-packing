#!/usr/bin/env python3
"""Blind checker for equilateral-triangle branch-and-bound certificates.

Written only from TREE_FORMAT (Moki&Julio, 2026-09-26). Python 3.8+ standard
library. No binary floats are formed; JSON non-integers are kept as text and
are rejected in every decision field.

Ambiguities — stricter reading taken (see the note after this file):
1. Op with a nonempty strip and an empty other region: maxS is undefined, so
   the op is rejected (not treated as a vacuous success).
2. Discard "some region is empty": only the cited pair counts. An empty third
   region does not excuse a pair whose maxS is not < Tc.
3. One or more positional files. Coverage is the union. Header `part` is
   checked as a shape only and is not a proof of completeness.
4. End record is required. `nodes` must equal the node-record count. `end`
   must be "PROVED" or "PROVED_MOD_CLUSTERS" matching the run, and is never
   enough to pass.
5. Any key not listed in the spec is rejected. Listed informational keys are
   ignored for the geometry and are still compared across part files.
6. An accept leaf whose framed region has no vertex is rejected.
7. A witness that meets the spec's integer tests falsifies a refutation, even
   if some other root would still cover that combo's orbit.
8. Parent-before-child in the same file is required. lo_id before hi_id is not.
9. If header.tiles (the count) is present it must equal k^2.
10. Cluster rho must be >= 0. Cluster points must be nonnegative and sum to G.
11. A pigeonhole file still has to satisfy k^2*m_t^2 > 1, the tile list, and Tc.
"""

import gzip
import itertools
import json
import sys


class CertError(Exception):
    def __init__(self, reason):
        super().__init__(reason)
        self.reason = reason


class CTX:
    path = ""
    line = 0


class RawFloat:
    """A JSON number that was not an integer. Never coerced to float."""

    __slots__ = ("text",)

    def __init__(self, text):
        self.text = text


def fail(msg):
    text = " ".join(str(msg).split())
    if CTX.path and CTX.line:
        loc = f"{CTX.path}:{CTX.line}: "
    elif CTX.path:
        loc = f"{CTX.path}: "
    else:
        loc = ""
    raise CertError(loc + text)


def loads(line):
    def pairs(ps):
        obj = {}
        for key, val in ps:
            if type(key) is not str:
                fail("non-string JSON key")
            if key in obj:
                fail(f"duplicate key {key!r}")
            obj[key] = val
        return obj

    def floats(token):
        return RawFloat(token)

    def constants(token):
        fail(f"non-finite number {token}")

    try:
        return json.loads(
            line,
            object_pairs_hook=pairs,
            parse_float=floats,
            parse_constant=constants,
        )
    except CertError:
        raise
    except json.JSONDecodeError as exc:
        fail(f"bad JSON ({exc.msg})")


def check_keys(obj, required, optional, what):
    if type(obj) is not dict:
        fail(f"{what} is not an object")
    keys = set(obj)
    missing = required - keys
    if missing:
        fail(f"{what} missing {sorted(missing)}")
    extra = keys - required - optional
    if extra:
        fail(f"{what} unknown fields {sorted(extra)}")


def recompute_Tc(p, q, G):
    """ceil(2 * (p/q)^2 * G^2) by integer arithmetic.

    For an integer S, S < 2*m_t^2*G^2 iff S < Tc, and S >= 2*m_t^2*G^2 iff
    S >= Tc. Every too-close test is strict < Tc; witnesses use >= Tc.
    """
    num = 2 * p * p * G * G
    den = q * q
    return (num + den - 1) // den


def parse_mt(text):
    if type(text) is not str:
        fail("m_t is not a string p/q")
    parts = text.split("/")
    if len(parts) != 2:
        fail("m_t is not p/q")
    num_s, den_s = parts
    if (
        not num_s
        or not den_s
        or any(ch < "0" or ch > "9" for ch in num_s)
        or any(ch < "0" or ch > "9" for ch in den_s)
    ):
        fail("m_t is not a positive p/q")
    p = int(num_s)
    q = int(den_s)
    if p <= 0 or q <= 0:
        fail("m_t is not positive")
    return p, q


def binom(n, k):
    if k < 0 or n < 0 or k > n:
        return 0
    k = min(k, n - k)
    c = 1
    for i in range(k):
        c = c * (n - i) // (i + 1)
    return c


def region_vertices(R, G):
    """Integer vertices of the tri-box, per the 12-candidate enumeration."""
    A0, A1, B0, B1, C0, C1 = R
    found = set()
    add = found.add
    for A in (A0, A1):
        for B in (B0, B1):
            C = G - A - B
            if C0 <= C <= C1:
                add((A, B, C))
    for A in (A0, A1):
        for C in (C0, C1):
            B = G - A - C
            if B0 <= B <= B1:
                add((A, B, C))
    for B in (B0, B1):
        for C in (C0, C1):
            A = G - B - C
            if A0 <= A <= A1:
                add((A, B, C))
    return found


def pairs_lt(vs, ws, Tc):
    """True iff both regions are nonempty and every vertex pair has S < Tc."""
    if not vs or not ws:
        fail("internal: maxS on an empty region")
    for A, B, C in vs:
        for Ap, Bp, Cp in ws:
            dA = A - Ap
            dB = B - Bp
            dC = C - Cp
            if dA * dA + dB * dB + dC * dC >= Tc:
                return False
    return True


def max_pair_s(vs):
    pts = list(vs)
    best = 0
    n = len(pts)
    for i in range(n):
        A, B, C = pts[i]
        for j in range(i, n):
            D, E, F = pts[j]
            dA = A - D
            dB = B - E
            dC = C - F
            s = dA * dA + dB * dB + dC * dC
            if s > best:
                best = s
    return best


def side_indices(side):
    """Return (coordinate, strip bound index, effect bound index).

    Low (even) strip is q <= t, so the strip replaces the upper bound and the
    surviving closed set is q >= t (lower bound moves to t, boundary kept).
    High (odd) is the mirror. Moving to t+1 would drop points that are not in
    the strip and can still be m_t-feasible.
    """
    q = side // 2
    lo = 2 * q
    hi = lo + 1
    if side % 2 == 0:
        return q, hi, lo
    return q, lo, hi


def reduced_region(R, side, t, get_other, G, Tc, where):
    _q, strip_at, effect_at = side_indices(side)
    lo = R[2 * _q]
    hi = R[2 * _q + 1]
    if t < lo or t > hi:
        fail(where + "t outside current bounds")
    strip = R[:]
    strip[strip_at] = t
    sv = region_vertices(strip, G)
    if sv:
        other = get_other()
        if not other:
            # Ambiguity 1: do not treat maxS over an empty region as success.
            fail(where + "strip nonempty and other region empty")
        if not pairs_lt(sv, other, Tc):
            fail(where + "strip is not strictly within m_t")
    out = R[:]
    out[effect_at] = t
    if out[0] > out[1] or out[2] > out[3] or out[4] > out[5]:
        fail(where + "inverted bounds")
    return out


def tile_region(i, j, l, g):
    return [i * g, (i + 1) * g, j * g, (j + 1) * g, l * g, (l + 1) * g]


def predicted_vertices(i, j, l, k, g):
    s = i + j + l
    if s == k - 1:
        return (
            ((i + 1) * g, j * g, l * g),
            (i * g, (j + 1) * g, l * g),
            (i * g, j * g, (l + 1) * g),
        )
    if s == k - 2:
        return (
            ((i + 1) * g, (j + 1) * g, l * g),
            ((i + 1) * g, j * g, (l + 1) * g),
            (i * g, (j + 1) * g, (l + 1) * g),
        )
    fail("internal: tile sum")


def expected_tiles(k):
    exp = []
    for i in range(k):
        for j in range(k):
            for l in range(k):
                s = i + j + l
                if s == k - 1 or s == k - 2:
                    exp.append((i, j, l))
    exp.sort()
    upright = k * (k + 1) // 2
    inverted = k * (k - 1) // 2
    if len(exp) != upright + inverted or upright + inverted != k * k:
        fail("internal: tile count != k^2")
    return exp


def check_grid_vertex_cover(tiles, index, k, g, G):
    """Frac-part 0 of the tile lemma: each k-grid vertex lies in a lowered upright tile."""
    for i in range(k + 1):
        for j in range(k + 1 - i):
            l = k - i - j
            A = i * g
            B = j * g
            C = l * g
            if A + B + C != G:
                fail("internal: grid vertex off the plane")
            lowered = []
            if i:
                lowered.append((i - 1, j, l))
            if j:
                lowered.append((i, j - 1, l))
            if l:
                lowered.append((i, j, l - 1))
            if not lowered:
                fail("grid vertex has no positive index to lower")
            for tri in lowered:
                if tri not in index:
                    fail("lowered triple is not a tile")
                if tri[0] + tri[1] + tri[2] != k - 1:
                    fail("lowered tile is not upright")
                R = tile_region(tri[0], tri[1], tri[2], g)
                if not (
                    R[0] <= A <= R[1]
                    and R[2] <= B <= R[3]
                    and R[4] <= C <= R[5]
                ):
                    fail("grid vertex lies outside its lowered tile")


def validate_tiles(raw, k, G, Tc):
    if type(raw) is not list:
        fail("tiles is not a list")
    exp = expected_tiles(k)
    parsed = []
    for tri in raw:
        if type(tri) is not list or len(tri) != 3:
            fail("tile is not a triple")
        if type(tri[0]) is not int or type(tri[1]) is not int or type(tri[2]) is not int:
            fail("tile coordinate is not an integer")
        parsed.append((tri[0], tri[1], tri[2]))
    if parsed != exp:
        fail("tile list does not match the specification")
    if G % k != 0:
        fail("G is not divisible by k")
    g = G // k
    if g <= 0:
        fail("tile scale g is not positive")
    cap = 2 * g * g
    # Diameter 1/k means max pairwise S == 2*g^2, and k^2*m_t^2 > 1 means that S < Tc.
    if cap >= Tc:
        fail("tile diameter is not strictly below m_t")
    index = {tri: i for i, tri in enumerate(exp)}
    for i, j, l in exp:
        R = tile_region(i, j, l, g)
        vs = region_vertices(R, G)
        pred = predicted_vertices(i, j, l, k, g)
        if vs != set(pred):
            fail(f"tile {(i, j, l)} vertices do not match the specification")
        if max_pair_s(vs) != cap:
            fail(f"tile {(i, j, l)} diameter is not 1/k")
    check_grid_vertex_cover(exp, index, k, g, G)
    return exp, g


# All 6 coordinate permutations. Composition used by the self-test is
# r[i] = p[q[i]], matching "permute by p, then by q" on a triple.
PERMS = (
    (0, 1, 2),
    (0, 2, 1),
    (1, 0, 2),
    (1, 2, 0),
    (2, 0, 1),
    (2, 1, 0),
)


def check_perm_group():
    if len(set(PERMS)) != 6:
        fail("internal: permutation list")
    have = set(PERMS)
    for p in PERMS:
        if sorted(p) != [0, 1, 2]:
            fail("internal: permutation list")
    for p in PERMS:
        for q in PERMS:
            r = (p[q[0]], p[q[1]], p[q[2]])
            if r not in have:
                fail("internal: permutations are not the full symmetry group")


def symmetry_maps(tiles):
    check_perm_group()
    index = {tri: i for i, tri in enumerate(tiles)}
    maps = []
    n = len(tiles)
    for p in PERMS:
        mp = [0] * n
        for i, tri in enumerate(tiles):
            image = (tri[p[0]], tri[p[1]], tri[p[2]])
            j = index.get(image)
            if j is None:
                fail("symmetry sends a tile outside the tile set")
            mp[i] = j
        if sorted(mp) != list(range(n)):
            fail("symmetry action is not a bijection on tiles")
        maps.append(mp)
    return maps


def apply_map(mask, mp, bit):
    out = 0
    while mask:
        lowest = mask & -mask
        b = lowest.bit_length() - 1
        out |= bit[mp[b]]
        mask &= mask - 1
    return out


def check_coverage(n_tiles, n, proved, maps):
    """Section 7: every N-subset has some coordinate perm whose image is a proved root.

    Masks are a faithful encoding of subsets (tile indices are small integers).
    Identity is maps[0], so a stored root is recognized on the first probe.
    """
    total = binom(n_tiles, n)
    bit = [1 << i for i in range(n_tiles)]
    proved_masks = set()
    for comb in proved:
        mask = 0
        for i in comb:
            mask |= bit[i]
        proved_masks.add(mask)
    seen = 0
    for comb in itertools.combinations(range(n_tiles), n):
        mask = 0
        for i in comb:
            mask |= bit[i]
        hit = False
        for mp in maps:
            if apply_map(mask, mp, bit) in proved_masks:
                hit = True
                break
        if not hit:
            return False, f"coverage {seen}/{total}", comb
        seen += 1
    if seen != total:
        fail("internal: coverage count")
    return True, f"coverage {total}/{total}", None


def parse_combo(raw, n, n_tiles, what):
    if type(raw) is not list or len(raw) != n:
        fail(f"{what} combo length")
    prev = -1
    out = []
    for x in raw:
        if type(x) is not int or x <= prev or x >= n_tiles:
            fail(f"{what} combo is not a strictly increasing tile subset")
        prev = x
        out.append(x)
    return tuple(out)


def parse_point(raw, G, what):
    if type(raw) is not list or len(raw) != 3:
        fail(f"{what} is not a triple")
    if type(raw[0]) is not int or type(raw[1]) is not int or type(raw[2]) is not int:
        fail(f"{what} coordinate is not an integer")
    if raw[0] < 0 or raw[1] < 0 or raw[2] < 0:
        fail(f"{what} has a negative coordinate")
    if raw[0] + raw[1] + raw[2] != G:
        fail(f"{what} does not sum to G")
    return (raw[0], raw[1], raw[2])


def freeze(value):
    """Structural identity for cross-file compares. Float tokens stay text."""
    if value is None or type(value) is bool or type(value) is int or type(value) is str:
        return value
    if type(value) is RawFloat:
        return ("float", value.text)
    if type(value) is list:
        return tuple(freeze(item) for item in value)
    if type(value) is dict:
        return tuple(sorted((key, freeze(item)) for key, item in value.items()))
    fail("unsupported JSON value")


HEADER_REQUIRED = {"n", "m_t", "k", "G", "Tc"}
HEADER_OPTIONAL = {"m_star", "eps", "tiles", "rule", "rounds", "cap", "part", "priority"}
CLUSTER_REQUIRED = {"id", "combo", "points", "frame", "rho"}
CLUSTER_OPTIONAL = {"value_float", "rattlers_by_cage_slack", "slack_max"}


def parse_part(text):
    if type(text) is not str:
        fail("part is not a string")
    pieces = text.split("/")
    if len(pieces) != 2 or not pieces[0].isdigit() or not pieces[1].isdigit():
        fail("part is not i/p")
    index = int(pieces[0])
    parts = int(pieces[1])
    if parts < 1 or index < 0 or index >= parts:
        fail("part is out of range")


def parse_clusters(raw, n, n_tiles, G):
    if type(raw) is not list:
        fail("accept is not a list")
    by_id = {}
    order = []
    for obj in raw:
        check_keys(obj, CLUSTER_REQUIRED, CLUSTER_OPTIONAL, "cluster")
        cid = obj["id"]
        if type(cid) is not int:
            fail("cluster id is not an integer")
        if cid in by_id:
            fail(f"duplicate cluster id {cid}")
        combo = parse_combo(obj["combo"], n, n_tiles, "cluster")
        points_raw = obj["points"]
        frame = obj["frame"]
        rho = obj["rho"]
        if type(points_raw) is not list or len(points_raw) != n:
            fail("cluster points length")
        if type(frame) is not list or len(frame) != n:
            fail("cluster frame length")
        if type(rho) is not int or rho < 0:
            fail("cluster rho")
        points = tuple(parse_point(p, G, "cluster point") for p in points_raw)
        flags = []
        for flag in frame:
            if type(flag) is not bool:
                fail("cluster frame entry is not a boolean")
            flags.append(flag)
        by_id[cid] = {
            "combo": combo,
            "points": points,
            "frame": tuple(flags),
            "rho": rho,
        }
        order.append(cid)
    return by_id, order


def vertices_in_box(vs, pt, rho):
    pa, pb, pc = pt
    for A, B, C in vs:
        if abs(A - pa) > rho or abs(B - pb) > rho or abs(C - pc) > rho:
            return False
    return True


def check_witness(raw, n, G, Tc):
    if type(raw) is not list or len(raw) != n:
        fail("witness length")
    pts = [parse_point(p, G, "witness point") for p in raw]
    for a in range(n):
        pa = pts[a]
        for b in range(a + 1, n):
            pb = pts[b]
            dA = pa[0] - pb[0]
            dB = pa[1] - pb[1]
            dC = pa[2] - pb[2]
            if dA * dA + dB * dB + dC * dC < Tc:
                fail("witness pair is closer than m_t")


def new_stats():
    return {
        "nodes": 0,
        "ops": 0,
        "discards": 0,
        "empties": 0,
        "splits": 0,
        "accepts": 0,
        "witnesses": 0,
        "undecided": 0,
        "roots": 0,
        "proved_roots": 0,
        "proved_combos": 0,
        "clusters": 0,
    }


def format_report(stats):
    if "n" not in stats:
        return ""
    lines = [
        (
            f"n={stats['n']} k={stats['k']} G={stats['G']} Tc={stats['Tc']} "
            f"m_t={stats['mt']} mode={stats['mode']} files={stats.get('files', 0)} "
            f"clusters={stats.get('clusters', 0)}"
        ),
        (
            f"nodes={stats['nodes']} ops={stats['ops']} discards={stats['discards']} "
            f"empties={stats['empties']} splits={stats['splits']} accepts={stats['accepts']} "
            f"witnesses={stats['witnesses']} undecided={stats['undecided']} "
            f"roots={stats['roots']} proved_roots={stats['proved_roots']} "
            f"proved_combos={stats['proved_combos']}"
        ),
    ]
    if stats.get("pigeon"):
        lines.append("pigeonhole=yes")
    if "coverage" in stats:
        lines.append(stats["coverage"])
    if "claim" in stats:
        lines.append("claim=" + stats["claim"])
    return "\n".join(lines)


class Block:
    def __init__(self, combo, root_id, regions):
        self.combo = combo
        self.ok = True
        self.open = 1
        self.pending = {root_id: regions}


def apply_fate(tag, fate, regions, getv, run, block, where):
    n = run.n
    stats = run.stats
    if tag == "D":
        if len(fate) != 3:
            fail(where + "discard arity")
        i, j = fate[1], fate[2]
        if type(i) is not int or type(j) is not int or not (0 <= i < n and 0 <= j < n) or i == j:
            fail(where + "discard pair")
        vi = getv(i)
        vj = getv(j)
        # Ambiguity 2: a third empty region does not satisfy this test.
        if vi and vj and not pairs_lt(vi, vj, run.Tc):
            fail(where + "discard maxS is not < Tc")
        stats["discards"] += 1
        return
    if tag == "E":
        if len(fate) != 2 or type(fate[1]) is not int or not (0 <= fate[1] < n):
            fail(where + "empty fate")
        if getv(fate[1]):
            fail(where + "empty-fate region is nonempty")
        stats["empties"] += 1
        return
    if tag == "S":
        if len(fate) != 6:
            fail(where + "split arity")
        a, q, c, lo_id, hi_id = fate[1], fate[2], fate[3], fate[4], fate[5]
        if (
            type(a) is not int
            or type(q) is not int
            or type(c) is not int
            or type(lo_id) is not int
            or type(hi_id) is not int
        ):
            fail(where + "split field type")
        if not (0 <= a < n) or q not in (0, 1, 2):
            fail(where + "split index")
        if lo_id == hi_id:
            fail(where + "split children are the same id")
        R = regions[a]
        lo_b = R[2 * q]
        hi_b = R[2 * q + 1]
        if c < lo_b or c > hi_b:
            fail(where + "split plane outside the region")
        if lo_id in run.used or hi_id in run.used:
            fail(where + "split id already used")
        lo_regs = [row[:] for row in regions]
        hi_regs = [row[:] for row in regions]
        lo_regs[a][2 * q + 1] = c
        hi_regs[a][2 * q] = c
        run.used.add(lo_id)
        run.used.add(hi_id)
        block.pending[lo_id] = lo_regs
        block.pending[hi_id] = hi_regs
        block.open += 2
        stats["splits"] += 1
        return
    if tag == "A":
        if run.mode != "cluster":
            fail(where + "accept fate in a refutation")
        if len(fate) != 2 or type(fate[1]) is not int:
            fail(where + "accept arity")
        cid = fate[1]
        cluster = run.clusters.get(cid)
        if cluster is None:
            fail(where + f"unknown cluster {cid}")
        if cluster["combo"] != block.combo:
            fail(where + "accept combo does not equal the root combo")
        rho = cluster["rho"]
        for i in range(n):
            if not cluster["frame"][i]:
                continue
            vs = getv(i)
            # Ambiguity 6.
            if not vs:
                fail(where + f"accept frame region {i} is empty")
            if not vertices_in_box(vs, cluster["points"][i], rho):
                fail(where + f"accept region {i} leaves the cluster box")
        stats["accepts"] += 1
        return
    if tag == "F":
        if len(fate) != 2:
            fail(where + "witness arity")
        check_witness(fate[1], n, run.G, run.Tc)
        stats["witnesses"] += 1
        block.ok = False
        # A checked-out witness is an m_t-feasible configuration.
        if run.mode == "proved":
            fail(where + "witness shows an m_t-feasible configuration")
        return
    if tag == "U":
        if len(fate) != 1:
            fail(where + "undecided arity")
        stats["undecided"] += 1
        block.ok = False
        return
    fail(where + f"unknown fate {tag!r}")


def apply_node(obj, block, run):
    nid = obj["id"]
    if type(nid) is not int:
        fail("node id is not an integer")
    if nid not in block.pending:
        if nid in run.used:
            fail(f"duplicate id {nid}")
        fail(f"unknown node {nid}")
    regions = block.pending.pop(nid)
    block.open -= 1
    if block.open != len(block.pending):
        fail("internal: open count does not match pending nodes")
    n = run.n
    cache = [None] * n

    def getv(i):
        cached = cache[i]
        if cached is None:
            cached = region_vertices(regions[i], run.G)
            cache[i] = cached
        return cached

    ops = obj["ops"]
    if type(ops) is not list:
        fail("ops is not a list")
    run.stats["ops"] += len(ops)
    for oi, op in enumerate(ops):
        where = f"op[{oi}] "
        if type(op) is not list or len(op) != 4:
            fail(where + "must be [i, side, t, j]")
        i, side, t, j = op
        if (
            type(i) is not int
            or type(side) is not int
            or type(t) is not int
            or type(j) is not int
        ):
            fail(where + "field type")
        if not (0 <= i < n and 0 <= j < n) or i == j:
            fail(where + "bad point index")
        if side < 0 or side > 5:
            fail(where + "side out of 0..5")

        def get_other(j=j):
            return getv(j)

        regions[i] = reduced_region(regions[i], side, t, get_other, run.G, run.Tc, where)
        cache[i] = None
    fate = obj["fate"]
    if type(fate) is not list or not fate or type(fate[0]) is not str:
        fail("bad fate")
    apply_fate(fate[0], fate, regions, getv, run, block, "fate ")
    if block.open != len(block.pending):
        fail("internal: open count does not match pending nodes after fate")
    run.stats["nodes"] += 1


def close_if_done(holder, run):
    block = holder["block"]
    if block is not None and block.open == 0:
        if block.pending:
            fail("internal: pending nodes with open == 0")
        if block.ok:
            run.proved.add(block.combo)
            run.stats["proved_roots"] += 1
            run.stats["proved_combos"] = len(run.proved)
        holder["block"] = None


def missing_nodes(block):
    pending = list(block.pending)
    shown = pending[:8]
    extra = "" if len(pending) <= 8 else f" (+{len(pending) - 8} more)"
    return f"missing nodes {shown}{extra}"


def read_header_line(obj, run, expect, seen_header):
    check_keys(obj, {"header", "tiles"}, {"accept"}, "header line")
    header = obj["header"]
    check_keys(header, HEADER_REQUIRED, HEADER_OPTIONAL, "header")
    n = header["n"]
    k = header["k"]
    G = header["G"]
    Tc = header["Tc"]
    if type(n) is not int or n < 1:
        fail("bad n")
    if type(k) is not int or k < 1:
        fail("bad k")
    if type(G) is not int or G <= 0:
        fail("bad G")
    if type(Tc) is not int or Tc <= 0:
        fail("bad Tc")
    if "tiles" in header:
        if type(header["tiles"]) is not int or header["tiles"] != k * k:
            fail("header.tiles count != k^2")
    if "part" in header:
        parse_part(header["part"])
    p, q = parse_mt(header["m_t"])
    got = recompute_Tc(p, q, G)
    if got != Tc:
        fail(f"Tc {Tc} != recomputed {got}")
    if k * k * p * p <= q * q:
        fail("k^2 m_t^2 is not > 1")
    mode = "cluster" if "accept" in obj else "proved"
    if expect == "proved" and mode != "proved":
        fail("expected a refutation, file has accept")
    if expect == "cluster" and mode != "cluster":
        fail("expected a cluster run, file has no accept")
    tiles, g = validate_tiles(obj["tiles"], k, G, Tc)
    clusters = {}
    order = []
    frozen_accept = None
    if mode == "cluster":
        clusters, order = parse_clusters(obj["accept"], n, k * k, G)
        frozen_accept = freeze(obj["accept"])
    if not seen_header:
        run.n = n
        run.k = k
        run.G = G
        run.g = g
        run.Tc = Tc
        run.mt = header["m_t"]
        run.p = p
        run.q = q
        run.mode = mode
        run.tiles = tiles
        run.n_tiles = k * k
        run.maps = symmetry_maps(tiles)
        run.clusters = clusters
        run.cluster_order = order
        run.accept_freeze = frozen_accept
        run.stats["n"] = n
        run.stats["k"] = k
        run.stats["G"] = G
        run.stats["Tc"] = Tc
        run.stats["mt"] = header["m_t"]
        run.stats["mode"] = mode
        run.stats["clusters"] = len(clusters)
        return
    if (
        n != run.n
        or k != run.k
        or G != run.G
        or Tc != run.Tc
        or header["m_t"] != run.mt
        or tiles != run.tiles
        or mode != run.mode
    ):
        fail("part files disagree on n, m_t, k, G, Tc, tiles, or run kind")
    if mode == "cluster" and frozen_accept != run.accept_freeze:
        fail("part files disagree on the cluster list")


def read_file(path, run, expect, seen_header):
    CTX.path = path
    CTX.line = 0
    file_nodes = 0
    holder = {"block": None}
    ended = False
    try:
        fh = gzip.open(path, "rt", encoding="utf-8")
    except OSError as exc:
        fail(str(exc))
    with fh:
        first = True
        for lineno, line in enumerate(fh, 1):
            CTX.line = lineno
            text = line.strip()
            if not text:
                fail("blank line")
            obj = loads(text)
            if type(obj) is not dict:
                fail("record is not an object")
            if first:
                first = False
                read_header_line(obj, run, expect, seen_header)
                seen_header = True
                continue
            if ended:
                fail("trailing record after end")
            keys = set(obj)
            if "end" in keys or "nodes" in keys:
                check_keys(obj, {"end", "nodes"}, set(), "end")
                if type(obj["end"]) is not str:
                    fail("end status is not a string")
                # LEAD'S PATCH (2026-09-26): TREE_FORMAT.md §3 says the end status is informational and must not be trusted;
                # a cluster-run part file whose roots all end without accept leaves legitimately says "PROVED". Only the
                # status VALUE set is checked here; the claim itself is decided from the replayed nodes and coverage below.
                allowed_status = {"PROVED"} if run.mode == "proved" else {"PROVED", "PROVED_MOD_CLUSTERS"}
                if obj["end"] not in allowed_status:
                    fail(f"end status {obj['end']!r} is not in {sorted(allowed_status)}")
                if type(obj["nodes"]) is not int or obj["nodes"] != file_nodes:
                    fail(f"end nodes {obj['nodes']} != node records {file_nodes}")
                block = holder["block"]
                if block is not None:
                    fail(missing_nodes(block))
                ended = True
                continue
            if "pigeonhole" in keys:
                check_keys(obj, {"pigeonhole"}, set(), "pigeonhole")
                if run.mode != "proved":
                    fail("pigeonhole record in a cluster run")
                if run.pigeon or run.stats["roots"] or holder["block"] is not None:
                    fail("pigeonhole mixed with roots or repeated")
                ph = obj["pigeonhole"]
                if (
                    type(ph) is not list
                    or len(ph) != 2
                    or type(ph[0]) is not int
                    or type(ph[1]) is not int
                ):
                    fail("pigeonhole is not [k^2, n]")
                if ph[0] != run.k * run.k or ph[1] != run.n:
                    fail("pigeonhole value is not [k^2, n]")
                if not (run.k * run.k < run.n):
                    fail("pigeonhole claimed but k^2 >= n")
                run.pigeon = True
                run.stats["pigeon"] = True
                continue
            if "root" in keys:
                check_keys(obj, {"id", "root"}, set(), "root")
                if holder["block"] is not None:
                    fail(missing_nodes(holder["block"]))
                if run.pigeon:
                    fail("root after a pigeonhole record")
                rid = obj["id"]
                if type(rid) is not int:
                    fail("root id is not an integer")
                if rid in run.used:
                    fail(f"duplicate id {rid}")
                combo = parse_combo(obj["root"], run.n, run.n_tiles, "root")
                regions = [
                    tile_region(run.tiles[s][0], run.tiles[s][1], run.tiles[s][2], run.g)
                    for s in combo
                ]
                run.used.add(rid)
                holder["block"] = Block(combo, rid, regions)
                run.stats["roots"] += 1
                continue
            if "ops" in keys or "fate" in keys or "id" in keys:
                check_keys(obj, {"id", "ops", "fate"}, set(), "node")
                block = holder["block"]
                if block is None:
                    fail("node outside a root block")
                apply_node(obj, block, run)
                file_nodes += 1
                close_if_done(holder, run)
                continue
            fail(f"unknown record keys {sorted(keys)}")
        if first:
            fail("missing header")
        if not ended:
            fail("missing end record")
    CTX.path = ""
    CTX.line = 0
    return seen_header


def verify(paths, expect, stats):
    run = Run(stats)
    seen = False
    for path in paths:
        # used-ids restart in each part file; proved roots accumulate.
        run.used = set()
        seen = read_file(path, run, expect, seen)
    if run.mode is None:
        fail("no header")
    if holder_open(run):
        fail("internal: block left open")
    if run.k * run.k < run.n:
        if not run.pigeon:
            fail("k^2 < n requires a pigeonhole record")
        if run.stats["roots"]:
            fail("pigeonhole file contains roots")
        stats["coverage"] = "coverage=pigeonhole"
        stats["claim"] = "no m_t-feasible configuration (pigeonhole)"
        return format_report(stats)
    if run.pigeon:
        fail("pigeonhole record but k^2 >= n")
    ok, msg, missing = check_coverage(run.n_tiles, run.n, run.proved, run.maps)
    stats["coverage"] = msg
    if not ok:
        tail = "" if missing is None else f" missing {list(missing)}"
        fail("coverage incomplete" + tail)
    if run.mode == "proved":
        stats["claim"] = "no m_t-feasible configuration"
    else:
        stats["claim"] = (
            "every m_t-feasible configuration meets a cluster box after a symmetry"
        )
    return format_report(stats)


def holder_open(run):
    return False


class Run:
    def __init__(self, stats):
        self.stats = stats
        self.proved = set()
        self.used = set()
        self.mode = None
        self.n = None
        self.k = None
        self.G = None
        self.g = None
        self.Tc = None
        self.mt = None
        self.p = None
        self.q = None
        self.tiles = None
        self.maps = None
        self.n_tiles = None
        self.clusters = {}
        self.cluster_order = []
        self.accept_freeze = None
        self.pigeon = False


def self_test():
    if recompute_Tc(1, 2, 3) != 5:
        fail("self-test: ceil(4.5)")
    if recompute_Tc(1, 2, 4) != 8:
        fail("self-test: exact 8")
    if recompute_Tc(1, 3, 3) != 2:
        fail("self-test: exact 2")
    published = (
        (115481600843309, 200000000000000, 33554432, 750750065388662),
        (144193229730109, 250000000000000, 33554432, 749099488619231),
        (216248892036713, 1000000000000000, 41943040, 164534731528431),
        (6757034594909, 31250000000000, 41943040, 164498539497368),
    )
    for p, q, G, Tc in published:
        got = recompute_Tc(p, q, G)
        if got != Tc:
            fail(f"self-test: Tc recomputed {got} != {Tc}")
    p, q, G, Tc = published[0]
    if 4 * p * p <= q * q:
        fail("self-test: ex4 separation")
    if binom(25, 16) != 2042975 or binom(4, 4) != 1 or binom(9, 9) != 1 or binom(4, 2) != 6:
        fail("self-test: binomial")
    tiles = expected_tiles(2)
    if tiles != [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0)]:
        fail("self-test: k=2 tiles")
    if len(expected_tiles(3)) != 9 or len(expected_tiles(5)) != 25:
        fail("self-test: tile counts")
    _exp, g = validate_tiles([[a, b, c] for a, b, c in tiles], 2, G, Tc)
    if g != 16777216:
        fail("self-test: g")
    if side_indices(0) != (0, 1, 0) or side_indices(1) != (0, 0, 1):
        fail("self-test: side 0/1")
    if side_indices(2) != (1, 3, 2) or side_indices(5) != (2, 4, 5):
        fail("self-test: side 2/5")
    t = 4428131
    central = [0, g, 0, g, 0, g]
    strip = [0, t, 0, g, 0, g]
    vs = region_vertices(strip, G)
    oracle = {(0, g, g), (t, g, G - t - g), (t, G - t - g, g)}
    if vs != oracle:
        fail("self-test: ex4 strip vertices")
    tile1 = region_vertices(tile_region(0, 0, 1, g), G)
    if tile1 != {(0, 0, G), (0, g, g), (g, 0, g)}:
        fail("self-test: tile (0,0,1) vertices")
    if not pairs_lt(vs, tile1, Tc):
        fail("self-test: ex4 first strip is not < Tc")
    out = reduced_region(central, 0, t, lambda: tile1, G, Tc, "self-test ")
    if out != [t, g, 0, g, 0, g] or central != [0, g, 0, g, 0, g]:
        fail("self-test: low-side effect")
    if pairs_lt([(G, 0, 0)], [(0, 0, G)], Tc) or 2 * G * G < Tc:
        fail("self-test: unit-side pair")
    if not pairs_lt([(0, 0, G)], [(0, 0, G)], Tc):
        fail("self-test: zero distance")
    maps = symmetry_maps(tiles)
    if maps[1] != [0, 2, 1, 3]:
        fail("self-test: B/C swap")
    bit = [1 << i for i in range(4)]
    if apply_map((1 << 1) | (1 << 3), maps[1], bit) != ((1 << 2) | (1 << 3)):
        fail("self-test: apply_map")
    ok, _msg, missing = check_coverage(4, 4, {(0, 1, 2, 3)}, maps)
    if not ok:
        fail("self-test: full cover")
    ok, _msg, missing = check_coverage(4, 4, set(), maps)
    if ok or list(missing) != [0, 1, 2, 3]:
        fail("self-test: coverage gap")
    lo, hi = None, None
    regions = [[0, 10, 0, 10, 0, 10]]
    lo = [row[:] for row in regions]
    hi = [row[:] for row in regions]
    lo[0][3] = 5
    hi[0][2] = 5
    if lo != [[0, 10, 0, 5, 0, 10]] or hi != [[0, 10, 5, 10, 0, 10]]:
        fail("self-test: split halves")
    if region_vertices([0, 0, 0, 0, 0, 0], 2) != set():
        fail("self-test: empty region")
    if region_vertices([0, 0, 0, 0, 2, 2], 2) != {(0, 0, 2)}:
        fail("self-test: point region")


def parse_args(argv):
    files = []
    expect = None
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--expect" or arg.startswith("--expect="):
            if expect is not None:
                fail("duplicate --expect")
            if arg == "--expect":
                if i + 1 >= len(argv):
                    fail("missing --expect value")
                expect = argv[i + 1]
                i += 2
            else:
                expect = arg.split("=", 1)[1]
                i += 1
            continue
        if arg.startswith("-"):
            fail(f"unknown argument {arg}")
        files.append(arg)
        i += 1
    if not files:
        fail(
            "usage: python verify_tri_tree.py <tree.jsonl.gz> "
            "[more part files...] [--expect proved|cluster]"
        )
    if expect not in (None, "proved", "cluster"):
        fail("--expect must be proved or cluster")
    if len(set(files)) != len(files):
        fail("duplicate file path")
    return files, expect


def main(argv):
    stats = new_stats()
    try:
        self_test()
        paths, expect = parse_args(argv)
        stats["files"] = len(paths)
        report = verify(paths, expect, stats)
    except CertError as exc:
        partial = format_report(stats)
        if partial:
            print(partial)
        print("VERDICT: FAIL " + exc.reason)
        return 1
    print(report)
    print("VERDICT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
