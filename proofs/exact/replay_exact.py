#!/usr/bin/env python3
"""Independent replay checker for the exact optimality certificates (Moki&Julio).
Status: reviewed (three adversarial rounds, 2026-09-25): PUBLISHABLE; published in v1.7. Mutation test: mutate_check.py.

Reads <case>_cert.json and <case>_tree.json.gz from this directory and re-verifies every step with exact arithmetic only:
fractions.Fraction, and pairs (a, b) meaning a + b*sqrt2 with a, b Fractions (exact field operations and an exact sign test).
Standard library only. It does NOT import exact_prove.py or grok/prove_small.py, and it takes no number in the files on
trust: every inequality that the proof uses is recomputed here.

Lemma 1 data (certificate): sqrt2 enclosure; c* feasible at r*; tight set T = the exact zero set; 0 < slack_lo <= slack for
every other constraint; stress lambda in Q(sqrt2)
with equilibrium and normalisation; rational lower bounds lambda_lo > 0; Lambda_hi; basis B, M recomputed from c*, M Minv = I
and Minv M = I in Q(sqrt2); ||Minv||_inf <= normMinv_hi (exact absolute values); rho by the stated formula; r_t < r*;
snapped images and the snap bound; rho' + snap_err_hi <= rho.
Lemma 2 data (tree): root = [0,G]^(2N) in the grid chart (walls), cells cover the root and each has diameter < d, every
injective assignment of cells accounted for; every reduction step valid (non-strict symmetry cuts, strip removals with a
pair range-maximum < d^2); every split covers its parent; every discard has an exact witness; every acceptance is the exact
subset test; the tree is finite with no undecided leaf. Positive control: every representative of c* that satisfies the
cuts lies in an accepted leaf.

usage: py -3.11 replay_exact.py [case ...]        (default: crc_800_4 crc_600_6)
exit code 0 iff every requested case prints PASS.
"""
import sys, os, json, gzip, hashlib, itertools
from fractions import Fraction as Fr
from math import lcm

DIR = os.path.dirname(os.path.abspath(__file__))


class Fail(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise Fail(msg)


# ------------------------------------------------------------------ Q(sqrt2): x = (a, b) = a + b sqrt2, a, b Fractions
def qp(v):
    return (Fr(v[0]), Fr(v[1]))


def qr(a):
    return (Fr(a), Fr(0))


def qadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def qsub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def qmul(x, y):
    return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def qneg(x):
    return (-x[0], -x[1])


def qsgn(x):
    a, b = x
    if a == 0 and b == 0:
        return 0
    if a >= 0 and b >= 0:
        return 1
    if a <= 0 and b <= 0:
        return -1
    t = a * a - 2 * b * b          # a, b of strictly opposite signs; t != 0 because sqrt2 is irrational
    need(t != 0, "impossible: a^2 = 2 b^2 with rational a, b != 0")
    if a > 0:                      # a + b sqrt2 > 0  <=>  a > |b| sqrt2  <=>  a^2 > 2 b^2
        return 1 if t > 0 else -1
    return 1 if t < 0 else -1      # a < 0 < b: positive <=> 2 b^2 > a^2


def qabs(x):
    return x if qsgn(x) >= 0 else qneg(x)


def qsum(xs):
    s = (Fr(0), Fr(0))
    for x in xs:
        s = qadd(s, x)
    return s


def qzero(x):
    return x[0] == 0 and x[1] == 0


# ----------------------------------------------------------------------------------------------- certificate (Lemma 1)
def all_constraints(c, r, h):
    """(descriptor, value g, gradient row, w) in the order walls L,R,B,T per circle, then pairs i<j"""
    n = len(c); out = []
    half = qr(Fr(1, 2)); hh = qr(h / 2)
    for i, (x, y) in enumerate(c):
        for side, g, v, s in (("L", qsub(qadd(x, half), r), 2 * i, 1), ("R", qsub(qsub(half, x), r), 2 * i, -1),
                              ("B", qsub(qadd(y, hh), r), 2 * i + 1, 1), ("T", qsub(qsub(hh, y), r), 2 * i + 1, -1)):
            row = [qr(0)] * (2 * n); row[v] = qr(s)
            out.append(({"kind": "wall", "i": i, "side": side}, g, row, qr(1)))
    for i in range(n):
        for j in range(i + 1, n):
            dx, dy = qsub(c[i][0], c[j][0]), qsub(c[i][1], c[j][1])
            g = qsub(qadd(qmul(dx, dx), qmul(dy, dy)), qmul(qr(4), qmul(r, r)))
            row = [qr(0)] * (2 * n)
            row[2 * i], row[2 * i + 1] = qmul(qr(2), dx), qmul(qr(2), dy)
            row[2 * j], row[2 * j + 1] = qmul(qr(-2), dx), qmul(qr(-2), dy)
            out.append(({"kind": "pair", "i": i, "j": j}, g, row, qmul(qr(8), r)))
    return out


def check_cert(cert, case, log):
    need(cert.get("case") == case, "certificate case name")
    n = cert["N"]; h = Fr(cert["h"]); need(n >= 2 and h > 0 and Fr(cert["width"]) == 1, "N, h, width")
    S2_LO, S2_HI = Fr(cert["S2_LO"]), Fr(cert["S2_HI"])
    need(0 < S2_LO < S2_HI and S2_LO * S2_LO < 2 < S2_HI * S2_HI, "sqrt2 enclosure S2_LO^2 < 2 < S2_HI^2")
    p, q = Fr(cert["r_star"]["p"]), Fr(cert["r_star"]["q"])
    r = (p, -q)
    need(q > 0 and qsgn(r) > 0, "r* = p - q sqrt2 with q > 0 and r* > 0")
    c = [(qp(x), qp(y)) for x, y in cert["c_star"]]
    need(len(c) == n, "c* has N centres")
    cons = all_constraints(c, r, h)
    for desc, g, _, _ in cons:
        need(qsgn(g) >= 0, f"c* violates {desc} at r*")
    T = [k for k in cons if qzero(k[1])]
    need([k[0] for k in T] == cert["T"], "tight set T is not the exact zero set of the constraints at (c*, r*)")
    slack_lo = cert["slack_lo"]                              # display bounds quoted in the write-up: verify them too
    need(len(slack_lo) == len(cons) - len(T), "slack_lo must list every non-tight constraint")
    for desc, g, _, _ in cons:
        if qzero(g):
            continue
        name = (f"wall {desc['side']} of circle {desc['i']}" if desc["kind"] == "wall" else f"pair ({desc['i']},{desc['j']})")
        need(name in slack_lo, f"slack_lo missing {name}")
        lo = Fr(slack_lo[name])
        need(lo > 0 and qsgn(qsub(g, qr(lo))) >= 0, f"slack_lo of {name} is not a positive lower bound of the exact slack")
    log(f"  c* in F(r*): all {len(cons)} constraints >= 0 exactly; tight set T = {len(T)} "
        f"({sum(1 for k in T if k[0]['kind'] == 'pair')} pairs, {sum(1 for k in T if k[0]['kind'] == 'wall')} walls), "
        f"the other {len(cons) - len(T)} strictly positive, each >= its certificate slack_lo > 0")
    m = len(T)
    lam = [qp(v) for v in cert["lambda"]]; need(len(lam) == m, "lambda length")
    for v in range(2 * n):
        need(qzero(qsum(qmul(lam[k], T[k][2][v]) for k in range(m))), f"equilibrium fails in coordinate {v}")
    need(qsum(qmul(lam[k], T[k][3]) for k in range(m)) == qr(1), "normalisation sum lambda_k w_k = 1")
    lam_lo = [Fr(v) for v in cert["lambda_lo"]]; need(len(lam_lo) == m, "lambda_lo length")
    for k in range(m):
        need(lam_lo[k] > 0, f"lambda_lo[{k}] > 0")
        need(qsgn(qsub(lam[k], qr(lam_lo[k]))) >= 0, f"lambda[{k}] >= lambda_lo[{k}]")
    Lam_hi = Fr(cert["Lambda_hi"])
    need(qsgn(qsub(qr(Lam_hi), qsum(lam[k] for k in range(m) if T[k][0]["kind"] == "pair"))) >= 0, "Lambda_hi >= sum pair lambda")
    log(f"  stress: equilibrium (2N = {2 * n} coordinates) and normalisation exact in Q(sqrt2); every lambda_k >= lambda_lo_k > 0; "
        f"Lambda <= Lambda_hi")
    B = cert["B"]
    need(len(B) == 2 * n and len(set(B)) == 2 * n and all(0 <= k < m for k in B), "basis B: 2N distinct indices of T")
    M = [T[k][2] for k in B]
    need([[list(map(str, x)) for x in row] for row in M] == cert["M"], "M in the certificate differs from the rows of B")
    Mi = [[qp(x) for x in row] for row in cert["Minv"]]
    need(len(Mi) == 2 * n and all(len(row) == 2 * n for row in Mi), "Minv shape")
    for X, Y, what in ((M, Mi, "M Minv"), (Mi, M, "Minv M")):
        for a in range(2 * n):
            for b in range(2 * n):
                need(qsum(qmul(X[a][k], Y[k][b]) for k in range(2 * n)) == qr(1 if a == b else 0), f"{what} != I at ({a},{b})")
    normM_hi = Fr(cert["normMinv_hi"])
    for a in range(2 * n):
        need(qsgn(qsub(qr(normM_hi), qsum(qabs(x) for x in Mi[a]))) >= 0, f"row {a}: sum |Minv| > normMinv_hi")
    rho = Fr(cert["rho"])
    need(rho == min(lam_lo[k] for k in B) / (2 * normM_hi) / (16 * Lam_hi), "rho formula")
    need(rho > 0, "rho > 0")
    log(f"  basis: M (from c*) times Minv = I and Minv times M = I exactly; ||Minv||_inf <= {float(normM_hi):.6g}; "
        f"rho = min_B lambda_lo / (2 normMinv_hi) / (16 Lambda_hi) = {float(rho):.9g} (display)")
    r_t = Fr(cert["r_t"])
    need(r_t == p - q * S2_HI - Fr(1, 10 ** 40), "r_t formula")
    need(r_t > 0 and qsgn(qsub(r, qr(r_t))) > 0, "0 < r_t < r*")
    snap = Fr(cert["snap_err_hi"]); rho_p = Fr(cert["rho_prime"])
    bmax = max(abs(v[1]) for pt in c for v in pt)
    need(snap >= bmax * (S2_HI - S2_LO), "snap_err_hi >= max|b| (S2_HI - S2_LO)")
    need(rho_p > 0 and rho_p + snap <= rho, "rho' > 0 and rho' + snap_err_hi <= rho")
    imgs = cert["images_point_frame"]
    need(sorted(tuple(e["mirror"]) for e in imgs) == sorted([(1, 1), (-1, 1), (1, -1), (-1, -1)]), "the four mirrors")
    images, exact_images = [], []
    for e in imgs:
        sx, sy = e["mirror"]; pts, ex = [], []
        for (x, y), (Xs, Ys) in zip(c, e["points"]):
            X = sx * (x[0] + x[1] * S2_LO) + Fr(1, 2) - r_t
            Y = sy * (y[0] + y[1] * S2_LO) + h / 2 - r_t
            need(Fr(Xs) == X and Fr(Ys) == Y, "snapped image point differs from its definition")
            xe = qadd(qmul(qr(sx), x), qr(Fr(1, 2) - r_t)); ye = qadd(qmul(qr(sy), y), qr(h / 2 - r_t))
            need(qsgn(qsub(qr(snap), qabs(qsub(qr(X), xe)))) > 0 and qsgn(qsub(qr(snap), qabs(qsub(qr(Y), ye)))) > 0,
                 "exact snap error exceeds snap_err_hi")
            pts.append((X, Y)); ex.append((xe, ye))
        images.append(pts); exact_images.append(ex)
    log(f"  r_t = p - q S2_HI - 1e-40 < r* (exact); snapped images match; |snap - exact| < snap_err_hi = {float(snap):.3g}; "
        f"rho' + snap_err_hi <= rho")
    return dict(n=n, h=h, r=r, c=c, r_t=r_t, rho_p=rho_p, images=images, mirrors=[tuple(e["mirror"]) for e in imgs],
                exact_images=exact_images, G=cert["grid_G"])


# --------------------------------------------------------------------------------------------------------- tree (Lemma 2)
def far(a0, a1, b0, b1):
    """max |u - v| over u in [a0, a1], v in [b0, b1] (both non-empty): the larger of the two end-point differences"""
    return max(a1 - b0, b1 - a0)


def check_tree(tree, cert_blob, ctx, case, log):
    n, G, r_t, h = ctx["n"], ctx["G"], ctx["r_t"], ctx["h"]
    need(tree.get("case") == case and tree["N"] == n and Fr(tree["h"]) == h and Fr(tree["r_t"]) == r_t, "tree header vs cert")
    need(tree["cert_sha256"] == hashlib.sha256(cert_blob).hexdigest(), "tree is not bound to this certificate (sha256)")
    need(tree["G"] == G and G >= 2 and G % 2 == 0, "grid G even (the y-cut Y_0 <= G/2 must be exact)")
    W, H, d = 1 - 2 * r_t, h - 2 * r_t, 2 * r_t
    need(W > 0 and H > 0 and d > 0, "W, H, d > 0")
    need(Fr(tree["W"]) == W and Fr(tree["H"]) == H and Fr(tree["d"]) == d, "W = 1 - 2 r_t, H = h - 2 r_t, d = 2 r_t")
    sW, sH = W / G, H / G
    # pair test: (fx W/G)^2 + (fy H/G)^2 < d^2, over one common positive denominator
    a2, b2, d2 = sW * sW, sH * sH, d * d
    L = lcm(a2.denominator, b2.denominator, d2.denominator)
    KA, KB, KC = a2.numerator * (L // a2.denominator), b2.numerator * (L // b2.denominator), d2.numerator * (L // d2.denominator)

    def lt_d(fx, fy):
        return fx * fx * KA + fy * fy * KB < KC

    def box_pair_lt(P, Q):          # max over P x Q of |p - q|^2 < d^2
        return lt_d(far(P[0], P[1], Q[0], Q[1]), far(P[2], P[3], Q[2], Q[3]))

    rho_p, images = ctx["rho_p"], ctx["images"]
    stats = {"ops": 0, "discard": {}, "accept": 0, "split": {}, "start_discard": {}}

    def discard_ok(S, w):
        t = w[0]
        if t == "empty":
            i, ax = w[1], w[2]; need(0 <= i < n and ax in (0, 1), "empty witness index")
            need(S[i][2 * ax] > S[i][2 * ax + 1], "empty witness: interval is not empty")
        elif t == "xorder":
            i = w[1]; need(0 <= i < n - 1, "xorder index")
            need(S[i][0] > S[i + 1][1], "xorder witness: X_i >= lo_i > hi_{i+1} >= X_{i+1} fails")
        elif t == "xsum":
            need(S[0][0] + S[n - 1][0] > G, "xsum witness: lo_1 + lo_N > G fails")
        elif t in ("pair", "forbid"):
            i, j = w[1], w[2]; need(0 <= i < n and 0 <= j < n and i != j, "pair witness indices")
            fx, fy = far(S[i][0], S[i][1], S[j][0], S[j][1]), far(S[i][2], S[i][3], S[j][2], S[j][3])
            if t == "pair" and len(w) > 3:
                need([fx, fy] == w[3:5], "pair witness: stated corner distances differ")
            need(lt_d(fx, fy), f"{t} witness: range-maximum of |p_i - p_j|^2 is not < d^2")
        else:
            raise Fail(f"unknown discard witness {w}")
        return t

    def apply_ops(S, ops):
        for op in ops:
            t = op[0]
            if t == "o_hi":                     # X_i <= X_{i+1} <= hi_{i+1}
                i = op[1]; need(0 <= i < n - 1, "o_hi index"); S[i][1] = min(S[i][1], S[i + 1][1])
            elif t == "o_lo":                   # X_{i+1} >= X_i >= lo_i
                i = op[1]; need(0 <= i < n - 1, "o_lo index"); S[i + 1][0] = max(S[i + 1][0], S[i][0])
            elif t == "s1":                     # X_1 <= G - X_N <= G - lo_N
                S[0][1] = min(S[0][1], G - S[n - 1][0])
            elif t == "sN":                     # X_N <= G - X_1 <= G - lo_1
                S[n - 1][1] = min(S[n - 1][1], G - S[0][0])
            elif t == "y1":                     # Y_1 <= G/2
                S[0][3] = min(S[0][3], G // 2)
            elif t == "sh":                     # strip of box i whose pair range-maximum with box j is < d^2
                _, i, j, side, v = op
                need(0 <= i < n and 0 <= j < n and i != j, "sh indices")
                x0, x1, y0, y1 = S[i]
                if side == "x0":
                    need(x0 <= v <= x1, "sh x0 range"); strip = (x0, v, y0, y1)
                elif side == "x1":
                    need(x0 <= v <= x1, "sh x1 range"); strip = (v, x1, y0, y1)
                elif side == "y0":
                    need(y0 <= v <= y1, "sh y0 range"); strip = (x0, x1, y0, v)
                elif side == "y1":
                    need(y0 <= v <= y1, "sh y1 range"); strip = (x0, x1, v, y1)
                else:
                    raise Fail("sh side")
                need(box_pair_lt(strip, S[j]), f"sh {side}: strip is not excluded by box {j}")
                S[i][("x0", "x1", "y0", "y1").index(side)] = v
            else:
                raise Fail(f"unknown op {op}")
            stats["ops"] += 1

    def accept_ok(S, mi, perm):
        need(0 <= mi < len(images) and sorted(perm) == list(range(n)), "accept: mirror / permutation")
        pts = images[mi]
        for i in range(n):
            px, py = pts[perm[i]]
            need(px - rho_p <= S[i][0] * sW and S[i][1] * sW <= px + rho_p and
                 py - rho_p <= S[i][2] * sH and S[i][3] * sH <= py + rho_p, f"accept: box {i} not inside the rho'-square")

    nodes = tree["nodes"]
    need(isinstance(nodes, list) and len(nodes) >= 1, "nodes")
    for k, rec in enumerate(nodes):
        need(rec["id"] == k, "node ids must be 0..K-1 in order")
    root = nodes[0]
    need(root["parent"] is None and root["box"] == [[0, G, 0, G]] * n and root["ops"] == [], "root = [0,G]^(2N)")
    expected = {}                                             # child id -> (parent id, box it must have)
    rf = root["fate"]
    if rf["t"] == "cells":
        cells = [tuple(cc) for cc in rf["cells"]]; m = len(cells)
        need(m == len(set(cells)) and m >= n, "cells distinct, at least N")
        xs = sorted(set((cc[0], cc[1]) for cc in cells)); ys = sorted(set((cc[2], cc[3]) for cc in cells))
        need(set(cells) == {(x[0], x[1], y[0], y[1]) for x in xs for y in ys}, "cells = product of x- and y-intervals")
        for iv in (xs, ys):
            reach = None
            for lo, hi in iv:
                need(lo <= hi, "cell interval")
                if reach is None:
                    need(lo <= 0, "cover starts at 0"); reach = hi
                else:
                    need(lo <= reach, "gap in the cell cover"); reach = max(reach, hi)
            need(reach >= G, "cover ends at G")
        for cc in cells:
            need(lt_d(cc[1] - cc[0], cc[3] - cc[2]), "a cell can hold two points (diameter >= d)")
        need(len(rf["children"]) == len(rf["assign"]), "children / assign lengths")
        table = {}
        for cid, asg in zip(rf["children"], rf["assign"]):
            key = tuple(asg); need(key not in table, "duplicate assignment")
            table[key] = ("child", cid)
        for asg, w in rf["discards"]:
            key = tuple(asg); need(key not in table, "duplicate assignment")
            table[key] = ("disc", w)
        count = 0
        for key in itertools.permutations(range(m), n):
            need(key in table, f"assignment {key} missing"); count += 1
            box = [list(cells[a]) for a in key]
            kind, val = table[key]
            if kind == "child":
                need(val not in expected, "child claimed twice"); expected[val] = (0, box)
            else:
                t = discard_ok(box, val); stats["start_discard"][t] = stats["start_discard"].get(t, 0) + 1
        need(count == len(table), "extra assignments in the log")
        log(f"  root: {m} cells cover [0,G]^2 (product cover), each of diameter < d; all {count} injective assignments "
            f"accounted for ({len(rf['children'])} children, {len(rf['discards'])} discarded with witnesses)")
    elif rf["t"] == "full":
        need(len(rf["children"]) == 1, "full: one child"); expected[rf["children"][0]] = (0, [[0, G, 0, G]] * n)
    else:
        raise Fail("root fate")
    tboxes = {}
    for rec in nodes[1:]:
        k = rec["id"]
        need(k in expected, f"node {k} is not a child of any processed node")
        par, box = expected.pop(k)
        need(rec["parent"] == par and par < k, f"node {k}: parent")
        need(rec["box"] == box, f"node {k}: box differs from what its parent's split produced")
        S = [list(b) for b in box]
        apply_ops(S, rec["ops"])
        fate = rec["fate"]; t = fate["t"]
        if "tbox" in rec:
            need(rec["tbox"] == S, f"node {k}: stored tightened box differs from the replay")
        if t == "discard":
            wt = discard_ok(S, fate["w"]); stats["discard"][wt] = stats["discard"].get(wt, 0) + 1
        elif t == "accept":
            accept_ok(S, fate["mirror"], fate["perm"]); stats["accept"] += 1; tboxes[k] = [tuple(b) for b in S]
        elif t == "split":
            how = fate["how"]; kids = []
            if how[0] == "bisect":
                _, i, ax, mid = how; need(0 <= i < n and ax in (0, 1), "bisect args")
                lo, hi = S[i][2 * ax], S[i][2 * ax + 1]; need(lo <= mid <= hi, "bisect: mid inside")
                for a, b in ((lo, mid), (mid, hi)):
                    bi = list(S[i]); bi[2 * ax], bi[2 * ax + 1] = a, b
                    kids.append([list(x) for x in S[:i]] + [bi] + [list(x) for x in S[i + 1:]])
            elif how[0] == "hole":
                _, i, j, hole = how; hx0, hx1, hy0, hy1 = hole
                need(0 <= i < n and 0 <= j < n and i != j, "hole indices")
                x0, x1, y0, y1 = S[i]
                need(x0 <= hx0 < hx1 <= x1 and y0 <= hy0 < hy1 <= y1, "hole: non-degenerate, inside box i")
                need(box_pair_lt((hx0, hx1, hy0, hy1), S[j]), "hole: not excluded by box j")
                parts = []                       # box i minus the hole is inside the union of these closed parts
                if x0 < hx0:
                    parts.append([x0, hx0, y0, y1])
                if hx1 < x1:
                    parts.append([hx1, x1, y0, y1])
                if y0 < hy0:
                    parts.append([hx0, hx1, y0, hy0])
                if hy1 < y1:
                    parts.append([hx0, hx1, hy1, y1])
                for part in parts:
                    kids.append([list(x) for x in S[:i]] + [part] + [list(x) for x in S[i + 1:]])
            else:
                raise Fail(f"node {k}: unknown split {how[0]}")
            need(len(kids) == len(fate["children"]) and len(kids) >= 1, f"node {k}: children count")
            for cid, kb in zip(fate["children"], kids):
                need(isinstance(cid, int) and cid > k and cid not in expected, f"node {k}: child id")
                expected[cid] = (k, kb)
            stats["split"][how[0]] = stats["split"].get(how[0], 0) + 1
        else:
            raise Fail(f"node {k}: fate '{t}' is neither discard, accept nor split (UNDECIDED leaf)")
    need(not expected, f"children never listed as nodes: {sorted(expected)[:5]}")
    nd = sum(stats["discard"].values())
    log(f"  tree: {len(nodes)} nodes, {stats['ops']} reduction steps replayed; splits {stats['split']}; leaves: "
        f"{stats['accept']} accepted + {nd} discarded {stats['discard']}; start assignments discarded {stats['start_discard']}; "
        f"0 undecided -> COMPLETE")
    # positive control: every labelling / mirror image of c* that satisfies the cuts must sit in an accepted leaf
    c, r_t = ctx["c"], ctx["r_t"]; reps = found = 0
    for (sx, sy), ex in zip(ctx["mirrors"], ctx["exact_images"]):
        for perm in itertools.permutations(range(n)):
            P = [ex[perm[i]] for i in range(n)]
            ok = all(qsgn(qsub(P[i + 1][0], P[i][0])) >= 0 for i in range(n - 1))
            ok = ok and qsgn(qsub(qr(W), qadd(P[0][0], P[n - 1][0]))) >= 0 and qsgn(qsub(qr(H / 2), P[0][1])) >= 0
            if not ok:
                continue
            reps += 1
            hit = any(all(qsgn(qsub(P[i][0], qr(b[i][0] * sW))) >= 0 and qsgn(qsub(qr(b[i][1] * sW), P[i][0])) >= 0 and
                          qsgn(qsub(P[i][1], qr(b[i][2] * sH))) >= 0 and qsgn(qsub(qr(b[i][3] * sH), P[i][1])) >= 0
                          for i in range(n)) for b in tboxes.values())
            need(hit, "positive control: a canonical representative of c* lies in no accepted leaf")
            found += 1
    need(reps >= 1, "positive control: no canonical representative of c*")
    log(f"  positive control: all {reps} canonical representatives of c* (exact, Q(sqrt2)) lie in accepted leaves")


def run(case):
    lines = []
    log = lines.append
    try:
        with open(os.path.join(DIR, f"{case}_cert.json"), "rb") as f:
            blob = f.read()
        cert = json.loads(blob)
        ctx = check_cert(cert, case, log)
        with gzip.open(os.path.join(DIR, f"{case}_tree.json.gz"), "rt", encoding="utf-8") as f:
            tree = json.load(f)
        check_tree(tree, blob, ctx, case, log)
        ok, why = True, ""
    except Fail as e:
        ok, why = False, str(e)
    except (KeyError, TypeError, ValueError, IndexError, OSError) as e:
        ok, why = False, f"malformed input: {type(e).__name__}: {e}"
    print(f"== {case}")
    for s in lines:
        print(s)
    print(f"PASS {case}" if ok else f"FAIL {case}: {why}", flush=True)
    return ok


if __name__ == "__main__":
    cases = sys.argv[1:] or ["crc_800_4", "crc_600_6"]
    res = [run(cs) for cs in cases]
    sys.exit(0 if all(res) else 1)
