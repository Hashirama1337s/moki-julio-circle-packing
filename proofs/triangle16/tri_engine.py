#!/usr/bin/env python3
"""Exact-integer tile-and-region branch-and-bound for N points in the unit-side equilateral triangle.  Moki&Julio.

Question answered: is there a placement of N points in the closed unit triangle with all pairwise distances >= m_t ?
PROVED means no (every tile combination refuted); FEASIBLE means an exactly checked placement exists; UNDECIDED means a
cap was hit.

Frame: barycentric grid A + B + C = G (A = G a, ...); squared distance of a difference = (dA^2 + dB^2 + dC^2) / (2 G^2).
A region is a closed tri-box {A in [A0,A1], B in [B0,B1], C in [C0,C1]} with integer bounds, kept TIGHT (each bound is
attained).  Every discard / reduction decision is an exact int64 comparison against Tc = ceil(2 m_t^2 G^2); int64 is
exact because G <= 2^26 (all squares and sums stay below 2^55).  Floats are never used in a decision.

Search:
  tiles   the k^2 closed triangles of side 1/k, k = floor(1/m_t) + 1 (diameter 1/k < m_t: <= 1 point per tile);
  combos  N-subsets of tiles, one per D3 orbit (D3 permutes the tile triples (i, j, l));
  node    tighten; pair test (max over the region pair of dS < Tc -> discard); reduction (for region i and region j,
          remove the largest strip at each of the 6 sides of i whose points are within distance < m_t of EVERY point
          of j -- exact binary search); witness test (one grid point per region, all pairs >= m_t exactly -> FEASIBLE);
          otherwise bisect one coordinate of one region.

Usage (examples):
  python tri_engine.py run --n 8 --m 0.3445 --eps 1e-4 --tree trees/n8_pos.jsonl.gz
  python tri_engine.py run --n 16 --m 0.216227269 --eps 1e-4 --sample 500 --cap 20000
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
import math
import os
import random
import sys
import time
from fractions import Fraction

import numpy as np
from numba import njit

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'solver', 'chx'))  # chx_prio: lowers process priority only
import chx_prio  # noqa: E402

PERMS = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]


# ----------------------------------------------------------------------------------------------------------------------
# problem set-up (exact, Python integers / Fractions)
# ----------------------------------------------------------------------------------------------------------------------
def frac_of(x: str) -> Fraction:
    return Fraction(x)


def target(m_star: Fraction, eps: Fraction, digits: int = 15) -> Fraction:
    """m_t = m*(1+eps) rounded AWAY from m* on a 10^-digits grid (up for eps > 0, down for eps < 0)."""
    raw = m_star * (1 + eps)
    scale = 10 ** digits
    if eps >= 0:
        return Fraction(-((-raw.numerator * scale) // raw.denominator), scale)
    return Fraction((raw.numerator * scale) // raw.denominator, scale)


def setup(n: int, mt: Fraction, gbits: int = 24):
    k = math.floor(1 / mt) + 1
    assert Fraction(1, k) < mt
    s = gbits
    G = k << s
    while G > (1 << 26):
        s -= 1
        G = k << s
    g = G // k
    Tc_frac = 2 * mt * mt * G * G
    Tc = -((-Tc_frac.numerator) // Tc_frac.denominator)  # ceil
    assert Tc < (1 << 60)
    tiles = []
    for i in range(k):
        for j in range(k - i):
            for l in (k - 1 - i - j, k - 2 - i - j):
                if l >= 0:
                    tiles.append((i, j, l))
    tiles.sort()
    assert len(tiles) == k * k
    tindex = {t: q for q, t in enumerate(tiles)}
    tperm = []
    for p in PERMS:
        tperm.append([tindex[(t[p[0]], t[p[1]], t[p[2]])] for t in tiles])
    tb = np.array([[t[0] * g, (t[0] + 1) * g, t[1] * g, (t[1] + 1) * g, t[2] * g, (t[2] + 1) * g] for t in tiles],
                  dtype=np.int64)
    return dict(n=n, mt=mt, k=k, G=G, g=g, Tc=Tc, tiles=tiles, tperm=tperm, tb=tb)


def canonical(combo, tperm):
    best = None
    for p in tperm:
        im = tuple(sorted(p[t] for t in combo))
        if best is None or im < best:
            best = im
    return best


def orbit_reps(T, n, tperm):
    for c in itertools.combinations(range(T), n):
        if canonical(c, tperm) == c:
            yield c


# ----------------------------------------------------------------------------------------------------------------------
# numba kernels (int64, exact)
# ----------------------------------------------------------------------------------------------------------------------
@njit(cache=True)
def _maxS(p0, p1, p2, p3, p4, p5, q0, q1, q2, q3, q4, q5):
    la = p0 - q1; ha = p1 - q0; lb = p2 - q3; hb = p3 - q2; lc = p4 - q5; hc = p5 - q4
    best = -1
    for xa in (la, ha):
        for yb in (lb, hb):
            z = -xa - yb
            if lc <= z <= hc:
                v = xa * xa + yb * yb + z * z
                if v > best:
                    best = v
        for zc in (lc, hc):
            y = -xa - zc
            if lb <= y <= hb:
                v = xa * xa + y * y + zc * zc
                if v > best:
                    best = v
    for yb in (lb, hb):
        for zc in (lc, hc):
            x = -yb - zc
            if la <= x <= ha:
                v = x * x + yb * yb + zc * zc
                if v > best:
                    best = v
    return best


@njit(cache=True)
def _far(p, q, Tc):
    """True if every pair of points of the two regions is at distance >= m_t (sound LOWER bound on the minimum)."""
    la = p[0] - q[1]; ha = p[1] - q[0]; lb = p[2] - q[3]; hb = p[3] - q[2]; lc = p[4] - q[5]; hc = p[5] - q[4]
    gmax = 0
    for lo, hi in ((la, ha), (lb, hb), (lc, hc)):
        gg = 0
        if lo > 0:
            gg = lo
        elif hi < 0:
            gg = -hi
        if gg > gmax:
            gmax = gg
    return 3 * gmax * gmax >= 2 * Tc


@njit(cache=True)
def _tighten(r, i, G):
    a0 = r[i, 0]; a1 = r[i, 1]; b0 = r[i, 2]; b1 = r[i, 3]; c0 = r[i, 4]; c1 = r[i, 5]
    na0 = max(a0, G - b1 - c1); na1 = min(a1, G - b0 - c0)
    nb0 = max(b0, G - a1 - c1); nb1 = min(b1, G - a0 - c0)
    nc0 = max(c0, G - a1 - b1); nc1 = min(c1, G - a0 - b0)
    if na0 > na1 or nb0 > nb1 or nc0 > nc1:
        return False
    r[i, 0] = na0; r[i, 1] = na1; r[i, 2] = nb0; r[i, 3] = nb1; r[i, 4] = nc0; r[i, 5] = nc1
    return True


@njit(cache=True)
def _strip_forbidden(r, i, q, low, t, j, G, Tc):
    """Is the closed strip of region i with coordinate q on the low (resp. high) side up to (from) t within distance
    < m_t of every point of region j?  Exact."""
    s = np.empty(6, np.int64)
    for u in range(6):
        s[u] = r[i, u]
    if low:
        s[2 * q + 1] = t
    else:
        s[2 * q] = t
    a0 = s[0]; a1 = s[1]; b0 = s[2]; b1 = s[3]; c0 = s[4]; c1 = s[5]
    na0 = max(a0, G - b1 - c1); na1 = min(a1, G - b0 - c0)
    nb0 = max(b0, G - a1 - c1); nb1 = min(b1, G - a0 - c0)
    nc0 = max(c0, G - a1 - b1); nc1 = min(c1, G - a0 - b0)
    if na0 > na1 or nb0 > nb1 or nc0 > nc1:
        return True  # empty strip: vacuous
    v = _maxS(na0, na1, nb0, nb1, nc0, nc1, r[j, 0], r[j, 1], r[j, 2], r[j, 3], r[j, 4], r[j, 5])
    return v < Tc


@njit(cache=True)
def _push(r, i, q, low, j, G, Tc):
    """Return the new bound (== old if nothing removable), or -1 if the whole region i is forbidden by j."""
    lo = r[i, 2 * q]; hi = r[i, 2 * q + 1]
    if low:
        if not _strip_forbidden(r, i, q, True, lo, j, G, Tc):
            return lo
        if _strip_forbidden(r, i, q, True, hi, j, G, Tc):
            return -1
        a = lo; b = hi  # invariant: strip(a) forbidden, strip(b) not
        while b - a > 1:
            m = (a + b) // 2
            if _strip_forbidden(r, i, q, True, m, j, G, Tc):
                a = m
            else:
                b = m
        return a
    else:
        if not _strip_forbidden(r, i, q, False, hi, j, G, Tc):
            return hi
        if _strip_forbidden(r, i, q, False, lo, j, G, Tc):
            return -1
        a = lo; b = hi  # invariant: strip(b) forbidden, strip(a) not
        while b - a > 1:
            m = (a + b) // 2
            if _strip_forbidden(r, i, q, False, m, j, G, Tc):
                b = m
            else:
                a = m
        return b


@njit(cache=True)
def process_node(r, n, G, Tc, maxrounds, rule, ops, wit):
    """In place on r (n x 6).  Returns (status, a, b, c, nops).
    status 0: discard, pair (a, b) too close everywhere;  1: split region a, coordinate b at c;
    2: feasible witness in wit (n x 3);  3: undecided (grid limit);  4: ops buffer full (treated as split)."""
    nops = 0
    maxops = ops.shape[0]
    for i in range(n):
        if not _tighten(r, i, G):
            return 5, i, 0, 0, nops
    # pair test
    for i in range(n):
        for j in range(i + 1, n):
            v = _maxS(r[i, 0], r[i, 1], r[i, 2], r[i, 3], r[i, 4], r[i, 5],
                      r[j, 0], r[j, 1], r[j, 2], r[j, 3], r[j, 4], r[j, 5])
            if v < Tc:
                return 0, i, j, 0, nops
    # reduction
    rounds = 0
    changed = True
    while changed and rounds < maxrounds:
        changed = False
        rounds += 1
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if _far(r[i], r[j], Tc):
                    continue
                for side in range(6):
                    q = side // 2
                    low = (side % 2) == 0
                    old = r[i, side]
                    t = _push(r, i, q, low, j, G, Tc)
                    if t == -1:
                        # whole region i forbidden by j: equivalent to the pair test on (i, j)
                        if i < j:
                            return 0, i, j, 0, nops
                        return 0, j, i, 0, nops
                    if t != old:
                        if nops >= maxops:
                            return 4, 0, 0, 0, nops
                        ops[nops, 0] = i; ops[nops, 1] = side; ops[nops, 2] = t; ops[nops, 3] = j
                        nops += 1
                        w = r[i, 2 * q + 1] - r[i, 2 * q]
                        r[i, side] = t
                        if not _tighten(r, i, G):
                            return 5, i, 0, 0, nops
                        if abs(t - old) * 32 >= w:
                            changed = True
    # witness: one grid point per region
    ok = True
    for i in range(n):
        A = (r[i, 0] + r[i, 1]) // 2
        blo = max(r[i, 2], G - A - r[i, 5]); bhi = min(r[i, 3], G - A - r[i, 4])
        B = (blo + bhi) // 2
        wit[i, 0] = A; wit[i, 1] = B; wit[i, 2] = G - A - B
    for i in range(n):
        for j in range(i + 1, n):
            da = wit[i, 0] - wit[j, 0]; db = wit[i, 1] - wit[j, 1]; dc = wit[i, 2] - wit[j, 2]
            if da * da + db * db + dc * dc < Tc:
                ok = False
                break
        if not ok:
            break
    if ok:
        return 2, 0, 0, 0, nops
    # split choice
    bi = -1; bq = -1; bw = 1
    if rule == 1:
        # most nearly violated conflicting pair; split the wider of its two regions
        best = -1.0; pi = -1; pj = -1
        for i in range(n):
            for j in range(i + 1, n):
                if _far(r[i], r[j], Tc):
                    continue
                v = _maxS(r[i, 0], r[i, 1], r[i, 2], r[i, 3], r[i, 4], r[i, 5],
                          r[j, 0], r[j, 1], r[j, 2], r[j, 3], r[j, 4], r[j, 5])
                wi = max(r[i, 1] - r[i, 0], r[i, 3] - r[i, 2], r[i, 5] - r[i, 4])
                wj = max(r[j, 1] - r[j, 0], r[j, 3] - r[j, 2], r[j, 5] - r[j, 4])
                if max(wi, wj) < 2:
                    continue
                ratio = (v - Tc) / Tc
                if pi == -1 or ratio < best:
                    best = ratio; pi = i; pj = j
        if pi >= 0:
            for i in (pi, pj):
                for q in range(3):
                    w = r[i, 2 * q + 1] - r[i, 2 * q]
                    if w > bw:
                        bw = w; bi = i; bq = q
    if bi == -1:
        for i in range(n):
            for q in range(3):
                w = r[i, 2 * q + 1] - r[i, 2 * q]
                if w > bw:
                    bw = w; bi = i; bq = q
    if bi == -1:
        return 3, 0, 0, 0, nops
    mid = (r[bi, 2 * bq] + r[bi, 2 * bq + 1]) // 2
    return 1, bi, bq, mid, nops


# ----------------------------------------------------------------------------------------------------------------------
# float witness heuristic (HEURISTIC ONLY: its output is accepted only after the exact integer check)
# ----------------------------------------------------------------------------------------------------------------------
def float_witness(P, combo, starts=3, seed=0):
    """Maximise the min distance with point i inside tile combo[i] (SLSQP), round to the grid, check exactly.
    Returns a list of integer (A, B, C) points or None."""
    from scipy.optimize import minimize
    n = P['n']; k = P['k']; G = P['G']; Tc = P['Tc']
    T = [P['tiles'][t] for t in combo]
    iu = np.triu_indices(n, 1)
    rng = np.random.default_rng(seed)
    lo = np.array([[t[0], t[1], t[2]] for t in T], float) / k
    hi = lo + 1.0 / k

    def cons(z):
        a = z[:n]; b = z[n:2 * n]; c = 1 - a - b; t = z[-1]
        da = a[iu[0]] - a[iu[1]]; db = b[iu[0]] - b[iu[1]]; dc = c[iu[0]] - c[iu[1]]
        pair = 0.5 * (da * da + db * db + dc * dc) - t
        return np.concatenate([pair, a - lo[:, 0], hi[:, 0] - a, b - lo[:, 1], hi[:, 1] - b, c - lo[:, 2], hi[:, 2] - c])

    cen = []
    for t in T:
        up = sum(t) == k - 1
        if up:
            v = np.array([[t[0] + 1, t[1], t[2]], [t[0], t[1] + 1, t[2]], [t[0], t[1], t[2] + 1]], float) / k
        else:
            v = np.array([[t[0] + 1, t[1] + 1, t[2]], [t[0] + 1, t[1], t[2] + 1], [t[0], t[1] + 1, t[2] + 1]], float) / k
        cen.append(v.mean(0))
    cen = np.array(cen)
    for s in range(starts):
        x0 = cen + (rng.random(cen.shape) - 0.5) * (0.3 / k if s else 0.0)
        z0 = np.concatenate([x0[:, 0], x0[:, 1], [0.0]])
        try:
            res = minimize(lambda z: -z[-1], z0, jac=lambda z: -np.eye(2 * n + 1)[-1],
                           constraints=[{'type': 'ineq', 'fun': cons}], method='SLSQP',
                           options={'maxiter': 300, 'ftol': 1e-14})
        except Exception:
            continue
        a = res.x[:n]; b = res.x[n:2 * n]
        pts = []
        for i in range(n):
            A = int(round(a[i] * G)); B = int(round(b[i] * G))
            A = min(max(A, 0), G); B = min(max(B, 0), G - A)
            pts.append((A, B, G - A - B))
        ok = all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2 >= Tc
                 for p, q in itertools.combinations(pts, 2))
        if ok:
            return [list(p) for p in pts]
    return None


# ----------------------------------------------------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------------------------------------------------
class TreeLog:
    def __init__(self, path):
        self.f = gzip.open(path, 'wt', compresslevel=6) if path else None

    def w(self, rec):
        if self.f:
            self.f.write(json.dumps(rec, separators=(',', ':')) + '\n')

    def close(self):
        if self.f:
            self.f.close()


def accepted_by(r, clusters):
    """Index of the first cluster whose box contains every frame region (rattler exempt), else -1."""
    for q, cl in enumerate(clusters):
        c = cl['pts']; rho = cl['rho']; ok = True
        for i in cl['frame_idx']:
            for u in range(3):
                if r[i, 2 * u] < c[i][u] - rho or r[i, 2 * u + 1] > c[i][u] + rho:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return q
    return -1


def solve_combo(P, combo, log: TreeLog, next_id, cap, maxrounds=30, rule=1, fwit=False, clusters=None):
    """Depth-first B&B on one tile combination.  Returns (status, nodes, next_id, witness, maxdepth).
    clusters (optional): boxes around known configurations of THIS combination; a node whose frame regions all lie in
    one of them is ACCEPTED (fate A) instead of being refined."""
    n = P['n']; G = P['G']; Tc = P['Tc']
    n_acc = 0
    nosplit = set(range(n))
    for cl in (clusters or []):
        nosplit -= set(cl['frame_idx'])
    root = P['tb'][list(combo)].copy()
    ops = np.zeros((4096, 4), np.int64)
    wit = np.zeros((n, 3), np.int64)
    rid = next_id; next_id += 1
    log.w({'id': rid, 'root': list(combo)})
    if fwit:
        w = float_witness(P, combo)
        if w is not None:
            log.w({'id': rid, 'ops': [], 'fate': ['F', w]})
            return 'FEASIBLE', 1, next_id, w, 0
    stack = [(root, rid, 0)]
    nodes = 0; maxdepth = 0; undecided = 0
    while stack:
        r, nid, depth = stack.pop()
        nodes += 1
        if depth > maxdepth:
            maxdepth = depth
        if nodes > cap:
            return 'CAP', nodes, next_id, None, maxdepth
        st, a, b, c, nops = process_node(r, n, G, Tc, maxrounds, rule, ops, wit)
        oplist = ops[:nops].tolist()
        if clusters and st in (1, 2, 3, 4):
            q = accepted_by(r, clusters)
            if q >= 0:
                log.w({'id': nid, 'ops': oplist, 'fate': ['A', clusters[q]['id']]})
                n_acc += 1
                continue
        if st == 0:
            log.w({'id': nid, 'ops': oplist, 'fate': ['D', int(a), int(b)]})
        elif st == 5:
            log.w({'id': nid, 'ops': oplist, 'fate': ['E', int(a)]})
        elif st == 2:
            log.w({'id': nid, 'ops': oplist, 'fate': ['F', wit.tolist()]})
            return 'FEASIBLE', nodes, next_id, wit.tolist(), maxdepth
        elif st == 3:
            log.w({'id': nid, 'ops': oplist, 'fate': ['U']})
            undecided += 1
        else:
            if clusters:
                # heuristic only: the cluster this node is heading to (smallest excess of its frame regions over the
                # box); never refine a region that cluster ignores (its rattler)
                bestq = None; bestx = None
                for cl in clusters:
                    c_ = cl['pts']; rho_ = cl['rho']; ex = 0
                    for i in cl['frame_idx']:
                        for u in range(3):
                            ex = max(ex, c_[i][u] - rho_ - r[i, 2 * u], r[i, 2 * u + 1] - c_[i][u] - rho_)
                    if bestx is None or ex < bestx:
                        bestx = ex; bestq = cl
                nosplit = set(range(n)) - set(bestq['frame_idx'])
            if clusters and nosplit and a in nosplit:  # never refine a region the acceptance test ignores (rattler)
                bw = 1; a2 = -1
                for i in range(n):
                    if i in nosplit:
                        continue
                    for q in range(3):
                        w = r[i, 2 * q + 1] - r[i, 2 * q]
                        if w > bw:
                            bw = w; a2 = i; b2 = q
                if a2 >= 0:
                    a = a2; b = b2; c = (r[a, 2 * b] + r[a, 2 * b + 1]) // 2
            if st == 4:  # ops buffer full: split the widest region
                bw = 1; a = -1
                for i in range(n):
                    for q in range(3):
                        w = r[i, 2 * q + 1] - r[i, 2 * q]
                        if w > bw:
                            bw = w; a = i; b = q
                if a == -1:
                    log.w({'id': nid, 'ops': oplist, 'fate': ['U']})
                    return 'UNDECIDED', nodes, next_id, None, maxdepth
                c = (r[a, 2 * b] + r[a, 2 * b + 1]) // 2
            lo_id = next_id; hi_id = next_id + 1; next_id += 2
            log.w({'id': nid, 'ops': oplist, 'fate': ['S', int(a), int(b), int(c), lo_id, hi_id]})
            rl = r.copy(); rl[a, 2 * b + 1] = c
            rh = r.copy(); rh[a, 2 * b] = c
            stack.append((rh, hi_id, depth + 1))
            stack.append((rl, lo_id, depth + 1))
    if undecided:
        return 'UNDECIDED', nodes, next_id, None, maxdepth
    if n_acc:
        return 'PROVED_MOD_CLUSTERS', nodes, next_id, n_acc, maxdepth
    return 'PROVED', nodes, next_id, None, maxdepth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['run'])
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--m', required=True, help='best-known m* (decimal string)')
    ap.add_argument('--eps', default='1e-4', help='relative offset; negative = control below m*')
    ap.add_argument('--tree', default='')
    ap.add_argument('--summary', default='')
    ap.add_argument('--cap', type=int, default=10 ** 7, help='node cap per combination')
    ap.add_argument('--sample', type=int, default=0, help='random sample of orbit reps (projection mode)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--part', default='0/1', help='i/p: process orbit reps with index %% p == i')
    ap.add_argument('--rule', type=int, default=0)
    ap.add_argument('--rounds', type=int, default=30)
    ap.add_argument('--gbits', type=int, default=24)
    ap.add_argument('--combos', default='', help='explicit combos, e.g. "0,1,2;3,4,5"')
    ap.add_argument('--timeout', type=float, default=1e9, help='wall seconds for the whole run')
    ap.add_argument('--fwit', action='store_true', help='float witness heuristic at each combination root')
    ap.add_argument('--accept', default='', help='JSON file of cluster boxes (acceptance mode)')
    a = ap.parse_args()
    prio = chx_prio.lower()
    m_star = Fraction(a.m)
    eps = Fraction(a.eps)
    mt = target(m_star, eps)
    P = setup(a.n, mt, a.gbits)
    T = P['k'] ** 2
    pi, pp = (int(x) for x in a.part.split('/'))
    t0 = time.time()
    summ = dict(n=a.n, m_star=a.m, eps=a.eps, m_t=str(mt), k=P['k'], G=P['G'], Tc=P['Tc'], tiles=T,
                rule=a.rule, rounds=a.rounds, cap=a.cap, part=a.part, priority=prio)
    log = TreeLog(a.tree)
    acc = json.load(open(a.accept)) if a.accept else None
    hdr = {'header': {k_: (v if not isinstance(v, Fraction) else str(v)) for k_, v in summ.items()}, 'tiles': P['tiles']}
    if acc:
        assert acc['G'] == P['G'], 'cluster file grid differs'
        hdr['accept'] = acc['clusters']
    log.w(hdr)
    by_combo = {}
    if acc:
        for cl in acc['clusters']:
            by_combo.setdefault(tuple(cl['combo']), []).append(
                {'id': cl['id'], 'pts': cl['points'], 'rho': cl['rho'],
                 'frame_idx': [i for i in range(a.n) if cl['frame'][i]]})
    if T < a.n:
        summ.update(status='PROVED', reason='pigeonhole', combos=0, nodes=0, seconds=0.0)
        log.w({'pigeonhole': [T, a.n]})
        log.close()
        print(json.dumps(summ))
        if a.summary:
            json.dump(summ, open(a.summary, 'w'), indent=1)
        return
    if a.combos:
        reps = [tuple(int(x) for x in s.split(',')) for s in a.combos.split(';')]
    else:
        reps = list(orbit_reps(T, a.n, P['tperm']))
    summ['raw_combos'] = math.comb(T, a.n)
    summ['orbit_reps'] = len(reps)
    if a.sample:
        rnd = random.Random(a.seed)
        reps = rnd.sample(reps, min(a.sample, len(reps)))
    reps = [c for q, c in enumerate(reps) if q % pp == pi]
    status = 'PROVED'
    next_id = 0
    per = []
    tot_nodes = 0
    witness = None
    for q, combo in enumerate(reps):
        tc = time.time()
        st, nodes, next_id, wit, md = solve_combo(P, combo, log, next_id, a.cap, a.rounds, a.rule, a.fwit,
                                                  by_combo.get(tuple(combo)))
        if st == 'PROVED_MOD_CLUSTERS':
            summ['accepted_leaves'] = summ.get('accepted_leaves', 0) + wit
            wit = None
            status = 'PROVED_MOD_CLUSTERS' if status == 'PROVED' else status
        dt = time.time() - tc
        tot_nodes += nodes
        per.append([list(combo), st, nodes, round(dt, 4), md])
        if st == 'FEASIBLE':
            status = 'FEASIBLE'; witness = {'combo': list(combo), 'grid_points': wit}
            break
        if st in ('UNDECIDED', 'CAP') and not a.sample:
            status = 'UNDECIDED'
        if time.time() - t0 > a.timeout:
            status = 'TIMEOUT'
            break
    summ.update(status=status, combos_done=len(per), nodes=tot_nodes, seconds=round(time.time() - t0, 3),
                witness=witness)
    if a.sample:
        arr = np.array([p[2] for p in per], float)
        tarr = np.array([p[3] for p in per], float)
        summ.update(sample_mean_nodes=float(arr.mean()), sample_max_nodes=float(arr.max()),
                    sample_mean_sec=float(tarr.mean()), sample_capped=sum(1 for p in per if p[1] == 'CAP'),
                    sample_proved=sum(1 for p in per if p[1] == 'PROVED'))
    summ['per_combo'] = per
    log.w({'end': status, 'nodes': tot_nodes})
    log.close()
    short = {k_: v for k_, v in summ.items() if k_ != 'per_combo'}
    print(json.dumps(short))
    if a.summary:
        json.dump(summ, open(a.summary, 'w'), indent=1)


if __name__ == '__main__':
    main()
