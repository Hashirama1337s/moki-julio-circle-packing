#!/usr/bin/env python3
"""Exact infeasibility proof: N equal circles in a rectangle of width 1 and height h.

Point form. Centers lie in [0, W] x [0, H], W = 1 - 2 r_t, H = h - 2 r_t,
r_t = r_known * (1 + eps), and every pair must have distance >= d = 2 r_t.
PROVED means that placement is impossible, so any feasible packing at r_known
is optimal within relative eps. This file does not re-check that r_known fits.

Search: one axis-aligned box per point, integer grid, decisions by integer
cross-multiplication of fractions.Fraction (no float accept/reject/split).
A node dies only if some pair has maximum distance strictly less than d,
symmetry bounds empty it, or an active-area cut removes the whole box.
Active area (Markot-Csendes): delete points whose distance to every point of
another box is < d. Boundary strips are cut exactly; a large interior forbidden
rectangle is branched out.

Symmetry, complete for identical points in a rectangle: x1 <= ... <= xN,
mirror so x1 + xN <= W, mirror so y1 <= H/2.

Grid cells are an outer cover. Hitting a cell-sized survivor is UNDECIDED,
not a proof. Packomania crc_800 N=4, crc_700 N=5, crc_600 N=6, crc_500 N=7,
crc_400 N=9, crc_300 N=12 are not certified here.

Usage:
  py -3.11 prove_small.py <h> <N> <r_known> <eps>
  py -3.11 prove_small.py --self-test
"""

from __future__ import annotations

import sys
from fractions import Fraction
from itertools import permutations
from math import gcd
from time import perf_counter_ns

MAX_NODES = 2_000_000
WAY_CAP = 20_000
G_CAP_BIT = 20


def format_seconds(ns: int) -> str:
    ms = ns // 1_000_000
    return f"{ms // 1000}.{ms % 1000:03d}"


def far(a0: int, a1: int, b0: int, b1: int) -> int:
    left = a1 - b0
    right = b1 - a0
    return left if left >= right else right


def gap(a0: int, a1: int, b0: int, b1: int) -> int:
    if a1 < b0:
        return b0 - a1
    if b1 < a0:
        return a0 - b1
    return 0


def parse_exact(text: str) -> Fraction:
    return Fraction(text.strip())


def choose_g(eps: Fraction) -> int:
    g = 1 << 12
    if eps <= 0:
        return g
    target = 1024
    while g < (1 << G_CAP_BIT) and (g * g) * eps < target:
        g <<= 1
    return g


class Engine:
    __slots__ = ("A", "B", "CG", "G", "n", "Wn", "Wd", "Hn", "Hd", "W", "H", "d")

    def __init__(self, W: Fraction, H: Fraction, d: Fraction, n: int, eps: Fraction) -> None:
        G = choose_g(eps)
        Wn, Wd = W.numerator, W.denominator
        Hn, Hd = H.numerator, H.denominator
        dn, dd = d.numerator, d.denominator
        A = Wn * Wn * Hd * Hd * dd * dd
        B = Hn * Hn * Wd * Wd * dd * dd
        C = dn * dn * Wd * Wd * Hd * Hd
        g = gcd(gcd(A, B), C)
        self.A = A // g
        self.B = B // g
        self.G = G
        self.CG = (C // g) * G * G
        self.n = n
        self.Wn = Wn
        self.Wd = Wd
        self.Hn = Hn
        self.Hd = Hd
        self.W = W
        self.H = H
        self.d = d

    def dist2_lt(self, dx: int, dy: int) -> bool:
        return self.A * dx * dx + self.B * dy * dy < self.CG

    def too_close(self, a: tuple, b: tuple) -> bool:
        return self.dist2_lt(
            far(a[0], a[1], b[0], b[1]),
            far(a[2], a[3], b[2], b[3]),
        )

    def always_far(self, a: tuple, b: tuple) -> bool:
        dx = gap(a[0], a[1], b[0], b[1])
        dy = gap(a[2], a[3], b[2], b[3])
        return self.A * dx * dx + self.B * dy * dy >= self.CG

    def can_forbid(self, other: tuple) -> bool:
        dx = other[1] - other[0]
        dy = other[3] - other[2]
        return self.A * dx * dx + self.B * dy * dy < (self.CG << 2)

    def point_forbidden(self, x: int, y: int, other: tuple) -> bool:
        dx = abs(x - other[0])
        alt = abs(x - other[1])
        if alt > dx:
            dx = alt
        dy = abs(y - other[2])
        alt = abs(y - other[3])
        if alt > dy:
            dy = alt
        return self.dist2_lt(dx, dy)

    def rect_forbidden(self, x0: int, x1: int, y0: int, y1: int, other: tuple) -> bool:
        return (
            self.point_forbidden(x0, y0, other)
            and self.point_forbidden(x0, y1, other)
            and self.point_forbidden(x1, y0, other)
            and self.point_forbidden(x1, y1, other)
        )

    def push_low_x(self, x0: int, x1: int, y0: int, y1: int, other: tuple) -> int:
        if not self.rect_forbidden(x0, x0, y0, y1, other):
            return x0
        lo, hi = x0, x1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if self.rect_forbidden(x0, mid, y0, y1, other):
                lo = mid
            else:
                hi = mid - 1
        return lo

    def push_high_x(self, x0: int, x1: int, y0: int, y1: int, other: tuple) -> int:
        if not self.rect_forbidden(x1, x1, y0, y1, other):
            return x1
        lo, hi = x0, x1
        while lo < hi:
            mid = (lo + hi) // 2
            if self.rect_forbidden(mid, x1, y0, y1, other):
                hi = mid
            else:
                lo = mid + 1
        return hi

    def push_low_y(self, x0: int, x1: int, y0: int, y1: int, other: tuple) -> int:
        if not self.rect_forbidden(x0, x1, y0, y0, other):
            return y0
        lo, hi = y0, y1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if self.rect_forbidden(x0, x1, y0, mid, other):
                lo = mid
            else:
                hi = mid - 1
        return lo

    def push_high_y(self, x0: int, x1: int, y0: int, y1: int, other: tuple) -> int:
        if not self.rect_forbidden(x0, x1, y1, y1, other):
            return y1
        lo, hi = y0, y1
        while lo < hi:
            mid = (lo + hi) // 2
            if self.rect_forbidden(x0, x1, mid, y1, other):
                hi = mid
            else:
                lo = mid + 1
        return hi

    def shrink_box(self, box: tuple, other: tuple):
        x0, x1, y0, y1 = box
        x0 = self.push_low_x(x0, x1, y0, y1, other)
        if x0 > x1:
            return None
        x1 = self.push_high_x(x0, x1, y0, y1, other)
        if x0 > x1:
            return None
        y0 = self.push_low_y(x0, x1, y0, y1, other)
        if y0 > y1:
            return None
        y1 = self.push_high_y(x0, x1, y0, y1, other)
        if y0 > y1:
            return None
        if self.rect_forbidden(x0, x1, y0, y1, other):
            return None
        return (x0, x1, y0, y1)

    def apply_symmetry(self, boxes: list) -> bool:
        n = len(boxes)
        G = self.G
        half = G // 2
        changed = True
        while changed:
            changed = False
            for i in range(n - 1):
                if boxes[i][1] > boxes[i + 1][1]:
                    boxes[i][1] = boxes[i + 1][1]
                    changed = True
                if boxes[i + 1][0] < boxes[i][0]:
                    boxes[i + 1][0] = boxes[i][0]
                    changed = True
                if boxes[i][0] > boxes[i][1] or boxes[i][2] > boxes[i][3]:
                    return False
            if boxes[-1][0] > boxes[-1][1] or boxes[-1][2] > boxes[-1][3]:
                return False
            cap = G - boxes[-1][0]
            if boxes[0][1] > cap:
                boxes[0][1] = cap
                changed = True
            capn = G - boxes[0][0]
            if boxes[-1][1] > capn:
                boxes[-1][1] = capn
                changed = True
            if boxes[0][0] > boxes[0][1] or boxes[-1][0] > boxes[-1][1]:
                return False
            if boxes[0][0] + boxes[-1][0] > G:
                return False
            if boxes[0][3] > half:
                boxes[0][3] = half
                changed = True
            if boxes[0][2] > boxes[0][3]:
                return False
        return True

    def tighten(self, state: tuple):
        boxes = [list(b) for b in state]
        n = self.n
        for _ in range(64):
            if not self.apply_symmetry(boxes):
                return None
            for i in range(n):
                bi = boxes[i]
                if bi[0] > bi[1] or bi[2] > bi[3]:
                    return None
                for j in range(i + 1, n):
                    if self.too_close(bi, boxes[j]):
                        return None
            changed = False
            for i in range(n):
                for j in range(n):
                    if i == j or not self.can_forbid(boxes[j]):
                        continue
                    nb = self.shrink_box(tuple(boxes[i]), tuple(boxes[j]))
                    if nb is None:
                        return None
                    if nb != tuple(boxes[i]):
                        boxes[i][:] = nb
                        changed = True
            if not changed:
                break
        return tuple(tuple(b) for b in boxes)

    def all_separated(self, state: tuple) -> bool:
        n = self.n
        for i in range(n):
            for j in range(i + 1, n):
                if not self.always_far(state[i], state[j]):
                    return False
        return True

    def max_side(self, state: tuple) -> int:
        m = 0
        for b in state:
            sx = b[1] - b[0]
            sy = b[3] - b[2]
            if sx > m:
                m = sx
            if sy > m:
                m = sy
        return m

    def x_impossible(self, state: tuple) -> bool:
        for i in range(self.n - 1):
            if state[i][0] > state[i + 1][1]:
                return True
        return False

    def axis_bounds(self, npart: int) -> list:
        G = self.G
        out = []
        for i in range(npart):
            a = (i * G) // npart
            b = G if i == npart - 1 else ((i + 1) * G + npart - 1) // npart
            if b < a:
                b = a
            out.append((a, b))
        return out

    def make_cells(self, nx: int, ny: int) -> list:
        xs = self.axis_bounds(nx)
        ys = self.axis_bounds(ny)
        cells = []
        for x0, x1 in xs:
            for y0, y1 in ys:
                cells.append((x0, x1, y0, y1))
        return cells

    def perm_count(self, m: int, n: int):
        ways = 1
        for k in range(m - n + 1, m + 1):
            ways *= k
            if ways > WAY_CAP:
                return None
        return ways

    def build_starts(self):
        """Return (kind, states, checked).

        kind 'pigeon' : fewer small-diameter cells than points.
        kind 'assign' : one point per chosen cell, filtered.
        kind 'full'   : one node, every point in the whole box.
        """
        n = self.n
        G = self.G
        limit = n + 1
        best = None
        for nx in range(1, limit + 1):
            for ny in range(1, limit + 1):
                cells = self.make_cells(nx, ny)
                if any(not self.too_close(c, c) for c in cells):
                    continue
                m = len(cells)
                if m < n:
                    return "pigeon", [], 1
                ways = self.perm_count(m, n)
                if ways is None:
                    continue
                if best is None or ways < best[0]:
                    best = (ways, cells)
        if best is None:
            full = tuple((0, G, 0, G) for _ in range(n))
            return "full", [full], 0
        cells = best[1]
        stack = []
        checked = 0
        for assign in permutations(range(len(cells)), n):
            checked += 1
            state = tuple(cells[i] for i in assign)
            if self.x_impossible(state):
                continue
            dead = False
            for i in range(n):
                for j in range(i + 1, n):
                    if self.too_close(state[i], state[j]):
                        dead = True
                        break
                if dead:
                    break
            if not dead:
                stack.append(state)
        return "assign", stack, checked

    def deepest_seed(self, box: tuple, other: tuple):
        x0, x1, y0, y1 = box
        pts = (
            (x0, y0),
            (x0, y1),
            (x1, y0),
            (x1, y1),
            ((x0 + x1) // 2, (y0 + y1) // 2),
            ((x0 + x1) // 2, y0),
            ((x0 + x1) // 2, y1),
            (x0, (y0 + y1) // 2),
            (x1, (y0 + y1) // 2),
            (
                min(max((other[0] + other[1]) // 2, x0), x1),
                min(max((other[2] + other[3]) // 2, y0), y1),
            ),
        )
        best = None
        best_score = 0
        for x, y in pts:
            if not self.point_forbidden(x, y, other):
                continue
            dx = abs(x - other[0])
            alt = abs(x - other[1])
            if alt > dx:
                dx = alt
            dy = abs(y - other[2])
            alt = abs(y - other[3])
            if alt > dy:
                dy = alt
            score = self.A * dx * dx + self.B * dy * dy
            if best is None or score < best_score:
                best = (x, y)
                best_score = score
        return best

    def expand_high_x(self, x0, x1, y0, y1, limit, other) -> int:
        if x1 >= limit or not self.rect_forbidden(x0, x1, y0, y1, other):
            return x1
        lo, hi = x1, limit
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if self.rect_forbidden(x0, mid, y0, y1, other):
                lo = mid
            else:
                hi = mid - 1
        return lo

    def expand_low_x(self, x0, x1, y0, y1, limit, other) -> int:
        if x0 <= limit or not self.rect_forbidden(x0, x1, y0, y1, other):
            return x0
        lo, hi = limit, x0
        while lo < hi:
            mid = (lo + hi) // 2
            if self.rect_forbidden(mid, x1, y0, y1, other):
                hi = mid
            else:
                lo = mid + 1
        return hi

    def expand_high_y(self, x0, x1, y0, y1, limit, other) -> int:
        if y1 >= limit or not self.rect_forbidden(x0, x1, y0, y1, other):
            return y1
        lo, hi = y1, limit
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if self.rect_forbidden(x0, x1, y0, mid, other):
                lo = mid
            else:
                hi = mid - 1
        return lo

    def expand_low_y(self, x0, x1, y0, y1, limit, other) -> int:
        if y0 <= limit or not self.rect_forbidden(x0, x1, y0, y1, other):
            return y0
        lo, hi = limit, y0
        while lo < hi:
            mid = (lo + hi) // 2
            if self.rect_forbidden(x0, x1, mid, y1, other):
                hi = mid
            else:
                lo = mid + 1
        return hi

    def grow_hole(self, box: tuple, other: tuple):
        seed = self.deepest_seed(box, other)
        if seed is None:
            return None
        x0 = x1 = seed[0]
        y0 = y1 = seed[1]
        for _ in range(3):
            x1 = self.expand_high_x(x0, x1, y0, y1, box[1], other)
            x0 = self.expand_low_x(x0, x1, y0, y1, box[0], other)
            y1 = self.expand_high_y(x0, x1, y0, y1, box[3], other)
            y0 = self.expand_low_y(x0, x1, y0, y1, box[2], other)
        if x1 <= x0 or y1 <= y0:
            return None
        if not self.rect_forbidden(x0, x1, y0, y1, other):
            return None
        return (x0, x1, y0, y1)

    def best_hole(self, state: tuple):
        best = None
        best_cut = 0
        n = self.n
        for i in range(n):
            bi = state[i]
            bw = bi[1] - bi[0]
            bh = bi[3] - bi[2]
            area = bw * bh
            if area <= 0:
                continue
            for j in range(n):
                if i == j or not self.can_forbid(state[j]):
                    continue
                hole = self.grow_hole(bi, state[j])
                if hole is None:
                    continue
                if hole[0] <= bi[0] and hole[1] >= bi[1] and hole[2] <= bi[2] and hole[3] >= bi[3]:
                    return i, hole, True
                cut = (hole[1] - hole[0]) * (hole[3] - hole[2])
                if cut * 8 >= area and cut > best_cut:
                    best_cut = cut
                    best = (i, hole, False)
        return best

    def hole_children(self, state: tuple, i: int, hole: tuple) -> list:
        x0, x1, y0, y1 = state[i]
        hx0, hx1, hy0, hy1 = hole
        parent = (x1 - x0) + (y1 - y0)
        parts = []
        if x0 < hx0:
            parts.append((x0, hx0, y0, y1))
        if hx1 < x1:
            parts.append((hx1, x1, y0, y1))
        if hx0 < hx1 and y0 < hy0:
            parts.append((hx0, hx1, y0, hy0))
        if hx0 < hx1 and hy1 < y1:
            parts.append((hx0, hx1, hy1, y1))
        out = []
        for p in parts:
            bulk = (p[1] - p[0]) + (p[3] - p[2])
            if bulk < parent:
                out.append(state[:i] + (p,) + state[i + 1 :])
        return out

    def physical_x_ge_y(self, b: tuple) -> bool:
        return (b[1] - b[0]) * self.Wn * self.Hd >= (b[3] - b[2]) * self.Hn * self.Wd

    def choose_bisect(self, state: tuple):
        n = self.n
        best = None
        best_slack = 0
        for i in range(n):
            for j in range(i + 1, n):
                dx = far(state[i][0], state[i][1], state[j][0], state[j][1])
                dy = far(state[i][2], state[i][3], state[j][2], state[j][3])
                slack = self.A * dx * dx + self.B * dy * dy
                if best is None or slack < best_slack:
                    best = (i, j, dx, dy)
                    best_slack = slack
        if best is None:
            i = 0
            axis = 0 if self.physical_x_ge_y(state[0]) else 1
            return i, axis
        i, j, dx, dy = best
        want_x = self.A * dx * dx >= self.B * dy * dy
        pick = i if (state[i][1] - state[i][0]) >= (state[j][1] - state[j][0]) else j
        if not want_x:
            pick = i if (state[i][3] - state[i][2]) >= (state[j][3] - state[j][2]) else j
        axis = 0 if want_x else 1
        b = state[pick]
        length = (b[1] - b[0]) if axis == 0 else (b[3] - b[2])
        if length < 2:
            # Fall back to the globally longest splittable side.
            pick = 0
            axis = 0
            best_len = -1
            for k in range(n):
                for ax, length_k in (
                    (0, state[k][1] - state[k][0]),
                    (1, state[k][3] - state[k][2]),
                ):
                    if length_k > best_len:
                        best_len = length_k
                        pick = k
                        axis = ax
        return pick, axis

    def bisect(self, state: tuple):
        i, axis = self.choose_bisect(state)
        b = state[i]
        if axis == 0:
            lo, hi = b[0], b[1]
            if hi - lo < 2:
                return []
            mid = (lo + hi) // 2
            c1 = (lo, mid, b[2], b[3])
            c2 = (mid, hi, b[2], b[3])
        else:
            lo, hi = b[2], b[3]
            if hi - lo < 2:
                return []
            mid = (lo + hi) // 2
            c1 = (b[0], b[1], lo, mid)
            c2 = (b[0], b[1], mid, hi)
        return [state[:i] + (c1,) + state[i + 1 :], state[:i] + (c2,) + state[i + 1 :]]

    def branch(self, state: tuple):
        found = self.best_hole(state)
        if found is not None:
            i, hole, covers = found
            if covers:
                return None
            kids = self.hole_children(state, i, hole)
            if kids:
                return kids
        kids = self.bisect(state)
        return kids if kids else None

    def score(self, state: tuple) -> int:
        total = 0
        n = self.n
        for i in range(n):
            for j in range(i + 1, n):
                dx = far(state[i][0], state[i][1], state[j][0], state[j][1])
                dy = far(state[i][2], state[i][3], state[j][2], state[j][3])
                total += self.A * dx * dx + self.B * dy * dy
        return total

    def real_boxes(self, state: tuple):
        G = self.G
        W = self.W
        H = self.H
        out = []
        for b in state:
            out.append(
                (
                    Fraction(b[0], G) * W,
                    Fraction(b[1], G) * W,
                    Fraction(b[2], G) * H,
                    Fraction(b[3], G) * H,
                )
            )
        return out

    def search(self, max_nodes: int):
        kind, stack, nodes = self.build_starts()
        if kind == "pigeon" or (kind == "assign" and not stack):
            return "PROVED", nodes, None, "tiling"
        while stack:
            state = stack.pop()
            nodes += 1
            if nodes > max_nodes:
                return "UNDECIDED", nodes, self.real_boxes(state), "node_limit"
            state = self.tighten(state)
            if state is None:
                continue
            if self.n >= 2 and self.all_separated(state):
                return "UNDECIDED", nodes, self.real_boxes(state), "separated"
            if self.max_side(state) <= 1:
                return "UNDECIDED", nodes, self.real_boxes(state), "tolerance"
            children = self.branch(state)
            if not children:
                continue
            scored = [(self.score(ch), idx, ch) for idx, ch in enumerate(children)]
            scored.sort()
            for _, _, ch in scored:
                stack.append(ch)
        return "PROVED", nodes, None, "search"


def prove(h, n: int, r_known, eps, max_nodes: int = MAX_NODES):
    """Return status, nodes, elapsed_ns, boxes_or_None, reason."""
    t0 = perf_counter_ns()
    h = Fraction(h)
    r_known = Fraction(r_known)
    eps = Fraction(eps)
    r_t = r_known * (1 + eps)
    W = 1 - 2 * r_t
    H = h - 2 * r_t
    d = 2 * r_t

    def done(status, nodes, boxes, reason):
        return status, nodes, perf_counter_ns() - t0, boxes, reason

    if n < 1 or h <= 0 or r_known <= 0 or eps < 0:
        raise ValueError("need h > 0, N >= 1, r_known > 0, eps >= 0")
    if W < 0 or H < 0:
        return done("PROVED", 1, None, "empty")
    if n == 1:
        return done("UNDECIDED", 1, [(Fraction(0), W, Fraction(0), H)], "single")
    if d <= 0:
        return done("UNDECIDED", 1, [(Fraction(0), W, Fraction(0), H)], "nonpositive_d")
    if W == 0 and H == 0:
        return done("PROVED", 1, None, "segment")
    if W == 0:
        if (n - 1) * d > H:
            return done("PROVED", 1, None, "segment")
        return done("UNDECIDED", 1, [(Fraction(0), W, Fraction(0), H)], "segment_fits")
    if H == 0:
        if (n - 1) * d > W:
            return done("PROVED", 1, None, "segment")
        return done("UNDECIDED", 1, [(Fraction(0), W, Fraction(0), H)], "segment_fits")
    if W * W + H * H < d * d:
        return done("PROVED", 1, None, "diagonal")

    engine = Engine(W, H, d, n, eps)
    status, nodes, boxes, reason = engine.search(max_nodes)
    return done(status, nodes, boxes, reason)


def n2_infeasible(h, r, eps) -> bool:
    r_t = r * (1 + eps)
    W = 1 - 2 * r_t
    H = h - 2 * r_t
    d = 2 * r_t
    if W < 0 or H < 0:
        return True
    return W * W + H * H < d * d


def self_test() -> int:
    ok = True

    def expect(h, n, r, eps, want, max_nodes=MAX_NODES):
        nonlocal ok
        status, nodes, ns, _boxes, reason = prove(h, n, r, eps, max_nodes=max_nodes)
        sec = format_seconds(ns)
        if status != want:
            ok = False
            print(
                f"FAIL want={want} got={status} N={n} h={h} r={r} eps={eps} "
                f"nodes={nodes} seconds={sec} reason={reason}",
                file=sys.stderr,
            )
            return
        print(f"PASS {status} N={n} h={h} r={r} eps={eps} nodes={nodes} seconds={sec} reason={reason}")

    # N=2: infeasible iff the center-box diagonal is shorter than d.
    n2_cases = [
        (1, Fraction(3, 10), Fraction(0)),
        (1, Fraction(1, 5), Fraction(1, 100)),
        (Fraction(4, 5), Fraction(1, 5), Fraction(1, 50)),
        (2, Fraction(1, 3), Fraction(1, 1000)),
        (Fraction(1, 2), Fraction(1, 5), Fraction(1, 100)),
        (1, Fraction(1, 4), Fraction(1, 1000)),
    ]
    for h, r, eps in n2_cases:
        expect(h, 2, r, eps, "PROVED" if n2_infeasible(h, r, eps) else "UNDECIDED", max_nodes=200_000)

    # Five circles of radius 1/4 cannot fit in the unit square: four quadrants,
    # each of diameter < d.
    expect(1, 5, Fraction(1, 4), Fraction(1, 1_000_000), "PROVED", max_nodes=10_000)

    # Four circles, r = 1/4. eps=1/10 dies at the quadrant assignment.
    # eps=1/100 and 1e-6 need active-area shrinkage plus branching.
    expect(1, 4, Fraction(1, 4), Fraction(1, 10), "PROVED", max_nodes=10_000)
    expect(1, 4, Fraction(1, 4), Fraction(1, 100), "PROVED", max_nodes=500_000)
    expect(1, 4, Fraction(1, 4), Fraction(1, 1_000_000), "PROVED", max_nodes=MAX_NODES)

    # Soundness: a single circle, and four circles known to fit, must not be proved impossible.
    expect(1, 1, Fraction(1, 4), Fraction(1, 100), "UNDECIDED", max_nodes=10)
    expect(1, 4, Fraction(1, 5), Fraction(1, 100), "UNDECIDED", max_nodes=100_000)

    if ok:
        print("SELF-TEST PASSED")
        return 0
    print("SELF-TEST FAILED", file=sys.stderr)
    return 1


def main(argv: list) -> int:
    if len(argv) == 1 or (len(argv) == 2 and argv[1] in ("--self-test", "self-test")):
        return self_test()
    if len(argv) != 5:
        print("usage: py -3.11 prove_small.py <h> <N> <r_known> <eps>", file=sys.stderr)
        print("       py -3.11 prove_small.py --self-test", file=sys.stderr)
        return 2
    try:
        h = parse_exact(argv[1])
        n = int(argv[2])
        r_known = parse_exact(argv[3])
        eps = parse_exact(argv[4])
        status, nodes, ns, boxes, reason = prove(h, n, r_known, eps)
    except (ValueError, ZeroDivisionError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    sec = format_seconds(ns)
    if status == "PROVED":
        print(f"PROVED nodes={nodes} seconds={sec}")
        return 0
    print(f"UNDECIDED nodes={nodes} seconds={sec} reason={reason}")
    if boxes:
        for i, b in enumerate(boxes):
            print(f"point {i} x=[{b[0]},{b[1]}] y=[{b[2]},{b[3]}]")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
