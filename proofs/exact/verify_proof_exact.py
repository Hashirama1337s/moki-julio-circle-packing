#!/usr/bin/env python3
"""Independent exact replay checker for the two circle-packing optimality proofs.

Usage:
    python verify_proof_exact.py <cert.json> <tree.json.gz>

Standard library only. Every decision is rational or in Q(sqrt(2)).
Exit status 0 if and only if the last line is PASS.
"""

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import math
import re
import sys
from fractions import Fraction


class Fail(Exception):
    pass


class Q:
    """a + b*sqrt(2), a and b rational. Sign is exact: sqrt(2) is irrational."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    def __repr__(self):
        return f"({self.a}+{self.b}*sqrt2)"

    def __eq__(self, other):
        return isinstance(other, Q) and self.a == other.a and self.b == other.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __neg__(self):
        return Q(-self.a, -self.b)

    def __add__(self, other):
        if not isinstance(other, Q):
            other = Q(other, 0)
        return Q(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __sub__(self, other):
        if not isinstance(other, Q):
            other = Q(other, 0)
        return Q(self.a - other.a, self.b - other.b)

    def __rsub__(self, other):
        return Q(other, 0) - self

    def __mul__(self, other):
        if isinstance(other, Q):
            return Q(
                self.a * other.a + 2 * self.b * other.b,
                self.a * other.b + self.b * other.a,
            )
        if type(other) is int or isinstance(other, Fraction):
            return Q(self.a * other, self.b * other)
        return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def sign(self):
        a, b = self.a, self.b
        if a == 0 and b == 0:
            return 0
        if a >= 0 and b >= 0:
            return 1
        if a <= 0 and b <= 0:
            return -1
        aa = a * a
        bb = 2 * b * b
        if aa > bb:
            return 1 if a > 0 else -1
        if aa < bb:
            return 1 if b > 0 else -1
        raise Fail("sign test hit a rational square root of 2")

    def abs(self):
        return -self if self.sign() < 0 else Q(self.a, self.b)


def need(obj, key):
    if not isinstance(obj, dict) or key not in obj:
        raise Fail(f"missing {key}")
    return obj[key]


def require_int(x, what):
    if type(x) is not int:
        raise Fail(f"{what} is not an int")
    return x


def parse_fraction(x, what):
    if type(x) is int:
        return Fraction(x)
    if type(x) is str:
        try:
            return Fraction(x.strip())
        except (ValueError, ZeroDivisionError):
            raise Fail(f"bad rational {what}")
    raise Fail(f"bad rational {what}")


def parse_q(pair, what):
    if not isinstance(pair, list) or len(pair) != 2:
        raise Fail(f"{what} is not a Q(sqrt(2)) pair")
    return Q(parse_fraction(pair[0], what), parse_fraction(pair[1], what))


def sqrt2_bounds(digits):
    """My own enclosure of sqrt(2), from the integer square root. No floats."""
    scale = 10 ** digits
    root = math.isqrt(2 * scale * scale)
    lo = Fraction(root, scale)
    hi = Fraction(root + 1, scale)
    if not (lo * lo <= 2 < hi * hi):
        raise Fail("internal sqrt(2) enclosure failed")
    return lo, hi


def rat_upper_abs(z, lo, hi):
    """Rational upper bound of |z| from a verified sqrt(2) enclosure."""
    if z.sign() < 0:
        z = -z
    if z.sign() == 0:
        return Fraction(0)
    if z.b >= 0:
        return z.a + z.b * hi
    return z.a + z.b * lo


def all_constraints(n):
    out = []
    for i in range(n):
        for side in "LRBT":
            out.append(("wall", i, side))
    for i in range(n):
        for j in range(i + 1, n):
            out.append(("pair", i, j))
    return out


def parse_constraint(entry, n):
    if not isinstance(entry, dict) or "kind" not in entry:
        raise Fail("bad constraint")
    kind = entry["kind"]
    if kind == "wall":
        i = require_int(entry.get("i"), "wall index")
        side = entry.get("side")
        if type(side) is not str or side not in ("L", "R", "B", "T"):
            raise Fail("bad wall")
        if not 0 <= i < n:
            raise Fail("wall index out of range")
        return ("wall", i, side)
    if kind == "pair":
        i = require_int(entry.get("i"), "pair i")
        j = require_int(entry.get("j"), "pair j")
        if i == j or not (0 <= i < n and 0 <= j < n):
            raise Fail("bad pair")
        if i > j:
            i, j = j, i
        return ("pair", i, j)
    raise Fail("unknown constraint kind")


def g_value(con, centres, radius, height):
    if con[0] == "wall":
        _, i, side = con
        x, y = centres[i]
        half = Q(Fraction(1, 2), 0)
        hh = Q(height / 2, 0)
        if side == "L":
            return x + half - radius
        if side == "R":
            return half - x - radius
        if side == "B":
            return y + hh - radius
        if side == "T":
            return hh - y - radius
        raise Fail("bad wall side")
    _, i, j = con
    dx = centres[i][0] - centres[j][0]
    dy = centres[i][1] - centres[j][1]
    return dx * dx + dy * dy - (Q(4, 0) * radius * radius)


def gradient(con, centres, n):
    """Rows of M: partial derivatives in centre order x0, y0, ..., x_{n-1}, y_{n-1}."""
    g = [Q(0, 0) for _ in range(2 * n)]
    if con[0] == "wall":
        _, i, side = con
        if side == "L":
            g[2 * i] = Q(1, 0)
        elif side == "R":
            g[2 * i] = Q(-1, 0)
        elif side == "B":
            g[2 * i + 1] = Q(1, 0)
        elif side == "T":
            g[2 * i + 1] = Q(-1, 0)
        else:
            raise Fail("bad wall side")
        return g
    _, i, j = con
    dx = centres[i][0] - centres[j][0]
    dy = centres[i][1] - centres[j][1]
    two = Q(2, 0)
    g[2 * i] = two * dx
    g[2 * i + 1] = two * dy
    g[2 * j] = -(two * dx)
    g[2 * j + 1] = -(two * dy)
    return g


_WALL_KEY = re.compile(r"wall\s+([LRBT])\s+of\s+circle\s+(\d+)\Z")
_PAIR_KEY = re.compile(r"pair\s+\(\s*(\d+)\s*,\s*(\d+)\s*\)\Z")


def parse_slack_key(key, n):
    if type(key) is not str:
        raise Fail("slack_lo key is not a string")
    wall = _WALL_KEY.fullmatch(key)
    if wall:
        i = int(wall.group(2))
        if not 0 <= i < n:
            raise Fail(f"bad slack_lo key {key}")
        return ("wall", i, wall.group(1))
    pair = _PAIR_KEY.fullmatch(key)
    if pair:
        i = int(pair.group(1))
        j = int(pair.group(2))
        if i > j:
            i, j = j, i
        if not 0 <= i < j < n:
            raise Fail(f"bad slack_lo key {key}")
        return ("pair", i, j)
    raise Fail(f"unrecognised slack_lo key {key}")


def parse_matrix(raw, n, what):
    if not isinstance(raw, list) or len(raw) != n:
        raise Fail(f"{what} is not {n} by {n}")
    matrix = []
    for r, row in enumerate(raw):
        if not isinstance(row, list) or len(row) != n:
            raise Fail(f"{what} row {r} has the wrong length")
        matrix.append([parse_q(entry, what) for entry in row])
    return matrix


def matmul(left, right):
    n = len(left)
    out = []
    for i in range(n):
        row = []
        for j in range(n):
            acc = Q(0, 0)
            for k in range(n):
                acc = acc + left[i][k] * right[k][j]
            row.append(acc)
        out.append(row)
    return out


def is_identity(matrix):
    n = len(matrix)
    one = Q(1, 0)
    zero = Q(0, 0)
    for i in range(n):
        for j in range(n):
            if matrix[i][j] != (one if i == j else zero):
                return False
    return True


def same_matrix(left, right):
    if len(left) != len(right):
        return False
    for a, b in zip(left, right):
        if len(a) != len(b):
            return False
        for x, y in zip(a, b):
            if x != y:
                return False
    return True


def snap_strict(coeff, snap_err, s2_lo):
    absolute = abs(coeff)
    if absolute == 0:
        return snap_err >= 0
    # |b| * (sqrt(2) - S2_LO) = -|b|*S2_LO + |b|*sqrt(2)
    err = Q(-(absolute * s2_lo), absolute)
    return (Q(snap_err, 0) - err).sign() > 0


def check_certificate(cert):
    n = require_int(need(cert, "N"), "N")
    if not 2 <= n <= 8:
        raise Fail("N is outside 2..8")
    height = parse_fraction(need(cert, "h"), "h")
    if height <= 0:
        raise Fail("h is not positive")
    if parse_fraction(need(cert, "width"), "width") != 1:
        raise Fail("width is not 1")

    r_star = need(cert, "r_star")
    if not isinstance(r_star, dict):
        raise Fail("r_star is not an object")
    p = parse_fraction(need(r_star, "p"), "p")
    q = parse_fraction(need(r_star, "q"), "q")
    if q <= 0:
        raise Fail("q must be positive, so r* decreases in sqrt(2)")
    radius = Q(p, -q)
    if radius.sign() <= 0:
        raise Fail("r* is not positive")

    s2_lo = parse_fraction(need(cert, "S2_LO"), "S2_LO")
    s2_hi = parse_fraction(need(cert, "S2_HI"), "S2_HI")
    unit = Fraction(1, 10**50)
    if s2_lo <= 0 or (s2_lo / unit).denominator != 1:
        raise Fail("S2_LO is not a positive multiple of 1e-50")
    if s2_hi - s2_lo != unit:
        raise Fail("S2_HI is not S2_LO + 1e-50")
    if not (s2_lo * s2_lo < 2 < s2_hi * s2_hi):
        raise Fail("S2 enclosure fails S2_LO^2 < 2 < S2_HI^2")

    r_t = parse_fraction(need(cert, "r_t"), "r_t")
    expected_rt = p - q * s2_hi - Fraction(1, 10**40)
    if r_t != expected_rt:
        raise Fail("r_t is not p - q*S2_HI - 1e-40")
    if r_t <= 0:
        raise Fail("r_t is not positive")
    if (radius - Q(r_t, 0)).sign() <= 0:
        raise Fail("r_t is not strictly below r*")

    width_pt = 1 - 2 * r_t
    height_pt = height - 2 * r_t
    diameter = 2 * r_t
    if width_pt <= 0 or height_pt <= 0 or diameter <= 0:
        raise Fail("W, H, or d is not positive")

    g_grid = require_int(need(cert, "grid_G"), "grid_G")
    if g_grid < 2 or g_grid % 2:
        raise Fail("G must be even and at least 2")

    raw_centres = need(cert, "c_star")
    if not isinstance(raw_centres, list) or len(raw_centres) != n:
        raise Fail("c* has the wrong length")
    centres = []
    for point in raw_centres:
        if not isinstance(point, list) or len(point) != 2:
            raise Fail("c* point is malformed")
        centres.append((parse_q(point[0], "c*"), parse_q(point[1], "c*")))

    raw_t = need(cert, "T")
    if not isinstance(raw_t, list):
        raise Fail("T is not a list")
    tight = [parse_constraint(entry, n) for entry in raw_t]
    if len(tight) != len(set(tight)):
        raise Fail("T repeats a constraint")

    zeros = []
    for con in all_constraints(n):
        sgn = g_value(con, centres, radius, height).sign()
        if sgn < 0:
            raise Fail(f"c* violates {con}")
        if sgn == 0:
            zeros.append(con)
    if set(tight) != set(zeros):
        raise Fail("T is not the exact zero set of the constraints at (c*, r*)")

    slack_raw = need(cert, "slack_lo")
    if not isinstance(slack_raw, dict):
        raise Fail("slack_lo is not an object")
    slack = {}
    for key, value in slack_raw.items():
        con = parse_slack_key(key, n)
        if con in slack:
            raise Fail(f"duplicate slack_lo for {con}")
        slack[con] = parse_fraction(value, "slack_lo")
    expected_slack = set(all_constraints(n)) - set(tight)
    if set(slack) != expected_slack:
        raise Fail("slack_lo does not list exactly the non-tight constraints")
    for con, lower in slack.items():
        if lower <= 0:
            raise Fail(f"slack_lo for {con} is not positive")
        gap = g_value(con, centres, radius, height)
        if (gap - Q(lower, 0)).sign() < 0:
            raise Fail(f"slack_lo exceeds the slack of {con}")

    raw_lambda = need(cert, "lambda")
    raw_lo = need(cert, "lambda_lo")
    if not isinstance(raw_lambda, list) or not isinstance(raw_lo, list):
        raise Fail("lambda or lambda_lo is not a list")
    multipliers = [parse_q(item, "lambda") for item in raw_lambda]
    lambda_lo = [parse_fraction(item, "lambda_lo") for item in raw_lo]
    if len(multipliers) != len(tight) or len(lambda_lo) != len(tight):
        raise Fail("lambda length does not match T")
    for index, (lam, lower) in enumerate(zip(multipliers, lambda_lo)):
        if lower <= 0:
            raise Fail(f"lambda_lo[{index}] is not positive")
        if (lam - Q(lower, 0)).sign() < 0:
            raise Fail(f"lambda[{index}] is below its stated lower bound")

    balance = [Q(0, 0) for _ in range(2 * n)]
    weight_sum = Q(0, 0)
    pair_sum = Q(0, 0)
    for lam, con in zip(multipliers, tight):
        for coord, partial in enumerate(gradient(con, centres, n)):
            balance[coord] = balance[coord] + lam * partial
        if con[0] == "wall":
            weight = Q(1, 0)
        else:
            weight = Q(8, 0) * radius
            pair_sum = pair_sum + lam
        weight_sum = weight_sum + lam * weight
    for coord, value in enumerate(balance):
        if value != Q(0, 0):
            raise Fail(f"equilibrium fails at coordinate {coord}")
    if weight_sum != Q(1, 0):
        raise Fail("normalisation sum lambda*w is not 1")

    lambda_hi = parse_fraction(need(cert, "Lambda_hi"), "Lambda_hi")
    if (Q(lambda_hi, 0) - pair_sum).sign() < 0:
        raise Fail("Lambda_hi is below the sum of the pair multipliers")

    raw_b = need(cert, "B")
    if not isinstance(raw_b, list):
        raise Fail("B is not a list")
    basis = [require_int(item, "B") for item in raw_b]
    if len(basis) != 2 * n or len(set(basis)) != len(basis):
        raise Fail("B is not 2N distinct tight constraints")
    if any(index < 0 or index >= len(tight) for index in basis):
        raise Fail("B index is outside T")

    built = [gradient(tight[index], centres, n) for index in basis]
    filed = parse_matrix(need(cert, "M"), 2 * n, "M")
    if not same_matrix(built, filed):
        raise Fail("M does not match the gradients of B at c*")
    inverse = parse_matrix(need(cert, "Minv"), 2 * n, "Minv")
    if not is_identity(matmul(built, inverse)) or not is_identity(matmul(inverse, built)):
        raise Fail("M and Minv are not inverses in both orders")

    norm_hi = parse_fraction(need(cert, "normMinv_hi"), "normMinv_hi")
    if norm_hi <= 0 or lambda_hi <= 0:
        raise Fail("normMinv_hi or Lambda_hi is not positive")
    enc_lo, enc_hi = sqrt2_bounds(80)
    for row_index, row in enumerate(inverse):
        upper = Fraction(0)
        for entry in row:
            bound = rat_upper_abs(entry, enc_lo, enc_hi)
            if (Q(bound, 0) - entry.abs()).sign() < 0:
                raise Fail("internal absolute-value bound failed")
            upper += bound
        if upper > norm_hi:
            raise Fail(f"normMinv_hi is below a rational row-sum bound (row {row_index})")

    rho = parse_fraction(need(cert, "rho"), "rho")
    rho_prime = parse_fraction(need(cert, "rho_prime"), "rho'")
    snap_err = parse_fraction(need(cert, "snap_err_hi"), "snap_err_hi")
    min_lo = min(lambda_lo[index] for index in basis)
    rho_calc = min_lo / (2 * norm_hi) / (16 * lambda_hi)
    if rho != rho_calc:
        raise Fail("rho is not min_B lambda_lo / (2*normMinv_hi) / (16*Lambda_hi)")
    if rho <= 0 or rho_prime <= 0 or snap_err < 0:
        raise Fail("rho, rho', or snap_err_hi has the wrong sign")
    if rho_prime + snap_err > rho:
        raise Fail("rho' + snap_err_hi exceeds rho")
    if (rho_prime * 10**30).denominator != 1:
        raise Fail("rho' is not on the 1e-30 grid")

    images = check_images(cert, centres, s2_lo, s2_hi, r_t, height, snap_err, n)
    return {
        "N": n,
        "h": height,
        "r": radius,
        "centres": centres,
        "r_t": r_t,
        "W": width_pt,
        "H": height_pt,
        "d": diameter,
        "G": g_grid,
        "rho_prime": rho_prime,
        "images": images,
    }


def check_images(cert, centres, s2_lo, s2_hi, r_t, height, snap_err, n):
    raw = need(cert, "images_point_frame")
    if not isinstance(raw, list) or len(raw) != 4:
        raise Fail("expected the four mirror images")
    shift_x = Fraction(1, 2) - r_t
    shift_y = height / 2 - r_t
    max_b = Fraction(0)
    for point in centres:
        max_b = max(max_b, abs(point[0].b), abs(point[1].b))
    if snap_err < max_b * (s2_hi - s2_lo):
        raise Fail("snap_err_hi is below max|b|*(S2_HI - S2_LO)")
    images = []
    seen = set()
    for image in raw:
        if not isinstance(image, dict):
            raise Fail("image is not an object")
        mirror = need(image, "mirror")
        if not isinstance(mirror, list) or len(mirror) != 2:
            raise Fail("mirror is malformed")
        sx = require_int(mirror[0], "sx")
        sy = require_int(mirror[1], "sy")
        if sx not in (1, -1) or sy not in (1, -1):
            raise Fail("mirror signs are not ±1")
        if (sx, sy) in seen:
            raise Fail("duplicate mirror")
        seen.add((sx, sy))
        raw_points = need(image, "points")
        if not isinstance(raw_points, list) or len(raw_points) != n:
            raise Fail("image has the wrong number of points")
        points = []
        for k, pair in enumerate(raw_points):
            if not isinstance(pair, list) or len(pair) != 2:
                raise Fail("image point is malformed")
            filed_x = parse_fraction(pair[0], "image x")
            filed_y = parse_fraction(pair[1], "image y")
            exact_x = sx * (centres[k][0].a + centres[k][0].b * s2_lo) + shift_x
            exact_y = sy * (centres[k][1].a + centres[k][1].b * s2_lo) + shift_y
            if filed_x != exact_x or filed_y != exact_y:
                raise Fail(f"mirror {(sx, sy)} point {k} is not the S2_LO snap")
            if not snap_strict(centres[k][0].b, snap_err, s2_lo):
                raise Fail("snap error exceeds snap_err_hi")
            if not snap_strict(centres[k][1].b, snap_err, s2_lo):
                raise Fail("snap error exceeds snap_err_hi")
            points.append((filed_x, filed_y))
        images.append(points)
    if seen != {(1, 1), (1, -1), (-1, 1), (-1, -1)}:
        raise Fail("mirrors are not exactly the four sign pairs")
    return images


def canon_rect(rect, g_grid):
    if not isinstance(rect, list) or len(rect) != 4:
        raise Fail("rectangle is malformed")
    vals = [require_int(v, "grid coordinate") for v in rect]
    for value in vals:
        if value < 0 or value > g_grid:
            raise Fail("grid coordinate is outside [0, G]")
    return vals


def canon_box(box, n, g_grid):
    if not isinstance(box, list) or len(box) != n:
        raise Fail("box does not have one rectangle per point")
    return [canon_rect(rect, g_grid) for rect in box]


def range_lt(far_x, far_y, width, height, diameter, g_grid):
    """Strict range maximum of |p-q|^2 over the two axis-aligned boxes, < d^2."""
    scale = Fraction(g_grid)
    value = (Fraction(far_x) * width / scale) ** 2 + (Fraction(far_y) * height / scale) ** 2
    return value < diameter * diameter


def separation(a, b, width, height, diameter, g_grid):
    """far on each axis, then whether the range maximum is strictly below d^2.

    Empty intervals (lo > hi) are not a pair witness. A line (lo == hi) is a
    real closed box. far(a0,a1,b0,b1) = max(a1-b0, b1-a0) is the max |u-v|.
    """
    if a[0] > a[1] or a[2] > a[3] or b[0] > b[1] or b[2] > b[3]:
        return None
    far_x = max(a[1] - b[0], b[1] - a[0])
    far_y = max(a[3] - b[2], b[3] - a[2])
    return far_x, far_y, range_lt(far_x, far_y, width, height, diameter, g_grid)


def check_witness(box, witness, n, g_grid, width, height, diameter):
    if not isinstance(witness, list) or not witness or type(witness[0]) is not str:
        raise Fail("witness is malformed")
    tag = witness[0]
    if tag == "pair":
        if len(witness) != 5:
            raise Fail("pair witness is malformed")
        i = require_int(witness[1], "pair i")
        j = require_int(witness[2], "pair j")
        far_x = require_int(witness[3], "fx")
        far_y = require_int(witness[4], "fy")
        if not (0 <= i < n and 0 <= j < n) or i == j:
            raise Fail("pair witness indices are bad")
        detail = separation(box[i], box[j], width, height, diameter, g_grid)
        if detail is None or detail[0] != far_x or detail[1] != far_y or not detail[2]:
            raise Fail("pair witness is false")
        return tag
    if tag == "forbid":
        if len(witness) != 3:
            raise Fail("forbid witness is malformed")
        i = require_int(witness[1], "forbid i")
        j = require_int(witness[2], "forbid j")
        if not (0 <= i < n and 0 <= j < n) or i == j:
            raise Fail("forbid witness indices are bad")
        detail = separation(box[i], box[j], width, height, diameter, g_grid)
        if detail is None or not detail[2]:
            raise Fail("forbid witness is false")
        return tag
    if tag == "empty":
        if len(witness) != 3:
            raise Fail("empty witness is malformed")
        i = require_int(witness[1], "empty i")
        axis = require_int(witness[2], "empty axis")
        if not 0 <= i < n or axis not in (0, 1):
            raise Fail("empty witness is out of range")
        x0, x1, y0, y1 = box[i]
        empty = x0 > x1 if axis == 0 else y0 > y1
        if not empty:
            raise Fail("empty witness is false")
        return tag
    if tag == "xorder":
        if len(witness) != 2:
            raise Fail("xorder witness is malformed")
        i = require_int(witness[1], "xorder")
        if not 0 <= i < n - 1:
            raise Fail("xorder index is out of range")
        if not box[i][0] > box[i + 1][1]:
            raise Fail("xorder witness is false")
        return tag
    if tag == "xsum":
        if len(witness) != 1:
            raise Fail("xsum witness is malformed")
        if not box[0][0] + box[n - 1][0] > g_grid:
            raise Fail("xsum witness is false")
        return tag
    raise Fail(f"unknown witness {tag}")


def apply_sh(box, op, n, g_grid, width, height, diameter):
    if len(op) != 5:
        raise Fail("sh is malformed")
    i = require_int(op[1], "sh i")
    j = require_int(op[2], "sh j")
    side = op[3]
    endpoint = require_int(op[4], "sh endpoint")
    if type(side) is not str or side not in ("x0", "x1", "y0", "y1"):
        raise Fail("sh side is not x0, x1, y0, or y1")
    if not (0 <= i < n and 0 <= j < n) or i == j:
        raise Fail("sh indices are bad")
    if not 0 <= endpoint <= g_grid:
        raise Fail("sh endpoint is outside [0, G]")
    x0, x1, y0, y1 = box[i]
    if side == "x0":
        if not x0 <= endpoint <= x1:
            raise Fail("sh does not move inward")
        strip = [x0, endpoint, y0, y1]
        updated = [endpoint, x1, y0, y1]
        moved = endpoint > x0
    elif side == "x1":
        if not x0 <= endpoint <= x1:
            raise Fail("sh does not move inward")
        strip = [endpoint, x1, y0, y1]
        updated = [x0, endpoint, y0, y1]
        moved = endpoint < x1
    elif side == "y0":
        if not y0 <= endpoint <= y1:
            raise Fail("sh does not move inward")
        strip = [x0, x1, y0, endpoint]
        updated = [x0, x1, endpoint, y1]
        moved = endpoint > y0
    else:
        if not y0 <= endpoint <= y1:
            raise Fail("sh does not move inward")
        strip = [x0, x1, endpoint, y1]
        updated = [x0, x1, y0, endpoint]
        moved = endpoint < y1
    if moved:
        detail = separation(strip, box[j], width, height, diameter, g_grid)
        if detail is None or not detail[2]:
            raise Fail("removed strip is not strictly within distance d")
    box[i] = updated


def apply_op(box, op, n, g_grid, width, height, diameter):
    """Clips are exactly the non-strict cuts. sh removes one boundary strip."""
    if not isinstance(op, list) or not op or type(op[0]) is not str:
        raise Fail("op is malformed")
    name = op[0]
    if name == "o_hi":
        if len(op) != 2:
            raise Fail("o_hi is malformed")
        i = require_int(op[1], "o_hi")
        if not 0 <= i < n - 1:
            raise Fail("o_hi index is out of range")
        box[i][1] = min(box[i][1], box[i + 1][1])
    elif name == "o_lo":
        if len(op) != 2:
            raise Fail("o_lo is malformed")
        i = require_int(op[1], "o_lo")
        if not 0 <= i < n - 1:
            raise Fail("o_lo index is out of range")
        box[i + 1][0] = max(box[i + 1][0], box[i][0])
    elif name == "s1":
        if len(op) != 1:
            raise Fail("s1 is malformed")
        box[0][1] = min(box[0][1], g_grid - box[n - 1][0])
    elif name == "sN":
        if len(op) != 1:
            raise Fail("sN is malformed")
        box[n - 1][1] = min(box[n - 1][1], g_grid - box[0][0])
    elif name == "y1":
        if len(op) != 1:
            raise Fail("y1 is malformed")
        box[0][3] = min(box[0][3], g_grid // 2)
    elif name == "sh":
        apply_sh(box, op, n, g_grid, width, height, diameter)
    else:
        raise Fail(f"unknown op {name}")
    for point in box:
        for value in point:
            if value < 0 or value > g_grid:
                raise Fail("op left the grid [0, G]")


def box_key(box):
    return tuple(tuple(point) for point in box)


def same_parts(boxes, parts):
    return sorted(box_key(b) for b in boxes) == sorted(box_key(b) for b in parts)


def child_start_boxes(ids, nodes, n, g_grid):
    boxes = []
    for cid in ids:
        if cid not in nodes:
            raise Fail(f"missing child {cid}")
        boxes.append(canon_box(nodes[cid]["box"], n, g_grid))
    return boxes


def bisect_parts(tbox, how, n):
    if len(how) != 4:
        raise Fail("bisect is malformed")
    i = require_int(how[1], "bisect point")
    axis = require_int(how[2], "bisect axis")
    mid = require_int(how[3], "bisect mid")
    if not 0 <= i < n or axis not in (0, 1):
        raise Fail("bisect is out of range")
    rect = tbox[i]
    lo, hi = (rect[0], rect[1]) if axis == 0 else (rect[2], rect[3])
    if mid != (lo + hi) // 2:
        raise Fail("bisect mid is not floor((lo+hi)/2)")
    if not lo < mid < hi:
        raise Fail("bisect does not split the side")

    def put(lo2, hi2):
        out = [row[:] for row in tbox]
        if axis == 0:
            out[i] = [lo2, hi2, rect[2], rect[3]]
        else:
            out[i] = [rect[0], rect[1], lo2, hi2]
        return out

    return [put(lo, mid), put(mid, hi)]


def hole_parts(tbox, how, n, g_grid, width, height, diameter):
    if len(how) != 4:
        raise Fail("hole is malformed")
    i = require_int(how[1], "hole point")
    j = require_int(how[2], "hole other")
    rect = how[3]
    if not isinstance(rect, list) or len(rect) != 4:
        raise Fail("hole rectangle is malformed")
    hx0, hx1, hy0, hy1 = [require_int(v, "hole") for v in rect]
    if not 0 <= i < n or not 0 <= j < n or i == j:
        raise Fail("hole indices are bad")
    x0, x1, y0, y1 = tbox[i]
    if not (x0 <= hx0 < hx1 <= x1 and y0 <= hy0 < hy1 <= y1):
        raise Fail("hole is not a positive sub-rectangle")
    detail = separation([hx0, hx1, hy0, hy1], tbox[j], width, height, diameter, g_grid)
    if detail is None or not detail[2]:
        raise Fail("hole is not forbidden")
    # Left, right, lower, upper closed remnants. Shared edges are allowed.
    cands = [
        [x0, hx0, y0, y1],
        [hx1, x1, y0, y1],
        [hx0, hx1, y0, hy0],
        [hx0, hx1, hy1, y1],
    ]

    def build(positive_only):
        parts = []
        for cand in cands:
            if positive_only:
                keep = cand[0] < cand[1] and cand[2] < cand[3]
            else:
                keep = cand[0] <= cand[1] and cand[2] <= cand[3]
            if keep:
                out = [row[:] for row in tbox]
                out[i] = cand[:]
                parts.append(out)
        return parts

    return build(True), build(False)


def abut(intervals, g_grid, name):
    if not intervals:
        raise Fail(f"{name} cover is empty")
    intervals = sorted(intervals)
    if intervals[0][0] != 0 or intervals[-1][1] != g_grid:
        raise Fail(f"{name} cells do not cover [0, G]")
    for lo, hi in intervals:
        if not lo < hi:
            raise Fail(f"{name} cell has no width")
    for (_, hi), (lo, _) in zip(intervals, intervals[1:]):
        # Closed cover: touch or overlap. A gap is next_lo > prev_hi.  (Grok's fix, 2026-09-25 reply-proofchecker2-fix.md)
        if lo > hi:
            raise Fail(f"{name} cell cover has a gap")


def cover_ok(cells, g_grid, width, height, diameter):
    tuples = []
    for cell in cells:
        if not (cell[0] < cell[1] and cell[2] < cell[3]):
            raise Fail("cell is degenerate")
        far_x = cell[1] - cell[0]
        far_y = cell[3] - cell[2]
        if not range_lt(far_x, far_y, width, height, diameter, g_grid):
            raise Fail("cell diameter is not strictly below d")
        tuples.append(tuple(cell))
    if len(tuples) != len(set(tuples)):
        raise Fail("duplicate cell")
    xs = sorted({(item[0], item[1]) for item in tuples})
    ys = sorted({(item[2], item[3]) for item in tuples})
    if len(tuples) != len(xs) * len(ys):
        raise Fail("cells are not a product")
    have = set(tuples)
    for x0, x1 in xs:
        for y0, y1 in ys:
            if (x0, x1, y0, y1) not in have:
                raise Fail("product cell is missing")
    abut(xs, g_grid, "x")
    abut(ys, g_grid, "y")


def as_injection(raw, n, n_cells):
    if not isinstance(raw, list) or len(raw) != n:
        raise Fail("assignment has the wrong length")
    assign = tuple(require_int(v, "assignment") for v in raw)
    if len(set(assign)) != n or any(v < 0 or v >= n_cells for v in assign):
        raise Fail("assignment is not an injection into the cells")
    return assign


def as_perm(raw, n):
    if not isinstance(raw, list) or len(raw) != n:
        raise Fail("perm has the wrong length")
    perm = [require_int(v, "perm") for v in raw]
    if sorted(perm) != list(range(n)):
        raise Fail("perm is not a permutation")
    return perm


def parse_discard(record):
    if isinstance(record, dict):
        assign = None
        for key in ("assign", "asg", "perm"):
            if key in record:
                assign = record[key]
                break
        witness = record.get("w", record.get("witness"))
        if assign is None or witness is None:
            raise Fail("discard record is missing an assignment or a witness")
        return assign, witness
    if (
        isinstance(record, list)
        and len(record) == 2
        and isinstance(record[0], list)
        and isinstance(record[1], list)
    ):
        return record[0], record[1]
    raise Fail("discard record is unrecognised")


def check_cells(fate, nodes, n, g_grid, width, height, diameter):
    raw_cells = fate.get("cells")
    assign_raw = fate.get("assign")
    children = fate.get("children")
    discards = fate.get("discards")
    if not isinstance(raw_cells, list) or not raw_cells:
        raise Fail("root cells are missing")
    if (
        not isinstance(assign_raw, list)
        or not isinstance(children, list)
        or not isinstance(discards, list)
    ):
        raise Fail("root is missing assign, children, or discards")
    cells = [canon_rect(cell, g_grid) for cell in raw_cells]
    cover_ok(cells, g_grid, width, height, diameter)
    n_cells = len(cells)
    if n_cells < n:
        raise Fail("fewer cells than points")
    injections = 1
    for k in range(n):
        injections *= n_cells - k
        if injections > 100000:
            raise Fail("too many cell assignments to enumerate")
    if len(assign_raw) != len(children):
        raise Fail("kept assignments and children differ in length")
    child_ids = [require_int(cid, "child") for cid in children]
    if len(child_ids) != len(set(child_ids)):
        raise Fail("duplicate kept start")
    kept = []
    for assign, cid in zip(assign_raw, child_ids):
        tup = as_injection(assign, n, n_cells)
        if cid not in nodes:
            raise Fail(f"missing start node {cid}")
        expected = [cells[index][:] for index in tup]
        if canon_box(nodes[cid]["box"], n, g_grid) != expected:
            raise Fail(f"start node {cid} is not its cell assignment")
        kept.append(tup)
    rejected = []
    for record in discards:
        assign, witness = parse_discard(record)
        tup = as_injection(assign, n, n_cells)
        box = [cells[index][:] for index in tup]
        check_witness(box, witness, n, g_grid, width, height, diameter)
        rejected.append(tup)
    if len(kept) != len(set(kept)):
        raise Fail("duplicate kept assignment")
    if len(rejected) != len(set(rejected)):
        raise Fail("duplicate rejected assignment")
    if set(kept) & set(rejected):
        raise Fail("an assignment is both kept and rejected")
    universe = set(itertools.permutations(range(n_cells), n))
    if set(kept) | set(rejected) != universe:
        raise Fail("an injective cell assignment is unaccounted for")
    return child_ids, len(kept), len(rejected), n_cells


def check_accept(tbox, fate, images, rho_prime, width, height, g_grid, n):
    mirror = require_int(fate.get("mirror"), "mirror")
    if not 0 <= mirror < len(images):
        raise Fail("mirror index is out of range")
    perm = as_perm(fate.get("perm"), n)
    scale = Fraction(g_grid)
    points = images[mirror]
    for i in range(n):
        x0, x1, y0, y1 = tbox[i]
        if x0 > x1 or y0 > y1:
            raise Fail("accepted box is empty")
        cx, cy = points[perm[i]]
        x_lo = Fraction(x0) * width / scale
        x_hi = Fraction(x1) * width / scale
        y_lo = Fraction(y0) * height / scale
        y_hi = Fraction(y1) * height / scale
        if not (cx - rho_prime <= x_lo and x_hi <= cx + rho_prime):
            raise Fail(f"accepted x box {i} is outside the rho' square")
        if not (cy - rho_prime <= y_lo and y_hi <= cy + rho_prime):
            raise Fail(f"accepted y box {i} is outside the rho' square")


def in_cut(xs, ys):
    for i in range(len(xs) - 1):
        if (xs[i + 1] - xs[i]).sign() < 0:
            return False
    if (xs[0] + xs[-1]).sign() > 0:
        return False
    return ys[0].sign() <= 0


def contains_config(tbox, xs, ys, shift_x, shift_y, width, height, g_grid):
    scale = Fraction(g_grid)
    for i in range(len(xs)):
        px = xs[i] + shift_x
        py = ys[i] + shift_y
        x0, x1, y0, y1 = tbox[i]
        if (px - Q(Fraction(x0) * width / scale, 0)).sign() < 0:
            return False
        if (Q(Fraction(x1) * width / scale, 0) - px).sign() < 0:
            return False
        if (py - Q(Fraction(y0) * height / scale, 0)).sign() < 0:
            return False
        if (Q(Fraction(y1) * height / scale, 0) - py).sign() < 0:
            return False
    return True


def check_representatives(centres, accepted, r_t, height, width, height_pt, g_grid):
    n = len(centres)
    shift_x = Q(Fraction(1, 2) - r_t, 0)
    shift_y = Q(height / 2 - r_t, 0)
    found = 0
    for sx in (1, -1):
        for sy in (1, -1):
            for perm in itertools.permutations(range(n)):
                xs = [Q(sx, 0) * centres[perm[i]][0] for i in range(n)]
                ys = [Q(sy, 0) * centres[perm[i]][1] for i in range(n)]
                if not in_cut(xs, ys):
                    continue
                found += 1
                hit = any(
                    contains_config(box, xs, ys, shift_x, shift_y, width, height_pt, g_grid)
                    for box in accepted
                )
                if not hit:
                    raise Fail("a D-representative of c* lies in no accepted leaf")
    if found == 0:
        raise Fail("no image of c* lies in the cut set D")
    return found


def walk(root, nodes, children_of):
    seen = set()
    stack = [root]
    while stack:
        nid = stack.pop()
        if nid in seen:
            raise Fail(f"node {nid} is reachable more than once")
        seen.add(nid)
        for cid in children_of[nid]:
            if cid not in nodes:
                raise Fail(f"node {nid} lists a missing child {cid}")
            parent = nodes[cid].get("parent", "missing")
            if parent != nid:
                raise Fail(f"node {cid} has a parent link that does not match the tree")
            stack.append(cid)
    if seen != set(nodes):
        orphan = min(set(nodes) - seen)
        raise Fail(f"orphan node {orphan}")


def tighten(box, ops, n, g_grid, width, height, diameter):
    current = [row[:] for row in box]
    if not isinstance(ops, list):
        raise Fail("ops is not a list")
    for index, op in enumerate(ops):
        try:
            apply_op(current, op, n, g_grid, width, height, diameter)
        except Fail as exc:
            raise Fail(f"op {index}: {exc}") from None
    return current


def process_node(node, ctx, nodes, state):
    n = ctx["N"]
    g_grid = ctx["G"]
    width, height, diameter = ctx["W"], ctx["H"], ctx["d"]
    nid = require_int(node["id"], "id")
    if "parent" not in node:
        raise Fail("parent is missing")
    parent = node["parent"]
    if parent is not None:
        parent = require_int(parent, "parent")
    box = canon_box(node["box"], n, g_grid)
    current = tighten(box, node.get("ops", None), n, g_grid, width, height, diameter)
    if "tbox" in node and canon_box(node["tbox"], n, g_grid) != current:
        raise Fail("tbox does not match the ops")
    fate = node.get("fate")
    if not isinstance(fate, dict) or type(fate.get("t")) is not str:
        raise Fail("fate is missing")
    kind = fate["t"]
    if kind == "cells":
        if parent is not None:
            raise Fail("cells fate is only valid at the root")
        if state["root"] is not None:
            raise Fail("more than one root")
        full = [0, g_grid, 0, g_grid]
        if any(point != full for point in box):
            raise Fail("root box is not [0, G]^2 for every point")
        ids, n_kept, n_rej, n_cells = check_cells(
            fate, nodes, n, g_grid, width, height, diameter
        )
        state["root"] = nid
        state["children"][nid] = ids
        state["n_kept"] = n_kept
        state["n_rej"] = n_rej
        state["n_cells"] = n_cells
        return
    if parent is None:
        raise Fail("root fate is not the cell cover")
    if kind == "split":
        if "tbox" not in node:
            raise Fail("split is missing tbox")
        how = fate.get("how")
        raw_ids = fate.get("children")
        if not isinstance(raw_ids, list):
            raise Fail("split is missing children")
        ids = [require_int(cid, "child") for cid in raw_ids]
        if len(ids) != len(set(ids)):
            raise Fail("split repeats a child")
        boxes = child_start_boxes(ids, nodes, n, g_grid)
        if not isinstance(how, list) or not how or type(how[0]) is not str:
            raise Fail("split how is malformed")
        if how[0] == "bisect":
            parts = bisect_parts(current, how, n)
            if not same_parts(boxes, parts):
                raise Fail("bisect children do not partition the tightened box")
        elif how[0] == "hole":
            positive, closed = hole_parts(current, how, n, g_grid, width, height, diameter)
            if not (same_parts(boxes, positive) or same_parts(boxes, closed)):
                raise Fail("hole children do not partition the tightened box")
        else:
            raise Fail(f"unknown split {how[0]}")
        state["children"][nid] = ids
        state["n_split"] += 1
        return
    if kind == "discard":
        if fate.get("children"):
            raise Fail("discard leaf has children")
        witness = fate.get("w", fate.get("witness"))
        if witness is None:
            raise Fail("discard is missing a witness")
        tag = check_witness(current, witness, n, g_grid, width, height, diameter)
        state["hist"][tag] = state["hist"].get(tag, 0) + 1
        state["children"][nid] = []
        state["n_dis"] += 1
        return
    if kind == "accept":
        if "tbox" not in node:
            raise Fail("accept is missing tbox")
        if fate.get("children"):
            raise Fail("accept leaf has children")
        check_accept(current, fate, ctx["images"], ctx["rho_prime"], width, height, g_grid, n)
        state["accepted"].append([row[:] for row in current])
        state["children"][nid] = []
        state["n_acc"] += 1
        return
    raise Fail(f"undecided or unknown fate {kind}")


def load_bytes(path, what):
    try:
        with open(path, "rb") as handle:
            return handle.read()
    except OSError as exc:
        raise Fail(f"cannot read {what}: {exc}") from None


def load_tree(path):
    raw = load_bytes(path, "tree")
    if raw.startswith(b"\x1f\x8b"):
        try:
            raw = gzip.decompress(raw)
        except Exception as exc:
            raise Fail(f"cannot gunzip tree: {exc}") from None
    try:
        tree = json.loads(raw.decode("utf-8"))
    except Exception as exc:
        raise Fail(f"tree is not json: {exc}") from None
    if not isinstance(tree, dict):
        raise Fail("tree is not an object")
    return tree


def verify(cert_path, tree_path, info):
    cert_bytes = load_bytes(cert_path, "certificate")
    try:
        cert = json.loads(cert_bytes.decode("utf-8"))
    except Exception as exc:
        raise Fail(f"certificate is not json: {exc}") from None
    if not isinstance(cert, dict):
        raise Fail("certificate is not an object")
    case = cert.get("case", "?")
    if type(case) is not str or not case or any(ch.isspace() for ch in case):
        raise Fail("bad case name")
    info["case"] = case

    ctx = check_certificate(cert)
    tree = load_tree(tree_path)
    if tree.get("case") != case:
        raise Fail("tree case does not match the certificate")
    digest = hashlib.sha256(cert_bytes).hexdigest()
    filed_hash = tree.get("cert_sha256")
    if type(filed_hash) is not str or len(filed_hash) != 64:
        raise Fail("cert_sha256 is not 64 hex digits")
    if any(ch not in "0123456789abcdefABCDEF" for ch in filed_hash):
        raise Fail("cert_sha256 is not hex")
    if filed_hash.lower() != digest:
        raise Fail("cert_sha256 does not match the certificate bytes")

    n = ctx["N"]
    if require_int(need(tree, "N"), "tree N") != n:
        raise Fail("tree N does not match")
    if parse_fraction(need(tree, "h"), "tree h") != ctx["h"]:
        raise Fail("tree h does not match")
    if parse_fraction(need(tree, "r_t"), "tree r_t") != ctx["r_t"]:
        raise Fail("tree r_t does not match the certificate")
    if require_int(need(tree, "G"), "tree G") != ctx["G"]:
        raise Fail("tree G does not match the certificate")
    if parse_fraction(need(tree, "W"), "W") != ctx["W"]:
        raise Fail("W is not 1 - 2*r_t")
    if parse_fraction(need(tree, "H"), "H") != ctx["H"]:
        raise Fail("H is not h - 2*r_t")
    if parse_fraction(need(tree, "d"), "d") != ctx["d"]:
        raise Fail("d is not 2*r_t")

    raw_nodes = need(tree, "nodes")
    if not isinstance(raw_nodes, list) or not raw_nodes:
        raise Fail("tree has no nodes")
    nodes = {}
    for node in raw_nodes:
        if not isinstance(node, dict) or type(node.get("id")) is not int:
            raise Fail("node id is missing or not an int")
        if node["id"] in nodes:
            raise Fail(f"duplicate node id {node['id']}")
        nodes[node["id"]] = node

    state = {
        "children": {},
        "accepted": [],
        "hist": {},
        "n_acc": 0,
        "n_dis": 0,
        "n_split": 0,
        "n_kept": 0,
        "n_rej": 0,
        "n_cells": 0,
        "root": None,
    }
    for node in raw_nodes:
        nid = node["id"]
        try:
            process_node(node, ctx, nodes, state)
        except Fail as exc:
            raise Fail(f"node {nid}: {exc}") from None
    if state["root"] is None:
        raise Fail("tree has no root")
    if set(state["children"]) != set(nodes):
        raise Fail("a node was not given a fate")
    walk(state["root"], nodes, state["children"])
    found = check_representatives(
        ctx["centres"],
        state["accepted"],
        ctx["r_t"],
        ctx["h"],
        ctx["W"],
        ctx["H"],
        ctx["G"],
    )
    if state["n_acc"] < 1:
        raise Fail("no accepted leaf")
    return [
        f"checked {case}: N={n} h={ctx['h']} G={ctx['G']}",
        "certificate: sqrt(2) enclosure, c* feasible, T the exact zero set, "
        "0 < slack_lo <= slack, lambda >= lambda_lo > 0, equilibrium, "
        "normalisation, Lambda <= Lambda_hi, M rebuilt from c*, "
        "Minv inverse both ways, rational norm bound, rho formula, "
        "r_t < r* via S2_HI, four snaps, rho' + snap_err_hi <= rho",
        "tree: sha256, even G, root [0,G]^(2N), product cell cover, "
        f"all injections ({state['n_kept']} kept, {state['n_rej']} rejected), "
        f"ops, witnesses, {state['n_split']} splits, {state['n_acc']} accepted, "
        f"{state['n_dis']} discarded, 0 undecided, each of {len(nodes)} nodes once, "
        f"{found} D-representatives of c* inside an accepted leaf",
    ]


def main(argv):
    info = {"case": "?"}
    try:
        if len(argv) != 3:
            raise Fail("usage: python verify_proof_exact.py <cert.json> <tree.json.gz>")
        notes = verify(argv[1], argv[2], info)
    except Fail as exc:
        text = str(exc).replace("\n", " ")
        print(f"FAIL {info['case']}: {text}")
        return 1
    except Exception as exc:
        text = str(exc).replace("\n", " ")
        print(f"FAIL {info['case']}: {type(exc).__name__}: {text}")
        return 1
    for line in notes:
        print(line)
    print(f"PASS {info['case']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
