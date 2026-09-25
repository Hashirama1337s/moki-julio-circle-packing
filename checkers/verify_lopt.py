#!/usr/bin/env python3
"""Local-optimality checker for equal-circle packings.

CLI: verify_lopt.py <container> <certificate.txt> <cert.json>
container is tri, rect:<h>, or quad.
Proves the rho / t0 / Delta certificate from positive multipliers,
feasibility of z*, and an exact-rational bound K on a float left inverse.
"""

from __future__ import annotations

import json
import math
import sys
import warnings
from fractions import Fraction

import numpy as np


def parse_decimal(s: str) -> Fraction:
    s = s.strip()
    if not s or len(s) > 20000:
        raise ValueError("bad decimal")
    sign = 1
    if s[0] == "+":
        s = s[1:]
    elif s[0] == "-":
        sign = -1
        s = s[1:]
    low = s.lower()
    exp = 0
    if "e" in low:
        mant, exp_s = low.split("e")
        exp = int(exp_s, 10)
        s = mant
    if not s or any(c in s.lower() for c in "naninf"):
        raise ValueError("bad decimal")
    if s.count(".") > 1:
        raise ValueError("bad decimal")
    if "." in s:
        ip, fp = s.split(".")
        if ip == "":
            ip = "0"
        if fp == "":
            fp = ""
        if not ip.isdigit() or not fp.isdigit():
            raise ValueError("bad decimal")
        digits = ip + fp
        scale = len(fp)
    else:
        if not s.isdigit():
            raise ValueError("bad decimal")
        digits = s
        scale = 0
    scale -= exp
    num = int(digits) if digits else 0
    if scale >= 0:
        return Fraction(sign * num, 10**scale)
    return Fraction(sign * num * 10 ** (-scale), 1)


def parse_number(s: str) -> Fraction:
    s = s.strip()
    if not s or len(s) > 20000:
        raise ValueError("bad number")
    if "e" in s.lower():
        return parse_decimal(s)
    if any(c in s.lower() for c in "naninf"):
        raise ValueError("non-finite number")
    return Fraction(s)


def parse_hex_float(s: str) -> Fraction:
    """Exact rational value of a Python float.hex() token, full written precision."""
    if not isinstance(s, str) or len(s) > 10000:
        raise ValueError("bad hex float")
    t = s.strip().lower()
    if not t:
        raise ValueError("bad hex float")
    sign = 1
    if t[0] == "+":
        t = t[1:]
    elif t[0] == "-":
        sign = -1
        t = t[1:]
    if not t.startswith("0x") or "p" not in t:
        raise ValueError("bad hex float")
    mant, exp_s = t[2:].split("p")
    exp = int(exp_s, 10)
    if mant.count(".") > 1:
        raise ValueError("bad hex float")
    if "." in mant:
        ip, fp = mant.split(".")
    else:
        ip, fp = mant, ""
    if ip == "":
        ip = "0"
    hexdigits = ip + fp
    if hexdigits == "" or any(c not in "0123456789abcdef" for c in hexdigits):
        raise ValueError("bad hex float")
    val = int(hexdigits, 16)
    shift = exp - 4 * len(fp)
    if shift >= 0:
        return Fraction(sign * val * (1 << shift), 1)
    return Fraction(sign * val, 1 << (-shift))


def _ge0_int(a: int, b: int) -> bool:
    """a + b*sqrt(2) >= 0 for integers a, b."""
    if a >= 0 and b >= 0:
        return True
    if a <= 0 and b <= 0:
        return False
    aa = a * a
    two = (b * b) << 1
    if a > 0:
        return aa >= two
    return two >= aa


def frac_ge0(a: Fraction, b: Fraction) -> bool:
    den = math.lcm(a.denominator, b.denominator)
    A = a.numerator * (den // a.denominator)
    B = b.numerator * (den // b.denominator)
    return _ge0_int(A, B)


def int_abs_sqrt2(a: int, b: int) -> tuple[int, int]:
    """|a + b*sqrt(2)| = c + d*sqrt(2) with c, d integers."""
    if a == 0 and b == 0:
        return 0, 0
    if _ge0_int(a, b):
        return a, b
    return -a, -b


def sqrt2_bounds(places: int) -> tuple[Fraction, Fraction]:
    """Rational bounds lo <= sqrt(2) < hi, hi - lo = 10^{-places}."""
    if places < 1:
        places = 1
    N = 10**places
    s = math.isqrt(2 * N * N)
    return Fraction(s, N), Fraction(s + 1, N)


class Alg:
    """Element a + b*sqrt(2) of Q(sqrt(2))."""

    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    def __add__(self, other):
        if not isinstance(other, Alg):
            other = Alg(other, 0)
        return Alg(self.a + other.a, self.b + other.b)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return Alg(-self.a, -self.b)

    def __sub__(self, other):
        if not isinstance(other, Alg):
            other = Alg(other, 0)
        return Alg(self.a - other.a, self.b - other.b)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if isinstance(other, Alg):
            return Alg(
                self.a * other.a + 2 * self.b * other.b,
                self.a * other.b + self.b * other.a,
            )
        f = other if isinstance(other, Fraction) else Fraction(other)
        return Alg(self.a * f, self.b * f)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        f = other if isinstance(other, Fraction) else Fraction(other)
        return Alg(self.a / f, self.b / f)

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def is_neg(self) -> bool:
        if self.is_zero():
            return False
        return not frac_ge0(self.a, self.b)

    def __lt__(self, other):
        if not isinstance(other, Alg):
            other = Alg(other, 0)
        return (self - other).is_neg()

    def abs(self) -> "Alg":
        if self.is_neg():
            return -self
        return Alg(self.a, self.b)


def parse_container(s: str):
    s = s.strip()
    if s == "tri":
        return "tri", None
    if s == "quad":
        return "quad", None
    if s == "semi":                     # PATCH (2026-09-24, checker B's spec): unit semicircle, y >= 0
        return "semi", None
    if s.startswith("rect:"):
        h = parse_number(s[5:])
        if h <= 0:
            raise ValueError("rectangle height must be positive")
        return "rect", h
    raise ValueError("container must be tri, rect:<h>, quad or semi")


def wall_ks(kind: str):
    if kind == "rect":
        return (0, 1, 2, 3)
    if kind == "semi":                  # PATCH (checker B's spec): k=1 y - r, k=2 arc (same polynomials as quad); k=0 is not a wall
        return (1, 2)
    return (0, 1, 2)


def kappa_of(kind: str, con) -> int:
    tag = con[0]
    if tag == "p":
        return 8
    k = con[2]
    if kind in ("quad", "semi") and k == 2:   # PATCH (checker B's spec): the semicircle arc has the quadrant arc's kappa
        return 2
    return 0


def g_value(kind, h, r, centres, con) -> Alg:
    tag, i, j = con
    if tag == "p":
        xi, yi = centres[i]
        xj, yj = centres[j]
        dx = xi - xj
        dy = yi - yj
        return Alg(dx * dx + dy * dy - 4 * r * r, 0)
    xi, yi = centres[i]
    k = j
    if kind == "tri":
        if k == 0:
            return Alg(xi - r, 0)
        if k == 1:
            return Alg(yi - r, 0)
        return Alg(1 - xi - yi, -r)
    if kind == "rect":
        if k == 0:
            return Alg(xi + Fraction(1, 2) - r, 0)
        if k == 1:
            return Alg(Fraction(1, 2) - xi - r, 0)
        if k == 2:
            return Alg(yi + h / 2 - r, 0)
        return Alg(h / 2 - yi - r, 0)
    if k == 0:
        return Alg(xi - r, 0)
    if k == 1:
        return Alg(yi - r, 0)
    return Alg((1 - r) * (1 - r) - xi * xi - yi * yi, 0)


def gradients(kind, h, r, centres, con):
    """Nonzero partial derivatives as (coord_key, Alg)."""
    tag, i, j = con
    if tag == "p":
        xi, yi = centres[i]
        xj, yj = centres[j]
        dx = xi - xj
        dy = yi - yj
        two = Fraction(2)
        return [
            (("x", i), Alg(two * dx, 0)),
            (("y", i), Alg(two * dy, 0)),
            (("x", j), Alg(-two * dx, 0)),
            (("y", j), Alg(-two * dy, 0)),
            (("r", None), Alg(-8 * r, 0)),
        ]
    xi, yi = centres[i]
    k = j
    if kind == "tri":
        if k == 0:
            return [(("x", i), Alg(1, 0)), (("r", None), Alg(-1, 0))]
        if k == 1:
            return [(("y", i), Alg(1, 0)), (("r", None), Alg(-1, 0))]
        return [
            (("x", i), Alg(-1, 0)),
            (("y", i), Alg(-1, 0)),
            (("r", None), Alg(0, -1)),
        ]
    if kind == "rect":
        if k == 0:
            return [(("x", i), Alg(1, 0)), (("r", None), Alg(-1, 0))]
        if k == 1:
            return [(("x", i), Alg(-1, 0)), (("r", None), Alg(-1, 0))]
        if k == 2:
            return [(("y", i), Alg(1, 0)), (("r", None), Alg(-1, 0))]
        return [(("y", i), Alg(-1, 0)), (("r", None), Alg(-1, 0))]
    if k == 0:
        return [(("x", i), Alg(1, 0)), (("r", None), Alg(-1, 0))]
    if k == 1:
        return [(("y", i), Alg(1, 0)), (("r", None), Alg(-1, 0))]
    return [
        (("x", i), Alg(-2 * xi, 0)),
        (("y", i), Alg(-2 * yi, 0)),
        (("r", None), Alg(-2 * (1 - r), 0)),
    ]


def parse_support(S, n: int, kind: str):
    if not isinstance(S, list) or len(S) == 0:
        raise ValueError("empty support")
    out = []
    allowed = wall_ks(kind)
    for raw in S:
        if not isinstance(raw, (list, tuple)) or len(raw) != 3:
            raise ValueError("bad support entry")
        tag, i, j = raw
        if tag not in ("p", "w"):
            raise ValueError("bad support tag")
        if type(i) is bool or type(j) is bool or not isinstance(i, int) or not isinstance(j, int):
            raise ValueError("support indices must be integers")
        if tag == "p":
            if not (0 <= i < n and 0 <= j < n) or i == j:
                raise ValueError("bad pair index")
            out.append(("p", i, j))
        else:
            if not (0 <= i < n) or j not in allowed:
                raise ValueError("bad wall index")
            out.append(("w", i, j))
    return out


def feasible_pairs(centres, r: Fraction):
    den = 1
    for x, y in centres:
        den = math.lcm(den, x.denominator)
        den = math.lcm(den, y.denominator)
    den = math.lcm(den, r.denominator)
    Xs = [x.numerator * (den // x.denominator) for x, _y in centres]
    Ys = [y.numerator * (den // y.denominator) for _x, y in centres]
    R = r.numerator * (den // r.denominator)
    gap = 4 * R * R
    n = len(centres)
    for i in range(n):
        xi = Xs[i]
        yi = Ys[i]
        for j in range(i + 1, n):
            dx = xi - Xs[j]
            dy = yi - Ys[j]
            if dx * dx + dy * dy < gap:
                return f"infeasible pair {i},{j}"
    return None


def feasible_walls(kind, h, r, centres):
    n = len(centres)
    for i in range(n):
        for k in wall_ks(kind):
            if g_value(kind, h, r, centres, ("w", i, k)).is_neg():
                return f"infeasible wall {i} k={k}"
    return None


def float_to_scaled(y: float, E: int) -> int:
    if y == 0.0:
        return 0
    n, d = y.as_integer_ratio()
    k = d.bit_length() - 1
    if d != (1 << k) or E < k:
        raise ValueError("float scaling failed")
    return n << (E - k)


def scaled_columns(Y: np.ndarray):
    """Y[i, a] = Ycols[a][i] / 2^E exactly."""
    nc, ns = Y.shape
    E = 0
    flat = Y.ravel()
    for val in flat:
        if val == 0.0:
            continue
        d = val.as_integer_ratio()[1]
        k = d.bit_length() - 1
        if k > E:
            E = k
    cols = [[0] * nc for _ in range(ns)]
    for a in range(ns):
        col = cols[a]
        base = a
        for i in range(nc):
            val = float(Y[i, base] if False else Y[i, a])
            if val != 0.0:
                col[i] = float_to_scaled(val, E)
    return cols, E


def infinity_norm_I_minus_YG(Y: np.ndarray, rows) -> Fraction:
    """Upper bound on |I - Y G|_inf.

    Y is float64, read as exact dyadics. Each row of G is a list of
    (column, Alg) with Alg = rat + coef*sqrt(2).
    """
    if Y.ndim != 2:
        raise ValueError("Y shape")
    nc, ns = Y.shape
    if len(rows) != ns:
        raise ValueError("G shape")
    if not np.isfinite(Y).all():
        raise ValueError("non-finite left inverse")
    Ycols, E = scaled_columns(Y)
    DG = 1
    any_coef = False
    for entries in rows:
        for _j, alg in entries:
            DG = math.lcm(DG, alg.a.denominator)
            DG = math.lcm(DG, alg.b.denominator)
            if alg.b != 0:
                any_coef = True
    rows_int = []
    for entries in rows:
        packed = []
        for j, alg in entries:
            gr = alg.a.numerator * (DG // alg.a.denominator)
            gc = alg.b.numerator * (DG // alg.b.denominator) if any_coef else 0
            if gr != 0 or gc != 0:
                packed.append((j, gr, gc))
        rows_int.append(packed)
    Nrat = [[0] * nc for _ in range(nc)]
    Ncoef = [[0] * nc for _ in range(nc)] if any_coef else None
    for a, packed in enumerate(rows_int):
        colY = Ycols[a]
        for j, gr, gc in packed:
            if gr:
                acc = Nrat[j]
                for i in range(nc):
                    yi = colY[i]
                    if yi:
                        acc[i] += yi * gr
            if gc:
                acc = Ncoef[j]
                for i in range(nc):
                    yi = colY[i]
                    if yi:
                        acc[i] += yi * gc
    Denom = (1 << E) * DG
    # 40 correct decimals of sqrt(2) make the enclosure gap ~1e-40.
    lo, hi = sqrt2_bounds(40)
    # PATCH (2026-09-23, reported to checker B's author): Fraction() reduces, so lo.denominator != 10**40 (s is even) and
    # lo.numerator was on a different scale -> every row with M < 0 was bounded with sqrt2 ~ lo/2^k (always an OVER-estimate,
    # never a false pass; crt 101: 1.25 instead of 1.27e-14). Use both bounds on the common denominator 10**40.
    N = 10**40
    hi_num = hi.numerator * (N // hi.denominator)
    lo_num = lo.numerator * (N // lo.denominator)
    assert Fraction(hi_num, N) == hi and Fraction(lo_num, N) == lo
    delta_up = Fraction(0)
    for i in range(nc):
        L = 0
        M = 0
        for j in range(nc):
            A = (Denom if i == j else 0) - Nrat[j][i]
            B = -(Ncoef[j][i] if any_coef else 0)
            if B == 0:
                L += A if A >= 0 else -A
            else:
                c, d = int_abs_sqrt2(A, B)
                L += c
                M += d
        if M >= 0:
            num = L * N + M * hi_num
        else:
            num = L * N + M * lo_num
        ub = Fraction(num, N * Denom)
        if ub > delta_up:
            delta_up = ub
    return delta_up


def y_inf_norm(Y: np.ndarray) -> Fraction:
    Ycols, E = scaled_columns(Y)
    nc = Y.shape[0]
    ns = Y.shape[1]
    twoE = 1 << E
    best = Fraction(0)
    for i in range(nc):
        s = 0
        for a in range(ns):
            yi = Ycols[a][i]
            s += -yi if yi < 0 else yi
        cand = Fraction(s, twoE)
        if cand > best:
            best = cand
    return best


def left_inverse(Gf: np.ndarray) -> np.ndarray:
    """Float left inverse of an (ns x nc) matrix, refined a few times."""
    ns, nc = Gf.shape
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Y = np.linalg.pinv(Gf, rcond=np.finfo(np.float64).eps)
    Y = np.asarray(Y, dtype=np.float64)
    if Y.shape != (nc, ns) or not np.isfinite(Y).all():
        raise ValueError("left inverse failed")
    best = Y.copy()
    best_rs = math.inf
    for _ in range(6):
        R = np.eye(nc, dtype=np.float64) - Y @ Gf
        rs = float(np.max(np.sum(np.abs(R), axis=1)))
        if not math.isfinite(rs):
            break
        if rs < best_rs:
            best_rs = rs
            best = Y.copy()
        if rs < 1e-15:
            break
        Yn = Y + R @ Y
        if not np.isfinite(Yn).all():
            break
        Y = Yn
    return best


def prove_K(rows) -> Fraction:
    """Rational K >= |Y|_inf / (1 - |I-YG|_inf) for the float Y read exactly."""
    ns = len(rows)
    if ns == 0:
        raise ValueError("empty jacobian")
    nc = max(j for entries in rows for j, _alg in entries) + 1 if any(rows) else 0
    # Column count is the coordinate count, which may exceed any stored index
    # when a backbone coordinate has an all-zero column. Caller passes nc.
    raise RuntimeError("internal: prove_K needs nc")


def prove_K_nc(rows, nc: int) -> Fraction:
    ns = len(rows)
    Gf = np.zeros((ns, nc), dtype=np.float64)
    s2 = float(math.sqrt(2.0))
    for a, entries in enumerate(rows):
        for j, alg in entries:
            Gf[a, j] = float(alg.a) + s2 * float(alg.b)
    if not np.isfinite(Gf).all():
        raise ValueError("jacobian does not fit in float64")
    Y = left_inverse(Gf)
    delta_up = infinity_norm_I_minus_YG(Y, rows)
    if delta_up >= 1:
        raise ValueError("K not proven")
    yinf = y_inf_norm(Y)
    return yinf / (1 - delta_up)


def rational_bound(alg: Alg, up: bool, places: int) -> Fraction:
    if alg.b == 0:
        return alg.a
    lo, hi = sqrt2_bounds(places)
    if up:
        root = hi if alg.b > 0 else lo
    else:
        root = lo if alg.b > 0 else hi
    return alg.a + alg.b * root


def _scale_q(fr: Fraction, digits: int, up: bool) -> int:
    """Integer q with q / 10^digits on the requested side of fr."""
    if fr == 0:
        return 0
    neg = fr < 0
    afr = -fr if neg else fr
    num = afr.numerator * (10**digits)
    den = afr.denominator
    floor_q = num // den
    ceil_q = floor_q if num % den == 0 else floor_q + 1
    if not neg:
        return ceil_q if up else floor_q
    # ceil(-a) = -floor(a); floor(-a) = -ceil(a)
    return -floor_q if up else -ceil_q


def format_scaled(q: int, digits: int) -> str:
    if q == 0:
        return "0"
    sign = "-" if q < 0 else ""
    s = str(abs(q))
    if digits <= 0:
        return sign + s
    if len(s) <= digits:
        s = s.zfill(digits + 1)
    return sign + s[:-digits] + "." + s[-digits:]


def format_alg(alg: Alg, digits: int, up: bool) -> str:
    if alg.is_zero():
        return "0"
    places = max(20, digits + 20)
    # Extra sqrt digits so the enclosure sits inside one last-place unit
    # whenever the sqrt coefficient is not enormous.
    mag = 0
    if alg.b != 0:
        mag = len(str(abs(alg.b.numerator))) - len(str(alg.b.denominator)) + 2
        if mag < 0:
            mag = 0
    bound = rational_bound(alg, up=up, places=places + mag)
    return format_scaled(_scale_q(bound, digits, up), digits)


def format_certificate(rho: Alg, t0: Alg, delta: Alg, r: Fraction):
    if t0.is_zero() and delta.is_zero() and not rho.is_neg() and not rho.is_zero():
        for digits in (20, 30, 40, 60, 80):
            rho_s = format_alg(rho, digits, up=False)
            if rho_s not in ("0", "-0") and not rho_s.startswith("-"):
                return [f"rho {rho_s}", "t0 0", "Delta 0 Delta/r 0"]
        return [
            f"rho {alg_exact(rho)}",
            "t0 0",
            "Delta 0 Delta/r 0",
        ]
    dr = delta / r
    for digits in (20, 30, 40, 60, 90, 120):
        rho_q = _scale_q(
            rational_bound(rho, up=False, places=digits + 30), digits, up=False
        )
        t0_q = _scale_q(
            rational_bound(t0, up=True, places=digits + 30), digits, up=True
        )
        if t0_q < rho_q:
            d_q = _scale_q(
                rational_bound(delta, up=True, places=digits + 30), digits, up=True
            )
            dr_q = _scale_q(
                rational_bound(dr, up=True, places=digits + 30), digits, up=True
            )
            return [
                f"rho {format_scaled(rho_q, digits)}",
                f"t0 {format_scaled(t0_q, digits)}",
                f"Delta {format_scaled(d_q, digits)} Delta/r {format_scaled(dr_q, digits)}",
            ]
    return [
        f"rho {alg_exact(rho)}",
        f"t0 {alg_exact(t0)}",
        f"Delta {alg_exact(delta)} Delta/r {alg_exact(dr)}",
    ]


def alg_exact(alg: Alg) -> str:
    if alg.b == 0:
        a = alg.a
        if a.denominator == 1:
            return str(a.numerator)
        return f"{a.numerator}/{a.denominator}"
    return f"({alg.a})+({alg.b})*sqrt(2)"


def certify(container: str, r: Fraction, centres, S_raw, lambdas) -> list[str]:
    kind, h = parse_container(container)
    n = len(centres)
    if r <= 0:
        return ["VERDICT: NOT CERTIFIED (r <= 0)"]
    S = parse_support(S_raw, n, kind)
    if len(lambdas) != len(S):
        return ["VERDICT: NOT CERTIFIED (lambda length)"]
    lams: list[Fraction] = []
    for t, lam in enumerate(lambdas):
        if not isinstance(lam, Fraction):
            raise ValueError("internal lambda type")
        if lam <= 0:
            return [f"VERDICT: NOT CERTIFIED (lambda {t} not strictly positive)"]
        lams.append(lam)
    bad = feasible_pairs(centres, r)
    if bad:
        return [f"VERDICT: NOT CERTIFIED ({bad})"]
    bad = feasible_walls(kind, h, r, centres)
    if bad:
        return [f"VERDICT: NOT CERTIFIED ({bad})"]

    backbone = set()
    for tag, i, j in S:
        backbone.add(i)
        if tag == "p":
            backbone.add(j)
    index = {}
    for i in sorted(backbone):
        index[("x", i)] = len(index)
        index[("y", i)] = len(index)
    r_index = len(index)
    index[("r", None)] = r_index
    nc = r_index + 1

    rows = []
    res = [Alg(0, 0) for _ in range(nc)]
    E = Alg(0, 0)
    Q = Fraction(0)
    max_eps = None
    max_kappa = 0
    for con, lam in zip(S, lams):
        eps = g_value(kind, h, r, centres, con)
        if eps.is_neg():
            return ["VERDICT: NOT CERTIFIED (infeasible support constraint)"]
        E = E + eps * lam
        if max_eps is None or max_eps < eps:
            max_eps = eps
        kap = kappa_of(kind, con)
        if kap > max_kappa:
            max_kappa = kap
        Q += lam * kap
        acc = {}
        for key, deriv in gradients(kind, h, r, centres, con):
            col = index[key]
            acc[col] = acc.get(col, Alg(0, 0)) + deriv
            res[col] = res[col] + deriv * lam
        rows.append([(col, acc[col]) for col in sorted(acc) if not acc[col].is_zero()])
    res[r_index] = res[r_index] + Alg(1, 0)

    lmin = min(lams)
    res_l1 = Alg(0, 0)
    for comp in res:
        res_l1 = res_l1 + comp.abs()
    A0 = E / lmin + max_eps
    A1 = res_l1 / lmin
    A2 = Q / lmin + max_kappa
    if A2 < 0:
        return ["VERDICT: NOT CERTIFIED (A2 < 0)"]

    try:
        K = prove_K_nc(rows, nc)
    except ValueError as exc:
        return [f"VERDICT: NOT CERTIFIED ({exc})"]
    if K <= 0:
        return ["VERDICT: NOT CERTIFIED (K not proven)"]

    KA1 = A1 * K
    half = Alg(Fraction(1, 2), 0)
    if not (KA1 - half).is_neg():
        return ["VERDICT: NOT CERTIFIED (K*A1 >= 1/2)"]
    if A2 == 0:
        rho = Alg(1, 0)
    else:
        rho = (half - KA1) / (K * A2)
    t0 = A0 * (K * 2)
    if not (t0 - rho).is_neg():
        return ["VERDICT: NOT CERTIFIED (t0 >= rho)"]
    Delta = E + res_l1 * t0 + (t0 * t0) * Q
    lines = format_certificate(rho, t0, Delta, r)
    lines.append("VERDICT: CERTIFIED")
    return lines


def load_certificate(path: str):
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    while lines and not lines[-1].strip():
        lines.pop()
    if not lines:
        raise ValueError("empty certificate")
    head = lines[0].split()
    if len(head) != 2 or head[0] != "r":
        raise ValueError('line 1 must be "r <decimal>"')
    r = parse_number(head[1])
    centres = []
    for line in lines[1:]:
        if not line.strip():
            raise ValueError("blank centre line")
        parts = line.split()
        if len(parts) != 2:
            raise ValueError("centre line must be 'x y'")
        centres.append((parse_number(parts[0]), parse_number(parts[1])))
    if not centres:
        raise ValueError("no centres")
    return r, centres


def load_json(path: str):
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict) or "S" not in data or "lambda" not in data:
        raise ValueError("json must contain S and lambda")
    raw_lams = data["lambda"]
    if not isinstance(raw_lams, list):
        raise ValueError("lambda must be a list")
    lams = []
    for item in raw_lams:
        if not isinstance(item, str):
            raise ValueError("lambda entries must be float.hex strings")
        lams.append(parse_hex_float(item))
    return data["S"], lams


def _self_test() -> None:
    if parse_hex_float("0x1.8p+0") != Fraction(3, 2):
        raise AssertionError("hex 1.5")
    if parse_hex_float("-0x1.0p+2") != -4:
        raise AssertionError("hex -4")
    if parse_hex_float("0x0.0p+0") != 0:
        raise AssertionError("hex 0")
    if parse_hex_float((0.1).hex()) != Fraction.from_float(0.1):
        raise AssertionError("hex roundtrip")
    if parse_decimal("0.1") != Fraction(1, 10):
        raise AssertionError("decimal 0.1")
    if parse_decimal("-1e-2") != Fraction(-1, 100):
        raise AssertionError("decimal sci")
    if parse_number("1/2") != Fraction(1, 2):
        raise AssertionError("fraction token")
    z = Alg(1, -1).abs()
    if z.a != -1 or z.b != 1:
        raise AssertionError("abs(1-sqrt(2))")
    if not Alg(Fraction(1, 4), Fraction(-1, 10)).abs().a:
        raise AssertionError("tri slack sign")
    if not Alg(Fraction(-17, 10), Fraction(-1, 10)).is_neg():
        raise AssertionError("outside hypotenuse")
    lo, hi = sqrt2_bounds(30)
    if not (lo * lo <= 2 < hi * hi):
        raise AssertionError("sqrt2 bounds")
    # One-entry residual: Y ≈ -1/sqrt(2), G = -sqrt(2), so YG ≈ 1.
    y = float(-math.sqrt(2.0) / 2.0)
    delta = infinity_norm_I_minus_YG(
        np.array([[y]], dtype=np.float64),
        [[(0, Alg(0, -1))]],
    )
    if delta >= Fraction(1, 1000):
        raise AssertionError(f"sqrt2 residual {delta}")
    if kappa_of("rect", ("p", 0, 1)) != 8:
        raise AssertionError("kappa pair")
    if kappa_of("quad", ("w", 0, 2)) != 2:
        raise AssertionError("kappa arc")
    if kappa_of("tri", ("w", 0, 2)) != 0 or kappa_of("quad", ("w", 0, 0)) != 0:
        raise AssertionError("kappa wall")

    # Two equal circles of radius 1/2 in a 1 x 2 rectangle, stacked.
    # Active set: both side walls of each circle, bottom of the lower,
    # top of the upper, and the pair. Multipliers are the exact KKT values.
    r = Fraction(1, 2)
    centres = [(Fraction(0), Fraction(-1, 2)), (Fraction(0), Fraction(1, 2))]
    S = [
        ["w", 0, 0],
        ["w", 0, 1],
        ["w", 0, 2],
        ["w", 1, 0],
        ["w", 1, 1],
        ["w", 1, 3],
        ["p", 0, 1],
    ]
    lam = [
        Fraction(1, 12),
        Fraction(1, 12),
        Fraction(1, 6),
        Fraction(1, 12),
        Fraction(1, 12),
        Fraction(1, 6),
        Fraction(1, 12),
    ]
    lines = certify("rect:2", r, centres, S, lam)
    if (
        lines[-1] != "VERDICT: CERTIFIED"
        or lines[1] != "t0 0"
        or lines[2] != "Delta 0 Delta/r 0"
    ):
        raise AssertionError("two-circle self-test failed: " + " | ".join(lines))
    bad = certify("rect:2", Fraction(3, 5), centres, S, lam)
    if "infeasible" not in bad[-1]:
        raise AssertionError(bad[-1])
    lam_neg = list(lam)
    lam_neg[0] = -lam_neg[0]
    neg = certify("rect:2", r, centres, S, lam_neg)
    if "lambda" not in neg[-1]:
        raise AssertionError(neg[-1])


def main(argv: list[str]) -> None:
    try:
        _self_test()
    except Exception as exc:
        print(f"VERDICT: NOT CERTIFIED (self-test failed: {exc})")
        return
    if len(argv) != 4:
        print(
            "VERDICT: NOT CERTIFIED "
            "(usage: verify_lopt.py <container> <certificate.txt> <cert.json>)"
        )
        return
    try:
        r, centres = load_certificate(argv[2])
        S, lams = load_json(argv[3])
        lines = certify(argv[1], r, centres, S, lams)
    except (OSError, json.JSONDecodeError, ValueError, KeyError, ZeroDivisionError) as exc:
        print(f"VERDICT: NOT CERTIFIED ({exc})")
        return
    for line in lines:
        print(line)


if __name__ == "__main__":
    main(sys.argv)
