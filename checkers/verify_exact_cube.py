#!/usr/bin/env python3
# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Exact checker for equal spheres in the unit cube (Packomania scu).

Cube side 1, centre at the origin. One shared radius. Decisions are integers
only: each decimal is n/10^d, then scaled by 2*10^D so 1/2 lies on the grid.

N <= 1500: every pair. N > 1500: hash grid of cell side exactly 2r.
Cell k is the half-open interval [k*L, (k+1)*L). If two indices differ by
2 or more on one axis, those integer coordinates differ by at least L+1.
L = 2*R is the scaled length 2r, so the gap is > 2r and the squared
distance is > (2r)^2. That pair cannot be an overlap. Every real
overlap sits in the same cell or the 3x3x3 neighbour block and is checked.
Tangency (squared distance == (2r)^2) is legal, including when a far
cell is not visited.

Usage:
    python verify_exact_cube.py <file> <N> <record_radius>
No arguments: self-tests.
"""

from __future__ import annotations

import os
import sys
import tempfile
from collections import defaultdict

ALL_PAIRS_MAX_N = 1500
_MAX_POW10 = 1_000_000


def pow10(k: int) -> int:
    if k < 0 or k > _MAX_POW10:
        raise ValueError("decimal scale out of range")
    return 10**k


def _ascii_digits(s: str) -> bool:
    return bool(s) and all("0" <= ch <= "9" for ch in s)


def _normalize(n: int, d: int) -> tuple[int, int]:
    if n == 0:
        return 0, 0
    while d > 0 and n % 10 == 0:
        n //= 10
        d -= 1
    return n, d


def parse_decimal(s: str) -> tuple[int, int]:
    """Parse a decimal or scientific string to n/10^d. No float."""
    t = s.strip()
    if not t:
        raise ValueError("empty number")
    sign = 1
    if t[0] == "+":
        t = t[1:]
    elif t[0] == "-":
        sign = -1
        t = t[1:]
    if not t:
        raise ValueError(f"bad number: {s!r}")
    exp = 0
    low = t.lower()
    if "e" in low:
        idx = low.index("e")
        mant = t[:idx]
        exp_s = t[idx + 1 :]
        if exp_s[:1] in "+-":
            exp_body = exp_s[1:]
        else:
            exp_body = exp_s
        if not _ascii_digits(exp_body):
            raise ValueError(f"bad number: {s!r}")
        exp = int(exp_s)
        t = mant
        if not t:
            raise ValueError(f"bad number: {s!r}")
    if t.count(".") > 1:
        raise ValueError(f"bad number: {s!r}")
    if "." in t:
        whole, frac = t.split(".", 1)
        if whole == "" and frac == "":
            raise ValueError(f"bad number: {s!r}")
        if whole == "":
            whole = "0"
        if not _ascii_digits(whole) or (frac != "" and not _ascii_digits(frac)):
            raise ValueError(f"bad number: {s!r}")
        digits = whole + frac
        scale = len(frac)
    else:
        if not _ascii_digits(t):
            raise ValueError(f"bad number: {s!r}")
        digits = t
        scale = 0
    n = sign * int(digits)
    shift = exp - scale
    if shift >= 0:
        return _normalize(n * pow10(shift), 0)
    return _normalize(n, -shift)


def compare(a: tuple[int, int], b: tuple[int, int]) -> int:
    """Compare n/10^d values. Negative if a < b, zero if equal, else positive."""
    an, ad = a
    bn, bd = b
    left = an * pow10(bd)
    right = bn * pow10(ad)
    if left < right:
        return -1
    if left > right:
        return 1
    return 0


def read_certificate(path: str, n_expected: int) -> tuple[tuple[int, int], list]:
    if n_expected < 0:
        raise ValueError("N is negative")
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    if not lines:
        raise ValueError("empty file")
    head = lines[0].strip().split()
    if len(head) != 2 or head[0] != "r":
        raise ValueError("first line must be 'r <radius>'")
    radius = parse_decimal(head[1])
    body = lines[1:]
    if len(body) != n_expected:
        raise ValueError(
            f"count mismatch: file has {len(body)} centers, expected {n_expected}"
        )
    pts = []
    for k, line in enumerate(body, start=2):
        parts = line.strip().split()
        if len(parts) != 3:
            raise ValueError(f"line {k} does not have 3 coordinates")
        pts.append(
            (parse_decimal(parts[0]), parse_decimal(parts[1]), parse_decimal(parts[2]))
        )
    return radius, pts


def _scale(num: int, den_pow: int, depth: int) -> int:
    # (num/10^den_pow) * 2 * 10^depth
    return num * 2 * pow10(depth - den_pow)


def _scaled(radius: tuple[int, int], pts: list) -> tuple[int, int, list]:
    depth = radius[1]
    for p in pts:
        for c in p:
            if c[1] > depth:
                depth = c[1]
    r_scaled = _scale(radius[0], radius[1], depth)
    half = pow10(depth)  # scaled 1/2
    out = []
    for x, y, z in pts:
        out.append(
            (
                _scale(x[0], x[1], depth),
                _scale(y[0], y[1], depth),
                _scale(z[0], z[1], depth),
            )
        )
    return r_scaled, half, out


def _outside(pts: list, r_scaled: int, half: int):
    lo = -half + r_scaled
    hi = half - r_scaled
    for i, (x, y, z) in enumerate(pts):
        if x < lo or x > hi:
            return i, "x"
        if y < lo or y > hi:
            return i, "y"
        if z < lo or z > hi:
            return i, "z"
    return None


def _conflict(pts: list, r_scaled: int, how: str):
    """Smallest (i, j) with squared distance strictly under (2r)^2, or None.

    how is 'all' or 'grid'. See the module docstring for the grid argument.
    """
    if r_scaled == 0:
        return None
    if r_scaled < 0:
        raise ValueError("negative scaled radius")
    limit = 4 * r_scaled * r_scaled
    n = len(pts)
    if how == "all":
        for i in range(n):
            x, y, z = pts[i]
            for j in range(i + 1, n):
                dx = x - pts[j][0]
                dy = y - pts[j][1]
                dz = z - pts[j][2]
                if dx * dx + dy * dy + dz * dz < limit:
                    return i, j
        return None
    if how != "grid":
        raise ValueError(f"unknown method: {how}")
    cell = 2 * r_scaled  # scaled 2r; world cell side is exactly 2r
    buckets: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    keys = []
    for i, (x, y, z) in enumerate(pts):
        key = (x // cell, y // cell, z // cell)
        keys.append(key)
        buckets[key].append(i)
    for i, key in enumerate(keys):
        x, y, z = pts[i]
        ix, iy, iz = key
        best = None
        for ox in (-1, 0, 1):
            for oy in (-1, 0, 1):
                for oz in (-1, 0, 1):
                    for j in buckets.get((ix + ox, iy + oy, iz + oz), ()):
                        if j <= i:
                            continue
                        if best is not None and j >= best:
                            continue
                        dx = x - pts[j][0]
                        dy = y - pts[j][1]
                        dz = z - pts[j][2]
                        if dx * dx + dy * dy + dz * dz < limit:
                            best = j
        if best is not None:
            return i, best
    return None


def resolve_method(n: int, method: str) -> str:
    if method == "auto":
        return "all" if n <= ALL_PAIRS_MAX_N else "grid"
    if method in ("all", "grid"):
        return method
    raise ValueError(f"unknown method: {method}")


def judge(radius: tuple[int, int], pts: list, record: tuple[int, int], method: str = "auto") -> str:
    if radius[0] < 0:
        return "VERDICT: INVALID negative radius"
    how = resolve_method(len(pts), method)
    r_scaled, half, spts = _scaled(radius, pts)
    bad = _outside(spts, r_scaled, half)
    if bad is not None:
        i, axis = bad
        return f"VERDICT: INVALID sphere {i} {axis} outside cube"
    hit = _conflict(spts, r_scaled, how)
    if hit is not None:
        i, j = hit
        return f"VERDICT: INVALID overlap between spheres {i} and {j}"
    if compare(radius, record) > 0:
        return "VERDICT: IMPROVES"
    return "VERDICT: VALID_NOT_BETTER"


def verdict_for(path: str, n: int, record: str, method: str = "auto") -> str:
    try:
        radius, pts = read_certificate(path, n)
        return judge(radius, pts, parse_decimal(record), method)
    except (ValueError, OSError) as exc:
        return f"VERDICT: INVALID {exc}"


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(
            "usage: python verify_exact_cube.py <file> <N> <record_radius>",
            file=sys.stderr,
        )
        return 2
    path, n_text, record_text = argv[1], argv[2], argv[3]
    try:
        if not _ascii_digits(n_text):
            raise ValueError(f"N is not a non-negative integer: {n_text}")
        n = int(n_text)
        record = parse_decimal(record_text)
        radius, pts = read_certificate(path, n)
        print(judge(radius, pts, record, "auto"))
    except (ValueError, OSError) as exc:
        print(f"VERDICT: INVALID {exc}")
    return 0


def dec_frac(num: int, places: int) -> str:
    """Exact decimal for num/10^places. Used by self-tests only."""
    if places < 0:
        raise ValueError("negative places")
    sign = "-" if num < 0 else ""
    v = -num if num < 0 else num
    if places == 0:
        return sign + str(v)
    digits = str(v)
    if len(digits) <= places:
        return sign + "0." + digits.zfill(places)
    return sign + digits[:-places] + "." + digits[-places:]


def _write_cert(directory: str, name: str, radius: str, rows: list[str]) -> str:
    path = os.path.join(directory, name)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"r {radius}\n")
        for row in rows:
            fh.write(row + "\n")
    return path


def _fail(note: str) -> None:
    raise SystemExit(f"SELF-TEST FAIL: {note}")


def _expect_eq(got: str, want: str, note: str) -> None:
    if got != want:
        _fail(f"{note}: got {got!r} want {want!r}")


def _expect_invalid(got: str, needle: str, note: str) -> None:
    if not got.startswith("VERDICT: INVALID ") or needle not in got:
        _fail(f"{note}: got {got!r}")


def _corner_rows() -> list[str]:
    rows = []
    for x in ("-0.25", "0.25"):
        for y in ("-0.25", "0.25"):
            for z in ("-0.25", "0.25"):
                rows.append(f"{x} {y} {z}")
    return rows


def _lattice_rows(overlap: bool) -> list[str]:
    """10x10x10 lattice on [-0.45, 0.45]^3, spacing 0.1.

    overlap=True moves the x-neighbour of the first point to distance
    2*(1e-6) - 1e-20 from it. Other gaps stay 0.1.
    """
    moved = dec_frac(-(449998 * 10**14 + 1), 20)
    rows = []
    for i in range(10):
        for j in range(10):
            for k in range(10):
                xs = dec_frac(10 * i - 45, 2)
                ys = dec_frac(10 * j - 45, 2)
                zs = dec_frac(10 * k - 45, 2)
                if overlap and (i, j, k) == (1, 0, 0):
                    xs = moved
                rows.append(f"{xs} {ys} {zs}")
    return rows


def run_self_tests() -> None:
    if dec_frac(-45, 2) != "-0.45" or dec_frac(45, 2) != "0.45" or dec_frac(0, 2) != "0.00":
        _fail("dec_frac lattice endpoints")
    if compare(parse_decimal("0.50"), parse_decimal("0.5")) != 0:
        _fail("0.50 vs 0.5")
    if compare(parse_decimal("-0.250"), parse_decimal("-0.25")) != 0:
        _fail("-0.250 vs -0.25")
    sci = "0." + ("0" * 16) + "15"  # 15/10^18 = 1.5e-17
    if compare(parse_decimal("1.5e-17"), parse_decimal(sci)) != 0:
        _fail("scientific parse")
    if parse_decimal("1e-6") != (1, 6):
        _fail("1e-6")

    r_lo = dec_frac(5 * 10**29 - 1, 30)  # 1/2 - 1e-30
    r_hi = dec_frac(5 * 10**29 + 1, 30)  # 1/2 + 1e-30
    if parse_decimal(r_lo) != (5 * 10**29 - 1, 30):
        _fail(f"r_lo parsed {parse_decimal(r_lo)!r}")
    if parse_decimal(r_hi) != (5 * 10**29 + 1, 30):
        _fail(f"r_hi parsed {parse_decimal(r_hi)!r}")
    if not compare(parse_decimal(r_lo), parse_decimal("0.5")) < 0:
        _fail("r_lo side of 1/2")
    if not compare(parse_decimal(r_hi), parse_decimal("0.5")) > 0:
        _fail("r_hi side of 1/2")

    quarter_hi = dec_frac(10**30 // 4 + 1, 30)  # 1/4 + 1e-30
    if parse_decimal(quarter_hi) != (10**30 // 4 + 1, 30):
        _fail("quarter_hi")
    if parse_decimal(dec_frac(10**30 // 4, 30)) != parse_decimal("0.25"):
        _fail("quarter normalises to 0.25")
    x_push = dec_frac(25 * 10**23 + 1, 25)  # 1/4 + 1e-25
    if parse_decimal(x_push) != (25 * 10**23 + 1, 25):
        _fail("x_push")
    if not compare(parse_decimal(x_push), parse_decimal("0.25")) > 0:
        _fail("x_push side of 1/4")
    moved = dec_frac(-(449998 * 10**14 + 1), 20)
    if parse_decimal(moved) != (-(449998 * 10**14 + 1), 20):
        _fail("moved x")

    with tempfile.TemporaryDirectory() as td:
        p1 = _write_cert(td, "n1.txt", r_lo, ["0 0 0"])
        _expect_eq(verdict_for(p1, 1, "0"), "VERDICT: IMPROVES", "n1 improves")
        _expect_eq(verdict_for(p1, 1, r_lo), "VERDICT: VALID_NOT_BETTER", "n1 equal record")
        _expect_eq(verdict_for(p1, 1, "0.5"), "VERDICT: VALID_NOT_BETTER", "n1 record above")
        code = main(["verify_exact_cube.py", p1, "1", "0"])
        if code != 0:
            _fail(f"main exit {code}")

        p1_bad = _write_cert(td, "n1_out.txt", r_hi, ["0 0 0"])
        _expect_invalid(verdict_for(p1_bad, 1, "0"), "outside", "n1 radius past half")

        corners = _corner_rows()
        if len(corners) != 8:
            _fail("corner count")
        p8 = _write_cert(td, "oct.txt", "0.25", corners)
        _expect_eq(verdict_for(p8, 8, "0.25", "all"), "VERDICT: VALID_NOT_BETTER", "oct all")
        _expect_eq(verdict_for(p8, 8, "0.25", "grid"), "VERDICT: VALID_NOT_BETTER", "oct grid")
        _expect_eq(verdict_for(p8, 8, "0", "auto"), "VERDICT: IMPROVES", "oct improves")

        p8_big = _write_cert(td, "oct_big.txt", quarter_hi, corners)
        big_all = verdict_for(p8_big, 8, "0", "all")
        big_grid = verdict_for(p8_big, 8, "0", "grid")
        _expect_invalid(big_all, "INVALID", "oct bigger all")
        _expect_invalid(big_grid, "INVALID", "oct bigger grid")
        if not big_all.startswith("VERDICT: INVALID") or not big_grid.startswith("VERDICT: INVALID"):
            _fail("oct bigger verdict")

        pushed_rows = []
        replaced = False
        for row in corners:
            if row == "0.25 0.25 0.25":
                pushed_rows.append(f"{x_push} 0.25 0.25")
                replaced = True
            else:
                pushed_rows.append(row)
        if not replaced:
            _fail("push target missing")
        p_wall = _write_cert(td, "wall.txt", "0.25", pushed_rows)
        _expect_invalid(verdict_for(p_wall, 8, "0", "all"), "outside", "wall all")
        _expect_invalid(verdict_for(p_wall, 8, "0", "grid"), "outside", "wall grid")

        _expect_invalid(verdict_for(p8, 7, "0"), "count", "short count")
        _expect_invalid(verdict_for(p1, 2, "0"), "count", "long count")

        rows_ok = _lattice_rows(False)
        rows_bad = _lattice_rows(True)
        if len(rows_ok) != 1000 or len(rows_bad) != 1000:
            _fail("lattice length")
        p_ok = _write_cert(td, "lat.txt", "1e-6", rows_ok)
        p_bad = _write_cert(td, "lat_hit.txt", "1e-6", rows_bad)
        ok_all = verdict_for(p_ok, 1000, "0", "all")
        ok_grid = verdict_for(p_ok, 1000, "0", "grid")
        ok_auto = verdict_for(p_ok, 1000, "0", "auto")
        _expect_eq(ok_all, "VERDICT: IMPROVES", "lattice all")
        _expect_eq(ok_grid, "VERDICT: IMPROVES", "lattice grid")
        _expect_eq(ok_auto, "VERDICT: IMPROVES", "lattice auto")
        bad_all = verdict_for(p_bad, 1000, "0", "all")
        bad_grid = verdict_for(p_bad, 1000, "0", "grid")
        _expect_invalid(bad_all, "overlap", "lattice hit all")
        _expect_invalid(bad_grid, "overlap", "lattice hit grid")
        _expect_eq(bad_all, bad_grid, "lattice hit same pair")

    print("SELF-TEST PASS")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        run_self_tests()
    else:
        sys.exit(main(sys.argv))
