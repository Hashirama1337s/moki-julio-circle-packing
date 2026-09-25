"""Exact checker for equal-circle packings. Python 3.11, standard library only.

Arithmetic that decides the verdict stays in the integers. Decimal strings are
scaled by one common power of ten. No float is used.
"""

from __future__ import annotations

import sys
import tempfile
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path


def _digits(text: str) -> bool:
    return bool(text) and all(c in "0123456789" for c in text)


def parse_decimal(token: str) -> tuple[int, int]:
    """Return (numerator, places) with value = numerator / 10**places."""
    s = token.strip()
    if not s:
        raise ValueError("format")
    sign = 1
    if s[0] == "+":
        s = s[1:]
    elif s[0] == "-":
        sign = -1
        s = s[1:]
    if not s or s == "." or s.count(".") > 1:
        raise ValueError("format")
    if "." in s:
        whole, frac = s.split(".")
        if whole == "":
            whole = "0"
        if frac == "":
            digits = whole
            places = 0
        else:
            digits = whole + frac
            places = len(frac)
    else:
        digits = s
        places = 0
    if not _digits(digits):
        raise ValueError("format")
    return sign * int(digits), places


def format_decimal(numerator: int, places: int) -> str:
    """Inverse of parse_decimal for an integer numerator and a power of ten."""
    if places < 0:
        raise ValueError("format")
    sign = "-" if numerator < 0 else ""
    magnitude = -numerator if numerator < 0 else numerator
    if places == 0:
        return sign + str(magnitude)
    scale = 10 ** places
    whole = magnitude // scale
    frac = magnitude % scale
    return f"{sign}{whole}.{frac:0{places}d}"


def ratio_gt(n1: int, p1: int, n2: int, p2: int) -> bool:
    """True iff n1/10**p1 > n2/10**p2. Both denominators are powers of ten."""
    return n1 * (10 ** p2) > n2 * (10 ** p1)


def parse_lines(lines: list[str]) -> tuple[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int]]]]:
    if not lines:
        raise ValueError("format")
    head = lines[0].strip().split()
    if len(head) != 2 or head[0] != "r":
        raise ValueError("format")
    radius = parse_decimal(head[1])
    points: list[tuple[tuple[int, int], tuple[int, int]]] = []
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) != 2:
            raise ValueError("format")
        points.append((parse_decimal(parts[0]), parse_decimal(parts[1])))
    return radius, points


def scale_all(
    radius: tuple[int, int],
    points: list[tuple[tuple[int, int], tuple[int, int]]],
) -> tuple[int, list[tuple[int, int]], int]:
    """Scale radius and centres by 10**max_places. Return (r, points, scale)."""
    max_places = radius[1]
    for (x, y) in points:
        if x[1] > max_places:
            max_places = x[1]
        if y[1] > max_places:
            max_places = y[1]
    factor: dict[int, int] = {}

    def convert(numerator: int, places: int) -> int:
        cached = factor.get(places)
        if cached is None:
            cached = 10 ** (max_places - places)
            factor[places] = cached
        return numerator * cached

    scaled = [(convert(x[0], x[1]), convert(y[0], y[1])) for (x, y) in points]
    return convert(radius[0], radius[1]), scaled, 10 ** max_places


def contained(container: str, pts: list[tuple[int, int]], r_scaled: int, scale: int) -> bool:
    """Square: |x| + r <= 1/2 and |y| + r <= 1/2.
    Disc: r <= 1 and x^2 + y^2 <= (1 - r)^2.
    """
    if container == "square":
        for x, y in pts:
            if 2 * (abs(x) + r_scaled) > scale:
                return False
            if 2 * (abs(y) + r_scaled) > scale:
                return False
        return True
    if container == "circle":
        if r_scaled > scale:
            return False
        inner = scale - r_scaled
        inner2 = inner * inner
        for x, y in pts:
            if x * x + y * y > inner2:
                return False
        return True
    return False


def pairs_separated(pts: list[tuple[int, int]], r_scaled: int) -> bool:
    """Return True when every pair has dx^2 + dy^2 >= 4 r^2.

    Cell side is L = 2 r in the scaled integers. Cell (i, j) is the set of
    centres with x // L == i and y // L == j. Python // is mathematical floor,
    so negative centres are indexed correctly. Geometrically that cell is the
    half-open square [iL, (i+1)L) x [jL, (j+1)L). On this integer lattice the
    x-coordinates in the cell are the integers iL <= x <= (i+1)L - 1.

    If two cells have x-indices at least 2 apart, take i2 >= i1 + 2. Then
    x1 <= (i1+1)*L - 1 and x2 >= i2*L, so x2 - x1 >= L + 1. The same bound
    holds for y. Whenever the larger index gap is at least 2,
    dx^2 + dy^2 >= (L+1)^2 > L^2 = 4 r^2, so separation already holds.
    In real coordinates the gap is already strict: a centre in cell i has
    x < (i+1)*L and a centre in cell i+2 has x >= (i+2)*L, so |dx| > L.

    A pair can overlap, or touch at distance exactly 2 r, only when its cells
    are the same cell or one of the eight edge and corner neighbours. Equality
    is not pushed into the exterior, because distance exactly L is not greater
    than L. The scan compares every pair in that 3x3 block once (higher index
    against lower index). Pairs outside the block are already strictly
    separated, so leaving them out cannot miss a violation or a tangency.
    """
    if r_scaled <= 0:
        return False
    length = r_scaled * 2
    minimum = length * length
    buckets: dict[tuple[int, int], list[tuple[int, int, int]]] = {}
    for index, (x, y) in enumerate(pts):
        key = (x // length, y // length)
        row = buckets.get(key)
        if row is None:
            buckets[key] = [(index, x, y)]
        else:
            row.append((index, x, y))
    for (ix, iy), members in list(buckets.items()):
        for ox in (-1, 0, 1):
            for oy in (-1, 0, 1):
                others = buckets.get((ix + ox, iy + oy))
                if others is None:
                    continue
                for i, xi, yi in members:
                    for j, xj, yj in others:
                        if j <= i:
                            continue
                        dx = xi - xj
                        dy = yi - yj
                        if dx * dx + dy * dy < minimum:
                            return False
    return True


def evaluate(
    container: str,
    radius: tuple[int, int],
    points: list[tuple[tuple[int, int], tuple[int, int]]],
    count: int,
    record: tuple[int, int],
) -> str:
    if container not in ("square", "circle"):
        return "VERDICT: INVALID arguments"
    if len(points) != count:
        return "VERDICT: INVALID count"
    if radius[0] <= 0:
        return "VERDICT: INVALID radius"
    r_scaled, scaled, scale = scale_all(radius, points)
    if r_scaled <= 0:
        return "VERDICT: INVALID radius"
    if not contained(container, scaled, r_scaled, scale):
        return "VERDICT: INVALID containment"
    if not pairs_separated(scaled, r_scaled):
        return "VERDICT: INVALID separation"
    if ratio_gt(radius[0], radius[1], record[0], record[1]):
        return "VERDICT: IMPROVES"
    return "VERDICT: VALID_NOT_BETTER"


def verdict_for(container: str, lines: list[str], count: int, record_token: str) -> str:
    try:
        record = parse_decimal(record_token)
    except ValueError:
        return "VERDICT: INVALID arguments"
    try:
        radius, points = parse_lines(lines)
    except ValueError:
        return "VERDICT: INVALID format"
    return evaluate(container, radius, points, count, record)


def verdict_from_file(container: str, path: str, count: int, record_token: str) -> str:
    try:
        record = parse_decimal(record_token)
    except ValueError:
        return "VERDICT: INVALID arguments"
    try:
        with open(path, "r", encoding="utf-8-sig") as handle:
            lines = handle.read().splitlines()
    except OSError:
        return "VERDICT: INVALID file"
    try:
        radius, points = parse_lines(lines)
    except ValueError:
        return "VERDICT: INVALID format"
    return evaluate(container, radius, points, count, record)


def _naive_separated(pts: list[tuple[int, int]], r_scaled: int) -> bool:
    minimum = (2 * r_scaled) * (2 * r_scaled)
    total = len(pts)
    for i in range(total):
        xi, yi = pts[i]
        for j in range(i + 1, total):
            dx = xi - pts[j][0]
            dy = yi - pts[j][1]
            if dx * dx + dy * dy < minimum:
                return False
    return True


def _facing(a: int, b: int, length: int) -> tuple[int, int]:
    """Closest lattice coordinates of cells a and b along one axis."""
    if a == b:
        low = a * length
        return low, low
    if b > a:
        return (a + 1) * length - 1, b * length
    return a * length, (b + 1) * length - 1


def _eq(failures: list[str], label: str, got: str, want: str) -> None:
    if got != want:
        failures.append(f"{label}: {got}")


def test_parse(failures: list[str]) -> None:
    samples = (
        (0, 0),
        (5, 0),
        (-5, 0),
        (25, 2),
        (-20, 2),
        (1, 40),
        (-1, 40),
        (25 * (10 ** 38) + 1, 40),
        (-(25 * (10 ** 38) + 1), 40),
        (25 * (10 ** 38) - 1, 40),
        (9 * (10 ** 39) + 1, 40),
        (9 * (10 ** 39) - 1, 40),
        (10 ** 39 + 1, 40),
    )
    for numerator, places in samples:
        token = format_decimal(numerator, places)
        got = parse_decimal(token)
        if got != (numerator, places):
            failures.append(f"roundtrip {token} -> {got}")
    expected = {
        "+0.5": (5, 1),
        "-.5": (-5, 1),
        "5.": (5, 0),
        "00.10": (10, 2),
        "-0.20": (-20, 2),
    }
    for token, want in expected.items():
        got = parse_decimal(token)
        if got != want:
            failures.append(f"parse {token} -> {got}")
    for token in ("1e-40", "", " ", ".", "-.", "+", "1.2.3", "--1", "1a", "1,2"):
        try:
            parse_decimal(token)
        except ValueError:
            continue
        failures.append(f"accepted {token!r}")


def test_compare(failures: list[str]) -> None:
    rows = (
        ((1, 1), (1, 1), False),
        ((1, 1), (10, 2), False),
        ((1001, 4), (1, 1), True),
        ((1, 1), (1001, 4), False),
        ((25, 2), (2499, 4), True),
        ((25, 2), (2500, 4), False),
        ((-1, 0), (0, 0), False),
        ((1, 0), (-5, 1), True),
    )
    for left, right, want in rows:
        got = ratio_gt(left[0], left[1], right[0], right[1])
        if got != want:
            failures.append(f"compare {left} > {right} -> {got}")


def test_verdicts(failures: list[str]) -> None:
    cases = (
        ("sq tangent =", "square", "r 0.25\n-0.25 0\n0.25 0\n", 2, "0.25", "VERDICT: VALID_NOT_BETTER"),
        ("sq tangent >", "square", "r 0.25\n-0.25 0\n0.25 0\n", 2, "0.2", "VERDICT: IMPROVES"),
        ("sq tangent <", "square", "r 0.25\n-0.25 0\n0.25 0\n", 2, "0.9", "VERDICT: VALID_NOT_BETTER"),
        ("sq spelling", "square", "r 0.25\n-0.25 0\n0.25 0\n", 2, "0.2500", "VERDICT: VALID_NOT_BETTER"),
        ("sq 0.2499", "square", "r 0.25\n-0.25 0\n0.25 0\n", 2, "0.2499", "VERDICT: IMPROVES"),
        ("sq corner", "square", "r 0.25\n0.25 0.25\n", 1, "0.25", "VERDICT: VALID_NOT_BETTER"),
        ("sq r big", "square", "r 0.6\n0 0\n", 1, "0.1", "VERDICT: INVALID containment"),
        ("tab header", "square", "r\t0.25\n-0.25   0\n0.25 0\n", 2, "0.25", "VERDICT: VALID_NOT_BETTER"),
        ("cci wall", "circle", "r 0.2\n-0.8 0\n-0.4 0\n", 2, "0.2", "VERDICT: VALID_NOT_BETTER"),
        ("cci 345", "circle", "r 0.2\n0.48 0.64\n", 1, "0.2", "VERDICT: VALID_NOT_BETTER"),
        ("cci overlap", "circle", "r 0.2\n0 0\n0.1 0\n", 2, "0.2", "VERDICT: INVALID separation"),
        ("cci r>1", "circle", "r 1.5\n0 0\n", 1, "0.1", "VERDICT: INVALID containment"),
        ("cci r=1", "circle", "r 1\n0 0\n", 1, "0.5", "VERDICT: IMPROVES"),
        ("cci r=1 off", "circle", "r 1\n0.1 0\n", 1, "0.5", "VERDICT: INVALID containment"),
        ("cci 1.000", "circle", "r 1.000\n0 0\n", 1, "1", "VERDICT: VALID_NOT_BETTER"),
        ("same cell", "square", "r 0.1\n0 0\n0.1 0\n", 2, "0.1", "VERDICT: INVALID separation"),
        ("adj +", "square", "r 0.1\n0.39 0\n0.40 0\n", 2, "0.1", "VERDICT: INVALID separation"),
        ("adj -", "square", "r 0.1\n-0.21 0\n-0.20 0\n", 2, "0.1", "VERDICT: INVALID separation"),
        ("adj y", "square", "r 0.1\n0 0.39\n0 0.40\n", 2, "0.1", "VERDICT: INVALID separation"),
        ("diag", "square", "r 0.05\n0.09 0.09\n0 0\n0.10 0.10\n", 3, "0.05", "VERDICT: INVALID separation"),
        ("neg tangent", "square", "r 0.1\n-0.2 0\n0 0\n", 2, "0.1", "VERDICT: VALID_NOT_BETTER"),
        ("far cells", "square", "r 0.05\n0 0\n0.2 0\n", 2, "0.05", "VERDICT: VALID_NOT_BETTER"),
        ("dup", "square", "r 0.1\n0 0\n0 0\n", 2, "0.1", "VERDICT: INVALID separation"),
        ("count", "square", "r 0.1\n0 0\n", 2, "0.1", "VERDICT: INVALID count"),
        ("n0", "square", "r 0.2\n", 0, "0.1", "VERDICT: IMPROVES"),
        ("header", "square", "R 0.1\n0 0\n", 1, "0.1", "VERDICT: INVALID format"),
        ("extra", "square", "r 0.1\n0 0 1\n", 1, "0.1", "VERDICT: INVALID format"),
        ("blank", "square", "r 0.1\n\n0 0\n", 1, "0.1", "VERDICT: INVALID format"),
        ("zero r", "square", "r 0\n0 0\n", 1, "0.1", "VERDICT: INVALID radius"),
        ("neg r", "square", "r -0.2\n0 0\n", 1, "0.1", "VERDICT: INVALID radius"),
        ("origin =", "square", "r 0.1\n0 0\n", 1, "0.1000", "VERDICT: VALID_NOT_BETTER"),
        ("origin >", "square", "r 0.1\n0 0\n", 1, "0.0999", "VERDICT: IMPROVES"),
        ("bad record", "square", "r 0.1\n0 0\n", 1, "1e-3", "VERDICT: INVALID arguments"),
        ("bad box", "hex", "r 0.1\n0 0\n", 1, "0.1", "VERDICT: INVALID arguments"),
    )
    for label, container, text, count, record, want in cases:
        _eq(failures, label, verdict_for(container, text.splitlines(), count, record), want)

    places = 40
    plus = 25 * (10 ** 38) + 1
    minus = 25 * (10 ** 38) - 1
    text = "r 0.25\n" + format_decimal(plus, places) + " 0\n"
    _eq(failures, "poke +x", verdict_for("square", text.splitlines(), 1, "0.1"), "VERDICT: INVALID containment")
    text = "r 0.25\n" + format_decimal(-plus, places) + " 0\n"
    _eq(failures, "poke -x", verdict_for("square", text.splitlines(), 1, "0.1"), "VERDICT: INVALID containment")
    text = "r 0.25\n0 " + format_decimal(plus, places) + "\n"
    _eq(failures, "poke +y", verdict_for("square", text.splitlines(), 1, "0.1"), "VERDICT: INVALID containment")
    text = "r " + format_decimal(plus, places) + "\n" + format_decimal(minus, places) + " 0\n"
    _eq(failures, "cancel square", verdict_for("square", text.splitlines(), 1, "0.25"), "VERDICT: IMPROVES")
    rim = 9 * (10 ** 39)
    text = "r 0.1\n" + format_decimal(rim + 1, places) + " 0\n"
    _eq(failures, "poke disc", verdict_for("circle", text.splitlines(), 1, "0.1"), "VERDICT: INVALID containment")
    text = "r " + format_decimal((10 ** 39) + 1, places) + "\n" + format_decimal(rim - 1, places) + " 0\n"
    _eq(failures, "cancel disc", verdict_for("circle", text.splitlines(), 1, "0.1"), "VERDICT: IMPROVES")
    text = "r 0.1\n" + format_decimal(rim, places) + " 0\n"
    _eq(failures, "exact rim", verdict_for("circle", text.splitlines(), 1, "0.1"), "VERDICT: VALID_NOT_BETTER")
    text = "r 0.25\n" + format_decimal(25 * (10 ** 38), places) + " 0\n"
    _eq(failures, "exact wall", verdict_for("square", text.splitlines(), 1, "0.25"), "VERDICT: VALID_NOT_BETTER")


def test_grid(failures: list[str]) -> None:
    """Overlaps in far cells are impossible; near misses and tangencies are kept.

    The closest lattice points of two cells at Chebyshev distance >= 2 are
    still strictly more than 2r apart, so an overlapping pair cannot sit there.
    An overlap shifted into a huge cell index is still reported. Pruned and
    all-pairs answers agree, including on a negative boundary.
    """
    for r_scaled in (1, 2, 5, 17):
        length = 2 * r_scaled
        minimum = length * length
        for base_x in (-4, -1, 0, 2):
            for base_y in (-3, 0, 1):
                for di in (-3, -2, 0, 2, 4):
                    for dj in (-2, 0, 2, 5):
                        if max(abs(di), abs(dj)) < 2:
                            continue
                        x1, x2 = _facing(base_x, base_x + di, length)
                        y1, y2 = _facing(base_y, base_y + dj, length)
                        if (x1 // length, y1 // length) != (base_x, base_y):
                            failures.append(f"cell1 {r_scaled} {base_x} {base_y}")
                            return
                        if (x2 // length, y2 // length) != (base_x + di, base_y + dj):
                            failures.append(f"cell2 {r_scaled} {di} {dj}")
                            return
                        gap_x = abs(x1 - x2)
                        gap_y = abs(y1 - y2)
                        if abs(di) >= 2 and gap_x < length + 1:
                            failures.append(f"gap x {di}")
                            return
                        if abs(dj) >= 2 and gap_y < length + 1:
                            failures.append(f"gap y {dj}")
                            return
                        dist2 = gap_x * gap_x + gap_y * gap_y
                        if dist2 <= minimum:
                            failures.append(f"far overlap possible {dist2} <= {minimum}")
                            return
                        pair = [(x1, y1), (x2, y2)]
                        if not pairs_separated(pair, r_scaled) or not _naive_separated(pair, r_scaled):
                            failures.append(f"far pair rejected r={r_scaled}")
                            return
        for index in (-10 ** 6, -5, -1, 0, 3, 10 ** 6):
            left = (index + 1) * length - 1
            overlap = left + (length - 1)
            touch = left + length
            if left // length != index or overlap // length != index + 1 or touch // length != index + 1:
                failures.append(f"adj cell {index} {left // length}")
                return
            if (overlap - left) * (overlap - left) >= minimum:
                failures.append(f"adj dist {index}")
                return
            if (touch - left) * (touch - left) != minimum:
                failures.append(f"touch dist {index}")
                return
            if pairs_separated([(left, 0), (overlap, 0)], r_scaled):
                failures.append(f"adj overlap missed {index}")
                return
            if _naive_separated([(left, 0), (overlap, 0)], r_scaled):
                failures.append(f"adj naive missed {index}")
                return
            if not pairs_separated([(left, 0), (touch, 0)], r_scaled):
                failures.append(f"lattice tangent rejected {index}")
                return
            if not _naive_separated([(left, 0), (touch, 0)], r_scaled):
                failures.append(f"lattice tangent naive {index}")
                return
        for ix, iy in ((-3, -4), (-1, 0), (0, 0), (2, -2), (10 ** 5, -10 ** 5)):
            x1 = (ix + 1) * length - 1
            y1 = (iy + 1) * length - 1
            x2 = (ix + 1) * length
            y2 = (iy + 1) * length
            if (x1 // length, y1 // length) != (ix, iy):
                failures.append(f"diag from {ix} {iy}")
                return
            if (x2 // length, y2 // length) != (ix + 1, iy + 1):
                failures.append(f"diag to {ix} {iy}")
                return
            if pairs_separated([(x1, y1), (x2, y2)], r_scaled):
                failures.append(f"diag overlap missed {ix} {iy}")
                return
        lattice = [(k * length, m * length) for k in range(-4, 5) for m in range(-4, 5)]
        if not pairs_separated(lattice, r_scaled) or not _naive_separated(lattice, r_scaled):
            failures.append(f"lattice tangent cloud r={r_scaled}")
            return
        moved: list[tuple[int, int]] = []
        changed = False
        for x, y in lattice:
            if not changed and x == -2 * length and y == 0:
                moved.append((x + 1, y))
                changed = True
            else:
                moved.append((x, y))
        if not changed:
            failures.append("move failed")
            return
        pruned = pairs_separated(moved, r_scaled)
        naive = _naive_separated(moved, r_scaled)
        if pruned or naive:
            failures.append(f"neg boundary overlap r={r_scaled} pruned={pruned} naive={naive}")
            return
        cloud = [
            (k * length + (k * 3 + 1) % length, m * length + (m * 5 + 2) % length)
            for k in range(-4, 5)
            for m in range(-4, 5)
        ]
        for shift in (0, -1, 8, -(10 ** 6)):
            sample = [(x + shift, y - 2 * shift) for x, y in cloud]
            if pairs_separated(sample, r_scaled) != _naive_separated(sample, r_scaled):
                failures.append(f"cloud mismatch r={r_scaled} shift={shift}")
                return
    length = 2
    big = [(k * length, m * length) for k in range(100) for m in range(100)]
    if not pairs_separated(big, 1):
        failures.append("10000 tangent lattice")


def test_main(failures: list[str]) -> None:
    def run(argv: list[str]) -> tuple[int, str]:
        buf = StringIO()
        with redirect_stdout(buf):
            code = main(argv)
        return code, buf.getvalue().strip()

    def check(label: str, argv: list[str], want: str) -> None:
        code, out = run(argv)
        if out != want:
            failures.append(f"{label}: {out}")
            return
        expect = 0 if want in ("VERDICT: IMPROVES", "VERDICT: VALID_NOT_BETTER") else 1
        if code != expect:
            failures.append(f"{label} rc {code}")

    check("argc", ["verify_exact_big.py"], "VERDICT: INVALID arguments")
    check("cont", ["verify_exact_big.py", "Square", "nope", "1", "0.1"], "VERDICT: INVALID arguments")
    check("nneg", ["verify_exact_big.py", "square", "nope", "-1", "0.1"], "VERDICT: INVALID arguments")
    check("rec", ["verify_exact_big.py", "square", "nope", "1", "1e-3"], "VERDICT: INVALID arguments")
    with tempfile.TemporaryDirectory() as directory:
        missing = str(Path(directory) / "no-such-pack.txt")
        check("missing", ["verify_exact_big.py", "square", missing, "1", "0.1"], "VERDICT: INVALID file")
        path = str(Path(directory) / "pack.txt")
        Path(path).write_text("r 0.25\n-0.25 0\n0.25 0\n", encoding="utf-8")
        check("main =", ["verify_exact_big.py", "square", path, "2", "0.25"], "VERDICT: VALID_NOT_BETTER")
        check("main >", ["verify_exact_big.py", "square", path, "2", "0.2"], "VERDICT: IMPROVES")
        check("main n", ["verify_exact_big.py", "square", path, "9", "0.2"], "VERDICT: INVALID count")
        Path(path).write_text("nope\n", encoding="utf-8")
        check("main fmt", ["verify_exact_big.py", "circle", path, "0", "0.1"], "VERDICT: INVALID format")
        Path(path).write_text("r 0.2\n-0.8 0\n-0.4 0\n", encoding="utf-8")
        check("main cci", ["verify_exact_big.py", "circle", path, "2", "0.2"], "VERDICT: VALID_NOT_BETTER")


def run_selftest() -> tuple[bool, str]:
    failures: list[str] = []
    try:
        test_parse(failures)
        test_compare(failures)
        test_verdicts(failures)
        test_grid(failures)
        test_main(failures)
    except Exception as exc:
        failures.append(f"exception {type(exc).__name__}: {exc}")
    if not failures:
        return True, "SELFTEST: PASS"
    shown = failures[:12]
    if len(failures) > 12:
        shown.append(f"(+{len(failures) - 12} more)")
    return False, "SELFTEST: FAIL " + "; ".join(shown)


def main(argv: list[str]) -> int:
    if "--selftest" in argv[1:]:
        ok, message = run_selftest()
        print(message)
        return 0 if ok else 1
    if len(argv) != 5:
        print("VERDICT: INVALID arguments")
        return 1
    container, path, n_text, record_token = argv[1], argv[2], argv[3], argv[4]
    if container not in ("square", "circle"):
        print("VERDICT: INVALID arguments")
        return 1
    if not _digits(n_text):
        print("VERDICT: INVALID arguments")
        return 1
    try:
        record = parse_decimal(record_token)
    except ValueError:
        print("VERDICT: INVALID arguments")
        return 1
    try:
        with open(path, "r", encoding="utf-8-sig") as handle:
            lines = handle.read().splitlines()
    except OSError:
        print("VERDICT: INVALID file")
        return 1
    try:
        radius, points = parse_lines(lines)
    except ValueError:
        print("VERDICT: INVALID format")
        return 1
    verdict = evaluate(container, radius, points, int(n_text), record)
    print(verdict)
    if verdict in ("VERDICT: IMPROVES", "VERDICT: VALID_NOT_BETTER"):
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
