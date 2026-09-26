# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Exact checker: N equal balls of radius r inside the unit ball in R^d.

Usage:
    python verify_exact_ball.py <d> <file> <N> <record_radius>

No float is used in any decision. Decimals are integers on a power-of-ten grid.
Run with no arguments to execute the self-tests.
"""

import io
import sys
import tempfile
import pathlib
from contextlib import redirect_stdout
from itertools import combinations

MAX_EXP = 100_000
MAX_DIGITS = 2_000


def ten_pow(n: int) -> int:
    if n < 0:
        raise ValueError("negative power")
    return 10 ** n


def normalize(mant: int, scale: int) -> tuple[int, int]:
    if mant == 0:
        return 0, 0
    while scale > 0 and mant % 10 == 0:
        mant //= 10
        scale -= 1
    return mant, scale


def parse_decimal(token: str) -> tuple[int, int]:
    """Value = mantissa / 10**scale, both integers, scale >= 0."""
    s = token.strip()
    if not s:
        raise ValueError("empty")
    sign = 1
    if s[0] in "+-":
        sign = -1 if s[0] == "-" else 1
        s = s[1:]
    if not s:
        raise ValueError("sign only")
    exp = 0
    e_at = -1
    for i, ch in enumerate(s):
        if ch in "eE":
            e_at = i
            break
    if e_at >= 0:
        exp_body = s[e_at + 1 :]
        s = s[:e_at]
        if not exp_body or not s:
            raise ValueError("exponent")
        exp_sign = 1
        if exp_body[0] in "+-":
            exp_sign = -1 if exp_body[0] == "-" else 1
            exp_body = exp_body[1:]
        if not exp_body or len(exp_body) > 6 or not exp_body.isdigit():
            raise ValueError("exponent")
        exp = exp_sign * int(exp_body)
    if abs(exp) > MAX_EXP:
        raise ValueError("exponent")
    if not s or s == "." or s.count(".") > 1:
        raise ValueError("digits")
    if "." in s:
        left, right = s.split(".", 1)
        if left == "":
            left = "0"
        if right == "":
            if not left.isdigit():
                raise ValueError("digits")
            digits = left
            frac = 0
        else:
            if not left.isdigit() or not right.isdigit():
                raise ValueError("digits")
            digits = left + right
            frac = len(right)
    else:
        if not s.isdigit():
            raise ValueError("digits")
        digits = s
        frac = 0
    if len(digits) > MAX_DIGITS:
        raise ValueError("digits")
    mant = sign * int(digits)
    pow10 = exp - frac
    if pow10 > MAX_EXP or pow10 < -MAX_EXP:
        raise ValueError("exponent")
    if pow10 >= 0:
        mant = mant * ten_pow(pow10)
        scale = 0
    else:
        scale = -pow10
    return normalize(mant, scale)


def greater(left: tuple[int, int], right: tuple[int, int]) -> bool:
    """Strict > for two mantissa/scale pairs. Powers of ten stay non-negative."""
    lm, ls = left
    rm, rs = right
    return lm * ten_pow(rs) > rm * ten_pow(ls)


def evaluate(d: int, text: str, n: int, record: str) -> str:
    if d < 1:
        return "VERDICT: INVALID dimension must be positive"
    if n < 0:
        return "VERDICT: INVALID count must be non-negative"
    try:
        record_pair = parse_decimal(record)
    except ValueError:
        return "VERDICT: INVALID bad record"
    if text.startswith("\ufeff"):
        text = text[1:]
    lines = text.splitlines()
    if not lines:
        return "VERDICT: INVALID empty file"
    head = lines[0].split()
    if len(head) != 2 or head[0] != "r":
        return "VERDICT: INVALID bad radius line"
    try:
        radius_pair = parse_decimal(head[1])
    except ValueError:
        return "VERDICT: INVALID bad radius"
    body = lines[1:]
    if len(body) != n:
        return f"VERDICT: INVALID expected {n} centers, found {len(body)}"
    parsed_rows: list[list[tuple[int, int]]] = []
    for line_no, line in enumerate(body, start=2):
        parts = line.split()
        if len(parts) != d:
            return (
                f"VERDICT: INVALID line {line_no} has {len(parts)} numbers, expected {d}"
            )
        row: list[tuple[int, int]] = []
        for tok in parts:
            try:
                row.append(parse_decimal(tok))
            except ValueError:
                return f"VERDICT: INVALID bad number on line {line_no}"
        parsed_rows.append(row)
    depth = radius_pair[1]
    for row in parsed_rows:
        for _mant, scale in row:
            if scale > depth:
                depth = scale
    ten = ten_pow(depth)
    radius = radius_pair[0] * ten_pow(depth - radius_pair[1])
    if radius <= 0:
        return "VERDICT: INVALID non-positive radius"
    # r > 1 still leaves (1-r)^2 positive, so the square test cannot catch it.
    if radius > ten:
        return "VERDICT: INVALID radius exceeds 1"
    centers: list[list[int]] = []
    for row in parsed_rows:
        centers.append([mant * ten_pow(depth - scale) for mant, scale in row])
    room = ten - radius
    limit = room * room
    for i, center in enumerate(centers, start=1):
        sq = 0
        for coord in center:
            sq += coord * coord
        if sq > limit:
            return f"VERDICT: INVALID center {i} outside"
    four_r2 = 4 * radius * radius
    count = len(centers)
    for i in range(count):
        ci = centers[i]
        for j in range(i + 1, count):
            cj = centers[j]
            dist = 0
            for a, b in zip(ci, cj):
                diff = a - b
                dist += diff * diff
            if dist < four_r2:
                return f"VERDICT: INVALID centers {i + 1} and {j + 1} overlap"
    if greater(radius_pair, record_pair):
        return "VERDICT: IMPROVES"
    return "VERDICT: VALID_NOT_BETTER"


def evaluate_file(d: int, path: str, n: int, record: str) -> str:
    try:
        text = pathlib.Path(path).read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        return "VERDICT: INVALID cannot read file"
    return evaluate(d, text, n, record)


def main(argv: list[str]) -> int:
    if len(argv) != 5:
        print("VERDICT: INVALID bad arguments")
        return 0
    d_s, path, n_s, record = argv[1], argv[2], argv[3], argv[4]
    if not d_s.isdigit() or not n_s.isdigit():
        print("VERDICT: INVALID bad arguments")
        return 0
    d = int(d_s)
    n = int(n_s)
    if d < 1:
        print("VERDICT: INVALID dimension must be positive")
        return 0
    print(evaluate_file(d, path, n, record))
    return 0


def format_int_scale(mant: int, scale: int) -> str:
    if scale < 0:
        raise ValueError("negative scale")
    sign = "-" if mant < 0 else ""
    n = abs(mant)
    if scale == 0:
        return sign + str(n)
    raw = str(n)
    if len(raw) <= scale:
        body = "0." + raw.zfill(scale)
    else:
        body = raw[: len(raw) - scale] + "." + raw[-scale:]
    return sign + body


def same_value(a: tuple[int, int], b: tuple[int, int]) -> bool:
    return (not greater(a, b)) and (not greater(b, a))


def cert(radius: str, rows: list[str]) -> str:
    return "r " + radius + "\n" + "".join(row + "\n" for row in rows)


def zeros(d: int) -> str:
    return " ".join(["0"] * d)


def axis(first: str, d: int) -> str:
    return " ".join([first] + ["0"] * (d - 1))


def twenty_four(comp: str) -> list[str]:
    neg = "-" + comp
    rows: list[str] = []
    for i, j in combinations(range(4), 2):
        for si in (comp, neg):
            for sj in (comp, neg):
                coords = ["0"] * 4
                coords[i] = si
                coords[j] = sj
                rows.append(" ".join(coords))
    return rows


def build_cell24() -> tuple[str, str, str]:
    """24 vectors +-e_i +-e_j with 1/sqrt(2) truncated inward, times (1-r).

    r is the 40-digit truncation of 1/3 - 10^(-12). Length is a hair under 1-r
    because the components are a hair under (1-r)/sqrt(2).
    """
    k = 30
    p = 40
    m = int("707106781186547524400844362104")
    ten_2k = ten_pow(2 * k)
    if not (2 * m * m <= ten_2k < 2 * ((m + 1) ** 2)):
        raise RuntimeError("1/sqrt(2) digits are not the inward truncation")
    r_mant = ((ten_pow(12) - 3) * ten_pow(p)) // (3 * ten_pow(12))
    one_minus = ten_pow(p) - r_mant
    if one_minus <= 0 or r_mant <= 0:
        raise RuntimeError("radius fixture")
    # (1-r)^2 * t^2 >= 2 r^2  is the closest-pair condition.
    if one_minus * one_minus * m * m < 2 * r_mant * r_mant * ten_2k:
        raise RuntimeError("24-point fixture overlaps")
    c_mant = m * one_minus
    scale = k + p
    limit = (ten_pow(k) * one_minus) ** 2
    if 2 * c_mant * c_mant > limit:
        raise RuntimeError("24-point fixture outside")
    moved_mant = c_mant + ten_pow(scale - 25)
    if moved_mant * moved_mant + c_mant * c_mant <= limit:
        raise RuntimeError("nudged centre still inside")
    comp = format_int_scale(c_mant, scale)
    moved = format_int_scale(moved_mant, scale)
    r_str = format_int_scale(r_mant, p)
    if not same_value(parse_decimal(comp), normalize(c_mant, scale)):
        raise RuntimeError("component round-trip")
    if not same_value(parse_decimal(moved), normalize(moved_mant, scale)):
        raise RuntimeError("moved round-trip")
    if not same_value(parse_decimal(r_str), normalize(r_mant, p)):
        raise RuntimeError("radius round-trip")
    return r_str, comp, moved


def run_case(d: int, text: str, n: int, record: str) -> str:
    with tempfile.TemporaryDirectory() as td:
        path = pathlib.Path(td) / "c.txt"
        path.write_text(text, encoding="utf-8")
        return evaluate_file(d, str(path), n, record)


def run_self_tests() -> int:
    failed = 0

    def expect(name: str, d: int, text: str, n: int, record: str, kind: str) -> None:
        nonlocal failed
        got = run_case(d, text, n, record)
        if kind == "INVALID":
            ok = got.startswith("VERDICT: INVALID ")
        else:
            ok = got == "VERDICT: " + kind
        print(f"{name}: {got} {'PASS' if ok else 'FAIL'}")
        if not ok:
            failed += 1

    try:
        r24, comp, moved = build_cell24()
    except (RuntimeError, ValueError) as exc:
        print(f"FIXTURE: {exc} FAIL")
        return 1
    rows24 = twenty_four(comp)
    if len(rows24) != 24:
        print(f"FIXTURE: expected 24 rows, found {len(rows24)} FAIL")
        return 1

    r_lo = "0." + ("9" * 30)
    r_hi = "1." + ("0" * 29) + "1"
    r_half_hi = "0.5" + ("0" * 28) + "1"
    record_lo_same = "9." + ("9" * 29) + "e-1"

    for d in (4, 5, 6):
        z = zeros(d)
        halves = [axis("0.5", d), axis("-0.5", d)]
        expect(f"d{d} centre r=1-1e-30", d, cert(r_lo, [z]), 1, "0", "IMPROVES")
        expect(
            f"d{d} centre equal record",
            d,
            cert(r_lo, [z]),
            1,
            record_lo_same,
            "VALID_NOT_BETTER",
        )
        expect(f"d{d} centre r=1+1e-30", d, cert(r_hi, [z]), 1, "0", "INVALID")
        expect(f"d{d} tangent halves", d, cert("0.5", halves), 2, "0", "IMPROVES")
        expect(
            f"d{d} tangent record equal",
            d,
            cert("0.5", halves),
            2,
            "5e-1",
            "VALID_NOT_BETTER",
        )
        expect(
            f"d{d} tangent record larger",
            d,
            cert("0.5", halves),
            2,
            "0.6",
            "VALID_NOT_BETTER",
        )
        expect(
            f"d{d} halves r=1/2+1e-30",
            d,
            cert(r_half_hi, halves),
            2,
            "0",
            "INVALID",
        )
        short = " ".join(["0"] * (d - 1))
        expect(f"d{d} short line", d, cert("0.1", [short]), 1, "0", "INVALID")
        expect(f"d{d} wrong count", d, cert("0.1", [z]), 2, "0", "INVALID")

    expect("d4 24-cell", 4, cert(r24, rows24), 24, "0", "IMPROVES")
    expect(
        "d4 24-cell equal record", 4, cert(r24, rows24), 24, r24, "VALID_NOT_BETTER"
    )
    nudged = rows24[:]
    parts = nudged[0].split(" ")
    parts[0] = moved
    nudged[0] = " ".join(parts)
    expect("d4 nudged outside", 4, cert(r24, nudged), 24, "0", "INVALID")
    expect(
        "d4 stacked overlap",
        4,
        cert("0.5", [zeros(4), zeros(4)]),
        2,
        "0",
        "INVALID",
    )
    expect("d4 long line", 4, cert("0.1", ["0 0 0 0 0"]), 1, "0", "INVALID")

    missing = evaluate_file(4, "this-file-does-not-exist-hsp.txt", 1, "0")
    ok_missing = missing.startswith("VERDICT: INVALID ")
    print(f"missing file: {missing} {'PASS' if ok_missing else 'FAIL'}")
    if not ok_missing:
        failed += 1

    with tempfile.TemporaryDirectory() as td:
        path = pathlib.Path(td) / "c.txt"
        path.write_text(cert(r_lo, [zeros(4)]), encoding="utf-8")
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = main(["verify_exact_ball.py", "4", str(path), "1", "0"])
        got = buf.getvalue().strip()
        ok_cli = rc == 0 and got == "VERDICT: IMPROVES"
        print(f"cli wiring: {got} {'PASS' if ok_cli else 'FAIL'}")
        if not ok_cli:
            failed += 1

    print("SELF-TEST: ALL PASS" if failed == 0 else f"SELF-TEST: {failed} FAILED")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] == "--self-test"):
        sys.exit(run_self_tests())
    sys.exit(main(sys.argv))
