# Checker B (independent implementation), saved verbatim 2026-09-24. USED ONLY FOR THE SEMICIRCLE (csc).
# NOTE: written without seeing verify_exact2.py; its rect:<h> branch assumes a [0,1] x [0,h] frame, NOT the centred
# frame of our certificates -> rectangles and the quadrant stay on verify_exact2.py.
# PATCH (2026-09-24, labelled below): Parser.expr/term looped past the end ('' in '+-' is True) -> every call raised
# IndexError (fail-safe: no verdict at all). Two guards added; nothing else changed. Original: verify_exact3_orig.py.
# verify_exact3.py
# py -3.11 verify_exact3.py <container> <file> <N> <record>
# py -3.11 verify_exact3.py --self-test
#
# container: rect:<h> | quad | semi
# file: "r <value>" then one "x y" per circle. Each value is one token:
#   finite decimal, scientific decimal, a/b, or + - * / and sqrt(integer)
#   in a single quadratic field. No spaces inside a token. 2*sqrt(3), not 2sqrt(3).
#
# rect:<h>  rectangle [0, 1] x [0, h]
# quad      unit quarter-disk, walls x >= r, y >= r, and the arc
# semi      unit upper half-disk: y >= r and x^2 + y^2 <= (1 - r)^2, with 1 - r >= 0
#
# stdout is exactly four lines:
#   count: <read>/<N>
#   fit: <inside>/<read>
#   pairs: <separated>/<pairs>
#   VERDICT: IMPROVES|TIES|WORSE|FAIL

import sys
from decimal import Decimal, InvalidOperation
from fractions import Fraction


class Q:
    """a + b*sqrt(d). b = 0 means a rational. Operations stay exact."""

    __slots__ = ("a", "b", "d")

    def __init__(self, a=0, b=0, d=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)
        self.d = int(d)
        if self.b == 0:
            self.d = 0
        elif self.d <= 0:
            raise ValueError("radicand")

    def __neg__(self):
        return Q(-self.a, -self.b, self.d)

    def __add__(self, other):
        if not isinstance(other, Q):
            other = Q(other)
        if self.b == 0 and other.b == 0:
            return Q(self.a + other.a)
        if self.b == 0:
            return Q(self.a + other.a, other.b, other.d)
        if other.b == 0:
            return Q(self.a + other.a, self.b, self.d)
        if self.d != other.d:
            raise ValueError("mixed square roots")
        return Q(self.a + other.a, self.b + other.b, self.d)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if not isinstance(other, Q):
            other = Q(other)
        return self + (-other)

    def __rsub__(self, other):
        return Q(other) - self

    def __mul__(self, other):
        if not isinstance(other, Q):
            other = Q(other)
        if self.b == 0 and other.b == 0:
            return Q(self.a * other.a)
        if self.b == 0:
            return Q(self.a * other.a, self.a * other.b, other.d)
        if other.b == 0:
            return Q(other.a * self.a, other.a * self.b, self.d)
        if self.d != other.d:
            raise ValueError("mixed square roots")
        return Q(
            self.a * other.a + self.b * other.b * self.d,
            self.a * other.b + self.b * other.a,
            self.d,
        )

    def __rmul__(self, other):
        return self * other

    def inv(self):
        if self.b == 0:
            if self.a == 0:
                raise ZeroDivisionError
            return Q(Fraction(1) / self.a)
        disc = self.a * self.a - self.b * self.b * self.d
        if disc == 0:
            raise ZeroDivisionError
        return Q(self.a / disc, -self.b / disc, self.d)

    def __truediv__(self, other):
        if not isinstance(other, Q):
            other = Q(other)
        return self * other.inv()

    def __rtruediv__(self, other):
        return Q(other) / self


def cmp0(q):
    """Sign of q as -1, 0, or +1."""
    if not isinstance(q, Q):
        q = Q(q)
    a, b, d = q.a, q.b, q.d
    if b == 0:
        return (a > 0) - (a < 0)
    if a == 0:
        return (b > 0) - (b < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    # Opposite signs: compare a^2 with b^2 * d.
    left = a * a
    right = b * b * d
    if a > 0:
        if left > right:
            return 1
        if left < right:
            return -1
        return 0
    if right > left:
        return 1
    if right < left:
        return -1
    return 0


def sqrt_int(n):
    if n < 0:
        raise ValueError("sqrt")
    if n == 0:
        return Q(0)
    a = 1
    d = n
    i = 2
    while i * i <= d:
        sq = i * i
        while d % sq == 0:
            d //= sq
            a *= i
        i += 1
    if d == 1:
        return Q(a)
    return Q(0, a, d)


class Parser:
    def __init__(self, s):
        self.s = s.strip()
        self.i = 0
        self.n = len(self.s)

    def parse(self):
        if self.n == 0:
            raise ValueError("empty")
        value = self.expr()
        if self.i != self.n:
            raise ValueError("trailing")
        return value

    def peek(self):
        return self.s[self.i] if self.i < self.n else ""

    def expr(self):
        value = self.term()
        while self.peek() != "" and self.peek() in "+-":   # PATCH (2026-09-24): "" in "+-" is True in Python
            op = self.s[self.i]
            self.i += 1
            rhs = self.term()
            value = value + rhs if op == "+" else value - rhs
        return value

    def term(self):
        value = self.factor()
        while self.peek() != "" and self.peek() in "*/":   # PATCH (2026-09-24): same end-of-string bug
            op = self.s[self.i]
            self.i += 1
            rhs = self.factor()
            value = value * rhs if op == "*" else value / rhs
        return value

    def factor(self):
        ch = self.peek()
        if ch == "+":
            self.i += 1
            return self.factor()
        if ch == "-":
            self.i += 1
            return -self.factor()
        if ch == "(":
            self.i += 1
            value = self.expr()
            if self.peek() != ")":
                raise ValueError("paren")
            self.i += 1
            return value
        if self.s.startswith("sqrt(", self.i):
            self.i += 5
            j = self.i
            if j < self.n and self.s[j] == "+":
                self.i += 1
                j = self.i
            if j >= self.n or not self.s[j].isdigit():
                raise ValueError("sqrt")
            while self.i < self.n and self.s[self.i].isdigit():
                self.i += 1
            n = int(self.s[j:self.i])
            if self.peek() != ")":
                raise ValueError("sqrt")
            self.i += 1
            return sqrt_int(n)
        return self.number()

    def number(self):
        j = self.i
        while self.i < self.n and (self.s[self.i].isdigit() or self.s[self.i] == "."):
            self.i += 1
        if self.i < self.n and self.s[self.i] in "eE":
            self.i += 1
            if self.i < self.n and self.s[self.i] in "+-":
                self.i += 1
            exp0 = self.i
            while self.i < self.n and self.s[self.i].isdigit():
                self.i += 1
            if self.i == exp0:
                raise ValueError("exponent")
        if self.i == j:
            raise ValueError("number")
        tok = self.s[j:self.i]
        try:
            return Q(Fraction(Decimal(tok)))
        except (InvalidOperation, ValueError) as exc:
            raise ValueError(tok) from exc


def parse_expr(text):
    return Parser(text).parse()


def emit(n_read, n_expect, n_fit, n_ok, n_pairs, verdict):
    # The 1,990-file diff hinges on these three labels.
    print(f"count: {n_read}/{n_expect}")
    print(f"fit: {n_fit}/{n_read}")
    print(f"pairs: {n_ok}/{n_pairs}")
    print(f"VERDICT: {verdict}")


def fits(kind, h, x, y, r):
    if cmp0(r) <= 0:
        return False
    if kind == "rect":
        if cmp0(x - r) < 0:
            return False
        if cmp0((Q(1) - r) - x) < 0:
            return False
        if cmp0(y - r) < 0:
            return False
        if cmp0((h - r) - y) < 0:
            return False
        return True
    one_m = Q(1) - r
    if cmp0(one_m) < 0:
        return False
    if kind == "quad":
        if cmp0(x - r) < 0:
            return False
    elif kind != "semi":
        raise ValueError(kind)
    if cmp0(y - r) < 0:
        return False
    arc = one_m * one_m - x * x - y * y
    return cmp0(arc) >= 0


def parse_container(spec):
    spec = spec.strip()
    if spec == "quad":
        return "quad", None
    if spec == "semi":
        return "semi", None
    prefix = "rect:"
    if spec.startswith(prefix):
        h = parse_expr(spec[len(prefix):])
        if cmp0(h) <= 0:
            raise ValueError("height")
        return "rect", h
    raise ValueError("container")


def split_r(head):
    if len(head) < 2 or head[0] != "r" or head[1] not in " \t=":
        raise ValueError("header")
    raw = head[2:].strip() if head[1] == "=" else head[1:].strip()
    if raw == "":
        raise ValueError("header")
    return raw


def load_packing(path):
    with open(path, encoding="utf-8-sig") as fh:
        lines = [ln.strip() for ln in fh.readlines()]
    lines = [ln for ln in lines if ln]
    if not lines:
        raise ValueError("empty file")
    r = parse_expr(split_r(lines[0]))
    pts = []
    for ln in lines[1:]:
        parts = ln.split()
        if len(parts) != 2:
            raise ValueError("coordinate")
        pts.append((parse_expr(parts[0]), parse_expr(parts[1])))
    return r, pts


def decide(kind, h, r, pts, n_expect, record):
    n_read = len(pts)
    n_fit = sum(1 for x, y in pts if fits(kind, h, x, y, r))
    n_pairs = n_read * (n_read - 1) // 2
    n_ok = 0
    four_rr = Q(4) * r * r
    for i in range(n_read):
        xi, yi = pts[i]
        for j in range(i + 1, n_read):
            dx = xi - pts[j][0]
            dy = yi - pts[j][1]
            if cmp0(dx * dx + dy * dy - four_rr) >= 0:
                n_ok += 1
    feasible = (
        n_expect >= 1
        and n_read == n_expect
        and n_fit == n_read
        and n_ok == n_pairs
    )
    if not feasible:
        verdict = "FAIL"
    else:
        c = cmp0(r - record)
        if c > 0:
            verdict = "IMPROVES"
        elif c == 0:
            verdict = "TIES"
        else:
            verdict = "WORSE"
    return n_read, n_fit, n_ok, n_pairs, verdict


def self_test():
    def packing(raw_r, coords):
        return parse_expr(raw_r), [(parse_expr(x), parse_expr(y)) for x, y in coords]

    record = parse_expr("1/3")
    r, pts = packing(
        "1/3",
        [("-sqrt(3)/3", "1/3"), ("sqrt(3)/3", "1/3"), ("0", "2/3")],
    )
    assert decide("semi", None, r, pts, 3, record)[-1] == "TIES"
    assert decide("semi", None, r, pts, 3, parse_expr("3/10"))[-1] == "IMPROVES"
    assert decide("semi", None, r, pts[:2], 3, record)[-1] == "FAIL"

    pushed = packing(
        "1/3",
        [("-sqrt(3)/3", "1/3"), ("sqrt(3)/3", "1/3"), ("0", "2/3+1e-30")],
    )
    assert decide("semi", None, pushed[0], pushed[1], 3, record)[-1] == "FAIL"

    bigger = parse_expr("1/3+1e-30")
    assert decide("semi", None, bigger, pts, 3, record)[-1] == "FAIL"

    overlap = packing("1/3", [("0", "1/3"), ("0", "1/3"), ("0", "2/3")])
    assert decide("semi", None, overlap[0], overlap[1], 3, record)[-1] == "FAIL"

    left = packing("1/5", [("-1/2", "1/2")])
    assert decide("semi", None, left[0], left[1], 1, parse_expr("1/5"))[-1] == "TIES"
    assert decide("quad", None, left[0], left[1], 1, parse_expr("1/5"))[-1] == "FAIL"

    half = packing("1/2", [("0", "1/2")])
    assert decide("semi", None, half[0], half[1], 1, parse_expr("1/2"))[-1] == "TIES"
    dec = packing("0.5", [("0", "0.5")])
    assert decide("semi", None, dec[0], dec[1], 1, parse_expr("1/2"))[-1] == "TIES"
    small = packing("2/5", [("0", "1/2")])
    assert decide("semi", None, small[0], small[1], 1, parse_expr("1/2"))[-1] == "WORSE"

    rq = parse_expr("sqrt(2)-1")
    assert decide("quad", None, rq, [(rq, rq)], 1, rq)[-1] == "TIES"

    box = packing("1/4", [("1/4", "1/4")])
    h = parse_expr("1")
    assert decide("rect", h, box[0], box[1], 1, parse_expr("1/4"))[-1] == "TIES"
    out = packing("1/4", [("1/10", "1/4")])
    assert decide("rect", h, out[0], out[1], 1, parse_expr("1/4"))[-1] == "FAIL"
    print("self-test: ok")
    return 0


def main(argv):
    if len(argv) == 2 and argv[1] == "--self-test":
        return self_test()
    if len(argv) != 5:
        sys.stderr.write(
            "usage: verify_exact3.py <container> <file> <N> <record>\n"
        )
        return 2
    try:
        n_expect = int(argv[3])
    except ValueError:
        emit(0, 0, 0, 0, 0, "FAIL")
        return 0
    try:
        kind, h = parse_container(argv[1])
        record = parse_expr(argv[4])
        r, pts = load_packing(argv[2])
        n_read, n_fit, n_ok, n_pairs, verdict = decide(
            kind, h, r, pts, n_expect, record
        )
    except (OSError, ValueError, ZeroDivisionError, InvalidOperation):
        emit(0, n_expect, 0, 0, 0, "FAIL")
        return 0
    emit(n_read, n_expect, n_fit, n_ok, n_pairs, verdict)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
