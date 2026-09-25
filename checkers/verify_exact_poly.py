#!/usr/bin/env python3
"""Exact checker for equal circles in a regular pentagon (Packomania cpt).

python verify_exact_poly.py 5 <file> <N> <record_radius>
  -> VERDICT: IMPROVES | VALID_NOT_BETTER | INVALID <why>

Frame: circumradius 1, centre at the origin, vertex at (0, 1), flat bottom.
Side lines are n·x = a with a = cos 36° = (1+√5)/4 and outward unit normals.
Decisions use Q(√5, σ), σ = √(10−2√5) > 0. No float is consulted.
"""

from __future__ import annotations

import os
import sys
import tempfile
from fractions import Fraction


class Q5:
    """u + v√5, u, v rational. The embedding takes √5 > 0."""

    __slots__ = ("u", "v")

    def __init__(self, u=0, v=0):
        if isinstance(u, float) or isinstance(v, float):
            raise TypeError("float")
        self.u = u if isinstance(u, Fraction) else Fraction(u)
        self.v = v if isinstance(v, Fraction) else Fraction(v)

    def __add__(self, other: "Q5") -> "Q5":
        return Q5(self.u + other.u, self.v + other.v)

    def __sub__(self, other: "Q5") -> "Q5":
        return Q5(self.u - other.u, self.v - other.v)

    def __neg__(self) -> "Q5":
        return Q5(-self.u, -self.v)

    def __mul__(self, other: "Q5") -> "Q5":
        return Q5(
            self.u * other.u + 5 * self.v * other.v,
            self.u * other.v + self.v * other.u,
        )

    def sign(self) -> int:
        u, v = self.u, self.v
        if u == 0 and v == 0:
            return 0
        if v == 0:
            return 1 if u > 0 else -1
        if u == 0:
            return 1 if v > 0 else -1
        uu = u * u
        five = 5 * v * v
        if uu < five:
            return 1 if v > 0 else -1
        if uu > five:
            return 1 if u > 0 else -1
        # u^2 = 5 v^2 is impossible for nonzero rationals. If it appeared,
        # opposite signs would be exact zero and equal signs would follow u.
        if (u > 0) == (v > 0):
            return 1 if u > 0 else -1
        return 0


class K:
    """P + σ Q, P, Q in Q(√5), σ = √(10−2√5) > 0.

    {1, σ} is a basis of K over Q(√5): σ is not in Q(√5), because
    Q(sin 36°, cos 36°) = Q(ζ5) has degree 4 and Q(cos 36°) = Q(√5).
    So P + σ Q = 0 only when P = Q = 0. Sign squares once against
    σ^2 = 10−2√5 and never evaluates a numeric square root.
    """

    __slots__ = ("p", "q")

    def __init__(self, p: Q5, q: Q5):
        if not isinstance(p, Q5) or not isinstance(q, Q5):
            raise TypeError("K expects Q5")
        self.p = p
        self.q = q

    @staticmethod
    def rat(x) -> "K":
        if isinstance(x, float):
            raise TypeError("float")
        fx = x if isinstance(x, Fraction) else Fraction(x)
        return K(Q5(fx, 0), Q5(0, 0))

    def __add__(self, other) -> "K":
        o = _as_k(other)
        return K(self.p + o.p, self.q + o.q)

    def __radd__(self, other) -> "K":
        return self + other

    def __sub__(self, other) -> "K":
        o = _as_k(other)
        return K(self.p - o.p, self.q - o.q)

    def __rsub__(self, other) -> "K":
        return _as_k(other) - self

    def __neg__(self) -> "K":
        return K(-self.p, -self.q)

    def __mul__(self, other) -> "K":
        o = _as_k(other)
        s2 = Q5(10, -2)  # σ^2 = 10 - 2√5
        return K(
            self.p * o.p + s2 * (self.q * o.q),
            self.p * o.q + self.q * o.p,
        )

    def __rmul__(self, other) -> "K":
        return self * other

    def sign(self) -> int:
        sq = self.q.sign()
        sp = self.p.sign()
        if sq == 0:
            return sp
        s2 = Q5(10, -2)
        diff = self.p * self.p - s2 * (self.q * self.q)
        sd = diff.sign()
        if sq > 0:
            # P + σ Q < 0 iff P < 0 and |P| > σ Q.
            if sp < 0 and sd > 0:
                return -1
            if sp < 0 and sd == 0:
                return 0
            return 1
        # Q < 0: P + σ Q > 0 iff P > 0 and |P| > σ |Q|.
        if sp > 0 and sd > 0:
            return 1
        if sp > 0 and sd == 0:
            return 0
        return -1


def _as_k(x) -> K:
    if isinstance(x, K):
        return x
    if isinstance(x, float):
        raise TypeError("float")
    if isinstance(x, Q5):
        return K(x, Q5(0, 0))
    return K.rat(x)


# a = cos 36° = (1+√5)/4
# cos 72° = (√5−1)/4
# sin 36° = σ/4
# sin 72° = σ(1+√5)/8 = sin 36° · (1+√5)/2
ALPHA = K(Q5(Fraction(1, 4), Fraction(1, 4)), Q5(0, 0))
GAMMA = K(Q5(Fraction(-1, 4), Fraction(1, 4)), Q5(0, 0))
BETA = K(Q5(0, 0), Q5(Fraction(1, 4), 0))
DELTA = K(Q5(0, 0), Q5(Fraction(1, 8), Fraction(1, 8)))
ZERO = K.rat(0)
ONE = K.rat(1)

# CCW edges from the top vertex. Each normal is outward and unit length,
# and both endpoints of the edge satisfy n·v = a.
NORMALS = (
    (-BETA, ALPHA),       # (0, 1) -- (-sin 72°, cos 72°)
    (-DELTA, -GAMMA),     # upper left -- bottom left
    (ZERO, -ONE),         # flat bottom, y = -a
    (DELTA, -GAMMA),      # bottom right -- upper right
    (BETA, ALPHA),        # upper right -- (0, 1)
)


def parse_number(tok: str) -> Fraction:
    if not isinstance(tok, str):
        raise ValueError("parse")
    try:
        return Fraction(tok.strip())
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("parse") from exc


def load(path: str):
    """Certificate: first `r <radius>`, then `x y` per centre.

    Also accepts Packomania `i x y` once a radius is known, a bare radius
    line, and a non-numeric author line. A coordinate file with no radius
    is an error: the radius is not in that format.
    """
    radius = None
    centres = []
    with open(path, "r", encoding="utf-8-sig") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line[0] in "#%":
                continue
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "r":
                if len(parts) != 2 or radius is not None:
                    raise ValueError("parse")
                radius = parse_number(parts[1])
                continue
            try:
                parse_number(parts[0])
            except ValueError:
                continue
            try:
                nums = [parse_number(p) for p in parts]
            except ValueError as exc:
                raise ValueError("parse") from exc
            if len(nums) == 1:
                if radius is not None or centres:
                    raise ValueError("parse")
                radius = nums[0]
            elif len(nums) == 2:
                centres.append((nums[0], nums[1]))
            elif len(nums) == 3:
                centres.append((nums[1], nums[2]))
            else:
                raise ValueError("parse")
    if radius is None:
        raise ValueError("parse")
    return radius, centres


def validate(centres, radius, n: int):
    if len(centres) != n:
        return False, "count"
    r = _as_k(radius)
    if r.sign() <= 0:
        return False, "radius"
    pts = [(_as_k(x), _as_k(y)) for x, y in centres]
    limit = r * r * 4
    for i, (x, y) in enumerate(pts):
        for w, (nx, ny) in enumerate(NORMALS, 1):
            # Distance to the side is a − n·c. Need that >= r.
            if (ALPHA - r - nx * x - ny * y).sign() < 0:
                return False, f"wall {w} circle {i + 1}"
        for j in range(i + 1, n):
            dx = x - pts[j][0]
            dy = y - pts[j][1]
            if (dx * dx + dy * dy - limit).sign() < 0:
                return False, f"overlap {i + 1} {j + 1}"
    return True, ""


def classify(centres, radius, n: int, record) -> str:
    ok, why = validate(centres, radius, n)
    if not ok:
        return f"VERDICT: INVALID {why}"
    if (_as_k(radius) - _as_k(record)).sign() > 0:
        return "VERDICT: IMPROVES"
    return "VERDICT: VALID_NOT_BETTER"


def evaluate(argv: list[str]) -> str:
    if len(argv) != 5:
        return "VERDICT: INVALID usage"
    if argv[1] != "5":
        return "VERDICT: INVALID sides"
    path, n_tok, rec_tok = argv[2], argv[3], argv[4]
    if n_tok.startswith("+"):
        n_tok = n_tok[1:]
    if not n_tok.isdigit():
        return "VERDICT: INVALID parse"
    try:
        record = parse_number(rec_tok)
        radius, centres = load(path)
    except (OSError, ValueError, TypeError):
        return "VERDICT: INVALID parse"
    return classify(centres, radius, int(n_tok), record)


def _dec_str(fr: Fraction) -> str:
    if fr < 0:
        return "-" + _dec_str(-fr)
    whole, rem = divmod(fr.numerator, fr.denominator)
    digs = []
    d = fr.denominator
    for _ in range(200):
        if rem == 0:
            break
        rem *= 10
        digs.append(str(rem // d))
        rem %= d
    else:
        raise AssertionError("nonterminating decimal")
    text = str(whole) if not digs else f"{whole}.{''.join(digs)}"
    if Fraction(text) != fr:
        raise AssertionError("decimal roundtrip")
    return text


def _need(got: str, want: str, label: str) -> None:
    if got != want:
        raise AssertionError(f"{label}: {got}")


def _file_verdict(text: str, n: int, record: str) -> str:
    with tempfile.TemporaryDirectory() as folder:
        path = os.path.join(folder, "pack.txt")
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        return evaluate(["verify_exact_poly.py", "5", path, str(n), record])


def self_test() -> None:
    try:
        _self_test_body()
    except AssertionError as exc:
        print(f"SELF-TEST: FAIL {exc}")
        raise SystemExit(1)
    print("SELF-TEST: PASS")


def _self_test_body() -> None:
    one = K.rat(1)
    verts = (
        (ZERO, ONE),
        (-DELTA, GAMMA),
        (-BETA, -ALPHA),
        (BETA, -ALPHA),
        (DELTA, GAMMA),
    )
    golden = K(Q5(Fraction(1, 2), Fraction(1, 2)), Q5(0, 0))
    if (DELTA - BETA * golden).sign() != 0:
        raise AssertionError("sin72 identity")
    if (K(Q5(1, 1), Q5(0, 0)) - ALPHA * 4).sign() != 0:
        raise AssertionError("cos36 identity")
    for nx, ny in NORMALS:
        if (nx * nx + ny * ny - one).sign() != 0:
            raise AssertionError("normal length")
    for i, (nx, ny) in enumerate(NORMALS):
        for step in (0, 1):
            vx, vy = verts[(i + step) % 5]
            if (nx * vx + ny * vy - ALPHA).sign() != 0:
                raise AssertionError(f"side {i + 1}")
    # 4/5 < cos 36° < 81/100, proved by squaring, locks the Q(√5) sign.
    if (ALPHA - Fraction(4, 5)).sign() <= 0:
        raise AssertionError("a > 4/5")
    if (ALPHA - Fraction(81, 100)).sign() >= 0:
        raise AssertionError("a < 81/100")

    eps30 = Fraction(1, 10**30)
    eps25 = Fraction(1, 10**25)
    origin = [(ZERO, ZERO)]
    _need(classify(origin, ALPHA - eps30, 1, 0), "VERDICT: IMPROVES", "n1 under")
    _need(
        classify(origin, ALPHA + eps30, 1, 0),
        "VERDICT: INVALID wall 1 circle 1",
        "n1 over",
    )
    _need(classify(origin, ALPHA, 1, ALPHA), "VERDICT: VALID_NOT_BETTER", "n1 touch-all")

    # Centre built on the bottom inner line y = -a + r. Distance to that
    # side is exactly r; the other four slacks stay positive for this r.
    r0 = Fraction(1, 5)
    touching = [(ZERO, -ALPHA + r0)]
    _need(classify(touching, r0, 1, 0), "VERDICT: IMPROVES", "side exact")
    _need(classify(touching, r0 - eps30, 1, 0), "VERDICT: IMPROVES", "side under")
    _need(
        classify(touching, r0 + eps30, 1, 0),
        "VERDICT: INVALID wall 3 circle 1",
        "side over",
    )

    r = Fraction(1, 10)
    pair = [(K.rat(-r), ZERO), (K.rat(r), ZERO)]
    _need(classify(pair, r, 2, 0), "VERDICT: IMPROVES", "pair exact")
    # Same centres, distance 2r. Claimed diameter 2r + 1e-25 overlaps by 1e-25.
    _need(
        classify(pair, r + eps25 / 2, 2, 0),
        "VERDICT: INVALID overlap 1 2",
        "pair overlap",
    )
    _need(
        classify([(K.rat(1), ZERO)], Fraction(1, 100), 1, 0),
        "VERDICT: INVALID wall 4 circle 1",
        "right wall",
    )

    _need(evaluate(["prog"]), "VERDICT: INVALID usage", "usage")
    _need(evaluate(["prog", "6", "x", "1", "0"]), "VERDICT: INVALID sides", "sides")
    _need(evaluate(["prog", "5", "x", "1.0", "0"]), "VERDICT: INVALID parse", "bad n")

    _need(
        _file_verdict("r 0.1\n0 0\n", 1, "0"),
        "VERDICT: IMPROVES",
        "file inside",
    )
    _need(
        _file_verdict("r 0.1\n0 0\n", 1, "0.1"),
        "VERDICT: VALID_NOT_BETTER",
        "file equal record",
    )
    _need(
        _file_verdict("r 0.1\n0 0\n", 1, "0.2"),
        "VERDICT: VALID_NOT_BETTER",
        "file worse",
    )
    tiny = _dec_str(Fraction(1, 2) + eps30)
    _need(
        _file_verdict(f"r {tiny}\n0 0\n", 1, "0.5"),
        "VERDICT: IMPROVES",
        "file 1e-30 above record",
    )
    _need(
        _file_verdict("r 0.01\n0 -1\n", 1, "0"),
        "VERDICT: INVALID wall 3 circle 1",
        "file below bottom",
    )
    _need(
        _file_verdict("r 0.1\n-0.1 0\n0.1 0\n", 2, "0.09"),
        "VERDICT: IMPROVES",
        "file pair exact",
    )
    gap = _dec_str(Fraction(1, 5) - eps25)
    _need(
        _file_verdict(f"r 0.1\n0 0\n{gap} 0\n", 2, "0"),
        "VERDICT: INVALID overlap 1 2",
        "file overlap",
    )
    _need(
        _file_verdict("r 0.1\n1 0 0\n2 0.05 0\n", 1, "0"),
        "VERDICT: INVALID count",
        "count",
    )
    _need(
        _file_verdict("0.1\nE. Specht\n# note\n\n1 0 0\n", 1, "0"),
        "VERDICT: IMPROVES",
        "pck shape",
    )
    _need(_file_verdict("1 0 0\n", 1, "0.1"), "VERDICT: INVALID parse", "no radius")
    _need(
        _file_verdict("r 0\n0 0\n", 1, "0"),
        "VERDICT: INVALID radius",
        "zero radius",
    )


if __name__ == "__main__":
    if len(sys.argv) == 1:
        self_test()
    else:
        print(evaluate(sys.argv))
