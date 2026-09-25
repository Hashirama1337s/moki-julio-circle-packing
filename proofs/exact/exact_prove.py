# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""EXACT optimality of two small rectangle entries (Moki&Julio).
Status: reviewed (three adversarial rounds, 2026-09-25): PUBLISHABLE; published in v1.7. Second checker: proofs/exact/verify_proof_exact.py (independent), v1.9.

Claim. N equal circles in the rectangle [-1/2, 1/2] x [-h/2, h/2] (THE frame: centred) cannot have radius > r* = p - q sqrt2,
and the exhibited configuration c* attains r*.
    crc_800_4 : N = 4, h = 4/5, r* = 7/20 - sqrt2/10
    crc_600_6 : N = 6, h = 3/5, r* = 1/5 - sqrt2/30

Lemma 1 (local, Q(sqrt2) stress certificate) gives a rational rho: the only configuration of radius >= r* within inf-norm rho of c*
(or of any mirror image / relabelling of c*) is that image itself.
Lemma 2 (global) runs the branch-and-bound engine (prove_small.py, class Engine; exact integer arithmetic on a 2^k grid) at a rational
r_t < r* through LoggedEngine below. LoggedEngine re-implements tighten / apply_symmetry / shrink_box / best_hole / bisect / the
start enumeration WITH A WITNESS LOG and CROSS-CHECKS every node against the unmodified Engine methods (identical output required).
Every node ends DISCARDED (exact witness) or ACCEPTED (all point boxes inside closed rho'-squares around the snapped points of one
mirror image / labelling of c*, rho' + snap error <= rho).

Arithmetic. Fractions and Q(sqrt2) = {a + b sqrt2 : a, b rational} (class Q2, exact field operations and an exact sign test).
sqrt2 enters numerically ONLY through the enclosure S2_LO < sqrt2 < S2_HI (S2_LO^2 < 2 < S2_HI^2 checked on integers):
    r_t := p - q S2_HI - 10^-40 < p - q sqrt2 = r*    (q > 0),
    snapped image coordinate a + b S2_LO, error |b| |sqrt2 - S2_LO| < |b| (S2_HI - S2_LO).
Floats appear ONLY in the LP that CHOOSES the free stress parameters and in ranking candidate bases; every chosen quantity is then
verified exactly. No sympy, no nsimplify.

Output (proof runs only): proofs/exact/<case>_cert.json and proofs/exact/<case>_tree.json.gz;
independent replay: py -3.11 proofs/exact/replay_exact.py

usage: py -3.11 exact_prove.py crc_800_4             proof run (writes certificate + tree)
       py -3.11 exact_prove.py crc_600_6 1/1000      NEGATIVE CONTROL: r_t = (p - q S2_HI)(1 - 1/1000); must end UNDECIDED;
                                                     writes nothing
"""
import sys, os, time, json, gzip, hashlib, itertools
from fractions import Fraction
from math import gcd, isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
PUBLISHED = os.path.basename(HERE) == "exact"              # published repo layout: proofs/exact/{exact_prove,prove_small}.py
sys.path.insert(0, HERE if PUBLISHED else os.path.join(HERE, "engine"))
import prove_small as PS                                   # the engine (read-only; subclassed below, never edited)

OUT = HERE if PUBLISHED else os.path.join(HERE, "proofs", "exact")
K50 = 10 ** 50
S2_LO = Fraction(isqrt(2 * K50 * K50), K50)                # floor(sqrt2 10^50) / 10^50, integer square root: no floats
S2_HI = S2_LO + Fraction(1, K50)
assert 0 < S2_LO < S2_HI and S2_LO * S2_LO < 2 < S2_HI * S2_HI
ROUND = 10 ** 30                                           # outward rounding grid for the printed rational bounds
MARGIN = Fraction(1, 10 ** 40)                             # r_t = p - q S2_HI - MARGIN
NODE_LIMIT = 5_000_000
CROSSCHECK = True                                          # compare every logged step with the unmodified Engine


def floor_q(x, D=ROUND):
    return Fraction((x.numerator * D) // x.denominator, D)


def dec(x, digits=30):
    """display only: x (Fraction) truncated toward -inf to `digits` decimals, exact integer arithmetic"""
    v = (x.numerator * 10 ** digits) // x.denominator; s = "-" if v < 0 else ""; v = abs(v)
    return f"{s}{v // 10 ** digits}.{v % 10 ** digits:0{digits}d}"


def ceil_q(x, D=ROUND):
    return Fraction(-((-x.numerator * D) // x.denominator), D)


class Q2:
    """a + b sqrt2, a and b Fractions: exact arithmetic in the field Q(sqrt2)."""
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    @staticmethod
    def of(v):
        return v if isinstance(v, Q2) else Q2(v)

    def __add__(s, o):
        o = Q2.of(o); return Q2(s.a + o.a, s.b + o.b)
    __radd__ = __add__

    def __sub__(s, o):
        o = Q2.of(o); return Q2(s.a - o.a, s.b - o.b)

    def __rsub__(s, o):
        return Q2.of(o) - s

    def __neg__(s):
        return Q2(-s.a, -s.b)

    def __mul__(s, o):
        o = Q2.of(o); return Q2(s.a * o.a + 2 * s.b * o.b, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__

    def inv(s):
        n = s.a * s.a - 2 * s.b * s.b                      # the field norm; 0 only for 0 because sqrt2 is irrational
        if n == 0:
            raise ZeroDivisionError("inverse of 0 in Q(sqrt2)")
        return Q2(s.a / n, -s.b / n)

    def __truediv__(s, o):
        return s * Q2.of(o).inv()

    def is_zero(s):
        return s.a == 0 and s.b == 0

    def __eq__(s, o):
        o = Q2.of(o); return s.a == o.a and s.b == o.b
    __hash__ = None

    def sign(s):
        """exact sign of a + b sqrt2"""
        a, b = s.a, s.b
        if a >= 0 and b >= 0:
            return 0 if (a == 0 and b == 0) else 1
        if a <= 0 and b <= 0:
            return -1
        t = a * a - 2 * b * b                              # opposite strict signs: compare |a| with |b| sqrt2; t != 0
        return 1 if (a > 0) == (t > 0) else -1

    def encl(s):
        """rigorous rational enclosure [lo, hi] (from S2_LO < sqrt2 < S2_HI)"""
        x, y = s.a + s.b * S2_LO, s.a + s.b * S2_HI
        return (x, y) if x <= y else (y, x)

    def lo(s):
        return s.encl()[0]

    def hi(s):
        return s.encl()[1]

    def abs_hi(s):
        x, y = s.encl(); return max(abs(x), abs(y))

    def absq(s):
        return s if s.sign() >= 0 else -s

    def js(s):
        return [str(s.a), str(s.b)]

    def __float__(s):                                      # display / LP ONLY
        return float(s.a) + float(s.b) * 1.4142135623730951

    def __repr__(s):
        return f"({s.a}) + ({s.b})*sqrt2"


def Pq(a, b):
    return Q2(Fraction(a), Fraction(b))


# c* exactly, centred frame (Packomania's packings identified in Q(sqrt2); every property used is re-verified exactly below)
CASES = {
    "crc_800_4": dict(N=4, h=Fraction(4, 5), p=Fraction(7, 20), q=Fraction(1, 10), c=[
        (Pq("-3/20", "-1/10"), Pq("-1/20", "-1/10")),
        (Pq("11/20", "-3/10"), Pq("-1/20", "-1/10")),
        (Pq("-11/20", "3/10"), Pq("1/20", "1/10")),
        (Pq("3/20", "1/10"), Pq("1/20", "1/10"))]),
    "crc_600_6": dict(N=6, h=Fraction(3, 5), p=Fraction(1, 5), q=Fraction(1, 30), c=[
        (Pq("-3/10", "-1/30"), Pq("-1/10", "-1/30")),
        (Pq("1/10", "-1/10"), Pq("-1/10", "-1/30")),
        (Pq("1/2", "-1/6"), Pq("-1/10", "-1/30")),
        (Pq("-1/2", "1/6"), Pq("1/10", "1/30")),
        (Pq("-1/10", "1/10"), Pq("1/10", "1/30")),
        (Pq("3/10", "1/30"), Pq("1/10", "1/30"))]),
}
MIRRORS = ((1, 1), (-1, 1), (1, -1), (-1, -1))
HALF = Fraction(1, 2)


# ----------------------------------------------------------------------------------------------------------------- Lemma 1
def constraints(c, r, h):
    """every constraint g_k(c, r) >= 0 of the problem, centred frame, fixed order (walls L R B T per circle, then pairs i<j).
    row = grad_c g_k (length 2N), w = -dg_k/dr (1 for walls, 8r for pairs)."""
    n = len(c); out = []
    for i, (x, y) in enumerate(c):
        for side, g, v, s in (("L", x + HALF - r, 2 * i, 1), ("R", HALF - x - r, 2 * i, -1),
                              ("B", y + h / 2 - r, 2 * i + 1, 1), ("T", h / 2 - y - r, 2 * i + 1, -1)):
            row = [Q2(0)] * (2 * n); row[v] = Q2(s)
            out.append(dict(kind="wall", i=i, side=side, g=g, row=row, w=Q2(1)))
    for i, j in itertools.combinations(range(n), 2):
        dx, dy = c[i][0] - c[j][0], c[i][1] - c[j][1]
        row = [Q2(0)] * (2 * n)
        row[2 * i], row[2 * i + 1], row[2 * j], row[2 * j + 1] = 2 * dx, 2 * dy, -2 * dx, -2 * dy
        out.append(dict(kind="pair", i=i, j=j, g=dx * dx + dy * dy - 4 * r * r, row=row, w=8 * r))
    return out


def cname(k):
    return f"wall {k['side']} of circle {k['i']}" if k["kind"] == "wall" else f"pair ({k['i']},{k['j']})"


def cdesc(k):
    return {"kind": "wall", "i": k["i"], "side": k["side"]} if k["kind"] == "wall" else {"kind": "pair", "i": k["i"], "j": k["j"]}


def rref(A):
    """reduced row echelon form over Q(sqrt2) of an augmented matrix (last column = right-hand side)"""
    A = [row[:] for row in A]; rows, cols = len(A), len(A[0]) - 1; piv = []; r = 0
    for col in range(cols):
        p = next((i for i in range(r, rows) if not A[i][col].is_zero()), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        iv = A[r][col].inv(); A[r] = [x * iv for x in A[r]]
        for i in range(rows):
            if i != r and not A[i][col].is_zero():
                f = A[i][col]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(col); r += 1
        if r == rows:
            break
    return A, piv


def inverse(M):
    """exact inverse over Q(sqrt2) by Gauss-Jordan, or None if singular"""
    n = len(M); A = [list(M[i]) + [Q2(1) if i == j else Q2(0) for j in range(n)] for i in range(n)]
    for col in range(n):
        p = next((i for i in range(col, n) if not A[i][col].is_zero()), None)
        if p is None:
            return None
        A[col], A[p] = A[p], A[col]
        iv = A[col][col].inv(); A[col] = [x * iv for x in A[col]]
        for i in range(n):
            if i != col and not A[i][col].is_zero():
                f = A[i][col]; A[i] = [x - f * y for x, y in zip(A[i], A[col])]
    return [row[n:] for row in A]


def matmul(X, Y):
    return [[sum((X[i][k] * Y[k][j] for k in range(len(Y))), Q2(0)) for j in range(len(Y[0]))] for i in range(len(X))]


def local_certificate(case):
    C = CASES[case]; n, h, p, q, c = C["N"], C["h"], C["p"], C["q"], C["c"]
    r = Q2(p, -q)
    assert q > 0 and r.sign() > 0
    cons = constraints(c, r, h)
    for k in cons:                                              # c* in F(r*): EXACT sign test for every constraint
        assert k["g"].sign() >= 0, ("c* infeasible at", cname(k))
    T = [k for k in cons if k["g"].is_zero()]                   # tight set: EXACT zero test in Q(sqrt2)
    slack = [k for k in cons if not k["g"].is_zero()]
    m = len(T)
    # stress system: sum_k lambda_k row_k = 0 (2N equations), sum_k lambda_k w_k = 1
    E = [[T[k]["row"][v] for k in range(m)] + [Q2(0)] for v in range(2 * n)] + [[T[k]["w"] for k in range(m)] + [Q2(1)]]
    R, piv = rref(E)
    for row in R[len(piv):]:
        assert all(x.is_zero() for x in row[:m]) and row[m].is_zero(), "no stress (inconsistent system)"
    free = [k for k in range(m) if k not in piv]
    lam_p = [Q2(0)] * m
    for ri, pc in enumerate(piv):
        lam_p[pc] = R[ri][m]
    null = []
    for f in free:
        v = [Q2(0)] * m; v[f] = Q2(1)
        for ri, pc in enumerate(piv):
            v[pc] = -R[ri][f]
        null.append(v)
    t = []
    if free:   # CHOOSE the free parameters: float LP maximising min lambda; the rational choice is verified exactly below
        import numpy as np
        from scipy.optimize import linprog
        coef = np.array([[float(null[f][k]) for f in range(len(free))] for k in range(m)])
        const = np.array([float(x) for x in lam_p])
        A = np.hstack([-coef, np.ones((m, 1))])
        res = linprog(np.r_[np.zeros(len(free)), -1.0], A_ub=A, b_ub=const,
                      bounds=[(None, None)] * len(free) + [(None, 1.0)], method="highs")
        assert res.status == 0 and res.x[-1] > 0, ("no strictly positive stress", res.x)
        t = [Fraction(float(x)).limit_denominator(10 ** 6) for x in res.x[:len(free)]]
    lam = [lam_p[k] + sum((t[f] * null[f][k] for f in range(len(free))), Q2(0)) for k in range(m)]
    for v in range(2 * n):                                      # equilibrium, EXACT
        assert sum((lam[k] * T[k]["row"][v] for k in range(m)), Q2(0)).is_zero()
    assert sum((lam[k] * T[k]["w"] for k in range(m)), Q2(0)) == Q2(1)   # normalisation, EXACT
    lam_lo = [floor_q(l.lo()) for l in lam]
    assert all(v > 0 for v in lam_lo), ("stress not strictly positive", [float(l) for l in lam])
    assert all((lam[k] - lam_lo[k]).sign() >= 0 for k in range(m))
    pairs = [k for k in range(m) if T[k]["kind"] == "pair"]
    Lam_hi = ceil_q(sum((lam[k].hi() for k in pairs), Fraction(0)))
    assert (Q2(Lam_hi) - sum((lam[k] for k in pairs), Q2(0))).sign() >= 0
    # basis B: 2N rows of T. Rank candidates with floats (CHOICE only), then take the best rigorous one among the top few.
    import numpy as np
    rowsf = np.array([[float(x) for x in k["row"]] for k in T]); lamf = np.array([float(x) for x in lam])
    cand = []
    for B in itertools.combinations(range(m), 2 * n):
        Mf = rowsf[list(B)]
        if abs(np.linalg.det(Mf)) < 1e-9:
            continue
        nf = np.abs(np.linalg.inv(Mf)).sum(axis=1).max()
        cand.append((-lamf[list(B)].min() / (2 * nf), B))
    cand.sort()
    best = None
    for _, B in cand[:8]:
        M = [T[k]["row"][:] for k in B]; Mi = inverse(M)
        if Mi is None:
            continue
        I = matmul(M, Mi)
        assert all(I[a][b] == Q2(1 if a == b else 0) for a in range(2 * n) for b in range(2 * n)), "M Minv != I"
        normM_hi = ceil_q(max(sum((x.abs_hi() for x in row), Fraction(0)) for row in Mi))
        eta_lo = min(lam_lo[k] for k in B) / (2 * normM_hi)
        if best is None or eta_lo > best[0]:
            best = (eta_lo, list(B), M, Mi, normM_hi)
    eta_lo, B, M, Mi, normM_hi = best
    for row in Mi:                                              # ||Minv||_inf <= normM_hi, EXACT (|x| exact in Q(sqrt2))
        assert (Q2(normM_hi) - sum((x.absq() for x in row), Q2(0))).sign() >= 0
    rho = min(lam_lo[k] for k in B) / (2 * normM_hi) / (16 * Lam_hi)
    return dict(n=n, h=h, p=p, q=q, r=r, c=c, T=T, slack=slack, lam=lam, lam_lo=lam_lo, free=free, t=t, Lam_hi=Lam_hi,
                B=B, M=M, Mi=Mi, normM_hi=normM_hi, rho=rho)


def snapped_images(c, h, r_t):
    """the 4 mirror images of c*, point frame at r_t (p = c + (1/2 - r_t, h/2 - r_t)), sqrt2 replaced by S2_LO"""
    out = []
    for sx, sy in MIRRORS:
        out.append([(sx * (x.a + x.b * S2_LO) + HALF - r_t, sy * (y.a + y.b * S2_LO) + h / 2 - r_t) for x, y in c])
    return out


# ------------------------------------------------------------------------------------------------ Lemma 2: logged engine
class LoggedEngine(PS.Engine):
    """The Engine with the grid G chosen by the caller, and witness-logging copies of the search steps.
    Nothing here changes a decision: every logged step is compared with the unmodified Engine method (CROSSCHECK)."""

    def __init__(self, W, H, d, n, G):
        super().__init__(W, H, d, n, Fraction(0))
        if G < 2 or G % 2:
            raise ValueError("G must be even: apply_symmetry uses half = G // 2 as the exact value G/2")
        self.G = G
        Wn, Wd, Hn, Hd, dn, dd = W.numerator, W.denominator, H.numerator, H.denominator, d.numerator, d.denominator
        A = Wn * Wn * Hd * Hd * dd * dd; B = Hn * Hn * Wd * Wd * dd * dd; C = dn * dn * Wd * Wd * Hd * Hd
        g = gcd(gcd(A, B), C)
        self.A, self.B, self.CG = A // g, B // g, (C // g) * G * G     # same integers as Engine.__init__ with this G

    def far_xy(self, a, b):
        return [PS.far(a[0], a[1], b[0], b[1]), PS.far(a[2], a[3], b[2], b[3])]

    def starts_log(self):
        """Engine.build_starts, returning the cell cover, the surviving assignments and a witness for every rejected one"""
        n, G = self.n, self.G
        best = None
        for nx in range(1, n + 2):
            for ny in range(1, n + 2):
                cells = self.make_cells(nx, ny)
                if any(not self.too_close(cc, cc) for cc in cells):
                    continue
                if len(cells) < n:
                    raise RuntimeError("pigeonhole start: nothing is feasible at r_t, but c* is -> data or engine error")
                ways = self.perm_count(len(cells), n)
                if ways is None:
                    continue
                if best is None or ways < best[0]:
                    best = (ways, cells)
        if best is None:
            return "full", None, [tuple((0, G, 0, G) for _ in range(n))], None, []
        cells = best[1]; kids, assigns, discards = [], [], []
        for assign in itertools.permutations(range(len(cells)), n):
            st = tuple(cells[i] for i in assign); w = None
            for i in range(n - 1):
                if st[i][0] > st[i + 1][1]:
                    w = ["xorder", i]; break
            if w is None:
                for i in range(n):
                    for j in range(i + 1, n):
                        if self.too_close(st[i], st[j]):
                            w = ["pair", i, j] + self.far_xy(st[i], st[j]); break
                    if w:
                        break
            if w is None:
                kids.append(st); assigns.append(list(assign))
            else:
                discards.append([list(assign), w])
        return "cells", [list(cc) for cc in cells], kids, assigns, discards

    def symmetry_log(self, boxes, ops):
        """Engine.apply_symmetry with every clip logged; returns a discard witness or None"""
        n, G = len(boxes), self.G; half = G // 2
        changed = True
        while changed:
            changed = False
            for i in range(n - 1):
                if boxes[i][1] > boxes[i + 1][1]:
                    boxes[i][1] = boxes[i + 1][1]; ops.append(["o_hi", i]); changed = True
                if boxes[i + 1][0] < boxes[i][0]:
                    boxes[i + 1][0] = boxes[i][0]; ops.append(["o_lo", i]); changed = True
                if boxes[i][0] > boxes[i][1]:
                    return ["empty", i, 0]
                if boxes[i][2] > boxes[i][3]:
                    return ["empty", i, 1]
            if boxes[-1][0] > boxes[-1][1]:
                return ["empty", n - 1, 0]
            if boxes[-1][2] > boxes[-1][3]:
                return ["empty", n - 1, 1]
            cap = G - boxes[-1][0]
            if boxes[0][1] > cap:
                boxes[0][1] = cap; ops.append(["s1"]); changed = True
            capn = G - boxes[0][0]
            if boxes[-1][1] > capn:
                boxes[-1][1] = capn; ops.append(["sN"]); changed = True
            if boxes[0][0] > boxes[0][1]:
                return ["empty", 0, 0]
            if boxes[-1][0] > boxes[-1][1]:
                return ["empty", n - 1, 0]
            if boxes[0][0] + boxes[-1][0] > G:
                return ["xsum"]
            if boxes[0][3] > half:
                boxes[0][3] = half; ops.append(["y1"]); changed = True
            if boxes[0][2] > boxes[0][3]:
                return ["empty", 0, 1]
        return None

    def shrink_log(self, i, j, boxes, ops):
        """Engine.shrink_box(boxes[i], boxes[j]) with every strip removal logged; returns (new box or None, witness)"""
        other = tuple(boxes[j]); x0, x1, y0, y1 = boxes[i]
        v = self.push_low_x(x0, x1, y0, y1, other)
        if v != x0:
            ops.append(["sh", i, j, "x0", v]); x0 = v
        if x0 > x1:
            return None, ["empty", i, 0]
        v = self.push_high_x(x0, x1, y0, y1, other)
        if v != x1:
            ops.append(["sh", i, j, "x1", v]); x1 = v
        if x0 > x1:
            return None, ["empty", i, 0]
        v = self.push_low_y(x0, x1, y0, y1, other)
        if v != y0:
            ops.append(["sh", i, j, "y0", v]); y0 = v
        if y0 > y1:
            return None, ["empty", i, 1]
        v = self.push_high_y(x0, x1, y0, y1, other)
        if v != y1:
            ops.append(["sh", i, j, "y1", v]); y1 = v
        if y0 > y1:
            return None, ["empty", i, 1]
        if self.rect_forbidden(x0, x1, y0, y1, other):
            return None, ["forbid", i, j]
        return (x0, x1, y0, y1), None

    def tighten_log(self, state):
        """Engine.tighten with a log; returns (tightened state or None, ops, discard witness or None)"""
        boxes = [list(b) for b in state]; n = self.n; ops = []
        for _ in range(64):
            w = self.symmetry_log(boxes, ops)
            if w:
                return None, ops, w
            for i in range(n):
                bi = boxes[i]
                if bi[0] > bi[1]:
                    return None, ops, ["empty", i, 0]
                if bi[2] > bi[3]:
                    return None, ops, ["empty", i, 1]
                for j in range(i + 1, n):
                    if self.too_close(bi, boxes[j]):
                        return None, ops, ["pair", i, j] + self.far_xy(bi, boxes[j])
            changed = False
            for i in range(n):
                for j in range(n):
                    if i == j or not self.can_forbid(boxes[j]):
                        continue
                    nb, w = self.shrink_log(i, j, boxes, ops)
                    if nb is None:
                        return None, ops, w
                    if nb != tuple(boxes[i]):
                        boxes[i][:] = nb; changed = True
            if not changed:
                break
        return tuple(tuple(b) for b in boxes), ops, None

    def best_hole_j(self, state):
        """Engine.best_hole, also returning the index j of the box that forbids the hole"""
        best = None; best_cut = 0; n = self.n
        for i in range(n):
            bi = state[i]; area = (bi[1] - bi[0]) * (bi[3] - bi[2])
            if area <= 0:
                continue
            for j in range(n):
                if i == j or not self.can_forbid(state[j]):
                    continue
                hole = self.grow_hole(bi, state[j])
                if hole is None:
                    continue
                if hole[0] <= bi[0] and hole[1] >= bi[1] and hole[2] <= bi[2] and hole[3] >= bi[3]:
                    return i, j, hole, True
                cut = (hole[1] - hole[0]) * (hole[3] - hole[2])
                if cut * 8 >= area and cut > best_cut:
                    best_cut = cut; best = (i, j, hole, False)
        return best

    def branch_log(self, state):
        """Engine.branch, returning ('discard', witness, None) | ('split', how, kids) | ('stuck', None, None)"""
        found = self.best_hole_j(state)
        if found is not None:
            i, j, hole, covers = found
            if covers:
                return "discard", ["forbid", i, j], None
            kids = self.hole_children(state, i, hole)
            if kids:
                return "split", ["hole", i, j, list(hole)], kids
        i, axis = self.choose_bisect(state)
        b = state[i]; lo, hi = (b[0], b[1]) if axis == 0 else (b[2], b[3])
        if hi - lo < 2:
            return "stuck", None, None          # Engine.branch would return None here: NOT a discard (see PROOFS_EXACT.md)
        kids = self.bisect(state)
        return "split", ["bisect", i, axis, (lo + hi) // 2], kids


# ------------------------------------------------------------------------------------------------------------------ main
def main(case, off=None):
    t_start = time.time()
    L = local_certificate(case)
    n, h, p, q, r, c = L["n"], L["h"], L["p"], L["q"], L["r"], L["c"]
    rho = L["rho"]
    bmax = max(abs(v.b) for pt in c for v in pt)
    snap_err_hi = bmax * (S2_HI - S2_LO)                   # |b| |sqrt2 - S2_LO| < |b| (S2_HI - S2_LO)
    rho_p = floor_q(rho - snap_err_hi)
    assert rho_p > 0 and rho_p + snap_err_hi <= rho
    print(f"{case}: N = {n}, h = {h}, r* = {p} - ({q}) sqrt2 in [{dec(r.lo())}, {dec(r.hi())}]  (Lemma 1 data: "
          f"{time.time() - t_start:.1f}s)", flush=True)
    print(f"  tight set T: {len(L['T'])} constraints ({sum(1 for k in L['T'] if k['kind'] == 'pair')} pairs): "
          + ", ".join(cname(k) for k in L["T"]), flush=True)
    print(f"  stress: free parameters {len(L['free'])} fixed at {[str(x) for x in L['t']]}; min lambda_lo = "
          f"{float(min(L['lam_lo'])):.6g}; Lambda_hi = {float(L['Lam_hi']):.6g}; ||Minv||_inf <= {float(L['normM_hi']):.6g}")
    print(f"  basis B = {[cname(L['T'][k]) for k in L['B']]}")
    print(f"  rho = {float(rho):.9g} (rational), snap error < {float(snap_err_hi):.3g}, rho' = {float(rho_p):.9g}", flush=True)

    if off is None:
        r_t = p - q * S2_HI - MARGIN
        assert r_t > 0 and (r - r_t).sign() > 0          # r_t < r*: q > 0 and S2_HI > sqrt2; also checked exactly in Q(sqrt2)
    else:
        r_t = (p - q * S2_HI) * (1 - Fraction(off))
        print(f"  NEGATIVE CONTROL: r_t = (p - q S2_HI)(1 - {off}) (rational); nothing is written", flush=True)
    W, H, d = 1 - 2 * r_t, h - 2 * r_t, 2 * r_t
    imgs = snapped_images(c, h, r_t)
    for (sx, sy), pts in zip(MIRRORS, imgs):              # snap bound, re-checked exactly in Q(sqrt2)
        for (X, Y), (x, y) in zip(pts, c):
            for S, e in ((X, sx * x + HALF - r_t), (Y, sy * y + h / 2 - r_t)):
                assert (Q2(snap_err_hi) - (Q2(S) - e).absq()).sign() > 0
    G = 2
    while G * rho_p < 16 * max(W, H):                     # at least ~16 grid cells per rho'
        G <<= 1
    eng = LoggedEngine(W, H, d, n, G)
    sW, sH = W / G, H / G
    print(f"  r_t = p - q S2_HI - 1e-40 (rational, < r*)" if off is None else "", f" grid G = {G}", flush=True)

    def accept(st):
        boxes = [(b[0] * sW, b[1] * sW, b[2] * sH, b[3] * sH) for b in st]
        for mi, pts in enumerate(imgs):
            cand = [[k for k, (px, py) in enumerate(pts)
                     if px - rho_p <= bx[0] and bx[1] <= px + rho_p and py - rho_p <= bx[2] and bx[3] <= py + rho_p]
                    for bx in boxes]
            if any(not cs for cs in cand):
                continue
            used, perm = set(), []

            def dfs(i):
                if i == n:
                    return True
                for k in cand[i]:
                    if k not in used:
                        used.add(k); perm.append(k)
                        if dfs(i + 1):
                            return True
                        used.discard(k); perm.pop()
                return False
            if dfs(0):
                return mi, perm
        return None

    kind, cells, starts, assigns, discards = eng.starts_log()
    if CROSSCHECK:
        k0, st0, _ = PS.Engine.build_starts(eng)
        assert {"assign": "cells", "full": "full"}[k0] == kind and list(st0) == list(starts), \
            "start enumeration differs from Engine.build_starts"
    print(f"  start: {kind}, {len(starts)} surviving assignments, {len(discards)} discarded "
          f"({time.time() - t_start:.1f}s)", flush=True)
    recs = {0: {"id": 0, "parent": None, "box": [[0, G, 0, G] for _ in range(n)], "ops": []}}
    root_children = []
    stack = []
    for sidx, st in enumerate(starts):
        nid = len(recs); recs[nid] = {"id": nid, "parent": 0, "box": [list(b) for b in st]}
        root_children.append(nid); stack.append((nid, st))
    if kind == "cells":
        recs[0]["fate"] = {"t": "cells", "cells": cells, "children": root_children, "assign": assigns, "discards": discards}
    else:
        recs[0]["fate"] = {"t": "full", "children": root_children}
    popped = acc = 0; wcount = {}
    for dsc in discards:
        wcount["start:" + dsc[1][0]] = wcount.get("start:" + dsc[1][0], 0) + 1
    status = "COMPLETE"
    while stack:
        nid, state = stack.pop(); popped += 1
        rec = recs[nid]
        if popped > NODE_LIMIT:
            status = "UNDECIDED node limit"; break
        tst, ops, w = eng.tighten_log(state)
        if CROSSCHECK:
            assert tst == PS.Engine.tighten(eng, state), f"tighten log differs from Engine.tighten at node {nid}"
        rec["ops"] = ops
        if tst is None:
            rec["fate"] = {"t": "discard", "w": w}; wcount[w[0]] = wcount.get(w[0], 0) + 1
            continue
        rec["tbox"] = [list(b) for b in tst]
        a = accept(tst)
        if a is not None:
            rec["fate"] = {"t": "accept", "mirror": a[0], "perm": a[1]}; acc += 1
            continue
        if eng.max_side(tst) <= 1:
            status = "UNDECIDED tolerance"; rec["fate"] = {"t": "undecided"}
            print("  undecided boxes:", [tuple(float(v) for v in bx) for bx in eng.real_boxes(tst)])
            break
        bk, how, kids = eng.branch_log(tst)
        if CROSSCHECK:
            orig = PS.Engine.branch(eng, tst)
            assert (orig is None and bk == "discard") or (bk == "split" and orig == kids), f"branch log differs at node {nid}"
        if bk == "stuck":
            status = "UNDECIDED stuck"; rec["fate"] = {"t": "undecided"}; break
        if bk == "discard":
            rec["fate"] = {"t": "discard", "w": how}; wcount["branch:" + how[0]] = wcount.get("branch:" + how[0], 0) + 1
            continue
        ids = [None] * len(kids)
        for _, idx, ch in sorted((eng.score(k), idx, k) for idx, k in enumerate(kids)):
            cid = len(recs); recs[cid] = {"id": cid, "parent": nid, "box": [list(b) for b in ch]}
            ids[idx] = cid; stack.append((cid, ch))
        rec["fate"] = {"t": "split", "how": how, "children": ids}
    elapsed = time.time() - t_start
    if status != "COMPLETE":
        print(f"{status}: {case}, r_t = {float(r_t):.12g}, nodes popped {popped}, {elapsed:.1f}s", flush=True)
        return status
    if off is not None:
        print(f"!! NEGATIVE CONTROL COMPLETED ({popped} nodes) -- this must not happen; investigate", flush=True)
        return "CONTROL FAILED"
    leaves = sum(1 for v in recs.values() if v["fate"]["t"] in ("discard", "accept"))
    print(f"EXACT PROOF COMPLETE: {case}: r* = {p} - ({q}) sqrt2 is optimal; tree nodes {len(recs)} (root + {popped} "
          f"processed), leaves {leaves}: accepted {acc}, discarded {leaves - acc} + {len(discards)} start assignments; "
          f"witness counts {dict(sorted(wcount.items()))}; {elapsed:.1f}s", flush=True)

    # --------------------------------------------------------------------------------------------- certificate + tree
    os.makedirs(OUT, exist_ok=True)
    T = L["T"]
    cert = {
        "format": "exact-optimality-certificate v1", "case": case, "authors": "Moki&Julio",
        "frame": "centred: rectangle [-1/2,1/2] x [-h/2,h/2]; a Q(sqrt2) number is [a, b] = a + b*sqrt2 (a, b rationals)",
        "N": n, "h": str(h), "width": "1",
        "r_star": {"p": str(p), "q": str(q), "meaning": "r* = p - q*sqrt2"},
        "c_star": [[x.js(), y.js()] for x, y in c],
        "constraint_order": "walls L,R,B,T per circle i = 0..N-1, then pairs (i,j), i<j lexicographic",
        "constraints": "L: x_i + 1/2 - r, R: 1/2 - x_i - r, B: y_i + h/2 - r, T: h/2 - y_i - r, pair: |c_i - c_j|^2 - 4 r^2; "
                       "all >= 0; w = -dg/dr = 1 (walls), 8r (pairs)",
        "T": [cdesc(k) for k in T],
        "slack_lo": {cname(k): str(floor_q(k["g"].lo())) for k in L["slack"]},
        "S2_LO": str(S2_LO), "S2_HI": str(S2_HI),
        "lambda": [l.js() for l in L["lam"]],
        "lambda_lo": [str(v) for v in L["lam_lo"]],
        "stress_free_params": [str(v) for v in L["t"]],
        "Lambda_hi": str(L["Lam_hi"]),
        "B": L["B"],
        "M": [[x.js() for x in row] for row in L["M"]],
        "Minv": [[x.js() for x in row] for row in L["Mi"]],
        "normMinv_hi": str(L["normM_hi"]),
        "rho": str(rho),
        "rho_formula": "rho = min_{k in B} lambda_lo_k / (2 normMinv_hi) / (16 Lambda_hi)",
        "snap_err_hi": str(snap_err_hi),
        "rho_prime": str(rho_p),
        "r_t": str(r_t), "r_t_formula": "r_t = p - q*S2_HI - 1/10^40",
        "grid_G": G,
        "images_point_frame": [{"mirror": list(m), "points": [[str(X), str(Y)] for X, Y in pts]} for m, pts in zip(MIRRORS, imgs)],
        "images_meaning": "mirror (sx, sy) of c*: point k = (sx*(a_x + b_x*S2_LO) + 1/2 - r_t, sy*(a_y + b_y*S2_LO) + h/2 - r_t)",
    }
    cpath = os.path.join(OUT, f"{case}_cert.json")
    blob = json.dumps(cert, indent=1).encode()
    with open(cpath, "wb") as f:
        f.write(blob)
    tree = {
        "format": "exact-bnb-tree v1", "case": case, "authors": "Moki&Julio", "engine": "prove_small.py Engine (exact integer branch and bound), logged",
        "cert_sha256": hashlib.sha256(blob).hexdigest(),
        "N": n, "h": str(h), "r_t": str(r_t), "G": G, "W": str(W), "H": str(H), "d": str(d),
        "chart": "point frame p = c + (1/2 - r_t, h/2 - r_t) in [0,W]x[0,H]; grid X = G p_x / W, Y = G p_y / H; "
                 "box [X0,X1,Y0,Y1] = closed real box [X0 W/G, X1 W/G] x [Y0 H/G, Y1 H/G]",
        "labels": "points 0..N-1 (the same labels as c_star in the certificate, up to the permutation of an acceptance)",
        "cuts": ["X_i <= X_{i+1}, i = 0..N-2 (centred: x_i <= x_{i+1})", "X_0 + X_{N-1} <= G (centred: x_0 + x_{N-1} <= 0)",
                 "Y_0 <= G/2 (centred: y_0 <= 0)"],
        "counts": {"nodes": len(recs), "processed": popped, "accepted": acc, "discarded_leaves": leaves - acc,
                   "start_assignments_discarded": len(discards), "witnesses": wcount},
        "nodes": [recs[k] for k in range(len(recs))],
    }
    tpath = os.path.join(OUT, f"{case}_tree.json.gz")
    with gzip.open(tpath, "wt", encoding="utf-8") as f:
        json.dump(tree, f, separators=(",", ":"))
    print(f"  wrote {cpath} and {tpath}", flush=True)
    return status


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in CASES:
        print(__doc__); sys.exit(2)
    control = len(sys.argv) > 2
    st = main(sys.argv[1], sys.argv[2] if control else None)
    ok = st.startswith("UNDECIDED") if control else st == "COMPLETE"   # a control must NOT complete
    sys.exit(0 if ok else 1)
