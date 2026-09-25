"""Checker A: exact checker for BIG packings (Packomania csq / cci frames), O(N): every number is read as an exact decimal and
scaled to a common integer grid (x = X / S), so every decision is integer arithmetic.
  square: the unit square centred at the origin (-1/2 <= x, y <= 1/2)   circle: the unit disc centred at the origin
Pairs: cells of side G = 2R (integers, floor division). If two centres lie in cells whose index differs by >= 2 in x (or y), then
|Xi - Xj| > G = 2R, so the pair is separated; hence checking every pair in the same or an adjacent cell (3 x 3 block) covers ALL
pairs. A half-open cell of side 2R holds at most 4 centres at mutual distance >= 2R, so > 4 in one cell is an immediate overlap.
usage: py -3.11 certify_big.py <square|circle> <file> <N> <record_radius>   -> "VERDICT: IMPROVES | VALID_NOT_BETTER | INVALID <why>"
file: line 1 "r <radius>", then N lines "x y" (plain decimals).
"""
import sys, collections
from fractions import Fraction as F

def dec(s):
    s = s.strip(); neg = s.startswith("-"); s = s.lstrip("+-")
    if "e" in s.lower() or not s.replace(".", "", 1).isdigit(): raise ValueError("not a plain decimal: " + s)
    a, _, b = s.partition("."); return (-1 if neg else 1) * int(a + b or "0"), len(b)

def check(container, path, n, rec):
    L = [l.split() for l in open(path) if l.strip()]
    if L[0][0] != "r" or len(L[0]) != 2: return "INVALID header"
    rows = L[1:]
    if len(rows) != n: return f"INVALID count {len(rows)} != {n}"
    if any(len(t) != 2 for t in rows): return "INVALID row format"
    vals = [dec(L[0][1])] + [dec(v) for t in rows for v in t]
    D = max(d for _, d in vals); S = 10 ** D
    I = [m * 10 ** (D - d) for m, d in vals]                   # exact integers: value = I / S
    R = I[0]; X = I[1::2]; Y = I[2::2]
    if R <= 0: return "INVALID r <= 0"
    for i in range(n):
        x, y = X[i], Y[i]
        if container == "square":
            if 2 * (abs(x) + R) > S or 2 * (abs(y) + R) > S: return f"INVALID circle {i + 1} outside the square"
        elif container == "circle":
            if R > S or x * x + y * y > (S - R) ** 2: return f"INVALID circle {i + 1} outside the disc"
        else: return "INVALID container"
    G = 2 * R; cells = collections.defaultdict(list)
    for i in range(n): cells[(X[i] // G, Y[i] // G)].append(i)
    if max(len(v) for v in cells.values()) > 4: return "INVALID more than 4 centres in one 2r cell (overlap)"
    R4 = 4 * R * R; pairs = 0
    for (cx, cy), mine in cells.items():
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                other = cells.get((cx + dx, cy + dy))
                if not other: continue
                for i in mine:
                    for j in other:
                        if j <= i: continue
                        pairs += 1
                        if (X[i] - X[j]) ** 2 + (Y[i] - Y[j]) ** 2 < R4: return f"INVALID overlap {i + 1} {j + 1}"
    return "IMPROVES" if F(R, S) > F(rec) else "VALID_NOT_BETTER"

def selftest():
    import os, tempfile
    def run(cont, r, pts, rec="0.0"):
        p = os.path.join(tempfile.gettempdir(), "certify_big_selftest.txt")
        open(p, "w").write(f"r {r}\n" + "".join(f"{x} {y}\n" for x, y in pts)); return check(cont, p, len(pts), rec)
    assert run("square", "0.25", [("-0.25", "-0.25"), ("0.25", "-0.25"), ("-0.25", "0.25"), ("0.25", "0.25")], "0.2") == "IMPROVES"
    assert run("square", "0.25", [("-0.25", "-0.25"), ("0.25", "-0.25"), ("-0.25", "0.25"), ("0.25", "0.25")], "0.25") == "VALID_NOT_BETTER"
    assert run("square", "0.25", [("-0.25", "-0.25"), ("0.2499999999999999999999999999999999999999", "-0.25")]).startswith("INVALID overlap")
    assert run("square", "0.1", [("0.4000000000000000000000000000000000000001", "0")]).startswith("INVALID circle")
    assert run("circle", "0.5", [("-0.5", "0"), ("0.5", "0")]) == "IMPROVES"                 # tangent to each other and to the disc
    assert run("circle", "0.5", [("-0.5", "0"), ("0.5000000000000000000000000000000000000001", "0")]).startswith("INVALID circle")
    far = [("-0.4", "-0.4"), ("0.4", "0.4"), ("-0.4", "-0.4000000000000000000000000000000000000001")]
    assert run("square", "0.01", far).startswith("INVALID"), "overlap far away in the list must be caught"
    print("selftest OK")

if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]: selftest()
    else: print("VERDICT:", check(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]))
