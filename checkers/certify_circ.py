# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Checker A: exact circle-form certificate for crc (rect:h), ccq (quad) and crt (tri). No floats in any decision.
File: line 1 "r <decimal>", then N lines "x y" (decimal strings). Record radius: 30-digit decimal from the shelf's radius.txt.
Fit rules (exact rationals; the only irrational terms are handled by sign checks + squaring):
  rect:h  -1/2 + r <= x <= 1/2 - r,  -h/2 + r <= y <= h/2 - r
  quad    x >= r, y >= r, 1 - r >= 0 and x^2 + y^2 <= (1 - r)^2
  semi    y >= r, 1 - r >= 0 and x^2 + y^2 <= (1 - r)^2          (unit semicircle, Packomania csc)
  tri     x >= r, y >= r, 1 - x - y >= 0 and (1 - x - y)^2 >= 2 r^2
Pairs: dx^2 + dy^2 >= 4 r^2. CLAIM RULE: IMPROVES requires r > r_rec (1 + 1e-10). Smaller gains are TIES (published values can be low by ~1e-24).
usage: py -3.11 certify_circ.py <container> <file> <record_radius_decimal>
"""
import sys
from fractions import Fraction as F

def fits(cont, x, y, r):
    if cont.startswith("rect:"):
        h = F(cont.split(":")[1])
        return -F(1, 2) + r <= x <= F(1, 2) - r and -h / 2 + r <= y <= h / 2 - r
    if cont == "quad":
        return x >= r and y >= r and 1 - r >= 0 and x * x + y * y <= (1 - r) ** 2
    if cont == "semi":
        return y >= r and 1 - r >= 0 and x * x + y * y <= (1 - r) ** 2
    if cont == "tri":
        s = 1 - x - y
        return x >= r and y >= r and s >= 0 and s * s >= 2 * r * r
    raise ValueError(cont)

def check(cont, path, rec, verbose=True):
    L = [l.split() for l in open(path) if l.strip()]
    assert L[0][0] == "r"; r = F(L[0][1]); C = [(F(a), F(b)) for a, b in L[1:]]; rec = F(rec)
    fit = r > 0 and all(fits(cont, x, y, r) for x, y in C)
    sep = all((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 >= 4 * r * r for i, a in enumerate(C) for b in C[i + 1:])
    if not (fit and sep): v = "INVALID"
    elif r - rec > rec * F(1, 10 ** 10): v = "IMPROVES"   # CLAIM RULE: relative gain > 1e-10
    elif r >= rec: v = "TIES"
    else: v = "WORSE"
    if verbose:
        print(f"{cont} N={len(C)}: fit {'PASS' if fit else 'FAIL'}, separation {'PASS' if sep else 'FAIL'}; r - rec = {float(r - rec):+.3e} (rel {float(r / rec - 1):+.3e})")
        print(f"VERDICT: {v}")
    return v, len(C), float(r / rec - 1)

if __name__ == "__main__":
    check(sys.argv[1], sys.argv[2], sys.argv[3])
