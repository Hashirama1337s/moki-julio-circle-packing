"""Claude's exact certificate (point form). No floats in any decision.
A candidate is a text file: line 1 "N <n>", then n lines "x y" (decimal strings) = points in the unit triangle
T = {x >= 0, y >= 0, x + y <= 1}. The packing value is d = min pairwise distance; circles: r = d / (2 + (2 + sqrt2) d), a strictly
increasing function of d, so "better than the record" <=> d^2 > d_rec^2 (compared exactly as rationals).
Record d_rec is the 30-digit decimal from Packomania's distance.txt, treated as exact; CLAIM RULE: IMPROVES only if the relative
gain exceeds 1e-10 — last-digit rounding / under-convergence of published values (seen at 1e-24) is a TIE, not a record.
usage: py -3.11 certify.py candidate.txt   ->   prints checks and VERDICT: IMPROVES / TIES / WORSE / INVALID
"""
import sys, os
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))

def record_d(n):
    for line in open(os.path.join(HERE, "data", "distance.txt")):
        p = line.split()
        if len(p) == 2 and int(p[0]) == n: return F(p[1])
    raise KeyError(n)

def load(path):
    L = [l.split() for l in open(path) if l.strip()]
    assert L[0][0] == "N"; n = int(L[0][1]); pts = [(F(a), F(b)) for a, b in L[1:]]
    return n, pts

def check(path, verbose=True):
    n, pts = load(path); out = {}
    out['count'] = len(pts) == n
    out['inside'] = all(x >= 0 and y >= 0 and x + y <= 1 for x, y in pts)
    d2 = min((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 for i, a in enumerate(pts) for b in pts[i + 1:])
    rec = record_d(n); rec2 = rec * rec
    margin2 = d2 - rec2
    if not (out['count'] and out['inside']): verdict = "INVALID"
    elif margin2 > 2 * rec2 * F(1, 10 ** 10): verdict = "IMPROVES"      # claim floor: relative gain > 1e-10 (see CLAIM RULE)
    elif margin2 >= 0: verdict = "TIES"
    else: verdict = "WORSE"
    if verbose:
        print(f"N={n}: count {'PASS' if out['count'] else 'FAIL'}, inside T {'PASS' if out['inside'] else 'FAIL'}")
        print(f"  exact min d^2 = {float(d2):.18e}  (d ~ {float(d2) ** 0.5:.18f})")
        print(f"  record d      = {float(rec):.18f}")
        print(f"  d^2 - rec^2   = {float(margin2):+.3e}   (relative d gain ~ {float(margin2 / (2 * rec2)):+.3e})")
        print(f"VERDICT: {verdict}")
    return verdict, n, float(margin2 / (2 * rec2))

if __name__ == "__main__":
    for p in sys.argv[1:]: check(p)
