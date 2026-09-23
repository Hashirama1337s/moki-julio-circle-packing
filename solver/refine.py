"""High-precision refinement of a locally optimal packing (point form) to ~40+ digits.
1. Active set from the float64 solution: pairs with d_ij < d (1 + tol), wall contacts (x=0, y=0, x+y=1) within tol*d.
2. Unknowns: coordinates of every point + d. Equations: |pi-pj|^2 = d^2 (active pairs); x_i = 0 / y_i = 0 / x_i + y_i = 1 (walls).
3. Mixed-precision Gauss-Newton: residuals in mpmath (60 digits), minimum-norm step from float64 lstsq on the Jacobian; each
   iteration gains ~12 digits. Stop when the max residual < 1e-45.
4. Wall-active coordinates are then set EXACTLY (0, or y = 1 - x) so rounding can never push a point outside.
5. Output: 45-digit decimal certificate (point form) + the refined d (60 digits). Every non-active constraint is re-checked.
usage: py -3.11 refine.py cert/cand_N.txt [--tol 1e-9]   (or import refine_points)
"""
import sys, numpy as np, mpmath as mp
from fractions import Fraction as F
mp.mp.dps = 60

def active_set(p, tol):
    n = len(p); D = np.sqrt(((p[:, None] - p[None]) ** 2).sum(-1)); iu = np.triu_indices(n, 1); d = D[iu].min()
    sel = D[iu] < d * (1 + tol); pairs = list(zip(iu[0][sel], iu[1][sel]))
    walls = [(i, 0) for i in range(n) if p[i, 0] < tol * d] + [(i, 1) for i in range(n) if p[i, 1] < tol * d] + \
            [(i, 2) for i in range(n) if 1 - p[i].sum() < tol * d]
    return pairs, walls, d

def refine_points(p, tol=1e-9, iters=8, verbose=False):
    p = np.asarray(p, dtype=np.float64); n = len(p)
    pairs, walls, d0 = active_set(p, tol)
    X = [mp.mpf(float(v)) for v in p.ravel()] + [mp.mpf(float(d0)) ** 2]      # unknowns: coords..., D = d^2
    k = len(X); m = len(pairs) + len(walls)
    for it in range(iters):
        R = []
        for i, j in pairs:
            R.append((X[2*i] - X[2*j]) ** 2 + (X[2*i+1] - X[2*j+1]) ** 2 - X[-1])
        for i, w in walls:
            R.append(X[2*i] if w == 0 else X[2*i+1] if w == 1 else X[2*i] + X[2*i+1] - 1)
        res = max(abs(r) for r in R)
        if verbose: print(f"  iter {it}: max residual {mp.nstr(res, 3)} (equations {m}, unknowns {k})")
        if res < mp.mpf(10) ** -52: break
        J = np.zeros((m, k)); xf = np.array([float(v) for v in X])
        for r, (i, j) in enumerate(pairs):
            dx, dy = xf[2*i] - xf[2*j], xf[2*i+1] - xf[2*j+1]
            J[r, 2*i], J[r, 2*i+1], J[r, 2*j], J[r, 2*j+1], J[r, -1] = 2*dx, 2*dy, -2*dx, -2*dy, -1
        for r, (i, w) in enumerate(walls, start=len(pairs)):
            if w == 0: J[r, 2*i] = 1
            elif w == 1: J[r, 2*i+1] = 1
            else: J[r, 2*i] = J[r, 2*i+1] = 1
        rf = np.array([float(v) for v in R])
        step = np.linalg.lstsq(J, -rf, rcond=None)[0]
        X = [X[q] + mp.mpf(float(step[q])) for q in range(k)]
    # exact wall placement
    P = [[X[2*i], X[2*i+1]] for i in range(n)]
    for i, w in walls:
        if w == 0: P[i][0] = mp.mpf(0)
        elif w == 1: P[i][1] = mp.mpf(0)
    for i, w in walls:
        if w == 2: P[i][1] = 1 - P[i][0]
    # corners (two walls at once) are set exactly; otherwise fixing one wall can break the other by ~1e-60 (found at N=88)
    W = {}
    for i, w in walls: W.setdefault(i, set()).add(w)
    for i, s in W.items():
        if {0, 1} <= s: P[i] = [mp.mpf(0), mp.mpf(0)]
        elif {1, 2} <= s: P[i] = [mp.mpf(1), mp.mpf(0)]
        elif {0, 2} <= s: P[i] = [mp.mpf(0), mp.mpf(1)]
    dmin = min(mp.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2) for q, a in enumerate(P) for b in P[q+1:])
    inside = all(x >= 0 and y >= 0 and x + y <= 1 for x, y in P)
    return P, dmin, res, inside, (len(pairs), len(walls), k)

def write_cert(n, P, path, digits=45):
    """Decimal strings; wall-active coordinates are exact (0 or 1-x). Returns the path."""
    with open(path, "w") as f:
        f.write(f"N {n}\n")
        for x, y in P:
            xs = mp.nstr(x, digits, strip_zeros=False, min_fixed=-mp.inf, max_fixed=mp.inf) if x != 0 else "0"
            ys = mp.nstr(y, digits, strip_zeros=False, min_fixed=-mp.inf, max_fixed=mp.inf) if y != 0 else "0"
            f.write(f"{xs} {ys}\n")
    # exact repair for hypotenuse contacts after rounding: y := 1 - x as an exact decimal
    L = [l.split() for l in open(path) if l.strip()]; out = [L[0]]
    for a, b in L[1:]:
        X, Y = F(a), F(b)
        if X + Y > 1: Y = 1 - X
        X = max(X, F(0)); Y = max(Y, F(0))
        out.append([a, f"{Y.numerator}/{Y.denominator}" if F(b) != Y else b])
    with open(path, "w") as f:
        f.write(" ".join(out[0]) + "\n"); [f.write(f"{a} {b}\n") for a, b in out[1:]]
    return path

if __name__ == "__main__":
    tol = float(sys.argv[sys.argv.index("--tol") + 1]) if "--tol" in sys.argv else 1e-9
    for path in [a for a in sys.argv[1:] if a.endswith(".txt")]:
        L = [l.split() for l in open(path) if l.strip()]; n = int(L[0][1])
        p = np.array([[float(F(a)), float(F(b))] for a, b in L[1:]])
        P, dmin, res, inside, shape = refine_points(p, tol, verbose=True)
        out = path.replace("cand_", "hp_")
        write_cert(n, P, out)
        print(f"N={n}: refined d = {mp.nstr(dmin, 40)} residual {mp.nstr(res, 3)} inside {inside} (pairs, walls, unknowns) {shape} -> {out}")
