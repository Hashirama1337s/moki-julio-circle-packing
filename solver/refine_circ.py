"""High-precision refinement, circle form (tri / rect:h / quad), mixed-precision Gauss-Newton on the contact graph.
Unknowns: centres + r. Equations: active pairs |ci-cj|^2 = 4 r^2; active walls:
  tri   x = r, y = r, 1 - x - y = sqrt2 r            rect  x + 1/2 = r, 1/2 - x = r, y + h/2 = r, h/2 - y = r
  quad  x = r, y = r, x^2 + y^2 = (1 - r)^2
Residuals at 60 digits, min-norm float64 steps (~12 digits per iteration). The certificate then takes the ROUNDED centres and the
largest radius they admit (computed at 80 digits, rounded down at 45) — so rounding can only lower r_c, never make it invalid.
"""
import numpy as np, mpmath as mp
mp.mp.dps = 60

def slacks_float(c, r, cont):
    x, y = c[:, 0], c[:, 1]
    if cont[0] == 'tri': return np.stack([x - r, y - r, (1 - x - y) / np.sqrt(2) - r], 1)
    if cont[0] == 'rect':
        h = cont[1]; return np.stack([x + 0.5 - r, 0.5 - x - r, y + h / 2 - r, h / 2 - y - r], 1)
    if cont[0] == 'quad': return np.stack([x - r, y - r, 1 - np.sqrt(x * x + y * y) - r], 1)
    if cont[0] == 'semi': return np.stack([np.full(len(x), 10.0), y - r, 1 - np.sqrt(x * x + y * y) - r], 1)   # slot 0 inert

def refine(c, cont, cont_s, tol=1e-9, iters=8):
    c = np.asarray(c, dtype=np.float64); n = len(c); iu = np.triu_indices(n, 1)
    D = np.sqrt(((c[:, None] - c[None]) ** 2).sum(-1)); d = D[iu]
    S = slacks_float(c, 0, cont); r0 = min(d.min() / 2, S.min())
    pairs = [(i, j) for i, j, dd in zip(iu[0], iu[1], d) if dd / 2 < r0 * (1 + tol)]
    Sr = slacks_float(c, r0, cont); walls = [(i, k) for i, k in zip(*np.nonzero(Sr < tol * r0))]
    h = mp.mpf(cont_s.split(":")[1]) if cont[0] == 'rect' else None; s2 = mp.sqrt(2)
    X = [mp.mpf(float(v)) for v in c.ravel()] + [mp.mpf(float(r0))]
    m, k = len(pairs) + len(walls), 2 * n + 1
    for it in range(iters):
        R = []; J = np.zeros((m, k)); rr = X[-1]; xf = np.array([float(v) for v in X])
        for q, (i, j) in enumerate(pairs):
            dx, dy = X[2*i] - X[2*j], X[2*i+1] - X[2*j+1]; R.append(dx * dx + dy * dy - 4 * rr * rr)
            fx, fy = float(dx), float(dy)
            J[q, 2*i], J[q, 2*i+1], J[q, 2*j], J[q, 2*j+1], J[q, -1] = 2*fx, 2*fy, -2*fx, -2*fy, -8 * xf[-1]
        for q, (i, w) in enumerate(walls, start=len(pairs)):
            x, y = X[2*i], X[2*i+1]
            if cont[0] == 'tri':
                e = [x - rr, y - rr, 1 - x - y - s2 * rr][w]; g = [(1, 0, -1), (0, 1, -1), (-1, -1, -float(s2))][w]
            elif cont[0] == 'rect':
                e = [x + mp.mpf(1)/2 - rr, mp.mpf(1)/2 - x - rr, y + h/2 - rr, h/2 - y - rr][w]
                g = [(1, 0, -1), (-1, 0, -1), (0, 1, -1), (0, -1, -1)][w]
            else:
                if w < 2: e = [x - rr, y - rr][w]; g = [(1, 0, -1), (0, 1, -1)][w]
                else: e = x * x + y * y - (1 - rr) ** 2; g = (2 * float(x), 2 * float(y), 2 * (1 - xf[-1]))
            R.append(e); J[q, 2*i], J[q, 2*i+1], J[q, -1] = g
        res = max(abs(v) for v in R)
        if res < mp.mpf(10) ** -52: break
        step = np.linalg.lstsq(J, -np.array([float(v) for v in R]), rcond=None)[0]
        X = [X[q] + mp.mpf(float(step[q])) for q in range(k)]
    C = [(X[2*i], X[2*i+1]) for i in range(n)]
    return C, X[-1], res, (len(pairs), len(walls), k)

def admissible_r(C, cont_s):
    mp.mp.dps = 80
    if cont_s.startswith("rect:"):
        h = mp.mpf(cont_s.split(":")[1]); w = min(min(x + mp.mpf(1)/2, mp.mpf(1)/2 - x, y + h/2, h/2 - y) for x, y in C)
    elif cont_s == "quad": w = min(min(x, y, 1 - mp.sqrt(x*x + y*y)) for x, y in C)
    elif cont_s == "semi": w = min(min(y, 1 - mp.sqrt(x*x + y*y)) for x, y in C)
    else: w = min(min(x, y, (1 - x - y) / mp.sqrt(2)) for x, y in C)
    pr = min(mp.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2) for i, a in enumerate(C) for b in C[i+1:]) / 2
    return min(w, pr)

def write_hp(C, cont_s, path, digits=45):
    """Round centres to `digits`, then r_c = admissible radius of the ROUNDED centres, rounded down at `digits`."""
    Cs = [(mp.nstr(x, digits, strip_zeros=False, min_fixed=-mp.inf, max_fixed=mp.inf),
           mp.nstr(y, digits, strip_zeros=False, min_fixed=-mp.inf, max_fixed=mp.inf)) for x, y in C]
    mp.mp.dps = 80
    rc = admissible_r([(mp.mpf(a), mp.mpf(b)) for a, b in Cs], cont_s)
    rc = mp.floor(rc * mp.mpf(10) ** digits) / mp.mpf(10) ** digits
    with open(path, "w") as f:
        f.write(f"r {mp.nstr(rc, digits, strip_zeros=False, min_fixed=-mp.inf, max_fixed=mp.inf)}\n")
        for a, b in Cs: f.write(f"{a} {b}\n")
    mp.mp.dps = 60
    return rc
