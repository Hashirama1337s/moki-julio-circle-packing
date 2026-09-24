"""Mixed-precision SLP: drive a packing to its EXACT local optimum (80-digit centres), then write a 45-digit certificate.
Why (09-23): lopt.py showed that many 'hp' certificates were Gauss-Newton solutions of a PARTIAL contact graph -- the float
search had stopped ~1e-9 short, and refine_circ tightened only the contacts it had found. Feasible and record-beating, but not
at the local peak.
Scheme: objective r(c) = min_a b_a(c), with b_a = half pair distance or distance to a wall (as slp_circ.rmin). Each step solves,
in float, the LP  max v  s.t.  (b_a - r)/D + grad b_a . u >= v  (|u| <= 1, step = D u), for every constraint within 4D of r,
with the slacks b_a - r computed in 80-digit arithmetic and SCALED by the trust radius D, so the float LP always works on O(1)
numbers; a tiny L1 penalty keeps free (flex) circles still. A step is kept only if the TRUE r (80 digits) grows.
Constraints farther than 4D from r cannot become the minimum in one step (|grad b . step| <= sqrt2 D).
"""
import numpy as np, mpmath as mp
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, hstack, vstack, csr_matrix, identity
mp.mp.dps = 80

class Box:
    def __init__(self, cont_s):
        self.kind = "rect" if cont_s.startswith("rect:") else cont_s
        self.h = mp.mpf(cont_s.split(":")[1]) if self.kind == "rect" else None
        self.hf = float(self.h) if self.h is not None else None
        self.nw = {"tri": 3, "rect": 4, "quad": 3, "semi": 3}[self.kind]

    def walls_f(self, cf):
        x, y = cf[:, 0], cf[:, 1]
        if self.kind == "tri": return np.stack([x, y, (1 - x - y) / np.sqrt(2)], 1)
        if self.kind == "rect": return np.stack([x + 0.5, 0.5 - x, y + self.hf / 2, self.hf / 2 - y], 1)
        if self.kind == "semi": return np.stack([np.full(len(x), 10.0), y, 1 - np.hypot(x, y)], 1)   # slot 0 inert
        return np.stack([x, y, 1 - np.hypot(x, y)], 1)

    def wall_mp(self, x, y, k):
        """(b, grad) in 80 digits / float for wall k of a centre (x, y)."""
        if self.kind == "tri":
            return [(x, (1.0, 0.0)), (y, (0.0, 1.0)), ((1 - x - y) / mp.sqrt(2), (-2 ** -0.5, -2 ** -0.5))][k]
        if self.kind == "rect":
            h = self.h
            return [(x + mp.mpf(1) / 2, (1.0, 0.0)), (mp.mpf(1) / 2 - x, (-1.0, 0.0)), (y + h / 2, (0.0, 1.0)), (h / 2 - y, (0.0, -1.0))][k]
        if k == 0: return (mp.mpf(10), (0.0, 0.0)) if self.kind == "semi" else (x, (1.0, 0.0))
        if k == 1: return y, (0.0, 1.0)
        rr = mp.sqrt(x * x + y * y); return 1 - rr, (-float(x / rr), -float(y / rr))

def rmin_mp(C, box, pairs=None):
    """Admissible radius of centres C (80 digits). pairs: candidate pairs (others are known to be far)."""
    n = len(C); w = min(box.wall_mp(x, y, k)[0] for x, y in C for k in range(box.nw))
    if pairs is None: pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    p = min((mp.sqrt((C[i][0] - C[j][0]) ** 2 + (C[i][1] - C[j][1]) ** 2) / 2 for i, j in pairs), default=mp.inf)
    return min(w, p)

def step_lp(C, box, r, D, eps=1e-6):
    n = len(C); cf = np.array([[float(x), float(y)] for x, y in C]); rf = float(r)
    reach = 4 * D + 1e-13 * rf
    Dm = np.sqrt(((cf[:, None] - cf[None]) ** 2).sum(-1)); iu = np.triu_indices(n, 1)
    selp = Dm[iu] / 2 - rf <= reach; I, J = iu[0][selp], iu[1][selp]
    Wf = box.walls_f(cf); wi, wk = np.nonzero(Wf - rf <= reach)
    rows, cols, vals, rhs, cons = [], [], [], [], []
    for q, (i, j) in enumerate(zip(I, J)):
        dx, dy = C[i][0] - C[j][0], C[i][1] - C[j][1]; d = mp.sqrt(dx * dx + dy * dy)
        ux, uy = float(dx / d), float(dy / d)
        rows += [q] * 4; cols += [2 * i, 2 * i + 1, 2 * j, 2 * j + 1]; vals += [ux / 2, uy / 2, -ux / 2, -uy / 2]
        rhs.append(float((d / 2 - r) / D)); cons.append(("p", int(i), int(j)))
    m0 = len(rhs)
    for q, (i, k) in enumerate(zip(wi, wk), start=m0):
        b, (gx, gy) = box.wall_mp(C[i][0], C[i][1], int(k))
        rows += [q, q]; cols += [2 * i, 2 * i + 1]; vals += [gx, gy]; rhs.append(float((b - r) / D)); cons.append(("w", int(i), int(k)))
    m = len(rhs); Gb = coo_matrix((vals, (rows, cols)), shape=(m, 2 * n)).tocsr()
    # LEXICOGRAPHIC (09-23: a fixed L1 penalty hid a real ascent of rate 8.6e-7 on ccq 202):
    # stage 1  max v  s.t.  -(G u) + v <= slack, |u| <= 1;   stage 2  min sum|u|  s.t. the same and  v >= v1 - 1e-9 |v1|
    opts = {"time_limit": 60.0, "primal_feasibility_tolerance": 1e-10, "dual_feasibility_tolerance": 1e-10}
    A1 = hstack([-Gb, csr_matrix(np.ones((m, 1)))]).tocsr(); b = np.array(rhs)
    r1 = linprog(np.r_[np.zeros(2 * n), -1.0], A_ub=A1, b_ub=b, bounds=[(-1, 1)] * (2 * n) + [(None, None)], method="highs", options=opts)
    if r1.status != 0 or r1.x is None or not np.all(np.isfinite(r1.x)): return None, None, list(zip(I.tolist(), J.tolist()))
    v1 = float(r1.x[-1])
    A2 = vstack([hstack([-Gb, Gb, csr_matrix(np.ones((m, 1)))]),
                 csr_matrix(np.r_[np.zeros(4 * n), -1.0][None, :])]).tocsr()
    r2 = linprog(np.r_[np.ones(4 * n), 0.0], A_ub=A2, b_ub=np.r_[b, -(v1 - 1e-9 * abs(v1))],
                 bounds=[(0, 1)] * (4 * n) + [(None, None)], method="highs", options=opts)
    if r2.status == 0 and r2.x is not None and np.all(np.isfinite(r2.x)):
        u = r2.x[:2 * n] - r2.x[2 * n:4 * n]; v = float(r2.x[-1])
    else:
        u = r1.x[:2 * n]; v = v1
    return u, v, list(zip(I.tolist(), J.tolist()))

def converge(C, cont_s, D0=1e-7, max_iter=600, log=None, floor_rel=1e-22, t_cap=None):
    """SLP phase: C list of (mpf, mpf) -> centres within ~1e-17 of the local optimum, where the contact set is identified.
    (Below that, the float LP cannot resolve the scaled slacks; refine_circ's Newton step on the identified contacts finishes.)
    Returns (C, r, info)."""
    box = Box(cont_s); C = [(mp.mpf(x), mp.mpf(y)) for x, y in C]
    mp.mp.dps = 80
    r = rmin_mp(C, box, near_pairs(C, box)); D = mp.mpf(D0); it = acc = 0; r0 = r; floor = mp.mpf(floor_rel) * r
    import time as _t; _t0 = _t.time()
    while it < max_iter and D > floor:
        if t_cap and _t.time() - _t0 > t_cap: break                  # 23:50: one ccq 551 job ran 30+ min (best-so-far kept)
        it += 1
        u, v, pairs = step_lp(C, box, r, D)
        if u is None or v <= 1e-9 or v * D < floor * mp.mpf(10) ** -4: D /= 16; continue   # nothing resolvable at this scale
        Cn = [(x + D * mp.mpf(float(u[2 * i])), y + D * mp.mpf(float(u[2 * i + 1]))) for i, (x, y) in enumerate(C)]
        rn = rmin_mp(Cn, box, pairs)
        if rn > r:
            acc += 1; step = max(abs(float(t)) for t in u) * D
            C, r = Cn, rn
            if log: log(f"it {it}: r gain {mp.nstr(rn - r0, 5)}, step {mp.nstr(step, 3)}, D {mp.nstr(D, 3)}, v {v:.3e}")
            D = 2 * D if max(abs(u)) > 0.999 else max(4 * step, floor)
        else:
            D /= 4
    rf = rmin_mp(C, box, near_pairs(C, box))
    return C, rf, {"iters": it, "accepted": acc, "gain": rf - r0, "final_D": D}

def near_pairs(C, box, rel=1e-6):
    cf = np.array([[float(x), float(y)] for x, y in C]); n = len(C)
    Dm = np.sqrt(((cf[:, None] - cf[None]) ** 2).sum(-1)); iu = np.triu_indices(n, 1)
    rf = min(float(box.walls_f(cf).min()), float(Dm[iu].min()) / 2)
    sel = Dm[iu] / 2 <= rf * (1 + rel)
    return list(zip(iu[0][sel].tolist(), iu[1][sel].tolist()))
