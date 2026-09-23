"""'Breathing radii' SLP (recon 09-23): radii are free in a band, then the band is annealed to zero.
Variables: centres c (2n), radii rho (n), t.  Constraints: |ci - cj| >= rho_i + rho_j (linearised: conservative, the distance is
convex), wall slack >= rho_i (straight walls exact; the quadrant arc is checked exactly after each step), t <= rho_i <= (1+beta) t.
Objective: maximise mean(rho).  beta = 0 is exactly the equal-circle problem (all rho_i = t).  A step is kept only if it is
TRULY feasible and the mean grows.  On each stage change the radii are clipped DOWN to (1+beta) min(rho): never an overlap.
Every run ends with the standard equal-radius polish (slp_circ.polish), so arms are compared on the same final objective.
"""
import time, numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack, hstack, csr_matrix
import slp_circ

SCHEDULE = (0.3, 0.15, 0.08, 0.04, 0.02, 0.01, 0.005, 0.0)

def feasible(c, rho, cont, tol=0.0):
    g, _ = slp_circ.walls(c, cont)
    if np.any(g - rho[:, None] < -tol): return False
    n = len(c); iu = np.triu_indices(n, 1)
    d = np.sqrt(((c[:, None] - c[None]) ** 2).sum(-1))[iu]
    return bool(np.all(d - rho[iu[0]] - rho[iu[1]] >= -tol))

def stage(c, rho, cont, beta, max_iter=40, t_cap=30.0):
    n = len(c); iu = np.triu_indices(n, 1); t0 = time.time()
    delta = 0.1 * rho.min(); cur = rho.mean()
    for it in range(max_iter):
        if time.time() - t0 > t_cap: break
        g, G = slp_circ.walls(c, cont)
        D = np.sqrt(((c[:, None] - c[None]) ** 2).sum(-1)); d = D[iu]
        sl = d - rho[iu[0]] - rho[iu[1]]
        sel = sl <= 4 * delta + 1e-12; I, J = iu[0][sel], iu[1][sel]; m = len(I)
        U = (c[I] - c[J]) / np.maximum(d[sel], 1e-300)[:, None]
        # vars: dc (2n) | drho (n) | t      pair:  -u.dci + u.dcj + drho_i + drho_j <= d - rho_i - rho_j
        nv = 3 * n + 1
        r1 = np.repeat(np.arange(m), 6)
        c1 = np.stack([2 * I, 2 * I + 1, 2 * J, 2 * J + 1, 2 * n + I, 2 * n + J], 1).ravel()
        v1 = np.stack([-U[:, 0], -U[:, 1], U[:, 0], U[:, 1], np.ones(m), np.ones(m)], 1).ravel()
        A1 = coo_matrix((v1, (r1, c1)), shape=(m, nv)); b1 = sl[sel]
        ws = g - rho[:, None]; wi, wk = np.nonzero(ws <= 4 * delta + 1e-12); mw = len(wi)
        r2 = np.repeat(np.arange(mw), 3)
        c2 = np.stack([2 * wi, 2 * wi + 1, 2 * n + wi], 1).ravel()
        v2 = np.stack([-G[wi, wk, 0], -G[wi, wk, 1], np.ones(mw)], 1).ravel()
        A2 = coo_matrix((v2, (r2, c2)), shape=(mw, nv)); b2 = ws[wi, wk]
        # band:  t - rho_i - drho_i <= 0  ->  -drho_i + t <= rho_i ;  rho_i + drho_i - (1+beta) t <= 0  ->  drho_i - (1+beta) t <= -rho_i
        ar = np.arange(n)
        A3 = coo_matrix((np.r_[-np.ones(n), np.ones(n)], (np.r_[ar, ar], np.r_[2 * n + ar, np.full(n, nv - 1)])), shape=(n, nv))
        A4 = coo_matrix((np.r_[np.ones(n), np.full(n, -(1 + beta))], (np.r_[ar, ar], np.r_[2 * n + ar, np.full(n, nv - 1)])), shape=(n, nv))
        A = vstack([A1, A2, A3, A4]).tocsr(); b = np.r_[b1, b2, rho, -rho]
        cost = np.r_[np.zeros(2 * n), -np.ones(n) / n, 0.0]
        res = linprog(cost, A_ub=A, b_ub=b, bounds=[(-delta, delta)] * (3 * n) + [(0, None)], method="highs",
                      options={"time_limit": 20.0})
        if res.status != 0 or res.x is None or not np.all(np.isfinite(res.x)):
            delta /= 2
            if delta < 1e-14: break
            continue
        cn = slp_circ.repair(c + res.x[:2 * n].reshape(n, 2), cont); rn = rho + res.x[2 * n:3 * n]
        # LP tolerance can leave ~1e-9 overlaps: shrink all radii by the measured excess, then re-clip the band
        gn, _ = slp_circ.walls(cn, cont); Dn = np.sqrt(((cn[:, None] - cn[None]) ** 2).sum(-1))[iu]
        exc = max(float((rn[:, None] - gn).max()), float((rn[iu[0]] + rn[iu[1]] - Dn).max()))
        if 0 < exc < 1e-7 * rn.mean(): rn = rn - (exc + 1e-16)
        rn = np.minimum(rn, (1 + beta) * rn.min())                              # exact band after rounding
        if rn.min() > 0 and feasible(cn, rn, cont) and rn.mean() > cur + 1e-15:
            gain = rn.mean() - cur; c, rho, cur = cn, rn, rn.mean()
            if gain < 1e-13 * cur and delta < 1e-10: break
            delta = min(delta * 1.5, 0.2 * rho.min())
        else:
            delta /= 4
            if delta < 1e-14: break
    return c, rho

def start_radii(c, cont):
    g, _ = slp_circ.walls(c, cont); n = len(c)
    D = np.sqrt(((c[:, None] - c[None]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf)
    rho = np.minimum(g.min(1), D.min(1) / 2) * 0.999
    return np.maximum(rho, 1e-9)

def breathe(c, cont, schedule=SCHEDULE, max_iter=40):
    c = slp_circ.unstick(slp_circ.repair(np.asarray(c, dtype=np.float64), cont), cont)
    rho = start_radii(c, cont)
    for beta in schedule:
        rho = np.minimum(rho, (1 + beta) * rho.min())
        c, rho = stage(c, rho, cont, beta, max_iter=max_iter)
    return slp_circ.polish(c, cont)
