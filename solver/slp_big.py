"""slp_circ.polish for BIG N (thousands of circles): the same trust-region SLP (HiGHS, a step is kept only if the TRUE radius
grows), but candidate pairs come from a KD-tree instead of an N x N matrix. Containers: ('rect', h) centred at the origin
(csq = ('rect', 1.0)) and ('circ',) = the unit disc centred at the origin (Packomania cci).
"""
import time, numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack
from scipy.spatial import cKDTree

def walls(c, cont):
    x, y = c[:, 0], c[:, 1]; n = len(c)
    if cont[0] == 'rect':
        h = cont[1]; g = np.stack([0.5 + x, 0.5 - x, h / 2 + y, h / 2 - y], 1)
        G = np.zeros((n, 4, 2)); G[:, 0, 0] = 1; G[:, 1, 0] = -1; G[:, 2, 1] = 1; G[:, 3, 1] = -1
    else:
        rr = np.sqrt(x * x + y * y) + 1e-300; g = (1 - rr)[:, None]
        G = np.zeros((n, 1, 2)); G[:, 0, 0] = -x / rr; G[:, 0, 1] = -y / rr
    return g, G

def rmin(c, cont):
    g, _ = walls(c, cont); d, _ = cKDTree(c).query(c, k=2); return float(min(g.min(), d[:, 1].min() / 2))

def repair(c, cont):
    c = c.copy()
    if cont[0] == 'rect': c[:, 0] = np.clip(c[:, 0], -0.5, 0.5); c[:, 1] = np.clip(c[:, 1], -cont[1] / 2, cont[1] / 2)
    else: rr = np.sqrt((c ** 2).sum(1)); m = rr > 1; c[m] /= rr[m, None]
    return c

def polish(c, cont, max_iter=300, t_cap=300.0, log=None, stop_at=None):   # 09-25: stop_at = stop once r >= it (claim bar + margin)
    c = np.asarray(c, dtype=np.float64); n = len(c); cur = rmin(c, cont); delta = 0.1 * cur; t0 = time.time()
    for it in range(max_iter):
        if time.time() - t0 > t_cap: break
        g, G = walls(c, cont)
        P = cKDTree(c).query_pairs(2 * (cur + 3 * delta) + 1e-12, output_type="ndarray"); I, J = P[:, 0], P[:, 1]; m = len(I)
        dv = c[I] - c[J]; d = np.sqrt((dv ** 2).sum(1)); U = dv / np.maximum(d, 1e-300)[:, None]
        r1 = np.repeat(np.arange(m), 5)
        c1 = np.stack([2 * I, 2 * I + 1, 2 * J, 2 * J + 1, np.full(m, 2 * n)], 1).ravel()
        v1 = np.stack([-U[:, 0] / 2, -U[:, 1] / 2, U[:, 0] / 2, U[:, 1] / 2, np.ones(m)], 1).ravel()
        A1 = coo_matrix((v1, (r1, c1)), shape=(m, 2 * n + 1)); b1 = d / 2
        wi, wk = np.nonzero(g <= cur + 3 * delta + 1e-12); mw = len(wi)
        r2 = np.repeat(np.arange(mw), 3)
        c2 = np.stack([2 * wi, 2 * wi + 1, np.full(mw, 2 * n)], 1).ravel()
        v2 = np.stack([-G[wi, wk, 0], -G[wi, wk, 1], np.ones(mw)], 1).ravel()
        A2 = coo_matrix((v2, (r2, c2)), shape=(mw, 2 * n + 1)); b2 = g[wi, wk]
        res = linprog(np.r_[np.zeros(2 * n), -1.0], A_ub=vstack([A1, A2]).tocsr(), b_ub=np.r_[b1, b2],
                      bounds=[(-delta, delta)] * (2 * n) + [(None, None)], method=('highs-ipm' if n >= 1500 else 'highs'), options={'time_limit': 120.0})   # 09-25: dual simplex took 166 s at N = 6336 (> the 60 s cap -> zero progress); IPM 18 s, same optimum
        if res.status != 0 or res.x is None or not np.all(np.isfinite(res.x)):
            delta /= 2
            if delta < 1e-16: break
            continue
        q = repair(c + res.x[:-1].reshape(n, 2), cont); new = rmin(q, cont)
        if new > cur + 1e-16:
            gain = new - cur; c, cur = q, new
            if log: log(it, cur, delta, time.time() - t0)
            if gain < 1e-15 and delta < 1e-12: break
            if stop_at is not None and cur >= stop_at: break
            delta = min(delta * 1.5, 0.2 * cur)
        else:
            delta /= 4
            if delta < 1e-16: break
    return c, float(cur)
