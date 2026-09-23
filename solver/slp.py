"""SLP local solver (numpy/scipy only, no torch) - see search.py for the GPU stages."""
import numpy as np

def polish(p, max_iter=300, delta0=None):
    """Sequential LP (HiGHS) with trust region, float64: maximise t s.t. d_ij + u_ij.(di - dj) >= t for near pairs, points stay in T.
    Accept a step only if the TRUE minimum distance grows; otherwise shrink the trust region. Converges fast at rigid optima."""
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    p = np.asarray(p, dtype=np.float64).copy(); n = len(p); iu = np.triu_indices(n, 1)
    def mind(q): return np.sqrt(((q[:, None] - q[None]) ** 2).sum(-1))[iu].min()
    p[:, 0] = np.clip(p[:, 0], 0, 1); p[:, 1] = np.clip(p[:, 1], 0, 1); o = p.sum(1) - 1; p[o > 0] -= o[o > 0, None] / 2
    d0 = np.sqrt(((p[:, None] - p[None]) ** 2).sum(-1))[iu]; bad = np.unique(iu[1][d0 < 1e-12])
    if len(bad):                  # coincident points (09-23 guard, as in slp_circ): move each a different distance toward (1/3,1/3)
        p[bad] += (np.array([1 / 3, 1 / 3]) - p[bad]) * np.linspace(1e-6, 1e-5, len(bad))[:, None]
    cur = mind(p); delta = delta0 or 0.05 * cur
    import time as _time; t_start = _time.time()
    for it in range(max_iter):
        if _time.time() - t_start > 90.0: break                   # wall-clock cap
        D = np.sqrt(((p[:, None] - p[None]) ** 2).sum(-1)); d = D[iu]
        sel = d <= cur + 4 * delta + 1e-12; I, J = iu[0][sel], iu[1][sel]; m = len(I)
        U = (p[I] - p[J]) / np.maximum(d[sel], 1e-300)[:, None]
        # rows: t - u.(di - dj) <= d_ij   (variables: dx0,dy0,...,t)
        r = np.repeat(np.arange(m), 5)
        c = np.stack([2 * I, 2 * I + 1, 2 * J, 2 * J + 1, np.full(m, 2 * n)], 1).ravel()
        v = np.stack([-U[:, 0], -U[:, 1], U[:, 0], U[:, 1], np.ones(m)], 1).ravel()
        A1 = coo_matrix((v, (r, c)), shape=(m, 2 * n + 1)); b1 = d[sel]
        # x+y + dx+dy <= 1
        A2 = coo_matrix((np.ones(2 * n), (np.repeat(np.arange(n), 2), np.arange(2 * n))), shape=(n, 2 * n + 1)); b2 = 1 - p.sum(1)
        from scipy.sparse import vstack
        lo = np.maximum(-delta, -p.ravel()); hi = np.full(2 * n, delta)
        res = linprog(np.r_[np.zeros(2 * n), -1.0], A_ub=vstack([A1, A2]).tocsr(), b_ub=np.r_[b1, b2],
                      bounds=list(zip(lo, hi)) + [(None, None)], method='highs', options={'time_limit': 20.0})
        if res.status != 0 or res.x is None or not np.all(np.isfinite(res.x)):
            delta /= 2
            if delta < 1e-16: break
            continue
        q = p + res.x[:-1].reshape(n, 2); q = np.clip(q, 0, None); o = q.sum(1) - 1; q[o > 0] -= o[o > 0, None] / 2
        new = mind(q)
        if new > cur + 1e-16:
            gain = new - cur; p, cur = q, new
            if gain < 1e-15 and delta < 1e-12: break
            delta = min(delta * 1.5, 0.1 * cur)
        else:
            delta /= 4
            if delta < 1e-16: break
    return p, float(cur)

