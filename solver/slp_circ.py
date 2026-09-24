"""Circle-form SLP for equal circles in a container: maximise r(c) = min( wall slacks of every centre, half of every pair distance ).
Containers:
  ('tri',)      Packomania crt: x >= r, y >= r, (1 - x - y)/sqrt2 >= r                (legs 1, right angle at the origin)
  ('rect', h)   Packomania crc: width 1, height h, centred at the origin: 1/2 -+ x >= r, h/2 -+ y >= r
  ('quad',)     Packomania ccq: x >= r, y >= r, 1 - |c| >= r                            (unit quarter disc)
  ('semi',)     Packomania csc: y >= r, 1 - |c| >= r            (unit semicircle, y >= 0; wall slot 0 is an inert dummy)
Sequential LP with a trust region (HiGHS); a step is accepted only if the TRUE r(c) grows. Same scheme as slp.py (point form).
"""
import time, numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack
S2 = np.sqrt(2.0)

def walls(c, cont):
    """(n, k) wall slacks and their gradients (n, k, 2)."""
    x, y = c[:, 0], c[:, 1]; n = len(c)
    if cont[0] == 'tri':
        g = np.stack([x, y, (1 - x - y) / S2], 1)
        G = np.zeros((n, 3, 2)); G[:, 0, 0] = 1; G[:, 1, 1] = 1; G[:, 2, :] = -1 / S2
    elif cont[0] == 'rect':
        h = cont[1]
        g = np.stack([0.5 + x, 0.5 - x, h / 2 + y, h / 2 - y], 1)
        G = np.zeros((n, 4, 2)); G[:, 0, 0] = 1; G[:, 1, 0] = -1; G[:, 2, 1] = 1; G[:, 3, 1] = -1
    elif cont[0] == 'quad':
        rr = np.sqrt(x * x + y * y) + 1e-300
        g = np.stack([x, y, 1 - rr], 1)
        G = np.zeros((n, 3, 2)); G[:, 0, 0] = 1; G[:, 1, 1] = 1; G[:, 2, 0] = -x / rr; G[:, 2, 1] = -y / rr
    elif cont[0] == 'semi':                 # slot 0 = inert dummy (slack 10, never active) so the arc keeps index 2 as in quad
        rr = np.sqrt(x * x + y * y) + 1e-300
        g = np.stack([np.full(n, 10.0), y, 1 - rr], 1)
        G = np.zeros((n, 3, 2)); G[:, 1, 1] = 1; G[:, 2, 0] = -x / rr; G[:, 2, 1] = -y / rr
    return g, G

def rmin(c, cont):
    g, _ = walls(c, cont); n = len(c); iu = np.triu_indices(n, 1)
    d = np.sqrt(((c[:, None] - c[None]) ** 2).sum(-1))[iu]
    return min(g.min(), d.min() / 2 if n > 1 else np.inf)

def repair(c, cont):
    """Keep centres inside the container (so wall slacks are >= 0)."""
    c = c.copy()
    if cont[0] == 'tri':
        c = np.clip(c, 0, 1); o = c.sum(1) - 1; c[o > 0] -= o[o > 0, None] / 2
    elif cont[0] == 'rect':
        c[:, 0] = np.clip(c[:, 0], -0.5, 0.5); c[:, 1] = np.clip(c[:, 1], -cont[1] / 2, cont[1] / 2)
    elif cont[0] == 'quad':
        c = np.clip(c, 0, None); rr = np.sqrt((c ** 2).sum(1)); m = rr > 1; c[m] /= rr[m, None]
    elif cont[0] == 'semi':
        c[:, 1] = np.clip(c[:, 1], 0, None); rr = np.sqrt((c ** 2).sum(1)); m = rr > 1; c[m] /= rr[m, None]
    return c

def unstick(c, cont, seed=0):
    """Separate exactly-coincident centres (repair() can clip several onto the same corner) and pull everything a hair inward,
    so the starting r(c) > 0 and no pair distance is 0."""
    anchor = np.asarray({'tri': (1 / 3, 1 / 3), 'rect': (0.0, 0.0), 'quad': (0.4, 0.4), 'semi': (0.0, 0.4)}[cont[0]])
    c = np.array(c, dtype=np.float64); n = len(c); iu = np.triu_indices(n, 1)
    d = np.sqrt(((c[:, None] - c[None]) ** 2).sum(-1))[iu]
    bad = np.unique(iu[1][d < 1e-12])
    if len(bad):                  # move each duplicate a DIFFERENT distance toward the interior anchor (convex: stays inside)
        u = np.linspace(1e-6, 1e-5, len(bad)); c[bad] += (anchor - c[bad]) * u[:, None]
    if rmin(c, cont) <= 0:        # a centre sitting exactly on a wall: pull everything a hair inward
        c = anchor + (c - anchor) * (1 - 1e-9)
    return c

def polish(c, cont, max_iter=300, delta0=None, t_cap=90.0):
    c = unstick(repair(np.asarray(c, dtype=np.float64), cont), cont); n = len(c); iu = np.triu_indices(n, 1)
    cur = rmin(c, cont); delta = delta0 or 0.1 * cur
    t_start = time.time()
    for it in range(max_iter):
        if time.time() - t_start > t_cap: break                    # wall-clock cap (09-23: 3 workers hung ~50 min)
        g, G = walls(c, cont)
        D = np.sqrt(((c[:, None] - c[None]) ** 2).sum(-1)); d = D[iu]
        sel = d / 2 <= cur + 3 * delta + 1e-12; I, J = iu[0][sel], iu[1][sel]; m = len(I)
        U = (c[I] - c[J]) / np.maximum(d[sel], 1e-300)[:, None]      # coincident centres can never divide by zero
        # pair rows: t - (u.(di - dj))/2 <= d/2
        r1 = np.repeat(np.arange(m), 5)
        c1 = np.stack([2 * I, 2 * I + 1, 2 * J, 2 * J + 1, np.full(m, 2 * n)], 1).ravel()
        v1 = np.stack([-U[:, 0] / 2, -U[:, 1] / 2, U[:, 0] / 2, U[:, 1] / 2, np.ones(m)], 1).ravel()
        A1 = coo_matrix((v1, (r1, c1)), shape=(m, 2 * n + 1)); b1 = d[sel] / 2
        # wall rows: t - grad.di <= g
        wi, wk = np.nonzero(g <= cur + 3 * delta + 1e-12); mw = len(wi)
        r2 = np.repeat(np.arange(mw), 3)
        c2 = np.stack([2 * wi, 2 * wi + 1, np.full(mw, 2 * n)], 1).ravel()
        v2 = np.stack([-G[wi, wk, 0], -G[wi, wk, 1], np.ones(mw)], 1).ravel()
        A2 = coo_matrix((v2, (r2, c2)), shape=(mw, 2 * n + 1)); b2 = g[wi, wk]
        res = linprog(np.r_[np.zeros(2 * n), -1.0], A_ub=vstack([A1, A2]).tocsr(), b_ub=np.r_[b1, b2],
                      bounds=[(-delta, delta)] * (2 * n) + [(None, None)], method='highs', options={'time_limit': 20.0})
        if res.status != 0 or res.x is None or not np.all(np.isfinite(res.x)):
            delta /= 2
            if delta < 1e-16: break
            continue
        q = repair(c + res.x[:-1].reshape(n, 2), cont); new = rmin(q, cont)
        if new > cur + 1e-16:
            gain = new - cur; c, cur = q, new
            if gain < 1e-15 and delta < 1e-12: break
            delta = min(delta * 1.5, 0.2 * cur)
        else:
            delta /= 4
            if delta < 1e-16: break
    return c, float(cur)
