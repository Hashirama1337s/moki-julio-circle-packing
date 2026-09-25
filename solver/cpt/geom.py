"""cpt probe: geometry + data loaders + SLP polish for equal circles in a REGULAR PENTAGON (Packomania cpt).

Container (derived from the cpt page and verified by check_geom.py):
  regular pentagon, circumradius R = 1, centred at the origin, one vertex at the top (0, 1), bottom side horizontal.
  vertices  v_k = (cos(90 + 72k deg), sin(90 + 72k deg)),  k = 0..4
  sides     outward unit normals n_k at angles 270 + 72k deg (270 = the bottom side), all at distance a = cos 36 deg
  a circle of radius t centred at c is inside  <=>  a - n_k . c >= t  for all k.
  N = 1 radius = a = cos 36 = 0.809016994375 (the inradius).
Float only. New file; nothing in the parent folder is edited.
"""
import os, re, time, numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # the folder holding data/refs/ and data/big/cpt/ (Packomania's cpt page and files)
PAGE = os.path.join(ROOT, 'data', 'refs', 'packomania_cpt_2026-09-25.html')
CDIR = os.path.join(ROOT, 'data', 'big', 'cpt')

A = np.cos(np.pi / 5)                                   # inradius (apothem) for circumradius 1
ANG = np.deg2rad(270.0 + 72.0 * np.arange(5))
NRM = np.stack([np.cos(ANG), np.sin(ANG)], 1)            # (5, 2) outward unit normals
VERT = np.stack([np.cos(np.deg2rad(90 + 72.0 * np.arange(5))), np.sin(np.deg2rad(90 + 72.0 * np.arange(5)))], 1)


def page_table():
    """{N: (radius_str, refs_tuple)} from the saved cpt page. ref '1' = Specht, program cpt, 2023; '2' = Amore & Morales 2023."""
    s = open(PAGE, encoding='utf-8', errors='replace').read()
    rows = re.findall(r'name="cpt(\d+)">.*?</a></td>\s*<td>([0-9.]+)</td>(.*?)</tr>', s, re.S)
    return {int(n): (r, tuple(re.findall(r'\[(\d+)\]', rest))) for n, r, rest in rows}


def load_coords(n):
    rows = []
    for line in open(os.path.join(CDIR, f'cpt{n}.txt')):
        p = line.split()
        if not p: continue
        rows.append((float(p[1]), float(p[2])) if len(p) >= 3 else (0.0, 0.0))   # cpt1.txt has an empty coordinate line
    c = np.array(rows, dtype=np.float64)
    assert len(c) == n, (n, len(c))
    return c


def walls(c, nrm=NRM, a=A):
    g = a - c @ nrm.T                                     # (n, 5) slacks
    G = np.broadcast_to(-nrm[None], (len(c), len(nrm), 2))  # d g / d c
    return g, G


def rmin(c, nrm=NRM, a=A):
    g, _ = walls(c, nrm, a)
    if len(c) == 1: return float(g.min())
    d, _ = cKDTree(c).query(c, k=2)
    return float(min(g.min(), d[:, 1].min() / 2))


def pair_min(c):
    d, _ = cKDTree(c).query(c, k=2); return float(d[:, 1].min() / 2)


def inside(p):
    return np.all(A - p @ NRM.T >= 0, axis=1)


def repair(c):
    """Project centres back into the pentagon (a few sequential half-plane projections)."""
    c = c.copy()
    for _ in range(3):
        v = c @ NRM.T - A
        if v.max() <= 0: break
        for k in range(5):
            m = v[:, k] > 0
            if m.any(): c[m] -= v[m, k, None] * NRM[k]
            v = c @ NRM.T - A
    return c


def polish(c, max_iter=400, t_cap=120.0, tol_rel=1e-15):
    """Trust-region SLP (HiGHS dual simplex): maximise t s.t. linearised |c_i - c_j| >= 2t and exact (linear) walls >= t.
    A step is kept only if the TRUE common radius grows. Returns (c, r, iterations)."""
    c = np.asarray(c, dtype=np.float64); n = len(c); cur = rmin(c); delta = 0.1 * max(cur, 1e-6); t0 = time.time()
    scale = max(cur, 1e-3); it = 0
    for it in range(max_iter):
        if time.time() - t0 > t_cap: break
        g, G = walls(c)
        P = cKDTree(c).query_pairs(2 * (cur + 3 * delta) + 1e-12, output_type='ndarray'); I, J = P[:, 0], P[:, 1]; m = len(I)
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
                      bounds=[(-delta, delta)] * (2 * n) + [(None, None)], method='highs', options={'time_limit': 60.0})
        if res.status != 0 or res.x is None or not np.all(np.isfinite(res.x)):
            delta /= 2
            if delta < 1e-16 * scale: break
            continue
        q = repair(c + res.x[:-1].reshape(n, 2)); new = rmin(q)
        if new > cur + 1e-17:
            gain = new - cur; c, cur = q, new; scale = max(cur, 1e-3)
            if gain < tol_rel * cur and delta < 1e-12 * cur: break
            delta = min(delta * 1.5, 0.2 * cur)
        else:
            delta /= 4
            if delta < 1e-15 * scale: break
    return c, float(cur), it + 1


def contacts(c, r, rel=1e-7):
    """Per-circle contact count (pairs within 2r(1+rel), walls within r(1+rel))."""
    n = len(c); cnt = np.zeros(n, int)
    P = cKDTree(c).query_pairs(2 * r * (1 + rel), output_type='ndarray')
    if len(P): np.add.at(cnt, P[:, 0], 1); np.add.at(cnt, P[:, 1], 1)
    g, _ = walls(c); cnt += (g <= r * (1 + rel)).sum(1)
    return cnt
