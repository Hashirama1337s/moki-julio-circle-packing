# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""k-gon tables (Packomania cxd = regular 16-gon, cpd = regular 15-gon): geometry + data loaders + SLP polish. Moki&Julio.
k-parametrised copy of ../cpt/geom.py (the pentagon).

FRAME (verified in float on every published cxd / cpd packing by check_geom.py):
  regular k-gon, circumradius 1, centred at the origin, one side horizontal at the bottom.
  side j = 0..k-1 has outward unit normal n_j = (cos t_j, sin t_j), t_j = 270 + 360 j / k degrees, at distance
  a = cos(180 / k deg) from the origin; vertices at angles 270 + 180/k + 360 j / k.
  a circle of radius t centred at c is inside  <=>  a - n_j . c - t >= 0  for all j.
  (k = 16: flat top as well; k = 15: vertex at (0, 1); k = 5 reproduces the cpt frame of ../cpt/geom.py.)

Page credits: on the cxd and cpd pages reference [1] is Amore (Phys. Fluids 35, 027130, 2023) and reference [2] is
E. Specht's program (cxd 2023, cpd 2020) -- the opposite numbering of the cpt page.  PROGRAM_REF holds that.
Float only.
"""
import os, re, time, numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # the folder holding data/refs/ and data/big/<table>/ (Packomania's cxd / cpd pages and files)
TABLES = {'cxd': 16, 'cpd': 15}                    # Packomania code -> number of sides
PROGRAM_REF = {'cxd': '2', 'cpd': '2', 'cpt': '1'}  # the page's reference number of Specht's program
PAGE_DATE = '2026-09-25'


class Frame:
    """Regular k-gon, circumradius 1, flat bottom side, centred at the origin."""

    def __init__(self, k):
        self.k = k
        self.A = float(np.cos(np.pi / k))                                  # apothem (inradius)
        self.ANG = np.deg2rad(270.0 + 360.0 * np.arange(k) / k)
        self.NRM = np.stack([np.cos(self.ANG), np.sin(self.ANG)], 1)     # (k, 2) outward unit normals
        va = np.deg2rad(270.0 + 180.0 / k + 360.0 * np.arange(k) / k)
        self.VERT = np.stack([np.cos(va), np.sin(va)], 1)

    # ------------------------------------------------------------------ basic geometry
    def walls(self, c):
        g = self.A - c @ self.NRM.T                                        # (n, k) slacks
        G = np.broadcast_to(-self.NRM[None], (len(c), self.k, 2))          # d g / d c
        return g, G

    def rmin(self, c):
        g, _ = self.walls(c)
        if len(c) == 1: return float(g.min())
        d, _ = cKDTree(c).query(c, k=2)
        return float(min(g.min(), d[:, 1].min() / 2))

    def dense_r(self, c):
        """Float min-radius from ALL pairs (O(N^2)) and all wall slacks: an independent recheck of rmin."""
        g, _ = self.walls(c); n = len(c)
        if n == 1: return float(g.min())
        iu = np.triu_indices(n, 1)
        d = np.sqrt(((c[:, None, :] - c[None, :, :]) ** 2).sum(-1))[iu]
        return float(min(g.min(), d.min() / 2))

    def inside(self, p):
        return np.all(self.A - p @ self.NRM.T >= 0, axis=1)

    def repair(self, c):
        """Project centres back into the polygon (a few sequential half-plane projections)."""
        c = c.copy()
        for _ in range(4):
            v = c @ self.NRM.T - self.A
            if v.max() <= 0: break
            for j in range(self.k):
                m = v[:, j] > 0
                if m.any(): c[m] -= v[m, j, None] * self.NRM[j]
                v = c @ self.NRM.T - self.A
        return c

    def contacts(self, c, r, rel=1e-7):
        """Per-circle contact count (pairs within 2r(1+rel), walls within r(1+rel))."""
        n = len(c); cnt = np.zeros(n, int)
        P = cKDTree(c).query_pairs(2 * r * (1 + rel), output_type='ndarray')
        if len(P): np.add.at(cnt, P[:, 0], 1); np.add.at(cnt, P[:, 1], 1)
        g, _ = self.walls(c); cnt += (g <= r * (1 + rel)).sum(1)
        return cnt

    # ------------------------------------------------------------------ SLP polish (same algorithm as ../cpt/geom.py)
    def polish(self, c, max_iter=400, t_cap=120.0, tol_rel=1e-15):
        """Trust-region SLP (HiGHS): maximise t s.t. linearised |c_i - c_j| >= 2t and exact (linear) walls >= t.
        A step is kept only if the TRUE common radius grows. Returns (c, r, iterations)."""
        k = self.k
        c = np.asarray(c, dtype=np.float64); n = len(c); cur = self.rmin(c); delta = 0.1 * max(cur, 1e-6); t0 = time.time()
        scale = max(cur, 1e-3); it = 0
        for it in range(max_iter):
            if time.time() - t0 > t_cap: break
            g, G = self.walls(c)
            P = cKDTree(c).query_pairs(2 * (cur + 3 * delta) + 1e-12, output_type='ndarray')
            if len(P) == 0: P = np.zeros((0, 2), int)
            I, J = P[:, 0], P[:, 1]; m = len(I)
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
            q = self.repair(c + res.x[:-1].reshape(n, 2)); new = self.rmin(q)
            if new > cur + 1e-17:
                gain = new - cur; c, cur = q, new; scale = max(cur, 1e-3)
                if gain < tol_rel * cur and delta < 1e-12 * cur: break
                delta = min(delta * 1.5, 0.2 * cur)
            else:
                delta /= 4
                if delta < 1e-15 * scale: break
        return c, float(cur), it + 1


FRAMES = {}


def frame(k):
    if k not in FRAMES: FRAMES[k] = Frame(k)
    return FRAMES[k]


# ---------------------------------------------------------------------- data loaders
def page_path(table):
    return os.path.join(ROOT, 'data', 'refs', f'packomania_{table}_{PAGE_DATE}.html')


def page_table(table):
    """{N: (radius_str, refs_tuple)} from the saved Packomania page of `table`."""
    s = open(page_path(table), encoding='utf-8', errors='replace').read()
    rows = re.findall(r'name="' + table + r'(\d+)">.*?</a></td>\s*<td>([0-9.]+)</td>(.*?)</tr>', s, re.S)
    return {int(n): (r, tuple(re.findall(r'\[(\d+)\]', rest))) for n, r, rest in rows}


def radius_file(table):
    """{N: radius_str} from data/big/<table>_radius.txt."""
    out = {}
    for line in open(os.path.join(ROOT, 'data', 'big', f'{table}_radius.txt')):
        p = line.split()
        if len(p) == 2: out[int(p[0])] = p[1]
    return out


def coord_strings(table, n):
    """Published centres as the printed decimal strings [(x, y)]; the N = 1 file prints no coordinates (circle at the centre)."""
    pts = []
    for line in open(os.path.join(ROOT, 'data', 'big', table, f'{table}{n}.txt')):
        p = line.split()
        if not p: continue
        if len(p) >= 3: pts.append((p[1], p[2]))
        elif len(p) == 1 and n == 1: pts.append(('0', '0'))
        else: raise ValueError(f'{table}{n}.txt: bad line {line!r}')
    assert len(pts) == n, (table, n, len(pts))
    return pts


def load_coords(table, n):
    return np.array([(float(x), float(y)) for x, y in coord_strings(table, n)], dtype=np.float64)


def program_only(table, lo=49, hi=None):
    """N in lo..hi whose page reference list is exactly [Specht's program]."""
    T = page_table(table); ref = PROGRAM_REF[table]
    hi = hi or max(T)
    return [n for n in range(lo, hi + 1) if n in T and T[n][1] == (ref,)]
