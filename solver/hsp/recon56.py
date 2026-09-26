# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Recover Packomania's hsp5 / hsp6 packings from their TRUNCATED coordinate files (2026-09-26; Moki&Julio). Float only.

Packomania serves hsp5-<N>.txt / hsp6-<N>.txt with only the first 4 of the 5 / 6 coordinates of every centre (measured over
all 299 + 249 files: |x_1..4|^2 <= (1 - r)^2 everywhere, to the 12-decimal rounding; hsp5 N = 2 is (0,0,0,0,+-1/2)). The missing
coordinates y_i (1 per ball in 5-D, 2 in 6-D) are pinned by the geometry: with the known part a_i and the printed radius r,
  wall : |y_i| <= m_i,        m_i^2 = (1 - r)^2 - |a_i|^2      (= for a ball touching the wall, i.e. almost every ball)
  pairs: |y_i - y_j| >= g_ij, g_ij^2 = 4 r^2 - |a_i - a_j|^2  (only pairs with g_ij^2 > 0 constrain anything)
(r lowered by 1e-9 to absorb the 12-decimal rounding). 5-D, all balls on the wall: y_i = s_i m_i and a pair with |m_i - m_j| <
g_ij forces s_i != s_j -> a 2-colouring of the conflict graph (BFS); an odd cycle marks an interior ball. 6-D: y_i = m_i (cos
phi_i, sin phi_i) from random angles. Then a smooth penalty (L-BFGS) on every wall and pair term with the first 4 columns fixed
finds y; the best of several starts is kept. The completed packing is then released (all coordinates free) to slpd.polish.
Success = the polished radius reaches the printed one (it is then Specht's packing or an equally good one); it is the seed for
the basin hopping that beat ~90 % of the hsp4 cells from Packomania's own coordinates.

usage: python recon56.py --test 5:30,6:100 [--polish 60] [--starts 12] [--workers 3]
"""
import os, sys, time, json, argparse
import numpy as np
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prio                                                      # noqa: E402
import geomd, slpd                                               # noqa: E402
OUT = os.path.join(geomd.OUT, 'recon56')


def trunc(d, n):
    """(printed radius float, known 4 columns (n, 4)) from Packomania's truncated file."""
    L = [l.split() for l in open(geomd.pub_path(d, n)) if l.strip()]
    A = np.array([[float(v) for v in t[1:]] for t in L[1:]], dtype=np.float64)
    assert A.shape == (n, 4), (d, n, A.shape)
    return float(L[0][0]), A


def needs(A, r, tol=1e-9, sq=1e-11):
    """m_i, and the constrained pairs (i, j, g_ij), at radius r - tol. sq = allowance on every SQUARED quantity for the files'
    12-decimal rounding (|a|^2 and |a_i - a_j|^2 are off by up to ~4e-12): without it a ball near the 'equator' (m ~ 1e-4) has
    its wall position wrong by ~1e-8 and cells like hsp5 N = 250 / 300 read as infeasible at the printed radius."""
    re_ = r - tol
    m = np.sqrt(np.maximum((1 - re_) ** 2 - (A * A).sum(1) + sq, 0.0))
    iu = np.triu_indices(len(A), 1)
    g2 = 4 * re_ * re_ - ((A[iu[0]] - A[iu[1]]) ** 2).sum(1) - sq; k = g2 > 0
    return m, iu[0][k], iu[1][k], np.sqrt(g2[k])


def colour5(m, I, J, G, rng):
    """5-D signs from the conflict graph (pairs where the same sign is impossible on the wall): BFS 2-colouring from a random
    root order; returns signs and the number of odd-cycle (unsatisfied) edges."""
    n = len(m); adj = [[] for _ in range(n)]
    for i, j, g in zip(I, J, G):
        if abs(m[i] - m[j]) < g: adj[i].append(j); adj[j].append(i)
    s = np.zeros(n, dtype=int)
    for root in rng.permutation(n):
        if s[root]: continue
        s[root] = 1 if rng.random() < 0.5 else -1; q = [root]
        while q:
            u = q.pop()
            for v in adj[u]:
                if not s[v]: s[v] = -s[u]; q.append(v)
    bad = sum(1 for u in range(n) for v in adj[u] if u < v and s[u] == s[v])
    return s, bad


def penalty(y, m, I, J, G, k):
    Y = y.reshape(-1, k); D = Y[I] - Y[J]; dist = np.sqrt((D * D).sum(1) + 1e-30)
    pv = np.maximum(G - dist, 0.0); nr = np.sqrt((Y * Y).sum(1) + 1e-30); wv = np.maximum(nr - m, 0.0)
    f = (pv * pv).sum() + (wv * wv).sum()
    gY = np.zeros_like(Y)
    coef = (-2 * pv / dist)[:, None] * D
    np.add.at(gY, I, coef); np.add.at(gY, J, -coef)
    gY += (2 * wv / nr)[:, None] * Y
    return f, gY.ravel()


def viol_of(Y, m, I, J, G):
    D = Y[I] - Y[J]
    return max(float(np.max(G - np.sqrt((D * D).sum(1)), initial=0.0)), float(np.max(np.sqrt((Y * Y).sum(1)) - m, initial=0.0)))


def complete5(A, r, time_limit=60.0):
    """5-D, exact combinatorics: the FEWEST interior balls z_i such that every other ball sits on the wall (x5 = +-m_i) with the
    conflict graph (same sign impossible: |m_i - m_j| < g_ij) 2-coloured. MILP (HiGHS): min sum z, per conflict edge
    s_i + s_j + z_i + z_j >= 1 and s_i + s_j - z_i - z_j <= 1. Interior balls are then placed one by one at the middle of the
    widest free gap of [-m_i, m_i] left by their neighbours, and everything is refined by L-BFGS-B (bounds |x5| <= m_i).
    Returns (C, violation, info)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix
    m, I, J, G = needs(A, r); n = len(m)
    e = np.abs(m[I] - m[J]) < G; Ie, Je = I[e], J[e]; k = len(Ie)
    rows = np.repeat(np.arange(2 * k).reshape(k, 2), 4, axis=1).ravel()       # edge t: row 2t (>= 1), row 2t + 1 (<= 1)
    cols = np.stack([Ie, Je, n + Ie, n + Je, Ie, Je, n + Ie, n + Je], 1).ravel()
    vals = np.tile([1, 1, 1, 1, 1, 1, -1, -1], k).astype(float)
    M = coo_matrix((vals, (rows, cols)), shape=(2 * k, 2 * n))
    lo = np.tile([1.0, -np.inf], k); hi = np.tile([np.inf, 1.0], k)
    res = milp(c=np.r_[np.zeros(n), np.ones(n)], constraints=LinearConstraint(M.tocsr(), lo, hi), integrality=np.ones(2 * n),
               bounds=Bounds(np.zeros(2 * n), np.ones(2 * n)), options=dict(time_limit=time_limit, disp=False))
    if res.x is None: return None, np.inf, dict(milp=res.status)
    s = np.round(res.x[:n]).astype(int); z = np.round(res.x[n:]).astype(int) == 1
    x = (2 * s - 1) * m; x[z] = 0.0; placed = ~z
    nb = [[] for _ in range(n)]
    for i, j, g in zip(I, J, G): nb[i].append((j, g)); nb[j].append((i, g))
    for i in np.flatnonzero(z)[np.argsort(-m[z])]:
        cuts = sorted((x[j] - g, x[j] + g) for j, g in nb[i] if placed[j])
        free, a = [], -m[i]
        for lo_, hi_ in cuts:
            if lo_ > a: free.append((a, min(lo_, m[i])))
            a = max(a, hi_)
        if a < m[i]: free.append((a, m[i]))
        free = [(u, v) for u, v in free if v > u]
        x[i] = (lambda t: (t[0] + t[1]) / 2)(max(free, key=lambda t: t[1] - t[0])) if free else 0.0
        placed[i] = True
    y0 = x.copy()
    f = lambda y: penalty(y, m, I, J, G, 1)
    res2 = minimize(f, y0, jac=True, method='L-BFGS-B', bounds=list(zip(-m, m)), options=dict(maxiter=5000, ftol=1e-30, gtol=1e-16))
    Y = res2.x[:, None]
    return np.hstack([A, Y]), viol_of(Y, m, I, J, G), dict(milp=int(res.status), interior=int(z.sum()), edges=k)


def complete5x(A, r, time_limit=120.0):
    """5-D, the full geometry as one MILP (HiGHS): x_i in [-m_i, m_i] continuous; s_i (sign), z_i (interior) binary with
    |x_i - (2 s_i - 1) m_i| <= 2 m_i z_i (a wall ball sits at +-m_i); for EVERY constrained pair an order binary o_p:
    x_i - x_j >= g_p - B_p (1 - o_p), x_j - x_i >= g_p - B_p o_p (B_p = m_i + m_j + g_p); min sum z. Its x is the completion
    (to the MIP tolerance), refined by L-BFGS-B. Returns (C, violation, info)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix
    m, I, J, G = needs(A, r); n = len(m); P = len(I)
    nv = 3 * n + P; X, S, Z, O = 0, n, 2 * n, 3 * n                  # variable blocks
    R, C, V, lo, hi = [], [], [], [], []
    def row(cols, vals, l, h):
        k = len(lo); R.extend([k] * len(cols)); C.extend(cols); V.extend(vals); lo.append(l); hi.append(h)
    for i in range(n):
        row([X + i, S + i, Z + i], [1.0, -2 * m[i], -2 * m[i]], -np.inf, -m[i])
        row([X + i, S + i, Z + i], [-1.0, 2 * m[i], -2 * m[i]], -np.inf, m[i])
    for p, (i, j, g) in enumerate(zip(I, J, G)):
        B = m[i] + m[j] + g
        row([X + i, X + j, O + p], [1.0, -1.0, -B], g - B, np.inf)
        row([X + i, X + j, O + p], [-1.0, 1.0, B], g, np.inf)
    M = coo_matrix((V, (R, C)), shape=(len(lo), nv)).tocsr()
    lb = np.r_[-m, np.zeros(2 * n + P)]; ub = np.r_[m, np.ones(2 * n + P)]; ub[S] = lb[S] = 1          # one sign fixed: x5 -> -x5
    integ = np.r_[np.zeros(n), np.ones(2 * n + P)]
    res = milp(c=np.r_[np.zeros(n), np.zeros(n), np.ones(n), np.zeros(P)], constraints=LinearConstraint(M, lo, hi),
               integrality=integ, bounds=Bounds(lb, ub), options=dict(time_limit=time_limit, disp=False, mip_rel_gap=0.0))
    if res.x is None: return None, np.inf, dict(milp=int(res.status), pairs=P)
    x = np.clip(res.x[:n], -m, m); z = int(np.round(res.x[Z:Z + n]).sum())
    res2 = minimize(lambda y: penalty(y, m, I, J, G, 1), x, jac=True, method='L-BFGS-B', bounds=list(zip(-m, m)),
                    options=dict(maxiter=5000, ftol=1e-30, gtol=1e-16))
    Y = res2.x[:, None]
    return np.hstack([A, Y]), viol_of(Y, m, I, J, G), dict(milp=int(res.status), interior=z, pairs=P, obj=float(res.fun))


def complete6x(A, r, rng, time_limit=120.0):
    """6-D, one MILP (HiGHS) on the ANGLES of the missing plane: a wall ball is y_i = m_i (cos phi_i, sin phi_i); a pair needs
    cos(phi_i - phi_j) <= kappa_ij = (m_i^2 + m_j^2 - g_ij^2) / (2 m_i m_j), i.e. with theta = arccos(kappa) the difference
    Delta = phi_i - phi_j (phi in [0, 2 pi]) lies in [theta, 2 pi - theta] (o_p = 1) or [-2 pi + theta, -theta] (o_p = 0);
    kappa >= 1: no constraint; kappa < -1: impossible on the wall -> z_i + z_j >= 1. An interior ball (z_i = 1) drops all its
    angle rows (big-M) and is placed afterwards at the best of 4,000 random points of its disc |y| <= m_i; then L-BFGS on the
    Cartesian penalty. min sum z; phi fixed to 0 for the ball with the largest m (rotation of the plane). Returns (C, viol, info)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix
    m, I, J, G = needs(A, r); n = len(m); TP = 2 * np.pi; BM = 4 * np.pi
    kap = (m[I] ** 2 + m[J] ** 2 - G ** 2) / np.maximum(2 * m[I] * m[J], 1e-300)
    act = kap < 1.0; imp = kap < -1.0
    Ia, Ja, th = I[act & ~imp], J[act & ~imp], np.arccos(np.clip(kap[act & ~imp], -1, 1))
    Ii, Ji = I[imp], J[imp]; P = len(Ia)
    PH, Z, O = 0, n, 2 * n; nv = 2 * n + P
    R, Cc, V, lo, hi = [], [], [], [], []
    def row(cols, vals, l, h):
        k = len(lo); R.extend([k] * len(cols)); Cc.extend(cols); V.extend(vals); lo.append(l); hi.append(h)
    for p, (i, j, t) in enumerate(zip(Ia, Ja, th)):
        # o = 1: t <= D <= 2pi - t ; o = 0: -2pi + t <= D <= -t ; each row relaxed by BM (z_i + z_j) and by BM for the other branch
        row([PH + i, PH + j, O + p, Z + i, Z + j], [1, -1, -BM, BM, BM], t - BM, np.inf)
        row([PH + i, PH + j, O + p, Z + i, Z + j], [1, -1, BM, -BM, -BM], -np.inf, TP - t + BM)
        row([PH + i, PH + j, O + p, Z + i, Z + j], [1, -1, BM, BM, BM], -TP + t, np.inf)
        row([PH + i, PH + j, O + p, Z + i, Z + j], [1, -1, -BM, -BM, -BM], -np.inf, -t)
    for i, j in zip(Ii, Ji): row([Z + i, Z + j], [1, 1], 1, np.inf)
    M = coo_matrix((V, (R, Cc)), shape=(len(lo), nv)).tocsr()
    lb = np.zeros(nv); ub = np.r_[np.full(n, TP), np.ones(n + P)]; k0 = int(np.argmax(m)); ub[PH + k0] = 0.0
    res = milp(c=np.r_[np.zeros(n), np.ones(n), np.zeros(P)], constraints=LinearConstraint(M, lo, hi),
               integrality=np.r_[np.zeros(n), np.ones(n + P)], bounds=Bounds(lb, ub), options=dict(time_limit=time_limit, disp=False))
    if res.x is None: return None, np.inf, dict(milp=int(res.status), pairs=P)
    ph = res.x[:n]; z = np.round(res.x[Z:Z + n]).astype(bool)
    Y = np.stack([m * np.cos(ph), m * np.sin(ph)], 1); placed = ~z
    nb = [[] for _ in range(n)]
    for i, j, g in zip(I, J, G): nb[i].append((j, g)); nb[j].append((i, g))
    for i in np.flatnonzero(z)[np.argsort(-m[z])]:
        u = rng.uniform(-1, 1, (4000, 2)); u = u[(u * u).sum(1) <= 1] * m[i]
        sl = np.full(len(u), np.inf)
        for j, g in nb[i]:
            if placed[j]: sl = np.minimum(sl, np.sqrt(((u - Y[j]) ** 2).sum(1)) - g)
        Y[i] = u[np.argmax(sl)] if len(u) else 0.0; placed[i] = True
    res2 = minimize(penalty, Y.ravel(), args=(m, I, J, G, 2), jac=True, method='L-BFGS-B', options=dict(maxiter=8000, ftol=1e-30, gtol=1e-16))
    Y = res2.x.reshape(-1, 2)
    return np.hstack([A, Y]), viol_of(Y, m, I, J, G), dict(milp=int(res.status), interior=int(z.sum()), pairs=P, impossible=len(Ii))


def minconf6(A, r, rng, sweeps=400, restarts=6, t_cap=120.0):
    """6-D min-conflicts on the angles of the wall balls: every ball in turn moves to the angle that maximises its worst angular
    slack min_j (|phi - phi_j|_circle - theta_ij) (evaluated exactly at every arc end phi_j +- theta_ij and the midpoints between
    them), random order, random restarts; balls whose best slack stays negative are declared interior and placed in their disc
    afterwards (as in complete6x). Returns (C, viol, info)."""
    m, I, J, G = needs(A, r); n = len(m); TP = 2 * np.pi
    kap = (m[I] ** 2 + m[J] ** 2 - G ** 2) / np.maximum(2 * m[I] * m[J], 1e-300)
    act = kap < 1.0; th = np.arccos(np.clip(kap, -1, 1)); th[kap < -1] = np.pi + 1e-3          # impossible pairs: always violated
    nb = [[] for _ in range(n)]
    for i, j, t, a in zip(I, J, th, act):
        if a: nb[i].append((j, t)); nb[j].append((i, t))
    nbj = [np.array([j for j, _ in L], dtype=int) for L in nb]; nbt = [np.array([t for _, t in L]) for L in nb]
    def worst(i, phi, ph):
        if not len(nbj[i]): return np.full(np.shape(phi), np.inf)
        dd = np.abs((np.asarray(phi)[..., None] - ph[nbj[i]] + np.pi) % TP - np.pi)
        return (dd - nbt[i]).min(-1)
    best = (None, -np.inf); t0 = time.process_time()
    for rs in range(restarts):
        ph = rng.uniform(0, TP, n)
        for sw in range(sweeps):
            moved = 0
            for i in rng.permutation(n):
                if not len(nbj[i]): continue
                e = np.concatenate([ph[nbj[i]] - nbt[i], ph[nbj[i]] + nbt[i]]) % TP; e.sort()
                mid = (e + np.r_[e[1:], e[0] + TP]) / 2 % TP
                cand = np.r_[e, mid, ph[i]]; w = worst(i, cand, ph); k = int(np.argmax(w))
                if w[k] > worst(i, ph[i], ph) + 1e-15: ph[i] = cand[k]; moved += 1
            sl = np.array([worst(i, ph[i], ph) for i in range(n)])
            if sl.min() >= -1e-12 or moved == 0 or time.process_time() - t0 > t_cap: break
            bad = np.flatnonzero(sl < -1e-12)                                    # kick: re-randomise the violated balls
            if sw % 20 == 19: ph[bad] = rng.uniform(0, TP, len(bad))
        sl = np.array([worst(i, ph[i], ph) for i in range(n)])
        score = -np.sum(sl < -1e-9) + 1e-3 * min(sl.min(), 0)
        if score > best[1]: best = (ph.copy(), score)
        if sl.min() >= -1e-12 or time.process_time() - t0 > t_cap: break
    ph = best[0]; sl = np.array([worst(i, ph[i], ph) for i in range(n)])
    z = sl < -1e-9
    Y = np.stack([m * np.cos(ph), m * np.sin(ph)], 1); placed = ~z
    nbc = [[] for _ in range(n)]
    for i, j, g in zip(I, J, G): nbc[i].append((j, g)); nbc[j].append((i, g))
    for i in np.flatnonzero(z)[np.argsort(-m[z])]:
        u = rng.uniform(-1, 1, (4000, 2)); u = u[(u * u).sum(1) <= 1] * m[i]
        s_ = np.full(len(u), np.inf)
        for j, g in nbc[i]:
            if placed[j]: s_ = np.minimum(s_, np.sqrt(((u - Y[j]) ** 2).sum(1)) - g)
        Y[i] = u[np.argmax(s_)]; placed[i] = True
    res2 = minimize(penalty, Y.ravel(), args=(m, I, J, G, 2), jac=True, method='L-BFGS-B', options=dict(maxiter=8000, ftol=1e-30, gtol=1e-16))
    Y = res2.x.reshape(-1, 2)
    return np.hstack([A, Y]), viol_of(Y, m, I, J, G), dict(method='minconf', interior=int(z.sum()), pairs=int(act.sum()))


def dfs6(A, r, rng, max_skips=4, tol=1e-7, node_cap=300000, t_cap=120.0):
    """6-D depth-first search on the angles of the wall balls. The ball with the most placed constraining neighbours goes next;
    its free set is the circle minus the open arcs (phi_j - theta_ij + tol, phi_j + theta_ij - tol); candidates = the free-arc
    ends, ranked by how many arc ends coincide there (a ball touching two placed balls sits where both arcs end), then the
    midpoints of long free arcs. An empty free set may be skipped (the ball is interior; at most max_skips) or backtracked.
    Skipped balls are placed in their disc afterwards; then L-BFGS on the Cartesian penalty. Returns (C, viol, info)."""
    m, I, J, G = needs(A, r); n = len(m); TP = 2 * np.pi
    kap = (m[I] ** 2 + m[J] ** 2 - G ** 2) / np.maximum(2 * m[I] * m[J], 1e-300)
    th = np.arccos(np.clip(kap, -1, 1)); act = kap < 1.0
    nb = [[] for _ in range(n)]
    for i, j, t, a in zip(I, J, th, act):
        if a: nb[i].append((j, t)); nb[j].append((i, t))
    ph = np.full(n, np.nan); state = np.zeros(n, dtype=int)          # 0 unplaced, 1 placed, 2 skipped (interior)
    k0 = int(np.argmax(m)); ph[k0] = 0.0; state[k0] = 1
    def cands(i):
        P = [(ph[j], t) for j, t in nb[i] if state[j] == 1]
        if not P: return [rng.uniform(0, TP)]
        ends = []
        for c, t in P: ends += [(c - t) % TP, (c + t) % TP]
        ends = np.array(ends)
        def free(x):
            return all(abs((x - c + np.pi) % TP - np.pi) >= t - tol for c, t in P)
        ok = [x for x in ends if free(x)]
        if not ok: return []
        mult = [int((np.abs((ends - x + np.pi) % TP - np.pi) < 1e-6).sum()) for x in ok]
        order = [x for _, x in sorted(zip(mult, ok), key=lambda t: -t[0])]
        es = np.sort(ends); mids = [((a + b) / 2) % TP for a, b in zip(es, np.r_[es[1:], es[0] + TP]) if b - a > 1e-3]
        return order + [x for x in mids if free(x)][:2]
    stack = []; nodes = 0; t0 = time.process_time(); best = (-1, None, None)
    def pick():
        U = np.flatnonzero(state == 0)
        if not len(U): return None
        sc = [sum(1 for j, _ in nb[i] if state[j] == 1) + 1e-3 * m[i] for i in U]
        return int(U[int(np.argmax(sc))])
    i = pick()
    while i is not None and nodes < node_cap and time.process_time() - t0 < t_cap:
        nodes += 1
        C_ = cands(i)
        if C_:
            stack.append((i, C_, 1)); ph[i] = C_[0]; state[i] = 1
        elif (state == 2).sum() < max_skips:
            stack.append((i, [], 0)); state[i] = 2
        else:                                                         # backtrack to the last ball with an untried candidate
            while stack:
                j, L, k = stack.pop(); state[j] = 0; ph[j] = np.nan
                if k < len(L): stack.append((j, L, k + 1)); ph[j] = L[k]; state[j] = 1; break
            if not stack: break
        placed = int((state == 1).sum())
        if placed > best[0]: best = (placed, ph.copy(), state.copy())
        i = pick()
    done = i is None
    if not done: ph, state = best[1], best[2]
    z = state != 1
    Y = np.zeros((n, 2)); Y[~z] = np.stack([m[~z] * np.cos(ph[~z]), m[~z] * np.sin(ph[~z])], 1); placed = ~z
    nbc = [[] for _ in range(n)]
    for a, b, g in zip(I, J, G): nbc[a].append((b, g)); nbc[b].append((a, g))
    for i in np.flatnonzero(z)[np.argsort(-m[z])]:
        u = rng.uniform(-1, 1, (4000, 2)); u = u[(u * u).sum(1) <= 1] * m[i]
        s_ = np.full(len(u), np.inf)
        for j, g in nbc[i]:
            if placed[j]: s_ = np.minimum(s_, np.sqrt(((u - Y[j]) ** 2).sum(1)) - g)
        Y[i] = u[np.argmax(s_)]; placed[i] = True
    res2 = minimize(penalty, Y.ravel(), args=(m, I, J, G, 2), jac=True, method='L-BFGS-B', options=dict(maxiter=8000, ftol=1e-30, gtol=1e-16))
    Y = res2.x.reshape(-1, 2)
    return np.hstack([A, Y]), viol_of(Y, m, I, J, G), dict(method='dfs', complete=done, nodes=nodes, interior=int(z.sum()))


def complete(d, A, r, rng, starts=12):
    """Best completion (n, d) over several starts; returns (C, worst violation of the constraints at r - 1e-9, start log)."""
    if d == 5:                                                    # exact MILP; if the file needs it, a slightly lower radius
        for f in (0.0, 1e-8, 1e-7, 1e-6, 1e-5):
            C, v, info = complete5x(A, r * (1 - f))
            if C is not None: return C, v, [dict(info, r_factor=f)]
    m, I, J, G = needs(A, r); k = d - 4; best = (None, np.inf); log = []
    for s in range(starts):
        if k == 1:
            sg, bad = colour5(m, I, J, G, rng); y0 = (sg * m)[:, None]
        else:
            ph = rng.uniform(0, 2 * np.pi, len(m)); y0 = np.stack([m * np.cos(ph), m * np.sin(ph)], 1); bad = None
        res = minimize(penalty, y0.ravel(), args=(m, I, J, G, k), jac=True, method='L-BFGS-B',
                       options=dict(maxiter=5000, ftol=1e-30, gtol=1e-16))
        Y = res.x.reshape(-1, k); D = Y[I] - Y[J]
        viol = max(float(np.max(G - np.sqrt((D * D).sum(1)), initial=0.0)), float(np.max(np.sqrt((Y * Y).sum(1)) - m, initial=0.0)))
        log.append(dict(start=s, odd_edges=bad, viol=viol))
        if viol < best[1]: best = (np.hstack([A, Y]), viol)
        if viol < 1e-12: break
    return best[0], best[1], log


def run(args):
    d, n, polish_cap, starts, seed = args
    prio.lower()
    rng = np.random.default_rng(seed); t0 = time.process_time()
    r_pr, A = trunc(d, n)
    C, viol, log = complete(d, A, r_pr, rng, starts)
    r_c = geomd.rmin(C); t1 = time.process_time()
    q, r_p = slpd.polish(C, t_cap=polish_cap, seed=seed)
    os.makedirs(OUT, exist_ok=True); np.save(os.path.join(OUT, f'hsp{d}_{n}.npy'), q)
    return dict(d=d, N=n, printed=r_pr, viol=viol, starts=len(log), r_completed=r_c, gap_completed=r_c - r_pr,
                r_polished=r_p, gap_polished=r_p - r_pr, t_complete=round(t1 - t0, 1), t_polish=round(time.process_time() - t1, 1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--test', required=True); ap.add_argument('--polish', type=float, default=60.0)
    ap.add_argument('--starts', type=int, default=12); ap.add_argument('--workers', type=int, default=3)
    a = ap.parse_args(); prio.lower()
    cells = [tuple(map(int, x.split(':'))) for x in a.test.split(',')]
    from multiprocessing import Pool
    with Pool(a.workers) as pool:
        for x in pool.imap_unordered(run, [(d, n, a.polish, a.starts, 20260926 + 1000 * d + n) for d, n in cells]):
            print(f"hsp{x['d']} N={x['N']:3d}: completion viol {x['viol']:.1e} ({x['starts']} starts), r_completed - printed {x['gap_completed']:+.3e}, "
                  f"polished - printed {x['gap_polished']:+.3e}  [{x['t_complete']} + {x['t_polish']} s cpu]", flush=True)
            with open(os.path.join(OUT, 'test.jsonl'), 'a') as f: f.write(json.dumps(x) + '\n')


if __name__ == '__main__':
    main()
