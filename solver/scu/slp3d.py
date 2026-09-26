# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""slp_big.polish in 3-D: equal spheres in the unit cube centred at the origin (Packomania scu: side 1, centres in
[-1/2 + r, 1/2 - r]^3). Same trust-region SLP: maximise t subject to linearised pair distances (U_ij . (c_i - c_j) / 2 >= t,
an inner approximation of |c_i - c_j| / 2 >= t) and the six walls; a step is kept only if the TRUE float min-radius grows.
Candidate pairs from a KD-tree. Float only; certification is exact and separate.

2026-09-25 instrument repair (polish; the first version is kept as polish_v1 so the probe numbers stay reproducible).
Measured: at N = 483 the first LP of polish_v1 (1,450 variables, ~2,600 rows) hit the 120 s time limit. The published
packing has 1,131 exact contacts, so the start is a vertex with > 1,000 tight rows: the LP is massively DEGENERATE (dual
simplex 500k iterations in 60 s, IPM + crossover also stalls). With a random right-hand-side perturbation it solves in
seconds, but HiGHS's default absolute primal feasibility tolerance (1e-7) is larger than the gains we chase, so the TRUE
radius of the step drops and the guard rejects every step. The repaired LP:
  (1) is SCALED: y = dx / delta in [-1, 1]^(3N), s = (t - cur) / delta, so every tolerance is relative to the step;
  (2) breaks degeneracy by TIGHTENING every row by an independent random amount in [pert/2, pert] (scaled units,
      default pert = 1e-9 > the feasibility tolerance), so the start vertex is no longer degenerate and a returned
      solution is truly feasible for the un-perturbed linearisation (the step's true radius is >= cur + delta s, up to
      float rounding in the evaluation);
  (3) passes tight HiGHS tolerances (primal / dual feasibility 1e-10) and a per-LP time limit (default 60 s);
      the LP is solved by HiGHS IPM (+ crossover) first, dual simplex (half the time limit) only if IPM fails.
      Measured under load: N = 866, delta = 0.1 r, 4,796 rows: IPM 3.3 s (393 it.) vs dual simplex 44-51 s (60-69k it.),
      same optimum; N = 738 published (degenerate) start: IPM 24 s, dual simplex > 30 s (time limit).
  (4) keeps the guard: a step is kept only if the true float min-radius (KD-tree nearest pair and six walls) grows.
An LP that fails (time limit, numerical trouble) quarters the trust region and is retried with a fresh perturbation.
STOP RULE: if the LP optimum is s <= 2 pert (no first-order improving direction) AND its step does not raise the true
radius, the point is first-order STATIONARY (a KKT point of max-min: the LP dual is a stress on the tight rows that bounds
t <= cur for every trust region), so the polish stops there instead of shrinking delta to 1e-16 (stats['stop'] =
'stationary'; stats['stress_rows'] = rows with a dual weight > 1e-9). Before stopping, the same LP is solved once more
WITHOUT the perturbation: its s = 0 step can still gain through the curvature of the pair distances (a second-order
gain the linearisation cannot see); if it gains, the walk goes on. Five accepted steps in a row with s <= 2 pert also
stop it ('near_stationary'). Stationary is NOT optimal: a second-order or combinatorial move may still gain.
"""
import time, numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack
from scipy.spatial import cKDTree


def walls(c):
    n = len(c); g = np.concatenate([0.5 + c, 0.5 - c], 1)                  # (n, 6): x-, y-, z- faces then +faces
    G = np.zeros((n, 6, 3))
    for a in range(3): G[:, a, a] = 1.0; G[:, 3 + a, a] = -1.0
    return g, G


def rmin(c):
    g, _ = walls(c); d, _ = cKDTree(c).query(c, k=2); return float(min(g.min(), d[:, 1].min() / 2))


def repair(c):
    return np.clip(c, -0.5, 0.5)


def _rows(c, cur, delta):
    """Linearised pair rows and wall rows around c (unscaled): A (dx, t) <= b."""
    n = len(c); g, G = walls(c)
    P = cKDTree(c).query_pairs(2 * (cur + 3 * delta) + 1e-12, output_type="ndarray"); I, J = P[:, 0], P[:, 1]; m = len(I)
    dv = c[I] - c[J]; d = np.sqrt((dv ** 2).sum(1)); U = dv / np.maximum(d, 1e-300)[:, None]
    wi, wk = np.nonzero(g <= cur + 3 * delta + 1e-12); mw = len(wi)
    return (I, J, U, d, m), (wi, wk, G, g, mw)


def polish(c, max_iter=400, t_cap=300.0, log=None, stop_at=None, lp_time=60.0, pert=1e-9, tol=1e-10, seed=0,
           stats=None, pretest=False, start_after_pretest=1e-3, method='highs-ipm', fallback='highs'):
    """Repaired SLP polish (see the module docstring). Returns (coordinates, true float min-radius).
    pretest=True first solves ONE LP at a tiny trust region (1e-6 cur: only the near-tight rows, fast even when the
    start is massively degenerate); if that LP is flat (s <= 2 pert) and its step does not gain, the start is
    first-order stationary for EVERY trust region (dual stress argument) and the polish returns at once
    (stats['stop'] = 'stationary_pretest'). Otherwise the trust region restarts at start_after_pretest * cur and grows
    x4 per kept step while it is below 0.01 cur and the LP was fast (< lp_time / 4), else x1.5 (degenerate published
    starts make the large-delta LPs time out; a small first trust region lets the contacts open up first)."""
    rng = np.random.default_rng(seed)
    c = np.asarray(c, dtype=np.float64).copy(); n = len(c); cur = rmin(c); delta = 0.1 * cur; t0 = time.time()
    if pretest:
        delta = 1e-6 * cur
    st = stats if stats is not None else {}
    for k in ('lp', 'lp_fail', 'lp_secs', 'kept', 'rejected'): st.setdefault(k, 0)
    st.pop('stop', None); nflat = 0; retry_flat = False
    for it in range(max_iter):
        if time.time() - t0 > t_cap: st['stop'] = 'time'; break
        (I, J, U, d, m), (wi, wk, G, g, mw) = _rows(c, cur, delta)
        # pair row  i<j :  -U/2 . y_i + U/2 . y_j + s <= (d/2 - cur) / delta - eps
        r1 = np.repeat(np.arange(m), 7)
        c1 = np.stack([3 * I, 3 * I + 1, 3 * I + 2, 3 * J, 3 * J + 1, 3 * J + 2, np.full(m, 3 * n)], 1).ravel()
        v1 = np.stack([-U[:, 0] / 2, -U[:, 1] / 2, -U[:, 2] / 2, U[:, 0] / 2, U[:, 1] / 2, U[:, 2] / 2, np.ones(m)], 1).ravel()
        A1 = coo_matrix((v1, (r1, c1)), shape=(m, 3 * n + 1)); b1 = (d / 2 - cur) / delta
        # wall row :  -G . y_i + s <= (g - cur) / delta - eps
        r2 = np.repeat(np.arange(mw), 4)
        c2 = np.stack([3 * wi, 3 * wi + 1, 3 * wi + 2, np.full(mw, 3 * n)], 1).ravel()
        v2 = np.stack([-G[wi, wk, 0], -G[wi, wk, 1], -G[wi, wk, 2], np.ones(mw)], 1).ravel()
        A2 = coo_matrix((v2, (r2, c2)), shape=(mw, 3 * n + 1)); b2 = (g[wi, wk] - cur) / delta
        b = np.r_[b1, b2]
        if not retry_flat: b = b - pert * (0.5 + 0.5 * rng.random(len(b)))
        tl = time.time(); A = vstack([A1, A2]).tocsr(); ok = False; res = None
        tries = ((method, min(lp_time, 10.0)),) if retry_flat else ((method, lp_time), (fallback, lp_time / 2))
        for meth, tlim in tries:
            if meth is None: continue
            try:
                res = linprog(np.r_[np.zeros(3 * n), -1.0], A_ub=A, b_ub=b, bounds=[(-1.0, 1.0)] * (3 * n) + [(None, None)],
                              method=meth, options={'time_limit': float(tlim), 'primal_feasibility_tolerance': tol,
                                                    'dual_feasibility_tolerance': tol})
                ok = res.status == 0 and res.x is not None and np.all(np.isfinite(res.x))
            except Exception:
                ok = False; res = None
            if ok: break
            st['lp_fallback'] = st.get('lp_fallback', 0) + (meth == method and fallback is not None)
        st['lp'] += 1; st['lp_secs'] += time.time() - tl
        if not ok:
            st['lp_fail'] += 1; delta /= 4
            if delta < 1e-16: st['stop'] = 'delta'; break
            continue
        s = float(res.x[-1]); st['last_s'] = s; st['last_delta'] = delta; flat = s <= 2 * pert
        if flat:
            try:                                           # the LP dual = a stress on the rows (weights sum to 1)
                lam = -np.asarray(res.ineqlin.marginals); st['stress_rows'] = int((lam > 1e-9).sum())
                st['stress_sum'] = float(lam.sum())
            except Exception:
                pass
        q = repair(c + delta * res.x[:-1].reshape(n, 3)); new = rmin(q)
        if new > cur + 1e-16:
            gain = new - cur; c, cur = q, new; st['kept'] += 1; retry_flat = False
            if log: log(it, cur, delta, time.time() - t0)
            nflat = nflat + 1 if flat else 0
            if nflat >= 5: st['stop'] = 'near_stationary'; break
            if gain < 1e-15 and delta < 1e-12: st['stop'] = 'tiny'; break
            if stop_at is not None and cur >= stop_at: st['stop'] = 'stop_at'; break
            grow = 4.0 if (delta < 0.01 * cur and time.time() - tl < lp_time / 4) else 1.5
            delta = min(delta * grow, 0.2 * cur)
            if pretest and it == 0: delta = start_after_pretest * cur   # the pretest step gained: ramp up from here
        else:
            st['rejected'] += 1
            if pretest and it == 0 and not flat:
                delta = start_after_pretest * cur; continue  # pretest inconclusive: ramp up from here
            if flat and not retry_flat and not (pretest and it == 0):
                # one retry of the same LP WITHOUT the perturbation: its s = 0 vertex step may still gain through the
                # curvature of the pair distances (second order; the linearisation is an inner approximation)
                retry_flat = True; continue
            if flat:
                # LP optimum s = 0 (up to the perturbation) and the step does not gain: no first-order improving
                # direction. The dual stress on the tight rows bounds t <= cur for EVERY trust region, so shrinking
                # delta cannot help: stop.
                st['stop'] = 'stationary_pretest' if (pretest and it == 0) else 'stationary'; break
            delta /= 4
            if delta < 1e-16: st['stop'] = 'delta'; break
    else:
        st['stop'] = 'max_iter'
    st.setdefault('stop', 'time')
    return c, float(cur)


def polish_v1(c, max_iter=400, t_cap=300.0, log=None, stop_at=None):
    """The 2026-09-25 probe version (unscaled LP, default HiGHS tolerances, 120 s per LP). Kept for the record only."""
    c = np.asarray(c, dtype=np.float64); n = len(c); cur = rmin(c); delta = 0.1 * cur; t0 = time.time()
    for it in range(max_iter):
        if time.time() - t0 > t_cap: break
        g, G = walls(c)
        P = cKDTree(c).query_pairs(2 * (cur + 3 * delta) + 1e-12, output_type="ndarray"); I, J = P[:, 0], P[:, 1]; m = len(I)
        dv = c[I] - c[J]; d = np.sqrt((dv ** 2).sum(1)); U = dv / np.maximum(d, 1e-300)[:, None]
        r1 = np.repeat(np.arange(m), 7)
        c1 = np.stack([3 * I, 3 * I + 1, 3 * I + 2, 3 * J, 3 * J + 1, 3 * J + 2, np.full(m, 3 * n)], 1).ravel()
        v1 = np.stack([-U[:, 0] / 2, -U[:, 1] / 2, -U[:, 2] / 2, U[:, 0] / 2, U[:, 1] / 2, U[:, 2] / 2, np.ones(m)], 1).ravel()
        A1 = coo_matrix((v1, (r1, c1)), shape=(m, 3 * n + 1)); b1 = d / 2
        wi, wk = np.nonzero(g <= cur + 3 * delta + 1e-12); mw = len(wi)
        r2 = np.repeat(np.arange(mw), 4)
        c2 = np.stack([3 * wi, 3 * wi + 1, 3 * wi + 2, np.full(mw, 3 * n)], 1).ravel()
        v2 = np.stack([-G[wi, wk, 0], -G[wi, wk, 1], -G[wi, wk, 2], np.ones(mw)], 1).ravel()
        A2 = coo_matrix((v2, (r2, c2)), shape=(mw, 3 * n + 1)); b2 = g[wi, wk]
        res = linprog(np.r_[np.zeros(3 * n), -1.0], A_ub=vstack([A1, A2]).tocsr(), b_ub=np.r_[b1, b2],
                      bounds=[(-delta, delta)] * (3 * n) + [(None, None)], method=('highs-ipm' if n >= 1500 else 'highs'),
                      options={'time_limit': 120.0})
        if res.status != 0 or res.x is None or not np.all(np.isfinite(res.x)):
            delta /= 2
            if delta < 1e-16: break
            continue
        q = repair(c + res.x[:-1].reshape(n, 3)); new = rmin(q)
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
