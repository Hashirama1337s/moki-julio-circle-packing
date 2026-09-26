# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Trust-region SLP polish for equal balls in the unit d-ball (2026-09-25; Moki&Julio). The d-dimensional, ball-container
port of the cube solver's slp3d.polish (../scu/slp3d.py) (the repaired version), with every degeneracy fix of that instrument kept:

LP (maximise s) around the current centres c (n x d), current true min-radius cur, trust region delta:
  variables  y = dx / delta in [-1, 1]^(n d)  (scaled step),  s = (t - cur) / delta  (free);
  pair rows  (i < j, |c_i - c_j| / 2 <= cur + margin delta):   s - U/2 . y_i + U/2 . y_j <= (d_ij / 2 - cur) / delta
             U = (c_i - c_j) / d_ij. INNER linearisation (|a + h| >= U . (a + h), Cauchy-Schwarz): a feasible step keeps
             every listed pair's half distance >= t.
  wall rows  (1 - |c_i| <= cur + margin delta):                s + u_i . y_i <= (1 - |c_i| - cur) / delta
             u_i = c_i / |c_i|. OUTER linearisation of the convex constraint |c_i| <= 1 - t (u . (c + h) <= |c + h|), so a
             step can overshoot the wall by O(delta^2) (|h_perp|^2 / 2|c|).
  margin = sqrt(d) + 1.
Fixes carried over from slp3d (see its docstring for the measurements behind them):
  (1) SCALED variables, so every solver tolerance is relative to the step;
  (2) every row TIGHTENED by an independent random amount in [pert/2, pert] (default 1e-9, scaled units) -> the start
      vertex of a massively degenerate (heavily contacted) packing is no longer degenerate;
  (3) tight HiGHS tolerances (primal / dual feasibility 1e-10), a per-LP time limit; interior point first, dual simplex
      (half the time) only if the IPM fails; a failed LP quarters the trust region;
  (4) GUARD: a step is kept only if the TRUE float min-radius (geomd.rmin: all pairs + the curved wall) grows.
  STOP RULE: LP optimum s <= 2 pert (no first-order improving direction) and the step does not gain -> first-order
  STATIONARY (the LP dual is a stress on the tight rows); before stopping, the same LP is re-solved once without the
  perturbation (its s = 0 step may still gain through curvature). Five kept steps in a row with s <= 2 pert ->
  'near_stationary'. pretest=True: one LP at delta = 1e-6 cur detects a stationary start; otherwise the trust region
  restarts at 1e-3 cur and grows x4 while below 0.01 cur and fast, else x1.5 (cap 0.2 cur).
Ball-specific REPAIR (the analogue of the cube's clip): after the LP step, a centre that overshoots the curved wall
(|q_i| > 1 - t_pred, t_pred = cur + delta s) is pulled back radially to |q_i| = 1 - t_pred. The pull is O(delta^2); the
guard (4) still decides on the TRUE min-radius of the repaired step.
Budgets are CPU seconds of this process (time.process_time) unless clock=time.time is passed.
"""
import time, numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack
import geomd


def repair(q, t_pred):
    nr = geomd.norms(q); lim = 1.0 - t_pred; bad = nr > lim
    if bad.any():
        q = q.copy(); q[bad] *= (lim / nr[bad])[:, None]
    return q


def _rows(c, cur, delta, margin):
    n, dim = c.shape
    D2 = geomd.pair_d2(c); lim = 2 * (cur + margin * delta) + 1e-12
    I, J = np.nonzero(np.triu(D2 <= lim * lim, 1))
    dv = c[I] - c[J]; dd = geomd.norms(dv); U = dv / np.maximum(dd, 1e-300)[:, None]
    nr = geomd.norms(c); g = 1.0 - nr
    wi = np.nonzero(g <= cur + margin * delta + 1e-12)[0]
    u = c[wi] / np.maximum(nr[wi], 1e-300)[:, None]
    return (I, J, U, dd), (wi, u, g[wi])


def _lp(c, cur, delta, margin, pert, rng, perturb):
    n, dim = c.shape; nv = dim * n + 1
    (I, J, U, dd), (wi, u, gw) = _rows(c, cur, delta, margin); m, mw = len(I), len(wi)
    ar = np.arange(dim)
    # pair rows
    r1 = np.repeat(np.arange(m), 2 * dim + 1)
    c1 = np.concatenate([dim * I[:, None] + ar[None, :], dim * J[:, None] + ar[None, :], np.full((m, 1), nv - 1)], 1).ravel()
    v1 = np.concatenate([-U / 2, U / 2, np.ones((m, 1))], 1).ravel()
    A1 = coo_matrix((v1, (r1, c1)), shape=(m, nv)); b1 = (dd / 2 - cur) / delta
    # wall rows
    r2 = np.repeat(np.arange(mw), dim + 1)
    c2 = np.concatenate([dim * wi[:, None] + ar[None, :], np.full((mw, 1), nv - 1)], 1).ravel()
    v2 = np.concatenate([u, np.ones((mw, 1))], 1).ravel()
    A2 = coo_matrix((v2, (r2, c2)), shape=(mw, nv)); b2 = (gw - cur) / delta
    b = np.r_[b1, b2]
    if perturb: b = b - pert * (0.5 + 0.5 * rng.random(len(b)))
    return vstack([A1, A2]).tocsr(), b, m, wi


def _solve(A, b, nv, tries, tol):
    """HiGHS: each (method, time limit) in turn until one returns a finite optimum. -> (x or None, n_fallbacks)."""
    nfb = 0
    for k, (meth, tlim) in enumerate(tries):
        if meth is None: continue
        try:
            res = linprog(np.r_[np.zeros(nv - 1), -1.0], A_ub=A, b_ub=b, bounds=[(-1.0, 1.0)] * (nv - 1) + [(None, None)],
                          method=meth, options={'time_limit': float(tlim), 'primal_feasibility_tolerance': tol,
                                                'dual_feasibility_tolerance': tol})
            if res.status == 0 and res.x is not None and np.all(np.isfinite(res.x)): return res, nfb
        except Exception:
            pass
        nfb += 1
    return None, nfb


def polish(c, max_iter=100000, t_cap=120.0, log=None, stop_at=None, lp_time=30.0, pert=1e-9, tol=1e-10, seed=0,
           stats=None, pretest=False, start_after_pretest=1e-3, method='highs-ipm', fallback='highs',
           clock=time.process_time, delta0=0.1, delta_max=0.2, soc=2):
    """Repaired SLP polish in the unit d-ball (module docstring). Returns (centres, TRUE float min-radius).
    soc = number of SECOND-ORDER CORRECTIONS of the wall rows per step: when the LP step overshoots the curved wall by
    eps_i > 5 % of the predicted gain, the wall rows' right-hand sides are lowered by the accumulated eps_i / delta and the
    LP is solved again (the standard SQP correction; measured at hsp4 N = 119: without it no step above delta = 1e-4 r is
    kept, with it delta = 1e-3 r steps gain 1.1e-8 each). Stationarity is decided on the FIRST (uncorrected) LP."""
    rng = np.random.default_rng(seed)
    c = np.asarray(c, dtype=np.float64).copy(); n, dim = c.shape; margin = np.sqrt(dim) + 1.0; nv = dim * n + 1
    cur = geomd.rmin(c); delta = delta0 * max(cur, 1e-6); t0 = clock()
    if pretest: delta = 1e-6 * cur
    st = stats if stats is not None else {}
    for k in ('lp', 'lp_fail', 'lp_secs', 'kept', 'rejected', 'soc'): st.setdefault(k, 0)
    st.pop('stop', None); nflat = 0; retry_flat = False
    for it in range(max_iter):
        if clock() - t0 > t_cap: st['stop'] = 'time'; break
        A, b, m, wi = _lp(c, cur, delta, margin, pert, rng, perturb=not retry_flat)
        st['rows'] = A.shape[0]
        tl = time.time()
        tries = ((method, min(lp_time, 10.0)),) if retry_flat else ((method, lp_time), (fallback, lp_time / 2))
        res, nfb = _solve(A, b, nv, tries, tol); st['lp'] += 1
        if nfb: st['lp_fallback'] = st.get('lp_fallback', 0) + nfb
        if res is None:
            st['lp_secs'] += time.time() - tl; st['lp_fail'] += 1; delta /= 4
            if delta < 1e-16: st['stop'] = 'delta'; break
            continue
        x = res.x; s0 = float(x[-1]); flat = s0 <= 2 * pert
        if flat:
            try:
                lam = -np.asarray(res.ineqlin.marginals); st['stress_rows'] = int((lam > 1e-9).sum())
            except Exception:
                pass
        nsoc = 0
        if not flat and soc and len(wi):
            acc = np.zeros(len(wi))
            for _ in range(soc):
                s = float(x[-1]); q0 = c + delta * x[:-1].reshape(n, dim)
                over = geomd.norms(q0[wi]) - (1.0 - (cur + delta * max(s, 0.0)))
                if over.max() <= max(0.05 * delta * max(s, 0.0), 1e-15): break
                acc += np.maximum(over, 0.0); b2 = b.copy(); b2[m:] -= acc / delta
                r2, _ = _solve(A, b2, nv, ((method, lp_time),), tol); st['lp'] += 1; st['soc'] += 1
                if r2 is None: break
                x = r2.x; nsoc += 1
        st['lp_secs'] += time.time() - tl
        s = float(x[-1]); st['last_s'] = s; st['last_delta'] = delta
        q = repair(c + delta * x[:-1].reshape(n, dim), cur + delta * max(s, 0.0)); new = geomd.rmin(q)
        if new > cur + 1e-16:
            gain = new - cur; c, cur = q, new; st['kept'] += 1; retry_flat = False
            if log: log(it, cur, delta, clock() - t0)
            nflat = nflat + 1 if flat else 0
            if nflat >= 5: st['stop'] = 'near_stationary'; break
            if gain < 1e-15 and delta < 1e-12: st['stop'] = 'tiny'; break
            if stop_at is not None and cur >= stop_at: st['stop'] = 'stop_at'; break
            # trust-region update by the ratio actual / predicted gain
            # (curvature: a step that needed a correction, or whose corrected s fell below half the first s, is at the
            # largest useful trust region: hold or shrink it instead of growing into a rejection)
            ratio = gain / (delta * s) if s > 2 * pert else 1.0
            if nsoc and (s < 0.5 * s0 or ratio < 0.5):
                grow = 0.5
            elif nsoc:
                grow = 1.0
            elif ratio > 0.75:
                grow = 4.0 if (delta < 0.01 * cur and ratio > 0.9 and time.time() - tl < lp_time / 4) else 1.5
            elif ratio > 0.25:
                grow = 1.0
            else:
                grow = 0.5
            delta = min(delta * grow, delta_max * cur)
            if pretest and it == 0: delta = start_after_pretest * cur
        else:
            st['rejected'] += 1
            if pretest and it == 0 and not flat:
                delta = start_after_pretest * cur; continue
            if flat and not retry_flat and not (pretest and it == 0):
                retry_flat = True; continue
            if flat:
                st['stop'] = 'stationary_pretest' if (pretest and it == 0) else 'stationary'; break
            delta /= (2 if nsoc else 4)
            if delta < 1e-16: st['stop'] = 'delta'; break
    else:
        st['stop'] = 'max_iter'
    st.setdefault('stop', 'time')
    return c, float(cur)
