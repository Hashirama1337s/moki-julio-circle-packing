#!/usr/bin/env python3
"""Local certificate for the N = 16 clusters (LOCAL.md §2), Moki&Julio.

Floats / mpmath only PROPOSE numbers (Newton-refined frame, LP stress, LP localisation vectors); every quantity that
enters the conclusion is recomputed in exact rational arithmetic (fractions.Fraction) from the rational data written
to out/local_cert16.json.  Skew coordinates (u, v) = (A/G, B/G); Q(du, dv) = du^2 + du dv + dv^2.

Usage: python local_cert.py [clusters.json] [rho_new_grid]
"""
import json
import os
import sys
import time
from fractions import Fraction as F

import mpmath as mp
import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'solver', 'chx'))  # chx_prio: lowers process priority only
import chx_prio  # noqa: E402

mp.mp.dps = 60
DEN = 2 ** 110


def Q(du, dv):
    return du * du + du * dv + dv * dv


def frame_setup(pts, frame):
    """pts: list of (u, v) floats for the 16 points.  Returns frame indices, tight pairs, walls (on frame indices)."""
    fr = [i for i in range(16) if frame[i]]
    n = len(fr)
    P = [pts[i] for i in fr]
    d2 = {}
    for a in range(n):
        for b in range(a + 1, n):
            d2[(a, b)] = Q(P[a][0] - P[b][0], P[a][1] - P[b][1])
    m2 = min(d2.values())
    pairs = sorted(k for k, v in d2.items() if v < m2 * (1 + 1e-5))
    walls = []
    for a in range(n):
        u, v = P[a]
        if u < 1e-6: walls.append((a, 'u'))
        if v < 1e-6: walls.append((a, 'v'))
        if 1 - u - v < 1e-6: walls.append((a, 'w'))
    return fr, P, pairs, walls


def parametrise(n, walls):
    """For each frame point: how (u, v) depend on free unknowns.  kind: 'free', 'u0' (u = 0, v free), 'v0' (v = 0,
    u free), 'w0' (v = 1 - u, u free), or a fixed corner."""
    ws = {a: set() for a in range(n)}
    for a, s in walls:
        ws[a].add(s)
    kinds = []
    for a in range(n):
        s = ws[a]
        if s == set():
            kinds.append('free')
        elif s == {'u'}:
            kinds.append('u0')
        elif s == {'v'}:
            kinds.append('v0')
        elif s == {'w'}:
            kinds.append('w0')
        elif s == {'u', 'v'}:
            kinds.append(('fix', 0, 0))
        elif s == {'u', 'w'}:
            kinds.append(('fix', 0, 1))
        elif s == {'v', 'w'}:
            kinds.append(('fix', 1, 0))
        else:
            raise ValueError(s)
    return kinds


def unpack(z, kinds, one):
    """z: list of free unknowns (numbers of any type); returns list of (u, v)."""
    out = []; p = 0
    for kd in kinds:
        if kd == 'free':
            out.append((z[p], z[p + 1])); p += 2
        elif kd == 'u0':
            out.append((0 * one, z[p])); p += 1
        elif kd == 'v0':
            out.append((z[p], 0 * one)); p += 1
        elif kd == 'w0':
            out.append((z[p], one - z[p])); p += 1
        else:
            out.append((kd[1] * one, kd[2] * one))
    return out


def pack(P, kinds):
    z = []
    for (u, v), kd in zip(P, kinds):
        if kd == 'free':
            z += [u, v]
        elif kd in ('u0',):
            z.append(v)
        elif kd in ('v0', 'w0'):
            z.append(u)
    return z


def newton(P, kinds, pairs, iters=12):
    z = [mp.mpf(x) for x in pack(P, kinds)]
    s = mp.mpf(min(Q(P[a][0] - P[b][0], P[a][1] - P[b][1]) for a, b in pairs))
    nz = len(z)
    one = mp.mpf(1)
    for it in range(iters):
        X = unpack(z, kinds, one)
        Fv = mp.matrix(len(pairs), 1)
        Jm = mp.matrix(len(pairs), nz + 1)
        # numeric-exact derivative via dual approach: analytic partials through the parametrisation
        for r, (a, b) in enumerate(pairs):
            du = X[a][0] - X[b][0]; dv = X[a][1] - X[b][1]
            Fv[r] = Q(du, dv) - s
            gu = 2 * du + dv; gv = du + 2 * dv
            for pt, sign in ((a, 1), (b, -1)):
                # column indices of point pt
                p = 0
                for q, kd in enumerate(kinds):
                    width = 2 if kd == 'free' else (1 if kd in ('u0', 'v0', 'w0') else 0)
                    if q == pt:
                        if kd == 'free':
                            Jm[r, p] += sign * gu; Jm[r, p + 1] += sign * gv
                        elif kd == 'u0':
                            Jm[r, p] += sign * gv
                        elif kd == 'v0':
                            Jm[r, p] += sign * gu
                        elif kd == 'w0':
                            Jm[r, p] += sign * (gu - gv)
                        break
                    p += width
            Jm[r, nz] = -1
        JT = Jm.T
        step = mp.lu_solve(JT * Jm, -(JT * Fv))
        for q in range(nz):
            z[q] += step[q]
        s += step[nz]
        res = max(abs(Fv[r]) for r in range(len(pairs)))
    X = unpack(z, kinds, one)
    Fv = [Q(X[a][0] - X[b][0], X[a][1] - X[b][1]) - s for a, b in pairs]
    return z, s, max(abs(f) for f in Fv)


def jac_rows(X, pairs, walls, n):
    """Exact (Fraction) gradient rows, length 2n; pairs first, then walls."""
    rows = []
    for a, b in pairs:
        du = X[a][0] - X[b][0]; dv = X[a][1] - X[b][1]
        g = [F(0)] * (2 * n)
        g[2 * a] = 2 * du + dv; g[2 * a + 1] = du + 2 * dv
        g[2 * b] = -(2 * du + dv); g[2 * b + 1] = -(du + 2 * dv)
        rows.append(g)
    for a, s in walls:
        g = [F(0)] * (2 * n)
        if s == 'u': g[2 * a] = F(1)
        elif s == 'v': g[2 * a + 1] = F(1)
        else: g[2 * a] = F(-1); g[2 * a + 1] = F(-1)
        rows.append(g)
    return rows


def wall_val(X, a, s):
    u, v = X[a]
    return u if s == 'u' else (v if s == 'v' else 1 - u - v)


def certify(cl, G, rho_new=None):
    t0 = time.time()
    pts16 = [(p[0] / G, p[1] / G) for p in cl['points']]
    fr, P, pairs, walls = frame_setup(pts16, cl['frame'])
    n = len(fr); npairs = len(pairs)
    kinds = parametrise(n, walls)
    z, s, res = newton(P, kinds, pairs)
    # rational frame
    zr = [F(int(mp.nint(x * DEN)), DEN) for x in z]
    X = unpack(zr, kinds, F(1))
    for a in range(n):
        assert X[a][0] >= 0 and X[a][1] >= 0 and X[a][0] + X[a][1] <= 1
    rows = jac_rows(X, pairs, walls, n)
    K = len(rows)
    Jf = np.array([[float(x) for x in r] for r in rows])
    gval = []  # g_k(c~; m) = base_k - m^2 * [k is pair]
    for a, b in pairs:
        gval.append(Q(X[a][0] - X[b][0], X[a][1] - X[b][1]))
    for a, sw in walls:
        gval.append(wall_val(X, a, sw))
    # 16-point configuration for the lower bound: frame + rattler (rational grid position)
    rat = [i for i in range(16) if not cl['frame'][i]]
    full = {}
    for q, i in enumerate(fr):
        full[i] = X[q]
    for i in rat:
        # rattler: maximise its min distance to the frame (float search), then rationalise; checked exactly below
        from scipy.optimize import minimize
        Xf = np.array([[float(a), float(b)] for a, b in X])
        def negmin(w):
            d = Xf - w
            return -np.min(d[:, 0] ** 2 + d[:, 0] * d[:, 1] + d[:, 1] ** 2)
        w0 = np.array([cl['points'][i][0] / G, cl['points'][i][1] / G])
        best = w0
        for trial in range(20):
            st = w0 + (np.random.default_rng(trial).random(2) - 0.5) * 0.01 if trial else w0
            rr = minimize(negmin, st, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-16, 'maxiter': 4000})
            if negmin(rr.x) < negmin(best) and rr.x[0] >= 0 and rr.x[1] >= 0 and rr.x.sum() <= 1:
                best = rr.x
        full[i] = (F(int(round(best[0] * 2 ** 40)), 2 ** 40), F(int(round(best[1] * 2 ** 40)), 2 ** 40))
    for i in range(16):
        u, v = full[i]
        assert u >= 0 and v >= 0 and u + v <= 1
    mlo2 = min(Q(full[i][0] - full[j][0], full[i][1] - full[j][1]) for i in range(16) for j in range(i + 1, 16))
    # stress LP
    A_eq = np.zeros((2 * n + 1, K + 1)); A_eq[:2 * n, :K] = Jf.T; A_eq[2 * n, :npairs] = 1
    b_eq = np.zeros(2 * n + 1); b_eq[-1] = 1
    A_ub = np.zeros((K, K + 1)); A_ub[:, :K] = -np.eye(K); A_ub[:, K] = 1
    lp = linprog(np.r_[np.zeros(K), -1], A_ub=A_ub, b_ub=np.zeros(K), A_eq=A_eq, b_eq=b_eq,
                 bounds=[(0, None)] * K + [(None, None)], method='highs')
    assert lp.success
    lam = [F(max(float(x), 0.0)) for x in lp.x[:K]]
    Lam = sum(lam[:npairs])
    r = [sum(lam[k] * rows[k][c] for k in range(K)) for c in range(2 * n)]
    r1 = sum(abs(x) for x in r)
    Phi = sum(lam[k] * gval[k] for k in range(K))
    # localisation LPs
    A_ = F(0); E_ = F(0); B_ = None; worst = None
    ys = []
    cpair = np.r_[np.ones(npairs), np.zeros(K - npairs)]
    for i in range(2 * n):
        for sg in (1, -1):
            e = np.zeros(2 * n); e[i] = -sg
            lpy = linprog(cpair, A_eq=Jf.T, b_eq=e, bounds=[(0, None)] * K, method='highs')
            assert lpy.success, (i, sg)
            y = [F(max(float(x), 0.0)) for x in lpy.x]
            eps = [sum(y[k] * rows[k][c] for k in range(K)) for c in range(2 * n)]
            eps[i] += sg
            e1 = sum(abs(x) for x in eps)
            Aq = 12 * sum(y[:npairs])
            Bq = sum(y[k] * (gval[k] - (mlo2 if k < npairs else 0)) for k in range(K))
            A_ = max(A_, Aq); E_ = max(E_, e1)
            B_ = Bq if B_ is None else max(B_, Bq)
            ys.append([str(x) for x in y])
    assert E_ < 1
    R = (1 - E_) / (2 * A_)
    delta = max(F(0), 2 * B_ / (1 - E_))
    mhi2 = (Phi + r1 * delta + 12 * Lam * delta * delta) / Lam
    # global box vs ball
    off = max(max(abs(F(cl['points'][i][0], G) - full[i][0]), abs(F(cl['points'][i][1], G) - full[i][1])) for i in fr)
    rho_g = F(cl['rho'], G)
    out = {
        'cluster': cl['id'], 'frame_points': n, 'tight_pairs': npairs, 'tight_walls': len(walls),
        'newton_residual': mp.nstr(res, 5), 'min_lambda_float': float(lp.x[K]),
        'A': float(A_), 'E': float(E_), 'R': float(R), 'B_mlo': float(B_), 'delta': float(delta),
        'Lambda': float(Lam), 'r1': float(r1),
        'm_lo': mp.nstr(mp.sqrt(mp.mpf(mlo2.numerator) / mlo2.denominator), 40),
        'm_hi': mp.nstr(mp.sqrt(mp.mpf(mhi2.numerator) / mhi2.denominator), 40),
        'm_hi_minus_m_lo': mp.nstr(mp.sqrt(mp.mpf(mhi2.numerator) / mhi2.denominator)
                                   - mp.sqrt(mp.mpf(mlo2.numerator) / mlo2.denominator), 5),
        'global_box_offset': float(off), 'global_box_rho': float(rho_g),
        'global_box_fits': bool(off + rho_g <= R),
        'rho_new_grid_max': int((R - off) * G) if R > off else 0,
        'seconds': round(time.time() - t0, 1),
    }
    data = {'frame_indices': fr, 'pairs': pairs, 'walls': walls, 'kinds': [str(k) for k in kinds],
            'frame_uv': [[str(a), str(b)] for a, b in X],
            'rattlers_uv': {str(i): [str(full[i][0]), str(full[i][1])] for i in rat},
            'm_lo_sq': str(mlo2), 'lambda': [str(x) for x in lam], 'y': ys,
            'A': str(A_), 'E': str(E_), 'R': str(R), 'B': str(B_), 'delta': str(delta), 'm_hi_sq': str(mhi2)}
    return out, data


def main():
    chx_prio.lower()
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'clusters16.json')
    cls = json.load(open(path))
    G = cls['G']
    res = []; datas = []
    for cl in cls['clusters']:
        o, d = certify(cl, G)
        print(json.dumps(o), flush=True)
        res.append(o); datas.append(d)
    json.dump({'summary': res, 'data': datas}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'local_cert16.json'), 'w'))


if __name__ == '__main__':
    main()
