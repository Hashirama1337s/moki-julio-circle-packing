#!/usr/bin/env python3
"""Exact re-check of the local certificate data (out/local_cert16.json) and of its link to a global cluster file.
Standard library only (fractions); no floats, no LP, no Newton: every claimed number is recomputed from the stored
rationals and every inequality is decided exactly.  Moki&Julio.

Checks per cluster (LOCAL.md §2):
  frame points and rattler(s) inside the triangle; the stored pair/wall lists are genuine constraints at the frame
  (the argument does not need them to be tight: any subset of valid constraints works);
  m_lo^2 = exact min of Q over all 120 pairs of the 16-point configuration;
  lambda >= 0, y^{i,s} >= 0; residuals r = J^T lambda, eps = J^T y + s e_i recomputed; A, E, B(m_lo) recomputed;
  R = (1 - E)/(2A), delta = 2 max(B, 0)/(1 - E), m_hi^2 = (Phi + |r|_1 delta + 12 Lambda delta^2)/Lambda recomputed;
  link: every frame point of the global box (centre P/G, radius rho/G in skew sup-norm) lies in the R-ball of c~,
  and the global target m_t <= m_lo.
Usage: python check_local.py out/local_cert16.json out/clusters16_rho4000.json <m_t as p/q>
"""
import json
import sys
from fractions import Fraction as F


def Q(du, dv):
    return du * du + du * dv + dv * dv


def main():
    cert = json.load(open(sys.argv[1]))
    clus = json.load(open(sys.argv[2]))
    mt = F(sys.argv[3])
    G = clus['G']
    ok_all = True
    lows = []; highs = []
    for q, (summ, d) in enumerate(zip(cert['summary'], cert['data'])):
        cl = clus['clusters'][q]
        fr = d['frame_indices']; n = len(fr)
        assert fr == [i for i in range(16) if cl['frame'][i]], 'frame indices differ from the cluster file'
        X = [(F(a), F(b)) for a, b in d['frame_uv']]
        pts = {fr[t]: X[t] for t in range(n)}
        for k_, uv in d['rattlers_uv'].items():
            pts[int(k_)] = (F(uv[0]), F(uv[1]))
        assert sorted(pts) == list(range(16))
        for u, v in pts.values():
            assert u >= 0 and v >= 0 and u + v <= 1, 'point outside the triangle'
        mlo2 = min(Q(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]) for i in range(16) for j in range(i + 1, 16))
        assert mlo2 == F(d['m_lo_sq']), 'm_lo^2 mismatch'
        pairs = [tuple(p) for p in d['pairs']]; walls = [tuple(w) for w in d['walls']]
        assert all(0 <= a < b < n for a, b in pairs) and len(set(pairs)) == len(pairs)
        assert all(0 <= a < n and s in ('u', 'v', 'w') for a, s in walls) and len(set(walls)) == len(walls)
        rows = []; base = []
        for a, b in pairs:
            du = X[a][0] - X[b][0]; dv = X[a][1] - X[b][1]
            g = [F(0)] * (2 * n)
            g[2 * a] = 2 * du + dv; g[2 * a + 1] = du + 2 * dv
            g[2 * b] = -(2 * du + dv); g[2 * b + 1] = -(du + 2 * dv)
            rows.append(g); base.append(Q(du, dv) - mlo2)
        for a, s in walls:
            g = [F(0)] * (2 * n)
            if s == 'u': g[2 * a] = F(1); w = X[a][0]
            elif s == 'v': g[2 * a + 1] = F(1); w = X[a][1]
            else: g[2 * a] = F(-1); g[2 * a + 1] = F(-1); w = 1 - X[a][0] - X[a][1]
            rows.append(g); base.append(w)
        K = len(rows); P = len(pairs)
        lam = [F(x) for x in d['lambda']]
        assert len(lam) == K and all(x >= 0 for x in lam)
        Lam = sum(lam[:P]); assert Lam > 0
        r1 = sum(abs(sum(lam[k] * rows[k][c] for k in range(K))) for c in range(2 * n))
        Phi = sum(lam[k] * (base[k] + (mlo2 if k < P else 0)) for k in range(K))
        A = F(0); E = F(0); B = None
        assert len(d['y']) == 4 * n
        idx = 0
        for i in range(2 * n):
            for sg in (1, -1):
                y = [F(x) for x in d['y'][idx]]; idx += 1
                assert len(y) == K and all(x >= 0 for x in y)
                eps = [sum(y[k] * rows[k][c] for k in range(K)) for c in range(2 * n)]
                eps[i] += sg
                E = max(E, sum(abs(x) for x in eps))
                A = max(A, 12 * sum(y[:P]))
                b = sum(y[k] * base[k] for k in range(K))
                B = b if B is None else max(B, b)
        assert E < 1
        R = (1 - E) / (2 * A)
        delta = 2 * max(B, F(0)) / (1 - E)
        mhi2 = (Phi + r1 * delta + 12 * Lam * delta * delta) / Lam
        same = (R == F(d['R']) and delta == F(d['delta']) and mhi2 == F(d['m_hi_sq']) and A == F(d['A'])
                and E == F(d['E']))
        # link to the global boxes
        rho = F(cl['rho'], G)
        off = max(max(abs(F(cl['points'][i][0], G) - pts[i][0]), abs(F(cl['points'][i][1], G) - pts[i][1])) for i in fr)
        fits = off + rho <= R
        below = mt * mt <= mlo2
        ok = same and fits and below
        ok_all &= ok
        lows.append(mlo2); highs.append(mhi2)
        print(f'cluster {q}: recomputed == stored: {same}; R = {float(R):.4e}; delta = {float(delta):.3e}; '
              f'box offset + rho = {float(off + rho):.4e} <= R: {fits}; m_t <= m_lo: {below}; '
              f'm_hi^2 - m_lo^2 = {float(mhi2 - mlo2):.3e}')
    print(('PASS' if ok_all else 'FAIL') + f': given the global stage, m16^2 in [{float(max(lows))!r}, {float(max(highs))!r}] '
          f'(exact rationals in the certificate), width {float(max(highs) - max(lows)):.3e}')
    sys.exit(0 if ok_all else 1)


if __name__ == '__main__':
    main()
