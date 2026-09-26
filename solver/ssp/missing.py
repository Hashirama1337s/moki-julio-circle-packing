# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Reconstruct the MISSING LAST SPHERE of Packomania's ssp files (2026-09-26; Moki&Julio). Float only.

Data-level finding (recon.py): every file ssp<N>.txt, N = 2..1000 (tarball txt/ssp_coords.tar.gz of 17-Jul-2026 AND the
online txt/ssp<N>.txt), lists N rows but row N is BLANK (index only, no coordinates): the files hold N - 1 centres.
The N-th centre is recovered here as the point of largest clearance
    t(p) = min( min_i |p - c_i| - 2 r,  1 - r - |p| )        (r = the printed radius)
among the N - 1 given centres: candidates = circumcentres of the Delaunay tetrahedra of the N - 1 centres (the Voronoi
vertices: every interior local maximum of t), + 60 000 random points in the outer shell (holes against the wall) and
20 000 in the whole ball; the best 40 by t are refined by a small trust-region LP (variables dp in R^3 and t; pair rows
inner-linearised, the wall outer-linearised + a radial pull-back; step kept only if the TRUE t grows) to 1e-16.
The best refined point is the reconstruction; its clearance t* is reported: t* >= -1e-11 means the full N-packing at
the printed radius is recovered to file precision (t* > 0 = the missing sphere is a rattler with room t*).
reconstruct(n) -> (centres (n, 3), t*, info). The second-best hole's clearance is reported too (uniqueness).
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sspgeo                                               # noqa: E402
import numpy as np
from scipy.spatial import Delaunay, cKDTree
from scipy.optimize import linprog
import geomd                                                # noqa: E402


def load_partial(n):
    """The N - 1 centres Packomania's file gives, and the printed radius (float)."""
    L = [l.split() for l in open(sspgeo.pub_path(n)) if l.strip()]
    rows = [t for t in L[1:] if not t[0].startswith('#')]
    full = [t for t in rows if len(t) >= 4]; blank = [t for t in rows if len(t) == 1]
    return np.array([[float(v) for v in t[1:4]] for t in full], dtype=np.float64), L[0][0], [int(t[0]) for t in blank]


def clear(P, C, r, tree=None):
    tree = tree or cKDTree(C)
    dd, _ = tree.query(P, k=1)
    return np.minimum(dd - 2 * r, 1 - r - geomd.norms(P))


def refine(p, C, r, tree, iters=400):
    p = p.copy(); t = float(clear(p[None], C, r, tree)[0]); delta = 0.05 * r
    for _ in range(iters):
        idx = tree.query_ball_point(p, 2 * r + max(t, 0) + 4 * delta + 1e-9)
        A = []; b = []
        for i in idx:
            v = p - C[i]; dv = np.linalg.norm(v); u = v / dv
            A.append([-u[0], -u[1], -u[2], 1.0]); b.append((dv - 2 * r) / delta)
        nr = np.linalg.norm(p); w = p / max(nr, 1e-300)
        A.append([w[0], w[1], w[2], 1.0]); b.append((1 - r - nr) / delta)
        # t in scaled units: t_new = delta * s
        res = linprog([0, 0, 0, -1.0], A_ub=np.array(A), b_ub=np.array(b), bounds=[(-1, 1)] * 3 + [(None, None)], method='highs')
        if res.status != 0:
            delta /= 4
            if delta < 1e-17: break
            continue
        q = p + delta * res.x[:3]; tp = delta * res.x[3]
        nq = np.linalg.norm(q)
        if nq > 1 - r - tp and nq > 0: q = q * ((1 - r - tp) / nq)
        tq = float(clear(q[None], C, r, tree)[0])
        if tq > t + 1e-17:
            p, t = q, tq; delta = min(delta * 1.5, 0.2 * r)
        else:
            delta /= 4
            if delta < 1e-17: break
    return p, t


def candidates(C, r, rng):
    n = len(C); out = []
    if n >= 4:
        try:
            D = Delaunay(C); S = C[D.simplices]                      # (m, 4, 3) circumcentres
            a = S[:, 1:] - S[:, :1]; rhs = 0.5 * (a ** 2).sum(-1)
            ok = np.abs(np.linalg.det(a)) > 1e-14
            cc = np.linalg.solve(a[ok], rhs[ok][..., None])[..., 0] + S[ok, 0]
            nr = geomd.norms(cc); lim = 1 - r
            cc[nr > lim] *= (lim / nr[nr > lim])[:, None]
            out.append(cc)
        except Exception:
            pass
    v = geomd.random_sphere(60000, 3, rng); out.append(v * (1 - r - r * rng.random(60000))[:, None])
    out.append(geomd.random_ball(20000, 3, 1 - r, rng))
    if n <= 3: out.append(-C.sum(0, keepdims=True))
    return np.vstack(out)


def reconstruct(n, seed=0, top=15):
    C, r_s, blank = load_partial(n); r = float(r_s); rng = np.random.default_rng(seed + n)
    if len(C) == n: return C, None, dict(n=n, missing=0)
    assert len(C) == n - 1 and blank == [n], (n, len(C), blank)
    if n == 2:
        p = -C[0]; t = float(clear(p[None], C, r)[0]); return np.vstack([C, p]), t, dict(n=n, missing=1, t=t, t2=None)
    tree = cKDTree(C); P = candidates(C, r, rng); tc = clear(P, C, r, tree)
    order = np.argsort(-tc)[:top * 5]
    # de-duplicate candidates within 0.5 r of a better one
    pick = []
    for k in order:
        if all(np.linalg.norm(P[k] - P[j]) > 0.5 * r for j in pick): pick.append(k)
        if len(pick) >= top: break
    ref = sorted((refine(P[k], C, r, tree) for k in pick), key=lambda x: -x[1])
    # distinct holes
    holes = []
    for p, t in ref:
        if all(np.linalg.norm(p - q) > 0.5 * r for q, _ in holes): holes.append((p, t))
    p, t = holes[0]
    c = np.vstack([C, p]); info = dict(n=n, missing=1, t=float(t), t2=float(holes[1][1]) if len(holes) > 1 else None,
                                       p=[float(x) for x in p], p_norm=float(np.linalg.norm(p)), rmin_full=geomd.rmin(c),
                                       gap_full_vs_printed=geomd.rmin(c) - r)
    return c, float(t), info


if __name__ == '__main__':
    import prio, json
    prio.lower()
    Ns = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 and sys.argv[1] != 'all' else list(range(2, 1001))
    outd = os.path.join(sspgeo.BIG, 'full'); os.makedirs(outd, exist_ok=True)
    rows_p = os.path.join(sspgeo.OUT, 'missing.jsonl'); os.makedirs(sspgeo.OUT, exist_ok=True)
    T0 = time.time()
    for n in Ns:
        t0 = time.time(); c, t, info = reconstruct(n); info['secs'] = round(time.time() - t0, 2)
        np.save(os.path.join(outd, f'ssp{n}.npy'), c)
        with open(rows_p, 'a') as f: f.write(json.dumps(info) + '\n')
        print(n, 't* %.3e' % info['t'], 't2 %s' % (None if info['t2'] is None else '%.3e' % info['t2']),
              '|p| %.4f' % info.get('p_norm', 0), 'gap_full %.3e' % info.get('gap_full_vs_printed', 0), '%.1fs' % info['secs'],
              '[%.0fs]' % (time.time() - T0), flush=True)
