# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Seeds for equal balls in the unit d-ball (2026-09-25; Moki&Julio). Float only.

(a) packomania(d, N)        Packomania's own coordinates (hsp4 only; the hsp5 / hsp6 files are truncated to 4 columns).
(b) code_seed(d, M, centre) a spherical code of M points (local Cohn-table file c<d>_<M>.txt, or cell24() / e6_roots())
                            on the shell |c| = 1 - r with r = s / (1 + s), s = sin(theta/2), theta = the code's minimal
                            angle; with centre=True one more ball at the origin (needs 1 - r >= 2 r, i.e. r <= 1/3; r is
                            capped at 1/3, where the shell code still fits because s > 1/2).
    shell_core(d, N, k, ...) code shell of N - k points (sphere_code: a penalty code at the Cohn-table angle for N - k
                            points) + k inner balls (a smaller penalty code on an inner shell, the rest random, or all
                            random), then inflate().
(c) cold(d, N, ...)         random points in the ball + inflate().
(d) insert / delete         neighbour seeds: insert k balls into the k largest holes of a packing of N - k balls, or
                            delete the k fewest-contact balls of a packing of N + k balls.
inflate(): L-BFGS on the penalty  sum_{i<j} max(0, (2 rho)^2 - |ci - cj|^2)^2 + sum_i max(0, |ci|^2 - (1 - rho)^2)^2
at rho = f * r_ref for a rising schedule f (a Lubachevsky-Stillinger-like inflation); the result is then polished by
slpd.polish (the only step that decides the radius).
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
import geomd


# ---------------------------------------------------------------- penalty relaxation in the ball
def _pen(x, n, d, rho):
    c = x.reshape(n, d); sq = (c * c).sum(1)
    D2 = sq[:, None] + sq[None, :] - 2 * (c @ c.T)
    H = 4 * rho * rho - D2; np.fill_diagonal(H, 0.0); np.maximum(H, 0.0, out=H)
    W = np.maximum(sq - (1 - rho) ** 2, 0.0)
    E = 0.5 * (H * H).sum() + (W * W).sum()
    g = -4 * (H.sum(1)[:, None] * c - H @ c) + 4 * W[:, None] * c
    k = 1.0 / rho ** 4
    return E * k, g.ravel() * k


def relax(c, rho, maxiter=600):
    c = np.asarray(c, dtype=np.float64); n, d = c.shape
    res = minimize(_pen, c.ravel(), args=(n, d, rho), jac=True, method='L-BFGS-B',
                   options=dict(maxiter=maxiter, gtol=1e-13, ftol=1e-16, maxcor=20))
    return res.x.reshape(n, d), float(res.fun)


def inflate(c, r_ref, stages=(0.8, 0.9, 0.96, 1.0, 1.004), maxiter=600):
    for f in stages: c, _ = relax(c, f * r_ref, maxiter)
    return c


# ---------------------------------------------------------------- spherical codes (penalty on the sphere)
def _sph(x, M, d, D2t):
    Y = x.reshape(M, d); ny = geomd.norms(Y); X = Y / ny[:, None]
    H = D2t - (2 - 2 * (X @ X.T)); np.fill_diagonal(H, 0.0); np.maximum(H, 0.0, out=H)
    E = 0.5 * (H * H).sum(); gX = 4 * (H @ X)
    gY = (gX - (gX * X).sum(1)[:, None] * X) / ny[:, None]
    return E, gY.ravel()


def sphere_code(M, d, rng, cos_target, X0=None, stages=(0.9, 0.97, 1.0, 1.003), maxiter=800):
    """M unit vectors in R^d with (about) the minimal angle arccos(cos_target): penalty on the chord, rising schedule."""
    X = geomd.random_sphere(M, d, rng) if X0 is None else np.asarray(X0, dtype=np.float64).copy()
    if M < 2: return X / geomd.norms(X)[:, None]
    ch = np.sqrt(2 - 2 * cos_target)
    for f in stages:
        res = minimize(_sph, X.ravel(), args=(M, d, (f * ch) ** 2), jac=True, method='L-BFGS-B',
                       options=dict(maxiter=maxiter, gtol=1e-14, ftol=1e-16, maxcor=20))
        X = res.x.reshape(M, d); X /= geomd.norms(X)[:, None]
    return X


_COHN = None


def cohn_cos(d, M):
    """Cohn-table max cosine for (d, M) (float), or an estimate for M outside the table (never a bar)."""
    global _COHN
    if _COHN is None: _COHN = geomd.cohn_table()
    if (d, M) in _COHN: return float(_COHN[(d, M)])
    if M <= d + 1: return -1.0 / d
    ks = sorted(m for (dd, m) in _COHN if dd == d)
    m0 = min(ks, key=lambda m: abs(m - M)); c0 = float(_COHN[(d, m0)])
    th = np.arccos(c0) * (m0 / M) ** (1.0 / (d - 1))                    # cap-area scaling
    return float(np.cos(th))


# ---------------------------------------------------------------- structured codes
def cell24():
    """The 24-cell: the 24 vectors (+-1, +-1, 0, 0) / sqrt 2 and permutations (D4 roots), minimal angle 60 degrees."""
    V = []
    for i in range(4):
        for j in range(i + 1, 4):
            for a in (1, -1):
                for b in (1, -1):
                    v = [0.0] * 4; v[i] = a; v[j] = b; V.append(v)
    V = np.array(V); return V / geomd.norms(V)[:, None]


def e6_roots():
    """The 72 roots of E6 in R^6 (the E6 kissing configuration, minimal angle 60 degrees): the E8 roots orthogonal to
    the A2 spanned by e1 - e2 and e2 - e3 (i.e. v1 = v2 = v3), written in an orthonormal basis of that 6-space."""
    R = []
    for i in range(8):
        for j in range(i + 1, 8):
            for a in (1, -1):
                for b in (1, -1):
                    v = [0.0] * 8; v[i] = a; v[j] = b; R.append(v)
    for m in range(256):
        v = [0.5 if (m >> k) & 1 == 0 else -0.5 for k in range(8)]
        if sum(x < 0 for x in v) % 2 == 0: R.append(v)
    R = np.array(R); R = R[(np.abs(R[:, 0] - R[:, 1]) < 1e-12) & (np.abs(R[:, 1] - R[:, 2]) < 1e-12)]
    assert len(R) == 72, len(R)
    A = np.array([[1, -1, 0, 0, 0, 0, 0, 0], [0, 1, -1, 0, 0, 0, 0, 0]], float)
    _, _, Vt = np.linalg.svd(A); B = Vt[2:]                     # orthonormal basis of the orthogonal complement (6 x 8)
    X = R @ B.T; assert np.allclose(geomd.norms(X), np.sqrt(2)), 'E6 roots not of norm sqrt 2'
    return X / np.sqrt(2)


def place_code(U, centre=False):
    """(centres, r) for unit vectors U on the shell |c| = 1 - r (+ one ball at the origin if centre)."""
    r = geomd.code_radius(geomd.max_cos(U))
    if centre: r = min(r, 1.0 / 3.0)
    c = U * (1 - r)
    if centre: c = np.vstack([c, np.zeros((1, U.shape[1]))])
    return c, r


def code_seed(d, M, centre=False):
    U = geomd.load_code(d, M)
    if U is None: return None
    return place_code(U, centre)


# ---------------------------------------------------------------- shell + core, cold
def shell_core(d, N, k, r_ref, rng, inner='code', stages=(0.8, 0.9, 0.96, 1.0, 1.004)):
    """Code shell of N - k balls + k inner balls, inflated at the reference radius r_ref."""
    M = N - k; rho = r_ref
    X = sphere_code(M, d, rng, cohn_cos(d, M))
    c = X * (1 - rho)
    if k > 0:
        R1 = 1 - rho - 1.8 * rho
        if inner == 'code' and R1 > rho:
            th = 2 * np.arcsin(min(1.0, rho / R1))                         # angle two touching balls subtend at R1
            cap = max(m for m in range(1, k + 1) if m == 1 or cohn_cos(d, m) <= np.cos(th) + 1e-12)
            k1 = min(k, cap)
            core = sphere_code(k1, d, rng, cohn_cos(d, k1)) * R1 if k1 >= 2 else R1 * geomd.random_sphere(1, d, rng)
            rest = k - k1
            if rest == 1:
                core = np.vstack([core, np.zeros((1, d))])
            elif rest > 1:
                core = np.vstack([core, geomd.random_ball(rest, d, max(R1 - 1.8 * rho, rho), rng)])
        else:
            core = geomd.random_ball(k, d, max(1 - 2.5 * rho, 0.05), rng)
        c = np.vstack([c, core])
    return inflate(c, r_ref, stages)


def cold(d, N, r_ref, rng, stages=(0.6, 0.8, 0.9, 0.96, 1.0, 1.004)):
    return inflate(geomd.random_ball(N, d, 1 - r_ref, rng), r_ref, stages)


# ---------------------------------------------------------------- holes, neighbour seeds
def clearance(p, cs, r):
    """Largest radius a ball at p could have, given balls of radius r at cs and the unit-ball wall."""
    dd, _ = cKDTree(cs).query(p, k=1)
    return np.minimum(dd - r, 1.0 - geomd.norms(p))


def find_holes(cs, r, k, rng, m=20000):
    """k greedy largest holes for centres cs (radius r): best of m random samples, refined by a pattern search."""
    d = cs.shape[1]
    p = geomd.random_ball(m, d, 1 - 0.5 * r, rng); rho = clearance(p, cs, r); holes = []
    D = np.vstack([np.eye(d), -np.eye(d), geomd.random_sphere(2 * d, d, rng)])
    for _ in range(k):
        ref = np.vstack([cs] + ([np.array(holes)] if holes else []))
        top = np.argsort(-rho)[:6]; bestp, bestv = None, -np.inf
        for t in top:
            x, v, s = p[t].copy(), rho[t], 0.5 * r
            while s > 1e-4 * r:
                cand = x + s * D; cv = clearance(cand, ref, r); j = int(np.argmax(cv))
                if cv[j] > v: x, v = cand[j], cv[j]
                else: s /= 2
            if v > bestv: bestp, bestv = x, v
        holes.append(bestp)
        rho = np.minimum(rho, geomd.norms(p - bestp) - r)
    return np.array(holes)


def insert(c, r, k, rng):
    h = find_holes(c, r, k, rng)
    return np.vstack([c, h + rng.uniform(-1e-3 * r, 1e-3 * r, h.shape)])


def delete(c, r, k, rng):
    key = geomd.contacts(c, r) + rng.random(len(c)) * 0.5
    return c[np.sort(np.argsort(key)[k:])]


def packomania(d, N):
    return geomd.load_pub(d, N)
