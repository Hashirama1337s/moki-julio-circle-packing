"""Equal circles in the unit square [-1/2, 1/2]^2. Two families.

stretched: affine triangular rows; in-row pitch and row gap scale separately
so the outer rows meet both pairs of walls. k rows of n, or alternating
n / n-1. Both axes. Largest lattice r with count >= N, then delete.

mixed: hex rows of length n or n-1 with the short rows on the top and bottom
walls, fitted by bisection; plus a phase sweep of affine hex lattices
(row lengths come out n or n-1, shorter parity can sit on both walls).
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def admitted_r(pts):
    pts = np.asarray(pts, np.float64)
    wall = float(np.min(0.5 - np.max(np.abs(pts), axis=1)))
    if len(pts) < 2:
        return wall
    dist, _ = cKDTree(pts).query(pts, k=2)
    return float(min(wall, 0.5 * float(np.min(dist[:, 1]))))


def delete_fewest(pts, r, N):
    """Drop N'-N circles, fewest contacts first; ties take the lowest index."""
    pts = np.asarray(pts, np.float64)
    m = int(pts.shape[0])
    if m <= N:
        return np.ascontiguousarray(pts, dtype=np.float64)
    deg = np.zeros(m, np.int32)
    adj = [[] for _ in range(m)]
    for a, b in cKDTree(pts).query_pairs(2.0 * r * (1.0 + 1e-7)):
        adj[a].append(b)
        adj[b].append(a)
        deg[a] += 1
        deg[b] += 1
    deg += (np.abs(np.abs(pts) - (0.5 - r)) <= 1e-8).sum(1).astype(np.int32)
    alive = np.ones(m, dtype=bool)
    for _ in range(m - N):
        i = int(np.argmin(np.where(alive, deg, 10**9)))
        alive[i] = False
        for j in adj[i]:
            if alive[j]:
                deg[j] -= 1
    return np.ascontiguousarray(pts[alive], dtype=np.float64)


def _refit(pts):
    """Uniform scale about the bbox centre. Keep the original if it does not gain."""
    base = admitted_r(pts)
    q = np.array(pts, np.float64, copy=True)
    q -= 0.5 * (q.min(0) + q.max(0))
    half = float(np.max(np.abs(q)))
    dist, _ = cKDTree(q).query(q, k=2)
    dmin = float(np.min(dist[:, 1]))
    if half <= 0.0 or dmin <= 0.0:
        return np.asarray(pts, np.float64), base
    q *= 0.5 / (half + 0.5 * dmin)
    rr = admitted_r(q)
    return (q, rr) if rr > base else (np.asarray(pts, np.float64), base)


def _r_lat(a, b):
    # Wall-tight triangular lattice. Constraints: in-row dx >= 2r,
    # diagonal hypot(dx/2, dy) >= 2r, same-phase gap 2*dy >= 2r.
    g = float(np.sqrt(0.25 / (a * a) + 1.0 / (b * b)))
    return min(g / (2.0 * (g + 1.0)), 1.0 / (2.0 * (a + 1.0)), 1.0 / (b + 2.0))


def _count_rows(k, n, alt):
    if not alt:
        return k * n
    return ((k + 1) // 2) * n + (k // 2) * (n - 1)


def _build_stretch(k, n, alt, r):
    a = (n - 1.0) if alt else (n - 0.5)
    W = 1.0 - 2.0 * r
    dx, dy = W / a, W / float(k - 1)
    rows = []
    for j in range(k):
        y = -0.5 + r + j * dy
        odd = j & 1
        if alt and odd:
            xs = -0.5 + r + 0.5 * dx + np.arange(n - 1) * dx
        else:
            off = 0.5 * dx if (odd and not alt) else 0.0
            xs = -0.5 + r + off + np.arange(n) * dx
        rows.append(np.column_stack((xs, np.full(xs.shape[0], y, np.float64))))
    return np.vstack(rows).astype(np.float64)


def _finish(pts, r, N):
    if len(pts) > N:
        pts = delete_fewest(pts, r, N)
    pts, rr = _refit(pts)
    flip = np.ascontiguousarray(pts[:, ::-1])
    flip, rr2 = _refit(flip)
    if rr2 > rr + 1e-15:
        return flip, rr2
    return np.ascontiguousarray(pts, dtype=np.float64), rr


def _stretched(N):
    k_hi = int(np.sqrt(N) * 3.0) + 8
    best = None
    pool = []
    for k in range(2, k_hi + 1):
        b = float(k - 1)
        for alt in (True, False):
            if alt:
                n = max(2, (N + k // 2 + k - 1) // k)
            else:
                n = max(1, (N + k - 1) // k)
            Np = _count_rows(k, n, alt)
            if Np < N:
                continue
            a = (n - 1.0) if alt else (n - 0.5)
            r = _r_lat(a, b)
            key = (r, -Np, -k, -n, 1 if alt else 0)
            pool.append((r, a / b))
            if best is None or key > best[0]:
                best = (key, k, n, alt, r)
    _, k, n, alt, r = best
    tag = "alt" if alt else "eq"
    # Both orientations: coordinate swap is tried inside _finish.
    pts = _build_stretch(k, n, alt, r)
    # Stretched spec: delete, do not rescale. Score the true radius.
    if len(pts) > N:
        pts = delete_fewest(pts, r, N)
    pts2 = np.ascontiguousarray(pts[:, ::-1])
    r1, r2 = admitted_r(pts), admitted_r(pts2)
    if r2 > r1 + 1e-15:
        pts, r1 = pts2, r2
    gammas = [0.50, 0.65, float(np.sqrt(3.0) / 2.0), 1.0, 1.25, 1.60, 2.0]
    for _, g in sorted(pool, reverse=True):
        if all(abs(g - h) > 1e-3 * max(abs(g), 1e-6) for h in gammas):
            gammas.append(float(g))
        if len(gammas) >= 10:
            break
    return f"stretched_{tag}_k{k}_n{n}", np.ascontiguousarray(pts, np.float64), float(r1), gammas


def _span_a(counts):
    a = 0.0
    for i, c in enumerate(counts):
        a = max(a, (c - 0.5) if (i & 1) else (c - 1.0))
    return a


def _make_counts(k, n, kind):
    if kind == 0:  # alternating, long rows on even indices
        c = [n if (i % 2 == 0) else n - 1 for i in range(k)]
    elif kind == 1:  # long rows on odd indices (shifted); ends short when k is even
        c = [n - 1 if (i % 2 == 0) else n for i in range(k)]
    elif kind == 2:  # alternating, then force both wall rows short
        c = [n if (i % 2 == 0) else n - 1 for i in range(k)]
        c[0] = n - 1
        c[-1] = n - 1
    else:  # short rows only on the two walls; interior stays long
        c = [n] * k
        c[0] = n - 1
        c[-1] = n - 1
    return c


def _ok_ab(a, b, r):
    if not (0.0 < r < 0.5) or a <= 0.0 or b <= 0.0:
        return False
    W = 1.0 - 2.0 * r
    dx, dy = W / a, W / b
    if dx + 1e-12 < 2.0 * r or dy + 1e-12 < r:
        return False
    return (0.5 * dx) ** 2 + dy * dy + 1e-12 >= 4.0 * r * r


def _bisect(ok, lo=1e-6, hi=0.49):
    if not ok(lo):
        return None
    if ok(hi):
        return hi
    for _ in range(52):
        mid = 0.5 * (lo + hi)
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo


def _place_counts(counts, r):
    a = _span_a(counts)
    k = len(counts)
    W = 1.0 - 2.0 * r
    dx = W / a
    dy = W / float(k - 1)
    rows = []
    for i, c in enumerate(counts):
        y = -0.5 + r + i * dy
        off = 0.5 * dx if (i & 1) else 0.0
        xs = -0.5 + r + off + np.arange(c, dtype=np.float64) * dx
        rows.append(np.column_stack((xs, np.full(c, y))))
    return np.vstack(rows).astype(np.float64)


def _bands(N, k_hi):
    found = []
    for k in range(2, k_hi + 1):
        b = float(k - 1)
        for kind in range(4):
            n = max(2, (N + k - 1) // k - 1)
            counts = _make_counts(k, n, kind)
            while sum(counts) < N and n < N:
                n += 1
                counts = _make_counts(k, n, kind)
            if sum(counts) < N:
                continue
            a = _span_a(counts)
            if a <= 0.0:
                continue
            closed = _r_lat(a, b)
            bis = _bisect(lambda t, a=a, b=b: _ok_ab(a, b, t))
            if bis is None:
                continue
            r = min(closed, bis)
            found.append((r, -(sum(counts) - N), -k, -kind, k, n, kind))
    found.sort(reverse=True)
    best = None
    for r, _, _, _, k, n, kind in found[:8]:
        pts = _place_counts(_make_counts(k, n, kind), r)
        if len(pts) < N:
            continue
        pts, rr = _finish(pts, r, N)
        if best is None or rr > best[0]:
            best = (rr, pts, f"mixed_band_k{k}_n{n}_c{kind}")
    return best


_FR = (0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875)


def _metric(r, gamma):
    fac = min(1.0, float(np.hypot(0.5, gamma)), 2.0 * gamma)
    dx = 2.0 * r / fac
    return dx, gamma * dx


def _row_index(r, dx, dy, fx, fy):
    lo, hi = -0.5 + r, 0.5 - r
    ox, oy = fx * dx, fy * dy
    j0 = int(np.ceil((lo - oy) / dy - 1e-8))
    j1 = int(np.floor((hi - oy) / dy + 1e-8))
    return lo, hi, ox, oy, j0, j1


def _lat_count(r, dx, dy, fx, fy):
    lo, hi, ox, oy, j0, j1 = _row_index(r, dx, dy, fx, fy)
    if j1 < j0:
        return 0
    js = np.arange(j0, j1 + 1, dtype=np.int64)
    off = ox + (js & 1) * (0.5 * dx)
    i0 = np.ceil((lo - off) / dx - 1e-8)
    i1 = np.floor((hi - off) / dx + 1e-8)
    return int(np.maximum(0.0, i1 - i0 + 1.0).sum())


def _lat_points(r, gamma, fx, fy):
    dx, dy = _metric(r, gamma)
    lo, hi, ox, oy, j0, j1 = _row_index(r, dx, dy, fx, fy)
    chunks = []
    for j in range(j0, j1 + 1):
        off = ox + (0.5 * dx if (j & 1) else 0.0)
        i0 = int(np.ceil((lo - off) / dx - 1e-8))
        i1 = int(np.floor((hi - off) / dx + 1e-8))
        if i1 < i0:
            continue
        xs = off + np.arange(i0, i1 + 1, dtype=np.float64) * dx
        chunks.append(np.column_stack((xs, np.full(xs.shape[0], oy + j * dy))))
    if not chunks:
        return np.zeros((0, 2), np.float64)
    pts = np.vstack(chunks).astype(np.float64)
    keep = np.max(np.abs(pts), axis=1) <= 0.5 - r + 1e-8
    return np.ascontiguousarray(pts[keep])


def _phase(N, gammas):
    scale = 1.0 / np.sqrt(N)
    lo, hi = 0.12 * scale, 0.57 * scale

    def ok(r):
        for gamma in gammas:
            dx, dy = _metric(r, gamma)
            for fy in _FR:
                for fx in _FR:
                    if _lat_count(r, dx, dy, fx, fy) >= N:
                        return True
        return False

    if not ok(lo):
        lo = 0.02 * scale
        if not ok(lo):
            return None
    if not ok(hi):
        for _ in range(46):
            mid = 0.5 * (lo + hi)
            if ok(mid):
                lo = mid
            else:
                hi = mid
    r = lo
    chosen = None
    for gi, gamma in enumerate(gammas):
        dx, dy = _metric(r, gamma)
        for fy in _FR:
            for fx in _FR:
                c = _lat_count(r, dx, dy, fx, fy)
                if c >= N and (chosen is None or (c - N, gi, fy, fx) < chosen[0]):
                    chosen = ((c - N, gi, fy, fx), gamma, fx, fy)
    if chosen is None:
        return None
    _, gamma, fx, fy = chosen
    pts = _lat_points(r, gamma, fx, fy)
    tries = 0
    while len(pts) < N and tries < 30:
        r *= 1.0 - 1e-4
        pts = _lat_points(r, gamma, fx, fy)
        tries += 1
    if len(pts) < N:
        return None
    pts, rr = _finish(pts, r, N)
    name = f"mixed_phase_g{gamma:.4f}_fx{fx:.3f}_fy{fy:.3f}"
    return rr, pts, name


def _mixed(N, gammas):
    k_hi = int(np.sqrt(N) * 3.0) + 8
    best = _bands(N, k_hi)
    phased = _phase(N, gammas)
    if phased is not None and (best is None or phased[0] > best[0]):
        best = phased
    if best is None:
        raise RuntimeError("mixed found nothing")
    rr, pts, name = best
    return name, np.ascontiguousarray(pts, np.float64), float(rr)


def candidates(N):
    """Yield (name, centres (N, 2) float64, admitted radius) for each family."""
    N = int(N)
    if N < 2:
        raise ValueError("N >= 2")
    name, pts, r, gammas = _stretched(N)
    if pts.shape != (N, 2):
        raise RuntimeError(name)
    yield name, pts, float(r)
    name, pts, r = _mixed(N, gammas)
    if pts.shape != (N, 2):
        raise RuntimeError(name)
    yield name, pts, float(r)


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    for name, pts, r in candidates(n):
        print(f"{name}  N={len(pts)}  r={r:.12e}  dens={n * np.pi * r * r:.6f}")
