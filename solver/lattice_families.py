# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
# recon2_families.py -- closed-form row families + monotone envelope vs Packomania crc tables (and our certified records).
# Standard float check only (margins found are >= 1e-6 relative); run the two exact checkers before any claim.
# py -3.11 vision/recon2_families.py
import math, re, itertools
import numpy as np
from scipy.spatial import cKDTree
ROOT = "./"   # run from the folder holding shelves/ and data/

def load(p):
    d = {}
    for l in open(p):
        n, r = l.split(); d[int(n)] = float(r)
    return d

def recs(name):  # our certified records (both checkers YES)
    d = {}
    try:
        for l in open(ROOT + "RECORDS_%s.md" % name, encoding="utf8"):
            m = re.match(r"\|\s*(\d+)\s*\|\s*([0-9.]+)\s*\|\s*([0-9.]+)\s*\|.*\|\s*(YES|NO)\s*\|", l)
            if m and m.group(4) == "YES": d[int(m.group(1))] = float(m.group(2))
    except FileNotFoundError:
        pass
    return d

def bisect_max(feas, lo, hi):
    if not feas(lo): return None
    if feas(hi): return hi
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if feas(mid) else (lo, mid)
    return lo

def family(W, H, k, n, kind):
    """k rows parallel to side W (length W), stacked across H. Returns (N, r) or None.
    alt  : rows of n, n-1, n, ... (short rows offset p/2), pitch p=(W-2r)/(n-1), spacing v=(H-2r)/(k-1)
    shift: k rows of n, odd rows offset by delta=sqrt(4r^2-v^2) (pitch 2r) or, if delta>r, stretched pitch p=2*delta
    sq   : k x n square grid"""
    if kind == "sq":
        return k * n, min(W / (2 * n), H / (2 * k))
    if k == 1:
        return n, min(W / (2 * n), H / 2)
    if kind == "alt":
        if n < 2: return None
        N = (k + 1) // 2 * n + k // 2 * (n - 1)
        def feas(r):
            p = (W - 2 * r) / (n - 1); v = (H - 2 * r) / (k - 1)
            if p < 2 * r or v < 0 or (k >= 3 and v < r): return False
            return (p / 2) ** 2 + v * v >= 4 * r * r
        r = bisect_max(feas, 1e-9, min(H / 2, W / (2 * n)))
        return (N, r) if r else None
    if kind == "shift":
        N = k * n
        def feas(r):
            v = (H - 2 * r) / (k - 1)
            if v < 0 or (k >= 3 and v < r): return False
            d2 = 4 * r * r - v * v
            if d2 > r * r:
                pp = 2 * math.sqrt(d2); return (n - 0.5) * pp + 2 * r <= W
            return 2 * r * n + math.sqrt(max(0.0, d2)) <= W
        r = bisect_max(feas, 1e-9, min(H / 2, W / 2))
        return (N, r) if r else None

def build(W, H, k, n, kind, r):
    pts = []
    if kind == "sq":
        return [(r + i * (W - 2 * r) / max(1, n - 1), r + j * (H - 2 * r) / max(1, k - 1)) for i in range(n) for j in range(k)]
    v = (H - 2 * r) / (k - 1) if k > 1 else 0
    if kind == "alt":
        p = (W - 2 * r) / (n - 1)
        for j in range(k):
            xs = [r + i * p for i in range(n)] if j % 2 == 0 else [r + p / 2 + i * p for i in range(n - 1)]
            pts += [(x, r + j * v) for x in xs]
    else:
        d2 = 4 * r * r - v * v
        pp, d = (2 * math.sqrt(d2), math.sqrt(d2)) if d2 > r * r else (2 * r, math.sqrt(max(0, d2)))
        for j in range(k):
            pts += [(r + (d if j % 2 else 0) + i * pp, r + j * v) for i in range(n)]
    return pts

def checked_radius(P, h):
    P = np.array(P); tr = cKDTree(P); dd, _ = tr.query(P, k=2)
    wall = min(P[:, 0].min(), (1 - P[:, 0]).min(), P[:, 1].min(), (h - P[:, 1]).min())
    return min(wall, dd[:, 1].min() / 2)

for hh in range(1, 9):
    h = hh / 10; name = "crc_%d00" % hh
    t = load(ROOT + "shelves/%s/radius.txt" % name); R = recs(name)
    best = {N: max(t[N], R.get(N, 0.0)) for N in t}
    # 1) monotone envelope: a larger-N packing with a larger radius, minus circles, beats entry N
    env, cur, src = {}, (0.0, None), {}
    for N in sorted(best, reverse=True):
        if best[N] > cur[0]: cur = (best[N], N)
        env[N] = cur
    for N in sorted(t):
        if env[N][0] > best[N] * (1 + 1e-10):
            print("%s N=%d  MONOTONE: delete %d circle(s) from N=%d -> r=%.13f (+%.2e)" % (name, N, env[N][1] - N, env[N][1], env[N][0], env[N][0] / best[N] - 1))
    # 2) closed-form families (both orientations), with deletion
    fams = {}
    for (W, H, orient) in [(1.0, h, "rows"), (h, 1.0, "cols")]:
        for k in range(1, 60):
            for n in range(1, 480):
                for kind in ("alt", "shift", "sq"):
                    res = family(W, H, k, n, kind)
                    if not res or res[0] > max(t) + 60: continue
                    N, r = res
                    if N not in fams or r > fams[N][0]: fams[N] = (r, (W, H, k, n, kind, orient))
    cur = (0.0, None)
    for N in sorted(fams, reverse=True):
        if fams[N][0] > cur[0]: cur = fams[N]
        fams[N] = cur
    for N in sorted(t):
        if N in fams and fams[N][0] > env[N][0] * (1 + 1e-10):
            r, (W, H, k, n, kind, orient) = fams[N]
            P = build(W, H, k, n, kind, r)
            if orient == "cols": P = [(y, x) for x, y in P]
            rc = checked_radius(P, h)
            print("%s N=%d  FAMILY %s k=%d n=%d %s (%d circles) r=%.13f checked=%.13f  (+%.2e vs best known)" % (name, N, orient, k, n, kind, len(P), r, rc, r / env[N][0] - 1))

# monotone check for the triangle and quadrant too (no row families there)
for name, path, rn in [("crt", ROOT + "data/radius.txt", "crt"), ("ccq", ROOT + "shelves/ccq/radius.txt", "ccq")]:
    t = load(path); R = recs(rn)
    if name == "crt" and not R: R = recs("")
    best = {N: max(t[N], R.get(N, 0.0)) for N in t}
    cur = (0.0, None)
    for N in sorted(best, reverse=True):
        if cur[0] > best[N] * (1 + 1e-10):
            print("%s N=%d  MONOTONE: delete %d circle(s) from N=%d -> r=%.13f (+%.2e)" % (name, N, cur[1] - N, cur[1], cur[0], cur[0] / best[N] - 1))
        if best[N] > cur[0]: cur = (best[N], N)
