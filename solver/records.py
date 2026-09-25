# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Packomania 'circles in an isosceles right triangle' (crt) records: legs of length 1, right angle at the origin.
Circle form: N circles of radius r, centres (x, y) with x >= r, y >= r, x + y <= 1 - sqrt(2) r, pairwise distance >= 2r.
Point form: N points in the unit triangle T = {x >= 0, y >= 0, x + y <= 1}, maximise the minimum pairwise distance d.
Map: inset the circle triangle by r on every side -> legs 1 - (2 + sqrt2) r, so d = 2r / (1 - (2 + sqrt2) r), r = d / (2 + (2 + sqrt2) d).
"""
import os, mpmath as mp
mp.mp.dps = 50
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")   # Packomania crt data
S2 = mp.sqrt(2)

def r_to_d(r): return 2 * r / (1 - (2 + S2) * r)
def d_to_r(d): return d / (2 + (2 + S2) * d)

def table(name):
    """{N: mpf} from radius.txt / distance.txt."""
    out = {}
    for line in open(os.path.join(D, name)):
        p = line.split()
        if len(p) == 2: out[int(p[0])] = mp.mpf(p[1])
    return out

def coords(n):
    """Circle centres of the published record for n, as mpf pairs (string-exact)."""
    pts = []
    for line in open(os.path.join(D, "coords", f"crt{n}.txt")):
        p = line.split()
        if len(p) == 3: pts.append((mp.mpf(p[1]), mp.mpf(p[2])))
    return pts

def circle_radius(pts):
    """Largest r such that circles of radius r at these centres fit (min of wall clearances and half pair distances)."""
    walls = min(min(x, y, (1 - x - y) / S2) for x, y in pts)
    pair = min(mp.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) for i, a in enumerate(pts) for b in pts[i + 1:]) / 2 if len(pts) > 1 else mp.inf
    return min(walls, pair)

def to_points(pts, r):
    """Circle centres -> points in the unit triangle (the inner triangle scaled to legs 1)."""
    L = 1 - (2 + S2) * r
    return [((x - r) / L, (y - r) / L) for x, y in pts]
