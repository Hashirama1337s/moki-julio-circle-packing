# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""chx = equal circles in a REGULAR HEXAGON (Packomania, page last updated 18-Dec-2020).  Shared loaders.  Moki&Julio 2026-09-26.

FRAME (checked exactly on every published file by the frame check): Packomania's chx coordinates are ALREADY in the k-gon
frame of certify_kgon.py with k = 6: circumradius 1, centred at the origin, flat side at the bottom (and top), vertices at
(+-1, 0) and (+-1/2, +-sqrt(3)/2); side normals at 270 + 60 j degrees; apothem a = sqrt(3)/2.  No rotation is applied.

Page (data/refs/chx/packomania_chx_2026-09-26.html): reference [1] = E. Specht, program chx (2020);
                                                     reference [2] = P. Amore et al., Phys. Fluids 35, 027130 (2023).
Coordinates: data/big/chx/chx<N>.txt (from txt/chx_coords.tar.gz, downloaded 2026-09-26); line = "i x y"; the N = 1 file
prints no coordinates (the centred circle).
The kgon float geometry (../kgon/geom.py: Frame + SLP polish) is imported by file path as module 'kgon_geom' (not edited).
"""
import os, re, sys, html, importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
K = 6
PAGE = os.path.join(ROOT, 'data', 'refs', 'chx', 'packomania_chx_2026-09-26.html')
BIG = os.path.join(ROOT, 'data', 'big', 'chx')
REFS = os.path.join(ROOT, 'data', 'refs', 'chx')


def _load(name, path):
    if name in sys.modules: return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m; spec.loader.exec_module(m); return m


kgeom = _load('kgon_geom', os.path.join(ROOT, 'kgon', 'geom.py'))
ck = _load('certify_kgon', os.path.join(os.path.dirname(ROOT), 'checkers', 'certify_kgon.py'))


def frame():
    return kgeom.frame(K)


def page_table():
    """{N: dict(radius, ratio, density, contacts, loose, boundary, symmetry, refs)} -- every column as printed."""
    s = open(PAGE, encoding='latin-1').read()
    out = {}
    for row in re.findall(r'<tr><td[^>]*><a href="chx\d+\.html"\s+name="chx\d+">.*?</tr>', s, re.S):
        n = int(re.search(r'name="chx(\d+)"', row).group(1))
        tds = re.findall(r'<td[^>]*>(.*?)</td>', row, re.S)
        txt = [html.unescape(re.sub(r'<[^>]+>', '', t)).strip() for t in tds]
        # txt: N, radius, ratio, density, contacts, loose, boundary, symmetry, reference
        out[n] = dict(radius=txt[1], ratio=txt[2], density=txt[3], contacts=txt[4], loose=txt[5], boundary=txt[6],
                      symmetry=txt[7], refs=re.findall(r'\[(\d+)\]', txt[8]))
    return out


def radius_txt():
    out = {}
    for line in open(os.path.join(REFS, 'chx_radius_2026-09-26.txt')):
        p = line.split()
        if len(p) == 2: out[int(p[0])] = p[1]
    return out


def coord_strings(n):
    pts = []
    for line in open(os.path.join(BIG, f'chx{n}.txt')):
        p = line.split()
        if not p: continue
        if len(p) == 3: pts.append((p[1], p[2]))
        elif len(p) == 1 and n == 1: pts.append(('0', '0'))
        else: raise ValueError(f'chx{n}.txt: bad line {line!r}')
    assert len(pts) == n, (n, len(pts))
    return pts


def load_coords(n):
    return np.array([(float(x), float(y)) for x, y in coord_strings(n)], dtype=np.float64)


def hex_family(nmax=1300):
    """Special lattice families on the page: H_k = 3k(k-1)+1, H_k - 1, H_k - k, H_k - 2k + 1 (k >= 2)."""
    s = {}
    for k in range(2, 40):
        H = 3 * k * (k - 1) + 1
        for v, name in ((H, 'H_k'), (H - 1, 'H_k-1'), (H - k, 'H_k-k'), (H - 2 * k + 1, 'H_k-2k+1')):
            if 1 <= v <= nmax: s.setdefault(v, f'{name} (k={k})')
    return s
