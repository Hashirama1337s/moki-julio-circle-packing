# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""ssp helpers (2026-09-26; Moki&Julio): equal spheres in a sphere (Packomania ssp, d = 3).

FRAME (verified by recon.py on the files and exactly by checkers/certify_ball.py on a few of them): container = the unit
sphere centred at the origin; the printed number is the radius r of the N equal spheres; a sphere centred at c is inside
iff |c| <= 1 - r; spheres do not overlap iff |c_i - c_j| >= 2 r. This is exactly the hsp frame at d = 3, so
../hsp/geomd.py (rmin, contacts, ...) and ../hsp/slpd.py (polish) are used unchanged (imported from ../hsp/ via sys.path).

Data:
  data/refs/ssp/packomania_ssp_2026-09-26.html   the table (Last update 13-Jul-2026; identical to the 2026-09-25 copy)
  data/big/ssp/coords/ssp<N>.txt                 Packomania's coordinates, N = 1..1000 (from txt/ssp_coords.tar.gz):
                                                 line 1 = r (12 decimals), '#' comment lines, then 'idx x y z r'
  data/big/ssp/ssp_<radius|author|...>_2026-09-26.txt   the per-column text tables linked from the page header
"""
import os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
HSP = os.path.join(ROOT, 'hsp')
if HSP not in sys.path: sys.path.insert(0, HSP)
REFS = os.path.join(ROOT, 'data', 'refs', 'ssp')
BIG = os.path.join(ROOT, 'data', 'big', 'ssp')
COORDS = os.path.join(BIG, 'coords')
PAGE = os.path.join(REFS, 'packomania_ssp_2026-09-26.html')
OUT = os.path.join(HERE, 'out')
NMAX = 1000
REFNAME = {'[1]': 'Fejes Toth 1943', '[12]': 'Zhou, Ren, He, Liu, Li (PESS, C&OR 164 (2024) 106522)',
           '[13]': 'Lai et al., private communication May 2024', '[31]': 'Specht, program ssp', '': '(no reference: trivial)'}
COLOR = {'#6495ed': '[31] Specht', '#3ff5a6': '[13] Lai et al.', '#c0e75f': '[12] Zhou et al.', '#ffd700': 'older / other'}


def _cells(line):
    return [re.sub(r'<[^>]+>', '', x).replace('&nbsp;', ' ').strip() for x in re.findall(r'<td[^>]*>(.*?)</td>', line)]


def table(all_rows=False):
    """{N: dict(r, distance, ratio, density, contacts, loose, boundary, core, ref, color, bold)} from the saved page.
    N <= 1000 rows have 10 cells (N, r, dist, ratio, density, contacts, loose, boundary, core, ref); the larger-N rows
    (the page's showcase packings) have 9 (no core column). all_rows=True keeps those too."""
    T = {}
    for line in open(PAGE, encoding='latin-1'):
        if not line.startswith('<tr><td'): continue
        c = _cells(line)
        if len(c) < 9 or not c[0].isdigit(): continue
        n = int(c[0]); col = re.search(r'bgcolor="(#[0-9a-fA-F]+)"', line)
        iv = lambda s: int(s) if s.strip().isdigit() else 0
        if len(c) >= 10:
            row = dict(r=c[1], distance=c[2], ratio=c[3], density=c[4], contacts=iv(c[5]), loose=iv(c[6]), boundary=iv(c[7]),
                       core=iv(c[8]), ref=c[9].strip())
        else:
            row = dict(r=c[1], distance=c[2], ratio=c[3], density=c[4], contacts=iv(c[5]), loose=iv(c[6]), boundary=iv(c[7]),
                       core=None, ref=c[8].strip())
        rad_html = re.findall(r'<td[^>]*>(.*?)</td>', line)[1]
        row['bold'] = bool(re.search(r'<(strong|b)>', rad_html))
        row['color'] = col.group(1).lower() if col else None
        if n <= NMAX or all_rows: T[n] = row
    return T


def pub_path(n):
    return os.path.join(COORDS, f'ssp{n}.txt')


def load_pub_strings(n):
    """(radius string of line 1, [(x, y, z) strings], [per-row radius strings], [indices of BLANK rows]) from
    Packomania's file. Every ssp file N >= 2 has row N blank (missing.py): the coordinate list is N - 1 long."""
    L = [l.split() for l in open(pub_path(n)) if l.strip()]
    r_s = L[0][0]; rows = [t for t in L[1:] if not t[0].startswith('#')]
    full = [t for t in rows if len(t) == 5]; blank = [int(t[0]) for t in rows if len(t) == 1]
    if len(rows) != n or len(full) + len(blank) != n:
        raise ValueError(f'ssp{n}.txt: {len(rows)} rows, widths {sorted({len(t) for t in rows})}')
    return r_s, [tuple(t[1:4]) for t in full], [t[4] for t in full], blank


def full_path(n):
    return os.path.join(BIG, 'full', f'ssp{n}.npy')


def load_full(n):
    """Packomania's N - 1 centres + the reconstructed N-th centre (missing.py), (n, 3) float, or None."""
    p = full_path(n)
    if not os.path.exists(p): return None
    c = np.load(p); return c if c.shape == (n, 3) else None
