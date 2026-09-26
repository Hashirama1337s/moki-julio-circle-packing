# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""scu shared helpers (2026-09-25; Moki&Julio): the Packomania scu table (HTML page of 2026-09-25), the cell classes,
the coordinate files (download with curl, 1 s between requests), and the float loaders. Float only; no decisions here.

Frame: cube of side 1 centred at the origin; a sphere (c, r) is inside iff -1/2 + r <= x, y, z <= 1/2 - r.

Cell classes in N 1-1008 (from the recon of the scu table):
  derivative : the printed radius string equals that of N + 1 (a larger root with spheres removed)  -> 376 in 201-1008
  root       : not a derivative, but N - 1 is a derivative with the same radius                     ->  27 in 201-1008
  free       : everything else                                                                       -> 405 in 201-1008
               (381 Specht program [31], 15 Cantrell [24], 9 Pack'n'tile [26])
Stale coordinate files (the file reproduces an older radius; the HTML is higher): 83 94 95 96 109 110 507 619.
"""
import os, re, time, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # the folder holding data/refs/ (Packomania's scu page) and data/big/scu/ (its files, Lai 2023, roc-climate N = 203)
REFS = os.path.join(ROOT, 'data', 'refs')
BIG = os.path.join(ROOT, 'data', 'big', 'scu')
LAI = os.path.join(BIG, 'lai2023_pesc')
ROC203 = os.path.join(BIG, 'roc_climate_2026-09-21', 'scu203_new.txt')
PAGE = os.path.join(REFS, 'packomania_scu_2026-09-25.html')
STALE = (83, 94, 95, 96, 109, 110, 507, 619)
URL = 'https://www.packomania.com/scu/txt/scu{}.txt'


def table():
    """{N: (radius_string, ref_string, pink)} from the HTML page (all 1,134 rows)."""
    T = {}
    for line in open(PAGE, encoding='utf-8', errors='replace'):
        m = re.search(r'name="scu(\d+)"', line)
        if not m or not line.startswith('<tr>'): continue
        cells = [re.sub(r'<[^>]+>', '', c).strip() for c in re.findall(r'<td[^>]*>(.*?)</td>', line)]
        T[int(m.group(1))] = (cells[1], cells[-1], 'ffbbff' in line)
    assert len(T) == 1134, len(T)
    return T


def classes(T=None):
    T = T or table()
    der = {n for n in range(1, 1008) if T[n][0] == T[n + 1][0]}
    roots = {n for n in range(2, 1009) if T[n - 1][0] == T[n][0] and n not in der}
    free = {n for n in range(1, 1009) if n not in der and n not in roots}
    return der, roots, free


def radius_txt():
    """The stale sidecar radius.txt (Nov 2012): the radius each stale coordinate file actually reproduces."""
    return {int(l.split()[0]): l.split()[1] for l in open(os.path.join(REFS, 'scu_radius_2026-09-25.txt')) if l.strip()}


def pub_path(n):
    return os.path.join(BIG, f'scu{n}.txt')


def fetch(n, sleep=1.0):
    """Download scu<N>.txt with curl (certificate checks ON; Python urllib fails TLS verification on this host)."""
    p = pub_path(n)
    if os.path.exists(p) and os.path.getsize(p) > 0: return p, False
    tmp = p + '.part'
    r = subprocess.run(['curl', '-s', '-f', '-o', tmp, URL.format(n)], capture_output=True)
    time.sleep(sleep)
    if r.returncode != 0 or not os.path.exists(tmp):
        if os.path.exists(tmp): os.remove(tmp)
        raise RuntimeError(f'curl failed for N={n}: rc={r.returncode}')
    rows = [l.split() for l in open(tmp) if l.strip()]
    if len(rows) != n or any(len(t) != 4 for t in rows):
        os.remove(tmp); raise RuntimeError(f'bad file for N={n}: {len(rows)} rows')
    os.replace(tmp, p)
    return p, True


def load_pub_strings(n):
    """[(x, y, z) strings] from scu<N>.txt (lines 'idx x y z')."""
    pts = []
    for line in open(pub_path(n)):
        t = line.split()
        if not t: continue
        if len(t) != 4: raise ValueError(f'scu{n}.txt: bad line {line!r}')
        pts.append((t[1], t[2], t[3]))
    if len(pts) != n: raise ValueError(f'scu{n}.txt: {len(pts)} rows != {n}')
    return pts


def load_pub(n):
    return np.array([[float(v) for v in p] for p in load_pub_strings(n)], dtype=np.float64)


def load_lai(n):
    """Lai-Hao-Xiao-Glover 2023 PESC file: line 1 'N L', then N centres of UNIT spheres in the cube [-L/2, L/2]^3.
    Returns (centres rescaled to the unit cube, L as float, raw strings)."""
    L = [l.split() for l in open(os.path.join(LAI, f'{n}_CubeSol.txt')) if l.strip()]
    assert int(L[0][0]) == n and len(L) == n + 1, (n, L[0], len(L))
    side = float(L[0][1]); raw = [tuple(t[:3]) for t in L[1:]]
    return np.array([[float(v) for v in t] for t in raw]) / side, side, (L[0][1], raw)


def load_roc203():
    pts = [l.split()[1:4] for l in open(ROC203) if l.strip()]
    assert len(pts) == 203
    return np.array([[float(v) for v in t] for t in pts]), pts


def float_r(c):
    """Float min-radius from ALL pairs (O(N^2), chunked) and the six wall slacks. Recheck only; never a decision."""
    c = np.asarray(c, dtype=np.float64); n = len(c)
    g = float(np.concatenate([0.5 + c, 0.5 - c], 1).min())
    if n == 1: return g
    best = np.inf
    for s in range(0, n, 256):
        blk = c[s:s + 256]
        d2 = ((blk[:, None, :] - c[None, :, :]) ** 2).sum(-1)
        for k in range(len(blk)): d2[k, :s + k + 1] = np.inf
        best = min(best, float(d2.min()))
    return float(min(g, np.sqrt(best) / 2))


def lower_priority():
    """BELOW_NORMAL for this process (Windows)."""
    try:
        import ctypes
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
    except Exception:
        pass
