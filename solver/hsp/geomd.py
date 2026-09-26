# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""hsp geometry (2026-09-25; Moki&Julio): equal d-dimensional balls in the unit d-ball (Packomania hsp4 / hsp5 / hsp6).

FRAME: container = the unit ball centred at the origin in R^d. A ball (c, r) is inside iff |c| <= 1 - r (and r <= 1);
balls i != j do not overlap iff |c_i - c_j| >= 2 r. The float min-radius of a set of centres is
    rmin(c) = min( min_i (1 - |c_i|),  min_{i<j} |c_i - c_j| / 2 ).
Float only; every decision about a record is exact and separate (checkers/certify_ball.py = checker A, checkers/verify_exact_ball.py =
checker B). N <= 300 here, so every pair computation is dense (O(N^2 d), a few ms).

Data (written by recon 7, the recon of the hsp tables):
  data/refs/packomania_hsp{4,5,6}_2026-09-25.html  the tables (radius, contacts, loose, boundary, core, reference)
  data/big/hsp4/coords/hsp4-<N>.txt                 Packomania's coordinates (line 1 = r, 12 decimals; 'idx x1..x4')
  data/big/hsp5|hsp6/coords/                         TRUNCATED to 4 columns -> unusable (load_pub refuses them)
  data/big/hsp_codes/c<d>_<M>.txt                    Cohn-table spherical codes (CSV, one point per line)
  data/refs/lit/cohn_codes/spherical_codes_plain_2026-09-25.txt   'd,M,max cosine[,minimal polynomial]' for d <= 32
"""
import os, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REFS = os.path.join(ROOT, 'data', 'refs')
BIG = os.path.join(ROOT, 'data', 'big')
CODES = os.path.join(BIG, 'hsp_codes')
COHN_PLAIN = os.path.join(REFS, 'lit', 'cohn_codes', 'spherical_codes_plain_2026-09-25.txt')
OUT = os.path.join(HERE, 'out')
NMAX = {4: 300, 5: 300, 6: 250}


# ---------------------------------------------------------------- float geometry
def norms(c):
    return np.sqrt((np.asarray(c) ** 2).sum(-1))


def pair_d2(c):
    """Dense squared distances (n, n), +inf on the diagonal (difference form: no Gram cancellation)."""
    c = np.asarray(c, dtype=np.float64)
    D2 = ((c[:, None, :] - c[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(D2, np.inf)
    return D2


def rmin_parts(c):
    """(min wall slack 1 - |c_i|, min half pair distance)."""
    c = np.asarray(c, dtype=np.float64)
    w = float((1.0 - norms(c)).min())
    if len(c) < 2: return w, np.inf
    return w, float(np.sqrt(pair_d2(c).min()) / 2)


def rmin(c):
    w, p = rmin_parts(c)
    return min(w, p)


float_r = rmin          # the dense recheck IS the evaluation here (N <= 300)


def contacts(c, r, rel=1e-7):
    """Contacts per ball at relative tolerance rel: pairs with |ci - cj| <= 2r(1 + rel) plus the wall if 1 - |ci| <= r(1 + rel)."""
    c = np.asarray(c, dtype=np.float64)
    cnt = (pair_d2(c) <= (2 * r * (1 + rel)) ** 2).sum(1)
    return cnt + ((1.0 - norms(c)) <= r * (1 + rel)).astype(int)


def contact_graph(c, r, rel=1e-7):
    """(pairs (m, 2) in contact, wall indices in contact)."""
    D2 = pair_d2(c); I, J = np.nonzero(np.triu(D2 <= (2 * r * (1 + rel)) ** 2, 1))
    return np.stack([I, J], 1), np.nonzero((1.0 - norms(c)) <= r * (1 + rel))[0]


def random_ball(n, d, R, rng):
    """n points uniform in the d-ball of radius R."""
    v = rng.normal(size=(n, d)); v /= norms(v)[:, None]
    return v * (R * rng.random(n) ** (1.0 / d))[:, None]


def random_sphere(n, d, rng):
    v = rng.normal(size=(n, d)); return v / norms(v)[:, None]


def lower_priority():
    """BELOW_NORMAL for this process (Windows)."""
    try:
        import ctypes
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
    except Exception:
        pass


# ---------------------------------------------------------------- Packomania tables
def page_path(d):
    return os.path.join(REFS, f'packomania_hsp{d}_2026-09-25.html')


def table(d):
    """{N: dict(r=<radius string>, contacts, loose, boundary, core, ref)} from the saved HTML page (all rows)."""
    T = {}
    for line in open(page_path(d), encoding='utf-8', errors='replace'):
        if not line.startswith('<tr><td'): continue
        cells = [re.sub(r'<[^>]+>', '', x).strip() for x in re.findall(r'<td[^>]*>(.*?)</td>', line)]
        if len(cells) < 10 or not cells[0].isdigit(): continue
        iv = lambda s: int(s) if s.strip().isdigit() else 0
        T[int(cells[0])] = dict(r=cells[1], contacts=iv(cells[5]), loose=iv(cells[6]), boundary=iv(cells[7]),
                                core=iv(cells[8]), ref=cells[9])
    assert len(T) == NMAX[d], (d, len(T))
    return T


def pub_path(d, n):
    return os.path.join(BIG, f'hsp{d}', 'coords', f'hsp{d}-{n}.txt')


def load_pub_strings(d, n):
    """(radius string, [tuple of d coordinate strings]) from Packomania's file. hsp5 / hsp6 files are truncated to 4
    columns (recon 7) and are refused."""
    L = [l.split() for l in open(pub_path(d, n)) if l.strip()]
    r_s = L[0][0]; pts = [tuple(t[1:]) for t in L[1:]]
    if len(pts) != n or any(len(p) != d for p in pts):
        raise ValueError(f'hsp{d}-{n}.txt: {len(pts)} rows, widths {sorted({len(p) for p in pts})} (expected {n} x {d})')
    return r_s, pts


def load_pub(d, n):
    return np.array([[float(v) for v in p] for p in load_pub_strings(d, n)[1]], dtype=np.float64)


# ---------------------------------------------------------------- spherical codes
def code_path(d, m):
    return os.path.join(CODES, f'c{d}_{m}.txt')


def load_code(d, m):
    """Unit vectors (m, d) of the local Cohn-table code c<d>_<m>.txt, or None if not on disk."""
    p = code_path(d, m)
    if not os.path.exists(p): return None
    X = np.array([[float(v) for v in l.replace(',', ' ').split()] for l in open(p) if l.strip()], dtype=np.float64)
    assert X.shape == (m, d), (p, X.shape)
    return X / norms(X)[:, None]


def cohn_table():
    """{(d, M): max-cosine string} from the plain-text Cohn table (d, M, cos rounded UP at 12 decimals, [polynomial])."""
    T = {}
    for line in open(COHN_PLAIN, encoding='utf-8', errors='replace'):
        t = line.strip().split(',')
        if len(t) >= 3 and t[0].isdigit() and t[1].isdigit(): T[(int(t[0]), int(t[1]))] = t[2]
    return T


def code_radius(maxcos):
    """Ball radius of a code with max cosine maxcos placed on the shell |c| = 1 - r: r = s / (1 + s), s = sin(theta/2)."""
    s = np.sqrt(max(0.0, (1.0 - maxcos) / 2)); return s / (1 + s)


def max_cos(U):
    G = U @ U.T; np.fill_diagonal(G, -np.inf); return float(G.max())
