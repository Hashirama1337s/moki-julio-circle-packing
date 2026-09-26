# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Exact checker A for equal spheres in a CUBE (Packomania scu). Authors: Moki&Julio.
Python 3.11, standard library only. No floating point number takes part in any decision.

FRAME: cube of side 1 centred at the origin; a sphere of radius r centred at c = (x, y, z) is inside iff
       -1/2 + r <= x, y, z <= 1/2 - r;  spheres i != j do not overlap iff |c_i - c_j|^2 >= (2r)^2.

FILE: line 1 "r <radius>", then exactly N lines "x y z" (blank lines are ignored). Every number must be a plain decimal,
optionally with an exponent (e.g. 1.5e-17); anything else (nan, inf, fractions, hex, non-ASCII digits) is refused. The
numbers are read EXACTLY: value = m * 10^e with integers m, e, and all of them are put on ONE integer grid
value = I / S with S = 10^D, D = the largest number of decimal places needed. Every test below is on these integers.

VALIDITY (all exact, all integers):
  count : exactly N centres; r > 0.
  walls : for every centre and axis,  2 X >= 2 R - S  and  2 X <= S - 2 R   (i.e. -1/2 + r <= x <= 1/2 - r, times 2S).
  pairs : (Xi - Xj)^2 + (Yi - Yj)^2 + (Zi - Zj)^2 >= (2R)^2 for EVERY pair i < j, checked COMPLETELY with a bucket grid:
          cell side C = 2R (an integer); sphere i goes to the cell (floor(Xi / C), floor(Yi / C), floor(Zi / C)) (exact
          integer floor division). Only pairs in the same or adjacent cells (27 cells, each unordered pair once) are
          tested. COMPLETENESS: if a pair overlaps, then (Xi - Xj)^2 < (2R)^2 so |Xi - Xj| < C; for Xi <= Xj this gives
          Xj / C < Xi / C + 1, hence floor(Xj / C) <= floor(Xi / C) + 1, and likewise for y and z; so every overlapping
          pair lies in the same or adjacent cells and is tested. Pairs at distance exactly 2R are valid (>=).
          --allpairs replaces the grid by the plain O(N^2) loop (used by the self-test to cross-check the grid).
VERDICT: IMPROVES iff valid at radius r and r > record_radius (both exact decimals); VALID_NOT_BETTER if valid and
         r <= record; else INVALID <why>. LOSS IN RADIUS: zero (the tests are exact, so VALID means the spheres fit).

usage: python checkers/certify_cube.py <file> <N> <record_radius> [--allpairs]   -> prints "VERDICT: ..."
       python checkers/certify_cube.py --selftest       -> unit tests + every published scu<N>.txt in data/big/scu/ (N <= 1008)
       python checkers/certify_cube.py --published      -> only the published-file test
"""
import sys, os, re, random
from fractions import Fraction as F

# ---------------------------------------------------------------- exact decimal parsing
DEC = re.compile(r'^([+-]?)([0-9]+)?(?:\.([0-9]*))?(?:[eE]([+-]?[0-9]+))?$')   # ASCII digits only (re.ASCII below)
DEC = re.compile(DEC.pattern, re.ASCII)


def dec(s):
    """Exact decimal string -> (integer mantissa m, exponent e) with value m * 10^e. Raises ValueError otherwise."""
    g = DEC.match(s.strip())
    if not g or (g.group(2) is None and not g.group(3)): raise ValueError('not a decimal: ' + s[:40])
    ip, fp = g.group(2) or '', g.group(3) or ''
    e = (int(g.group(4)) if g.group(4) else 0) - len(fp)
    if abs(e) > 5000 or len(ip) + len(fp) > 5000: raise ValueError('number too long: ' + s[:40])
    m = int(ip + fp or '0')
    return (-m if g.group(1) == '-' else m), e


def dec_fraction(s):
    m, e = dec(s)
    return F(m * 10 ** e) if e >= 0 else F(m, 10 ** -e)


# ---------------------------------------------------------------- the check
def first_overlap_grid(X, Y, Z, R):
    """First overlapping pair (1-based) or None; complete by the bucket argument in the module docstring."""
    C = 2 * R; R4 = 4 * R * R; cells = {}
    for i in range(len(X)):
        cells.setdefault((X[i] // C, Y[i] // C, Z[i] // C), []).append(i)
    bad = None
    for (a, b, c), L in cells.items():
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if (da, db, dc) < (0, 0, 0): continue                  # each unordered cell pair once
                    M = cells.get((a + da, b + db, c + dc))
                    if M is None: continue
                    same = (da, db, dc) == (0, 0, 0)
                    for s, i in enumerate(L):
                        xi, yi, zi = X[i], Y[i], Z[i]
                        for j in (L[s + 1:] if same else M):
                            dx = xi - X[j]; dy = yi - Y[j]; dz = zi - Z[j]
                            if dx * dx + dy * dy + dz * dz < R4:
                                p = (min(i, j) + 1, max(i, j) + 1)
                                if bad is None or p < bad: bad = p
    return bad


def first_overlap_all(X, Y, Z, R):
    R4 = 4 * R * R
    for i in range(len(X)):
        xi, yi, zi = X[i], Y[i], Z[i]
        for j in range(i + 1, len(X)):
            dx = xi - X[j]; dy = yi - Y[j]; dz = zi - Z[j]
            if dx * dx + dy * dy + dz * dz < R4: return (i + 1, j + 1)
    return None


def check_values(r_s, pts_s, n, rec_s, allpairs=False):
    """radius string, list of (x, y, z) strings, expected count, record radius string -> verdict string."""
    if len(pts_s) != n: return f'INVALID count {len(pts_s)} != {n}'
    try:
        vals = [dec(r_s)] + [dec(t) for p in pts_s for t in p]
        rec = dec_fraction(rec_s)
    except ValueError as ex:
        return f'INVALID {ex}'
    D = max(0, max(-e for _, e in vals)); S = 10 ** D
    I = [m * 10 ** (e + D) for m, e in vals]                     # exact: value = I / S
    R = I[0]; X = I[1::3]; Y = I[2::3]; Z = I[3::3]
    if R <= 0: return 'INVALID r <= 0'
    lo, hi = 2 * R - S, S - 2 * R                                 # 2X in [2R - S, S - 2R]
    if lo > hi: return 'INVALID r > 1/2'
    for i in range(n):
        for ax, v in (('x', X[i]), ('y', Y[i]), ('z', Z[i])):
            if 2 * v < lo: return f'INVALID sphere {i + 1} crosses the face {ax} = -1/2'
            if 2 * v > hi: return f'INVALID sphere {i + 1} crosses the face {ax} = +1/2'
    bad = (first_overlap_all if allpairs else first_overlap_grid)(X, Y, Z, R)
    if bad: return f'INVALID overlap {bad[0]} {bad[1]}'
    return 'IMPROVES' if F(R, S) > rec else 'VALID_NOT_BETTER'


def read_cert(path):
    L = [l.split() for l in open(path, encoding='ascii') if l.strip()]
    if not L or L[0][0] != 'r' or len(L[0]) != 2: raise ValueError('header must be "r <radius>"')
    if any(len(t) != 3 for t in L[1:]): raise ValueError('every centre line must be "x y z"')
    return L[0][1], [tuple(t) for t in L[1:]]


def check_file(path, n, rec_s, allpairs=False):
    try:
        r_s, pts = read_cert(path)
    except (ValueError, OSError, UnicodeDecodeError) as ex:
        return f'INVALID {ex}'
    return check_values(r_s, pts, n, rec_s, allpairs)


# ---------------------------------------------------------------- published scu packings
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PUB = os.path.join(ROOT, 'data', 'big', 'scu')
PAGE = os.path.join(ROOT, 'data', 'refs', 'packomania_scu_2026-09-25.html')
SIDE = os.path.join(ROOT, 'data', 'refs', 'scu_radius_2026-09-25.txt')
STALE = (83, 94, 95, 96, 109, 110, 507, 619)       # coordinate file reproduces the old (radius.txt) radius, HTML is higher


def page_radii():
    T = {}
    for line in open(PAGE, encoding='utf-8', errors='replace'):
        m = re.search(r'name="scu(\d+)"', line)
        if not m or not line.startswith('<tr>'): continue
        cells = [re.sub(r'<[^>]+>', '', c).strip() for c in re.findall(r'<td[^>]*>(.*?)</td>', line)]
        T[int(m.group(1))] = cells[1]
    return T


def side_radii():
    return {int(l.split()[0]): l.split()[1] for l in open(SIDE) if l.strip()}


def load_pub(n, d=PUB):
    pts = []
    for line in open(os.path.join(d, f'scu{n}.txt')):
        p = line.split()
        if not p: continue
        if len(p) != 4: raise ValueError(f'scu{n}.txt: bad line {line!r}')
        pts.append((p[1], p[2], p[3]))
    return pts


def fmt_dec(q, places=40):
    """Fraction -> decimal string rounded toward zero to `places` digits (for building test radii)."""
    neg = q < 0; q = abs(q); m = (q.numerator * 10 ** places) // q.denominator
    s = str(m).rjust(places + 1, '0'); s = s[:-places] + '.' + s[-places:]
    return ('-' if neg else '') + s


def published_check(d=PUB, verbose=True):
    """Every downloaded scu<N>.txt (N <= 1008). Non-stale: VALID_NOT_BETTER at printed - 1e-15 (record = printed) and
    INVALID at printed + 1e-9. Stale (83 94 95 96 109 110 507 619): the file must be INVALID at the HTML radius - 1e-15,
    VALID at its own old radius (radius.txt) - 1e-15 and INVALID at the old radius + 1e-9."""
    T, SR = page_radii(), side_radii()
    Ns = sorted(int(m.group(1)) for f in os.listdir(d) for m in [re.fullmatch(r'scu(\d+)\.txt', f)] if m and int(m.group(1)) <= 1008)
    lo_e, hi_e = F(1, 10 ** 15), F(1, 10 ** 9); fails = []; n_ok = 0; stale_ok = 0
    for n in Ns:
        pts = load_pub(n, d)
        if n in STALE:
            rp, ro = dec_fraction(T[n]), dec_fraction(SR[n])
            v1 = check_values(fmt_dec(rp - lo_e), pts, n, T[n])
            v2 = check_values(fmt_dec(ro - lo_e), pts, n, SR[n])
            v3 = check_values(fmt_dec(ro + hi_e), pts, n, SR[n])
            good = v1.startswith('INVALID') and v2 == 'VALID_NOT_BETTER' and v3.startswith('INVALID')
            stale_ok += good
            if not good: fails.append((n, 'stale', v1, v2, v3))
            continue
        rp = dec_fraction(T[n])
        v1 = check_values(fmt_dec(rp - lo_e), pts, n, T[n])
        v2 = check_values(fmt_dec(rp + hi_e), pts, n, T[n])
        if v1 == 'VALID_NOT_BETTER' and v2.startswith('INVALID'): n_ok += 1
        else: fails.append((n, v1, v2))
    ns = sum(n in STALE for n in Ns)
    if verbose:
        print(f'published scu files checked: {len(Ns)} (N = {Ns[0]}..{Ns[-1]}); non-stale {len(Ns) - ns}: {n_ok} VALID_NOT_BETTER at '
              f'printed - 1e-15 and INVALID at printed + 1e-9; stale {ns}: {stale_ok} INVALID at the HTML radius, VALID at '
              f'the old radius - 1e-15, INVALID at old + 1e-9; failures: {fails[:10] or "none"}')
    return Ns, fails


# ---------------------------------------------------------------- self-tests
def selftest():
    def run(r, pts, rec='0', allpairs=False):
        return check_values(r if isinstance(r, str) else fmt_dec(r, 70), [tuple(p) for p in pts], len(pts), rec, allpairs)
    e40 = F(1, 10 ** 40)
    # 1. one sphere: centred, r = 1/2 exactly valid, +1e-40 invalid; touching each face exactly, and 1e-40 beyond
    assert run('0.5', [('0', '0', '0')]) == 'IMPROVES'
    assert run(F(1, 2) + e40, [('0', '0', '0')]) == 'INVALID r > 1/2'
    assert run(F(1, 2) - e40, [('0', '0', '0')]) == 'IMPROVES'
    r = F(1, 10); edge = F(1, 2) - r                                   # centre coordinate at exact tangency
    for ax in range(3):
        for sg, face in ((1, '+1/2'), (-1, '-1/2')):
            p = ['0', '0', '0']; p[ax] = fmt_dec(sg * edge, 70)
            assert run('0.1', [tuple(p)]) == 'IMPROVES', (ax, sg)
            p[ax] = fmt_dec(sg * (edge + e40), 70)
            v = run('0.1', [tuple(p)]); assert v == f'INVALID sphere 1 crosses the face {"xyz"[ax]} = {face}', v
            p[ax] = fmt_dec(sg * (edge - e40), 70)
            assert run('0.1', [tuple(p)]) == 'IMPROVES'
    # 2. two spheres at exactly 2r: along an axis, along (1,2,2)/3 and (2,3,6)/7 (exact rational unit vectors); 1e-40 closer
    for (a, b, c), q in (((1, 0, 0), 1), ((1, 2, 2), 3), ((2, 3, 6), 7), ((0, 3, 4), 5)):
        for scale, want in ((F(1), 'IMPROVES'), (1 - e40, 'INVALID overlap 1 2')):
            h = F(1, 10) * scale                                        # half separation; r = 0.1, sep = 2r * scale
            p1 = tuple(fmt_dec(h * t / q, 70) for t in (a, b, c)); p2 = tuple(fmt_dec(-h * t / q, 70) for t in (a, b, c))
            exact = all(F(h * t / q) == dec_fraction(s) for t, s in zip((a, b, c), p1))
            if not exact: continue                                      # only exactly representable decimals
            assert run('0.1', [p1, p2]) == want, (a, b, c, scale, run('0.1', [p1, p2]))
    assert run('0.1', [('-0.1', '0', '0'), ('0.1', '0', '0')]) == 'IMPROVES'
    assert run('0.1', [('-0.1', '0', '0'), ('0.0999999999999999999999999999999999999999', '0', '0')]) == 'INVALID overlap 1 2'
    assert run('0.1', [('0.06', '0.08', '0'), ('-0.06', '-0.08', '0')]) == 'IMPROVES'
    assert run('0.1', [('0.06', '0.08', '0'), ('-0.06', '-0.0799999999999999999999999999999999999999', '0')]) == 'INVALID overlap 1 2'
    # 3. overlaps far apart in the list, across bucket boundaries, at negative coordinates (floor division), and the
    #    reported pair is the lexicographically first one
    far = [('-0.3', '-0.3', '-0.3'), ('0.3', '0.3', '0.3'), ('0', '0.3', '-0.3'), ('-0.3', '-0.29999999999999999999', '-0.3')]
    assert run('0.05', far) == 'INVALID overlap 1 4'
    assert run('0.05', far, allpairs=True) == 'INVALID overlap 1 4'
    # C = 0.1: points at -0.10000..01 and -0.0000..01 sit in cells -2 and -1 (adjacent); distance just below 2r
    assert run('0.05', [('-0.1000000000000000000001', '0', '0'), ('-0.0000000000000000000002', '0', '0')]) == 'INVALID overlap 1 2'
    assert run('0.05', [('-0.1000000000000000000001', '0', '0'), ('0.0000000000000000000000', '0', '0')]) == 'IMPROVES'
    # diagonal neighbours: cells (0, 0, 0) and (1, 1, 0), separation sqrt(2) * 0.03 < 2r
    assert run('0.05', [('0.0999999999999999999999', '0.0999999999999999999999', '0'), ('0.13', '0.13', '0')]) == 'INVALID overlap 1 2'
    assert run('0.05', [('-0.0999999999999999999999', '0.0999999999999999999999', '-0.05'), ('-0.13', '0.13', '-0.06')]) == 'INVALID overlap 1 2'
    # 4. record comparison and format guards
    assert run('0.1', [('0', '0', '0')], '0.1') == 'VALID_NOT_BETTER'
    assert run('0.1', [('0', '0', '0')], '0.0999999999999999999999999999999999999999') == 'IMPROVES'
    assert run('0.1', [('0', '0', '0')], '0.1000000000000000000000000000000000000001') == 'VALID_NOT_BETTER'
    assert run('0', [('0', '0', '0')]) == 'INVALID r <= 0'
    assert run('-0.1', [('0', '0', '0')]) == 'INVALID r <= 0'
    for junk in ('nan', 'inf', '.', '-', '1/3', '0x1', '٣', '1e', '--1', '1.2.3', '1,5', ' '):
        assert run('0.1', [('0', junk, '0')]).startswith('INVALID not a decimal'), junk
    assert check_values('0.1', [('0', '0', '0')], 2, '0') == 'INVALID count 1 != 2'
    assert run('1.0e-1', [('-1e-1', '0', '0'), ('0.1', '0E0', '4.1e-1')]).startswith('INVALID sphere 2 crosses the face z = +1/2')
    assert run('1.0e-1', [('-1e-1', '0', '0'), ('0.1', '0E0', '0e-5')]) == 'IMPROVES'        # exponents read exactly
    assert run('1.5e-1', [('-1e-1', '0', '0'), ('0.1', '0E0', '0')]) == 'INVALID overlap 1 2'
    assert run('0.1', [('-0.1', '0', '0'), ('0.1', '0', '0')], '1e-1') == 'VALID_NOT_BETTER'
    # 5. the grid agrees with the plain all-pairs loop on random packings made tight on purpose (centres in the middle
    #    half of the cube so the walls never decide; R = floor(dmin / 2) + {0, +1, -1, +k}: valid, just invalid, ...)
    from math import isqrt
    rng = random.Random(20260925); agree = 0; kinds = {}
    for trial in range(1500):
        n = rng.randrange(2, 50); S = 10 ** rng.randrange(3, 25); w = S // 4
        P = [tuple(rng.randrange(-w, w + 1) for _ in range(3)) for _ in range(n)]
        d2 = min(sum((a - b) ** 2 for a, b in zip(P[i], P[j])) for i in range(n) for j in range(i + 1, n))
        R0 = isqrt(d2) // 2                                               # (2 R0)^2 <= d2 < (2 R0 + 2)^2
        for R in {R0, R0 + 1, max(1, R0 - 1), R0 + rng.randrange(1, 3 + R0 // 3)}:
            if R <= 0 or 2 * R > S // 4: continue
            pts = [tuple(fmt_dec(F(c, S), 30) for c in p) for p in P]
            rs = fmt_dec(F(R, S), 30)
            a1 = check_values(rs, pts, n, '0'); a2 = check_values(rs, pts, n, '0', allpairs=True)
            assert a1 == a2, (trial, R, a1, a2); agree += 1; kinds[a1.split()[0]] = kinds.get(a1.split()[0], 0) + 1
    print(f'grid vs all-pairs on random tight packings: {agree} verdicts agree, 0 disagree ({kinds})')
    print('unit self-tests OK')
    Ns, fails = published_check()
    # 6. negative control on EVERY non-stale published file: take its closest pair (exact integer search on a grid of
    #    side 4R) and move the first sphere towards the second by 1e-12 of their separation -> must be INVALID overlap
    T = page_radii(); neg_bad = []; neg_n = 0
    for n in Ns:
        if n in STALE or n < 2: continue
        pts = load_pub(n); rp = dec_fraction(T[n]); r_s = fmt_dec(rp - F(1, 10 ** 15))
        vals = [dec(t) for p in pts for t in p]; D = max(0, max(-e for _, e in vals)); S = 10 ** D
        I = [m * 10 ** (e + D) for m, e in vals]; X, Y, Z = I[0::3], I[1::3], I[2::3]
        C = 4 * (rp.numerator * S // rp.denominator) + 1; cells = {}
        for k in range(n): cells.setdefault((X[k] // C, Y[k] // C, Z[k] // C), []).append(k)
        best = None
        for (a, b, c), L in cells.items():
            for i in L:
                for da in (-1, 0, 1):
                    for db in (-1, 0, 1):
                        for dc in (-1, 0, 1):
                            for j in cells.get((a + da, b + db, c + dc), ()):
                                if j <= i: continue
                                d2 = (X[i] - X[j]) ** 2 + (Y[i] - Y[j]) ** 2 + (Z[i] - Z[j]) ** 2
                                if best is None or d2 < best[0]: best = (d2, i, j)
        _, i, j = best; t = F(1, 10 ** 12); moved = list(pts)
        moved[i] = tuple(fmt_dec(F(u, S) + t * F(w - u, S), D + 14) for u, w in ((X[i], X[j]), (Y[i], Y[j]), (Z[i], Z[j])))
        v = check_values(r_s, moved, n, T[n]); neg_n += 1
        if not v.startswith('INVALID overlap'): neg_bad.append((n, v))
    print(f'negative control (closest pair pushed together by 1e-12 of its separation): {neg_n - len(neg_bad)}/{neg_n} '
          f'INVALID overlap; not: {neg_bad[:10] or "none"}')
    ok = not fails and not neg_bad
    print('SELFTEST ' + ('OK' if ok else f'FAILED: {fails[:10]} {neg_bad}'))
    return ok


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if a != '--allpairs']; ap = '--allpairs' in sys.argv[1:]
    if args[:1] == ['--selftest']:
        sys.exit(0 if selftest() else 1)
    elif args[:1] == ['--published']:
        sys.exit(0 if not published_check()[1] else 1)
    elif len(args) == 3:
        print('VERDICT:', check_file(args[0], int(args[1]), args[2], ap))
    else:
        print(__doc__); sys.exit(2)
