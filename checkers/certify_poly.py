# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Exact checker A for equal circles in a REGULAR k-GON (Packomania cpt = the pentagon, k = 5). Authors: Moki&Julio.
Standard library only. No floating point number takes part in any decision.

FRAME: circumradius 1, centred at the origin, one vertex at (0, 1) (for odd k the bottom side is flat); side j has outward
unit normal n_j at angle 90 + 360 (j + 1/2) / k degrees and lies at distance a = cos(180/k deg) from the origin.
For k = 5: a = (1 + sqrt5) / 4, normals at 126, 198, 270, 342, 54 degrees.

FILE: line 1 "r <radius>", then one "x y" line per centre. Every number must be a plain decimal (optionally with an
exponent, e.g. 1.5e-3); it is read EXACTLY. All numbers are put on one integer grid: value = I / S with S = 10^D.

VALIDITY (all exact):
  count : exactly N centres, r > 0.
  pairs : (Xi - Xj)^2 + (Yi - Yj)^2 >= (2R)^2 for ALL pairs (plain O(N^2) loop over every pair, integers only).
  walls : for every centre and every side, the signed slack  a - n_j . c - r >= 0.
          The normals and a are irrational. Each slack, times a positive constant (4S for k = 5), is written EXACTLY as
              E = c0 + c1 sqrt(d) - T sqrt(u + v sqrt(d))          (c0, c1, T integers, linear in S, R, X, Y)
          and its sign is decided exactly in the field Q(sqrt d, sqrt(u + v sqrt d)):
              sign(p + q sqrt d): trivial if p, q have the same sign, else compare p^2 with d q^2 (integers);
              E >= 0 with T != 0: decide sign(c0 + c1 sqrt d) exactly, then square once (both sides >= 0) and decide
              the sign of  (c0 + c1 sqrt d)^2 - T^2 (u + v sqrt d)  = p + q sqrt d  the same way.
          For k = 5 (d = 5):  4S.slack for the side at 270 deg = (S - 4R + 4Y) + S sqrt5
                              at 342 / 198 deg = (S - 4R - Y) + (S + Y) sqrt5 -/+ X sqrt(10 + 2 sqrt5)
                              at  54 / 126 deg = (S - 4R - Y) + (S - Y) sqrt5 -/+ X sqrt(10 - 2 sqrt5)
          (sin 72 = sqrt(10 + 2 sqrt5)/4, cos 72 = (sqrt5 - 1)/4, sin 36 = sqrt(10 - 2 sqrt5)/4, cos 36 = (1 + sqrt5)/4.)
          LOSS IN RADIUS: ZERO. The wall test is an exact sign test, not an inner polygon, so it is both sound and
          complete: a file is VALID iff the circles truly fit. (Exact tangency to a wall is impossible with decimal data,
          since a is irrational and sqrt(10 +- 2 sqrt5) is not in Q(sqrt5); the self-test shows the test resolves a wall
          to 1e-60 on both sides.)
  Supported k: 3, 4, 5, 6 (the regular polygons whose normals live in Q(sqrt d) or one square root above it). Others
  are refused.

VERDICT: IMPROVES iff valid at radius r and r > record_radius (both exact decimals); VALID_NOT_BETTER if valid and
r <= record; else INVALID <why>.

usage: py -3.11 certify_poly.py <k> <file> <N> <record_radius>        -> prints "VERDICT: ..."
       py -3.11 certify_poly.py --selftest                            -> unit tests + all 200 published cpt packings
       py -3.11 certify_poly.py --published [dir]                     -> only the 200 published cpt packings
"""
import sys, os, re, math, random
from fractions import Fraction as F

# ---------------------------------------------------------------- side tables
# Each side: (angle_deg, c0, c1, T, d, u, v), coefficient vectors over (S, R, X, Y); positive scale factor noted per k.
SIDES = {
    5: [  # 4S * slack, d = 5
        (270, (1, -4, 0, 4), (1, 0, 0, 0), (0, 0, 0, 0), 5, 0, 0),
        (342, (1, -4, 0, -1), (1, 0, 0, 1), (0, 0, 1, 0), 5, 10, 2),
        (54, (1, -4, 0, -1), (1, 0, 0, -1), (0, 0, 1, 0), 5, 10, -2),
        (126, (1, -4, 0, -1), (1, 0, 0, -1), (0, 0, -1, 0), 5, 10, -2),
        (198, (1, -4, 0, -1), (1, 0, 0, 1), (0, 0, -1, 0), 5, 10, 2),
    ],
    3: [  # 2S * slack, a = 1/2, d = 3
        (270, (1, -2, 0, 2), (0, 0, 0, 0), (0, 0, 0, 0), 3, 0, 0),
        (30, (1, -2, 0, -1), (0, 0, -1, 0), (0, 0, 0, 0), 3, 0, 0),
        (150, (1, -2, 0, -1), (0, 0, 1, 0), (0, 0, 0, 0), 3, 0, 0),
    ],
    4: [  # sqrt2 * S * slack, a = sqrt2/2, d = 2
        (45, (1, 0, -1, -1), (0, -1, 0, 0), (0, 0, 0, 0), 2, 0, 0),
        (135, (1, 0, 1, -1), (0, -1, 0, 0), (0, 0, 0, 0), 2, 0, 0),
        (225, (1, 0, 1, 1), (0, -1, 0, 0), (0, 0, 0, 0), 2, 0, 0),
        (315, (1, 0, -1, 1), (0, -1, 0, 0), (0, 0, 0, 0), 2, 0, 0),
    ],
    6: [  # 2S * slack, a = sqrt3/2, d = 3
        (0, (0, -2, -2, 0), (1, 0, 0, 0), (0, 0, 0, 0), 3, 0, 0),
        (60, (0, -2, -1, 0), (1, 0, 0, -1), (0, 0, 0, 0), 3, 0, 0),
        (120, (0, -2, 1, 0), (1, 0, 0, -1), (0, 0, 0, 0), 3, 0, 0),
        (180, (0, -2, 2, 0), (1, 0, 0, 0), (0, 0, 0, 0), 3, 0, 0),
        (240, (0, -2, 1, 0), (1, 0, 0, 1), (0, 0, 0, 0), 3, 0, 0),
        (300, (0, -2, -1, 0), (1, 0, 0, 1), (0, 0, 0, 0), 3, 0, 0),
    ],
}


def sgn(p, q, d):
    """Exact sign of p + q sqrt(d) for integers p, q and a positive non-square integer d."""
    if p >= 0 and q >= 0: return 1 if (p or q) else 0
    if p <= 0 and q <= 0: return -1
    t = p * p - d * q * q                     # != 0 because d is not a square and p, q != 0
    return (1 if t > 0 else -1) if p > 0 else (1 if t < 0 else -1)


def geq0(c0, c1, T, d, u, v):
    """Exact test  c0 + c1 sqrt(d) - T sqrt(u + v sqrt(d)) >= 0   (u + v sqrt d > 0 checked at import)."""
    e = sgn(c0, c1, d)
    if T == 0: return e >= 0
    P0, P1 = c0 * c0 + d * c1 * c1, 2 * c0 * c1          # (c0 + c1 sqrt d)^2
    W0, W1 = T * T * u, T * T * v                          # T^2 (u + v sqrt d)
    if T < 0:                                              # E + |T| sqrt(w) >= 0
        if e >= 0: return True
        return sgn(W0 - P0, W1 - P1, d) >= 0              # |T| sqrt w >= -E > 0  <=>  T^2 w >= E^2
    if e <= 0: return False                                # E <= 0 < T sqrt w
    return sgn(P0 - W0, P1 - W1, d) >= 0                   # E >= T sqrt w > 0  <=>  E^2 >= T^2 w


for _k, _t in SIDES.items():                               # structural guards on the tables
    for _s in _t:
        _d = _s[4]; assert math.isqrt(_d) ** 2 != _d, "d must be a non-square"
        if any(_s[3]): assert sgn(_s[5], _s[6], _d) > 0, "u + v sqrt d must be > 0"
    assert len(_t) == _k


def lin(cf, S, R, X, Y):
    return cf[0] * S + cf[1] * R + cf[2] * X + cf[3] * Y


def wall_ok(k, S, R, X, Y):
    """None if the circle (X/S, Y/S; R/S) is inside every side, else the angle (deg) of the first side it crosses."""
    for ang, a0, a1, aT, d, u, v in SIDES[k]:
        if not geq0(lin(a0, S, R, X, Y), lin(a1, S, R, X, Y), lin(aT, S, R, X, Y), d, u, v): return ang
    return None


# ---------------------------------------------------------------- exact decimal parsing
DEC = re.compile(r'^([+-]?)([0-9]+)?(?:\.([0-9]*))?(?:[eE]([+-]?[0-9]+))?$')   # ASCII digits only


def dec(s):
    """Exact decimal string -> (integer mantissa m, exponent e) with value m * 10^e. Raises ValueError otherwise."""
    g = DEC.match(s.strip())
    if not g or (g.group(2) is None and not g.group(3)): raise ValueError('not a decimal: ' + s)
    ip, fp = g.group(2) or '', g.group(3) or ''
    e = (int(g.group(4)) if g.group(4) else 0) - len(fp)
    if abs(e) > 5000 or len(ip) + len(fp) > 5000: raise ValueError('number too long: ' + s[:40])
    m = int(ip + fp or '0')
    return (-m if g.group(1) == '-' else m), e


def dec_fraction(s):
    m, e = dec(s)
    return F(m * 10 ** e) if e >= 0 else F(m, 10 ** -e)


# ---------------------------------------------------------------- the check
def check_values(k, r_s, pts_s, n, rec_s):
    """k-gon, radius string, list of (x, y) strings, expected count, record radius string -> verdict string."""
    if k not in SIDES: return f'INVALID unsupported k={k} (supported: {sorted(SIDES)})'
    if len(pts_s) != n: return f'INVALID count {len(pts_s)} != {n}'
    try:
        vals = [dec(r_s)] + [dec(t) for p in pts_s for t in p]
        rec = dec_fraction(rec_s)
    except ValueError as ex:
        return f'INVALID {ex}'
    D = max(0, max(-e for _, e in vals)); S = 10 ** D
    I = [m * 10 ** (e + D) for m, e in vals]                  # exact: value = I / S
    R = I[0]; X = I[1::2]; Y = I[2::2]
    if R <= 0: return 'INVALID r <= 0'
    for i in range(n):
        bad = wall_ok(k, S, R, X[i], Y[i])
        if bad is not None: return f'INVALID circle {i + 1} crosses the side with normal at {bad} deg'
    R4 = 4 * R * R
    for i in range(n):
        xi, yi = X[i], Y[i]
        for j in range(i + 1, n):
            dx = xi - X[j]; dy = yi - Y[j]
            if dx * dx + dy * dy < R4: return f'INVALID overlap {i + 1} {j + 1}'
    return 'IMPROVES' if F(R, S) > rec else 'VALID_NOT_BETTER'


def read_cert(path):
    L = [l.split() for l in open(path) if l.strip()]
    if not L or L[0][0] != 'r' or len(L[0]) != 2: raise ValueError('header must be "r <radius>"')
    if any(len(t) != 2 for t in L[1:]): raise ValueError('every centre line must be "x y"')
    return L[0][1], [tuple(t) for t in L[1:]]


def check_file(k, path, n, rec_s):
    try:
        r_s, pts = read_cert(path)
    except (ValueError, OSError) as ex:
        return f'INVALID {ex}'
    return check_values(k, r_s, pts, n, rec_s)


# ---------------------------------------------------------------- published cpt packings
ROOT = os.path.dirname(os.path.abspath(__file__))
PUB = os.path.join(ROOT, 'data', 'big', 'cpt')
PAGE = os.path.join(ROOT, 'data', 'refs', 'packomania_cpt_2026-09-25.html')


def page_radii():
    s = open(PAGE, encoding='utf-8', errors='replace').read()
    return {int(n): r for n, r in re.findall(r'name="cpt(\d+)">.*?</a></td>\s*<td>([0-9.]+)</td>', s, re.S)}


def load_pub(n, d=PUB):
    pts = []
    for line in open(os.path.join(d, f'cpt{n}.txt')):
        p = line.split()
        if not p: continue
        if len(p) == 3: pts.append((p[1], p[2]))
        elif len(p) == 1 and n == 1: pts.append(('0', '0'))     # cpt1.txt prints no coordinates: the centred circle
        else: raise ValueError(f'cpt{n}.txt: bad line {line!r}')
    return pts


def fmt_dec(q, places=40):
    """Fraction -> decimal string rounded toward zero to `places` digits (for building test radii)."""
    neg = q < 0; q = abs(q); m = (q.numerator * 10 ** places) // q.denominator
    s = str(m).rjust(places + 1, '0'); s = s[:-places] + '.' + s[-places:]
    return ('-' if neg else '') + s


def published_check(d=PUB, verbose=True):
    """Every published cpt packing, at printed radius - off, must be VALID_NOT_BETTER (record = printed radius)."""
    T = page_radii(); out = {}
    for off in ('0', '0.000000000001', '0.000000000002'):
        bad = []
        for n in range(1, 201):
            pts = load_pub(n, d); rp = T[n]
            r = fmt_dec(dec_fraction(rp) - dec_fraction(off), 15)
            v = check_values(5, r, pts, n, rp)
            if v != 'VALID_NOT_BETTER': bad.append((n, v))
        out[off] = bad
        if verbose:
            print(f'published cpt, r = printed - {off}: {200 - len(bad)}/200 VALID_NOT_BETTER; not: '
                  + (', '.join(f'N={n}: {v}' for n, v in bad[:12]) + (' ...' if len(bad) > 12 else '') if bad else 'none'))
    return out


# ---------------------------------------------------------------- self-tests
def selftest():
    from decimal import Decimal, getcontext
    getcontext().prec = 120
    P = 60                                                    # digits used to resolve walls in the tests
    s5 = math.isqrt(5 * 10 ** (2 * P))                        # floor(sqrt5 * 10^P), proved by the next line
    assert s5 * s5 <= 5 * 10 ** (2 * P) < (s5 + 1) ** 2
    a_lo = F((10 ** P + s5) // 4, 10 ** P)                    # a - 1e-60 < a_lo < a   (a irrational)
    a_hi = a_lo + F(1, 10 ** P)                               # a < a_hi
    e30 = F(1, 10 ** 30)

    def run(k, r, pts, rec='0'):
        return check_values(k, r if isinstance(r, str) else fmt_dec(r, 70), [tuple(p) for p in pts], len(pts), rec)

    # 1. N = 1 at the centre, r = a -/+ 1e-30, and the sharp version at 1e-60
    assert run(5, a_lo - e30, [('0', '0')]) == 'IMPROVES'
    assert run(5, a_hi + e30, [('0', '0')]).startswith('INVALID circle 1')
    assert run(5, a_lo, [('0', '0')]) == 'IMPROVES'                       # a - 1e-60 < r < a
    assert run(5, a_hi, [('0', '0')]).startswith('INVALID')               # a < r < a + 1e-60
    assert run(5, a_lo, [('0', '0')], fmt_dec(a_lo, 70)) == 'VALID_NOT_BETTER'   # r == record -> not better
    assert run(5, a_lo, [('0', '0')], '0.809016994375') == 'VALID_NOT_BETTER'  # printed N=1 radius is above a
    # 2. one circle touching a side (to 1e-60 on each side of tangency; exact tangency is impossible with decimals)
    r = F(1, 10)
    assert run(5, r, [('0', fmt_dec(-(a_lo - r), 70))]) == 'IMPROVES'       # y = r - a_lo > r - a : inside
    assert run(5, r, [('0', fmt_dec(-(a_hi - r), 70))]).startswith('INVALID circle 1 crosses the side with normal at 270')
    D5 = Decimal(5).sqrt(); aD = (1 + D5) / 4
    for ang, (nx, ny) in {342: ((10 + 2 * D5).sqrt() / 4, -(D5 - 1) / 4), 54: ((10 - 2 * D5).sqrt() / 4, aD),
                          126: (-(10 - 2 * D5).sqrt() / 4, aD), 198: (-(10 + 2 * D5).sqrt() / 4, -(D5 - 1) / 4)}.items():
        for sgn_, want in ((-1, 'IMPROVES'), (1, 'INVALID')):
            t = aD - Decimal('0.1') + sgn_ * Decimal('1e-60')              # centre at distance t along the normal
            x, y = (t * nx).quantize(Decimal('1e-66')), (t * ny).quantize(Decimal('1e-66'))   # rounding moves it < 2e-66
            v = run(5, '0.1', [(str(x), str(y))])
            assert v.startswith(want), (ang, sgn_, v)
            if want == 'INVALID': assert f'normal at {ang} deg' in v, (ang, v)
    # 3. two circles at exactly 2r (horizontal, and a 3-4-5 direction), and tiny overlaps
    assert run(5, '0.1', [('-0.1', '0'), ('0.1', '0')]) == 'IMPROVES'
    assert run(5, '0.1', [('0.06', '0.08'), ('-0.06', '-0.08')]) == 'IMPROVES'
    assert run(5, '0.1', [('-0.1', '0'), ('0.0999999999999999999999999999999999999999', '0')]) == 'INVALID overlap 1 2'
    assert run(5, '0.1', [('0.06', '0.08'), ('-0.06', '-0.0799999999999999999999999999999999999999')]) == 'INVALID overlap 1 2'
    far = [('-0.3', '-0.3'), ('0.3', '0.3'), ('0', '0.4'), ('-0.3', '-0.29999999999999999999999999999999999999999')]
    assert run(5, '0.05', far) == 'INVALID overlap 1 4', 'overlap far apart in the list must be caught'
    # 4. orientation (vertex up, flat bottom) and format guards
    assert run(5, '0.01', [('0', '0.9')]) == 'IMPROVES'                   # under the top vertex: inside
    assert run(5, '0.01', [('0', '-0.9')]).startswith('INVALID')          # below the flat bottom (y = -0.809)
    assert run(5, '0.05', [('0', '-0.759')]) == 'IMPROVES'                # 0.809017 - 0.759 = 0.050017 >= 0.05
    assert run(5, '0.05', [('0', '-0.759')], '0.05') == 'VALID_NOT_BETTER'
    assert run(5, '0.05', [('0', '-0.759')], '0.0500000000000000000000000000000000001') == 'VALID_NOT_BETTER'
    assert run(5, '0.05', [('0', '-0.7591')]).startswith('INVALID circle 1 crosses the side with normal at 270')
    assert run(5, '0', [('0', '0')]) == 'INVALID r <= 0'
    assert run(5, '0.1', [('0', 'abc')]).startswith('INVALID not a decimal')
    for junk in ('nan', 'inf', '.', '-', '1/3', '0x1', '٣', '1e', '--1', '1.2.3'):
        assert run(5, '0.1', [('0', junk)]).startswith('INVALID not a decimal'), junk
    assert run(7, '0.1', [('0', '0')]).startswith('INVALID unsupported k')
    assert check_values(5, '0.1', [('0', '0')], 2, '0') == 'INVALID count 1 != 2'
    assert run(5, '1.0e-1', [('-1e-1', '0'), ('0.1', '0E0')]) == 'IMPROVES'          # exponent form is read exactly
    assert run(5, '1.5e-1', [('-1e-1', '0'), ('0.1', '0E0')]) == 'INVALID overlap 1 2'
    # 5. float cross-check of every side table (the float is only the reference in this TEST, never in a decision):
    #    random rational circles; wherever the float slack is clearly nonzero, the exact sign must agree.
    rng = random.Random(20260925); agree = 0
    for k in SIDES:
        a = math.cos(math.pi / k)
        for _ in range(3000):
            x, y =F(rng.randrange(-10 ** 9, 10 ** 9), 10 ** 9), F(rng.randrange(-10 ** 9, 10 ** 9), 10 ** 9)
            rr = F(rng.randrange(1, 4 * 10 ** 8), 10 ** 9)
            S = 10 ** 9; X, Y, R = int(x * S), int(y * S), int(rr * S)
            for ang, a0, a1, aT, d, u, v in SIDES[k]:
                th = math.radians(ang); fs = a - float(rr) - (math.cos(th) * float(x) + math.sin(th) * float(y))
                if abs(fs) < 1e-9: continue
                ex = geq0(lin(a0, S, R, X, Y), lin(a1, S, R, X, Y), lin(aT, S, R, X, Y), d, u, v)
                assert ex == (fs > 0), (k, ang, x, y, rr, fs)
                agree += 1
        # the side normals must be exactly the k directions 90 + 360 (j + 1/2) / k
        assert sorted(s[0] % 360 for s in SIDES[k]) == sorted(round(90 + 360 * (j + 0.5) / k) % 360 for j in range(k)), k
    print(f'float cross-check of the side tables (k = 3, 4, 5, 6): {agree} side tests agree, 0 disagree')
    # 6. sgn / geq0 unit checks
    assert sgn(3, -1, 5) == 1 and sgn(2, -1, 5) == -1 and sgn(-3, 1, 5) == -1 and sgn(-2, 1, 5) == 1 and sgn(0, 0, 5) == 0
    assert geq0(0, 0, 0, 5, 0, 0) and not geq0(-1, 0, 0, 5, 0, 0)
    assert geq0(4, 0, 1, 5, 10, 2) and not geq0(3, 0, 1, 5, 10, 2)        # sqrt(10 + 2 sqrt5) = 3.804
    assert geq0(-3, 0, -1, 5, 10, 2) and not geq0(-4, 0, -1, 5, 10, 2)
    assert geq0(3, 0, 1, 5, 10, -2) and not geq0(2, 0, 1, 5, 10, -2)       # sqrt(10 - 2 sqrt5) = 2.351
    print('unit self-tests OK')
    res = published_check()
    bad = res['0.000000000002']
    # 7. negative control: the same packings mirrored y -> -y (= the pentagon rotated 36 deg, vertex down) must fail
    T = page_radii(); mir_ok = []
    for n in range(2, 201):
        pts = [(x, y[1:] if y.startswith('-') else '-' + y) for x, y in load_pub(n)]
        r = fmt_dec(dec_fraction(T[n]) - F(2, 10 ** 12), 15)
        if not check_values(5, r, pts, n, T[n]).startswith('INVALID'): mir_ok.append(n)
    print(f'negative control (mirrored frame, N = 2-200): {199 - len(mir_ok)}/199 INVALID; still valid: {mir_ok or "none"}')
    ok = not bad and not mir_ok
    print('SELFTEST ' + ('OK' if ok else f'FAILED: published not VALID_NOT_BETTER at printed - 2e-12: {bad}; mirrored still valid: {mir_ok}'))
    return ok


if __name__ == '__main__':
    if sys.argv[1:2] == ['--selftest']:
        sys.exit(0 if selftest() else 1)
    elif sys.argv[1:2] == ['--published']:
        r = published_check(sys.argv[2] if len(sys.argv) > 2 else PUB)
        sys.exit(0 if not r['0.000000000002'] else 1)
    elif len(sys.argv) == 5:
        print('VERDICT:', check_file(int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), sys.argv[4]))
    else:
        print(__doc__); sys.exit(2)
