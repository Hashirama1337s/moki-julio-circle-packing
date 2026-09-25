# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Exact checker A for equal circles in a REGULAR k-GON, any k >= 3 (Packomania cxd = 16-gon, cpd = 15-gon, cpt = pentagon).
Authors: Moki&Julio.  Python standard library only.  No floating point number takes part in any decision.

FRAME: circumradius 1, centred at the origin, one side horizontal at the bottom.  Side j = 0..k-1 has outward unit normal
    n_j = (cos t_j, sin t_j),  t_j = 270 + 360 j / k degrees = pi u_j with u_j = 3/2 + 2 j / k,
at distance a = cos(pi / k) from the origin.  A circle (c, r) is inside iff  a - n_j . c - r >= 0  for every j.
(k = 16: flat top too; k = 15: vertex at (0, 1); for odd k this is exactly the frame of certify_poly.py.)

FILE: line 1 "r <radius>", then one "x y" line per centre.  Every number must be a plain decimal (optionally with an
exponent, e.g. 1.5e-3); it is read EXACTLY.  All numbers are put on one integer grid: value = I / S with S = 10^D.

VALIDITY (all integer arithmetic):
  count : exactly N centres, r > 0.
  pairs : (Xi - Xj)^2 + (Yi - Yj)^2 >= (2R)^2 for ALL pairs (plain O(N^2) loop over every pair, integers only).
  walls : RIGOROUS INTERVAL ENCLOSURES with explicit remainder bounds, adaptive precision.  At precision level G the checker
          holds integers with
              A_lo <= a 10^G <= A_hi,   NX_lo <= cos(t_j) 10^G <= NX_hi,   NY_lo <= sin(t_j) 10^G <= NY_hi,
          and bounds the scaled slack  S 10^G (a - n_j . c - r) = a10^G S - X nx10^G - Y ny10^G - R 10^G  by
              lower = A_lo S - max(X NX_lo, X NX_hi) - max(Y NY_lo, Y NY_hi) - R 10^G
              upper = A_hi S - min(X NX_lo, X NX_hi) - min(Y NY_lo, Y NY_hi) - R 10^G .
          lower >= 0 -> the circle is inside that side (proved);  upper < 0 -> it crosses it (proved, INVALID);
          otherwise the level is raised: G = 42, 82, 142, 202.  Still undecided at G = 202 -> INVALID (conservative).
    How the enclosures are made (integers and exact rationals only):
      pi   : Machin, pi = 16 atan(1/5) - 4 atan(1/239).  atan(1/b) = sum (-1)^n / ((2n+1) b^(2n+1)); every term is floored
             (for the lower sum) or ceiled (for the upper sum) on the grid 10^-W, W = G + 20; the series is alternating with
             strictly decreasing terms, so the omitted tail is bounded by the first omitted term, which is added (ceiled) to
             both sides.  -> integers PI_LO <= pi 10^W <= PI_HI.
      cos(pi u), sin(pi u), u rational: EXACT reduction to u in [0, 1/4] by  (c, s)(u + 2) = (c, s)(u);
             (c, s)(u) = (-c, -s)(u - 1);  (c, s)(u) = (-c, s)(1 - u);  (c, s)(u) = (s, c)(1/2 - u).
             Exact values (Niven): u = 0 -> (1, 0);  u = 1/6 -> sin = 1/2.  Otherwise, with x~ = PI_LO u / 10^W (rational,
             0 < x~ < 1): the Taylor polynomial T of degree d is evaluated EXACTLY (one integer numerator over the common
             integer denominator D^d d!), the Lagrange remainder is <= x~^(d+1)/(d+1)! <= 1/(d+1)!  with (d+1)! >= 10^(G+10),
             and the Lipschitz term |pi u - x~| <= (PI_HI - PI_LO) u / 10^W is added (|sin'|, |cos'| <= 1).  The enclosure
             is rounded OUTWARD to the grid 10^-G (floor / ceil) and clipped to [-1, 1].
          Every enclosure is at most 3 units of 10^-G wide, so for |x|, |y| <= 1 the slack interval is at most 9e-G wide.
    LOSS: the wall test is SOUND (an invalid packing is never accepted: acceptance needs a proved lower bound >= 0).  It is
          complete up to 9e-202: a packing all of whose wall slacks are >= 9e-202 is always accepted, i.e. the radius loss
          is < 1e-201 absolute.  Slacks below that (exact tangency to an irrational side cannot happen with decimal data
          except through cancellation, e.g. k = 3, 6) are reported INVALID "undecided".  The pair test is exact (zero loss).
    Supported k: every integer k >= 3.

VERDICT: IMPROVES iff valid at radius r and r > record_radius (both exact decimals); VALID_NOT_BETTER if valid and
r <= record; else INVALID <why>.

usage: py -3.11 certify_kgon.py <k> <file> <N> <record_radius>   -> prints "VERDICT: ..."
       py -3.11 certify_kgon.py --selftest                       -> unit tests + all published cxd / cpd packings + k = 5
                                                                    cross-check against certify_poly.py on all 200 cpt packings
       py -3.11 certify_kgon.py --published                      -> only the published cxd / cpd packings
"""
import sys, os, re, math, random
from fractions import Fraction as F

LEVELS = (42, 82, 142, 202)                     # G: enclosure grid 10^-G; undecided at the last level -> INVALID


# ---------------------------------------------------------------- rigorous enclosures (integers only)
def _atan_inv(b, W):
    """Integers (lo, hi) with lo <= atan(1/b) 10^W <= hi, b >= 2 an integer."""
    one = 10 ** W; lo = hi = 0; n = 0; bp = b                      # bp = b^(2n+1)
    while True:
        den = (2 * n + 1) * bp
        tf, tc = one // den, -((-one) // den)                      # floor / ceil of the term on the grid
        if tc <= 1 and n > 0:                                       # term < 1 ulp: it is the first omitted term
            return lo - tc, hi + tc
        if n % 2 == 0: lo += tf; hi += tc
        else: lo -= tc; hi -= tf
        n += 1; bp *= b * b


_PI = {}


def pi_enclosure(W):
    """Integers (PI_LO, PI_HI) with PI_LO <= pi 10^W <= PI_HI (Machin)."""
    if W not in _PI:
        l5, h5 = _atan_inv(5, W); l239, h239 = _atan_inv(239, W)
        _PI[W] = (16 * l5 - 4 * h239, 16 * h5 - 4 * l239)
    return _PI[W]


def _ceil_div(a, b):
    return -((-a) // b)


def _taylor(u, G, which):
    """Integers (lo, hi) with lo <= f(pi u) 10^G <= hi, f = sin or cos, rational u in (0, 1/4]."""
    W = G + 20; PL, PH = pi_enclosure(W)
    p, D = PL * u.numerator, 10 ** W * u.denominator              # x~ = p / D = PI_LO u / 10^W
    assert 0 < p < D, 'x~ must lie in (0, 1)'
    d = 1
    while math.factorial(d + 1) < 10 ** (G + 10): d += 1         # Lagrange remainder <= 1/(d+1)! <= 10^-(G+10)
    par = 1 if which == 'sin' else 0
    if d % 2 != par: d += 1                                        # last kept power has the right parity
    fd = math.factorial(d)
    # T = num / den exactly, den = D^d d!,  num = sum over kept powers m of  sign_m (d!/m!) p^m D^(d-m)
    num = 0
    for m in range(par, d + 1, 2):
        num += (-1) ** ((m - par) // 2) * (fd // math.factorial(m)) * p ** m * D ** (d - m)
    den = D ** d * fd
    tf = (num * 10 ** G) // den; tc = _ceil_div(num * 10 ** G, den)
    err = F(10 ** G, math.factorial(d + 1)) + F((PH - PL) * u.numerator * 10 ** G, u.denominator * 10 ** W)
    e = _ceil_div(err.numerator, err.denominator)
    one = 10 ** G
    return max(tf - e, -one), min(tc + e, one)


def cos_sin(u, G):
    """((clo, chi), (slo, shi)) integers: clo <= cos(pi u) 10^G <= chi, slo <= sin(pi u) 10^G <= shi, u rational."""
    u = F(u); u = u - 2 * math.floor(u / 2)                       # [0, 2)
    one = 10 ** G
    neg = lambda iv: (-iv[1], -iv[0])
    if u >= 1:
        c, s = cos_sin(u - 1, G); return neg(c), neg(s)
    if u > F(1, 2):
        c, s = cos_sin(1 - u, G); return neg(c), s
    if u > F(1, 4):
        c, s = cos_sin(F(1, 2) - u, G); return s, c
    if u == 0: return (one, one), (0, 0)
    if u == F(1, 6): return _taylor(u, G, 'cos'), (one // 2, one // 2)     # sin(pi/6) = 1/2 exactly (G >= 1)
    return _taylor(u, G, 'cos'), _taylor(u, G, 'sin')


_TAB = {}


def table(k, G):
    """(A_lo, A_hi), [(NX_lo, NX_hi, NY_lo, NY_hi)] for the k sides, at grid 10^-G."""
    key = (k, G)
    if key not in _TAB:
        A, _ = cos_sin(F(1, k), G)
        sides = []
        for j in range(k):
            c, s = cos_sin(F(3, 2) + F(2 * j, k), G)
            sides.append((c[0], c[1], s[0], s[1]))
        _TAB[key] = (A, sides)
    return _TAB[key]


def wall_status(k, S, R, X, Y):
    """None if the circle (X/S, Y/S; R/S) is proved inside every side; else (j, 'crosses' | 'undecided')."""
    todo = range(k)
    for G in LEVELS:
        (Alo, Ahi), sides = table(k, G); Rg = R * 10 ** G
        pend = []
        for j in todo:
            xl, xh, yl, yh = sides[j]
            a1, a2, b1, b2 = X * xl, X * xh, Y * yl, Y * yh
            lower = Alo * S - (a1 if a1 > a2 else a2) - (b1 if b1 > b2 else b2) - Rg
            if lower >= 0: continue                                  # proved inside this side
            upper = Ahi * S - (a1 if a1 < a2 else a2) - (b1 if b1 < b2 else b2) - Rg
            if upper < 0: return j, 'crosses'                        # proved outside this side
            pend.append(j)                                           # undecided at this level
        if not pend: return None
        todo = pend                                                  # raise the precision for the undecided sides only
    return todo[0], 'undecided'


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
def angle_str(j, k):
    """Normal angle of side j, 270 + 360 j / k mod 360 degrees, as an exact decimal (or p/q when it does not terminate)."""
    ang = (F(270) + F(360 * j, k)) % 360
    for places in range(7):
        if (ang * 10 ** places).denominator == 1:
            return fmt_dec(ang, places).rstrip('.') if places else str(int(ang))
    return f'{ang.numerator}/{ang.denominator}'


def check_values(k, r_s, pts_s, n, rec_s):
    """k-gon, radius string, list of (x, y) strings, expected count, record radius string -> verdict string."""
    if not (isinstance(k, int) and k >= 3): return f'INVALID unsupported k={k} (need an integer k >= 3)'
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
        bad = wall_status(k, S, R, X[i], Y[i])
        if bad is not None:
            j, why = bad
            return (f'INVALID circle {i + 1} {why} the side with normal at {angle_str(j, k)} deg'
                    + (' (undecided at 1e-200: treated as invalid)' if why == 'undecided' else ''))
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


# ---------------------------------------------------------------- published packings
ROOT = os.path.dirname(os.path.abspath(__file__))
TABLES = {'cxd': 16, 'cpd': 15}


def page_radii(code):
    s = open(os.path.join(ROOT, 'data', 'refs', f'packomania_{code}_2026-09-25.html'), encoding='utf-8', errors='replace').read()
    return {int(n): r for n, r in re.findall(r'name="' + code + r'(\d+)">.*?</a></td>\s*<td>([0-9.]+)</td>', s, re.S)}


def load_pub(code, n):
    pts = []
    for line in open(os.path.join(ROOT, 'data', 'big', code, f'{code}{n}.txt')):
        p = line.split()
        if not p: continue
        if len(p) == 3: pts.append((p[1], p[2]))
        elif len(p) == 1 and n == 1: pts.append(('0', '0'))     # the N = 1 file prints no coordinates: the centred circle
        else: raise ValueError(f'{code}{n}.txt: bad line {line!r}')
    return pts


def fmt_dec(q, places=40):
    """Fraction -> decimal string rounded toward zero to `places` digits (for building test radii)."""
    neg = q < 0; q = abs(q); m = (q.numerator * 10 ** places) // q.denominator
    s = str(m).rjust(places + 1, '0'); s = s[:-places] + '.' + s[-places:]
    return ('-' if neg else '') + s


def write_cert_file(path, r_s, pts):
    with open(path, 'w', newline='\n') as f:
        f.write(f'r {r_s}\n')
        for x, y in pts: f.write(f'{x} {y}\n')


def published_check(verbose=True, via_files=None):
    """Every published cxd / cpd packing at printed radius - off; must be VALID_NOT_BETTER at off = 2e-12 (record = printed).
    via_files: a folder -> the packings are first written as certificate files and read back through check_file."""
    out = {}
    for code, k in TABLES.items():
        T = page_radii(code); Ns = sorted(T)
        for off in ('0', '0.000000000001', '0.000000000002'):
            bad = []
            for n in Ns:
                pts = load_pub(code, n); rp = T[n]
                r = fmt_dec(dec_fraction(rp) - dec_fraction(off), 15)
                if via_files and off == '0.000000000002':
                    path = os.path.join(via_files, f'{code}_{n}.txt'); write_cert_file(path, r, pts)
                    v = check_file(k, path, n, rp)
                else:
                    v = check_values(k, r, pts, n, rp)
                if v != 'VALID_NOT_BETTER': bad.append((n, v))
            out[(code, off)] = bad
            if verbose:
                print(f'published {code} (k={k}), r = printed - {off}: {len(Ns) - len(bad)}/{len(Ns)} VALID_NOT_BETTER; not: '
                      + ((', '.join(f'N={n}: {v}' for n, v in bad[:8]) + (' ...' if len(bad) > 8 else '')) if bad else 'none'),
                      flush=True)
    return out


# ---------------------------------------------------------------- independent Decimal trig (test reference only)
def _dpi(prec):
    """pi to `prec` digits by the decimal-module documentation recipe (independent of the checker's Machin enclosure)."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec + 10
        three = Decimal(3); lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
        while s != lasts:
            lasts = s; n, na = n + na, na + 8; d, da = d + da, da + 32; t = (t * n) / d; s += t
        ctx.prec = prec; return +s                       # rounded to prec digits, NOT to the caller context


def _dcos_sin(x, prec):
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec + 10
        i, lasts, s, fact, num, sign = 0, 0, Decimal(1), 1, Decimal(1), 1
        while s != lasts:                                       # cos
            lasts = s; i += 2; fact *= i * (i - 1); num *= x * x; sign *= -1; s += num / fact * sign
        c = s
        i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
        while s != lasts:                                       # sin
            lasts = s; i += 2; fact *= i * (i - 1); num *= x * x; sign *= -1; s += num / fact * sign
        ctx.prec = prec; return +c, +s


# ---------------------------------------------------------------- self-tests
def selftest():
    import time, tempfile
    from decimal import Decimal, getcontext, localcontext
    getcontext().prec = 120
    t0 = time.time()

    def run(k, r, pts, rec='0'):
        return check_values(k, r if isinstance(r, str) else fmt_dec(r, 90), [tuple(p) for p in pts], len(pts), rec)

    # 0. the enclosures themselves: pi and cos(pi/k) against the independent Decimal references, at every level
    PI = _dpi(110)
    for G in LEVELS:
        W = G + 20; PL, PH = pi_enclosure(W); sh = max(0, W - 100)
        ref = int((PI * 10 ** (W - sh)).to_integral_value(rounding='ROUND_FLOOR'))   # floor(pi 10^min(W,100))
        assert 0 <= PH - PL <= 10 ** 4 and PL // 10 ** sh <= ref + 1 and PH // 10 ** sh >= ref, (G, PH - PL)   # width <= 1e-(G+16)
    for k in (3, 4, 5, 6, 7, 8, 12, 15, 16, 17, 24):
        for G in LEVELS[:2]:
            (lo, hi), sides = table(k, G)
            aD, _ = _dcos_sin(PI / k, 100)
            ref = int((aD * 10 ** G).to_integral_value(rounding='ROUND_FLOOR'))
            assert lo <= ref + 1 and hi >= ref and hi - lo <= 3, (k, G, lo, hi, ref)
            for j, (xl, xh, yl, yh) in enumerate(sides):
                cD, sD = _dcos_sin(PI * (Decimal(3) / 2 + Decimal(2 * j) / k), 100)
                assert xl - 1 <= cD * 10 ** G <= xh + 1 and yl - 1 <= sD * 10 ** G <= yh + 1, (k, j)
                assert xh - xl <= 3 and yh - yl <= 3
    # exact special values (Niven): axis normals and sin(pi/6) = 1/2 are exact degenerate intervals
    assert cos_sin(F(3, 2), 42) == ((0, 0), (-10 ** 42, -10 ** 42))
    assert cos_sin(F(1, 6), 42)[1] == (5 * 10 ** 41, 5 * 10 ** 41) and cos_sin(F(1, 3), 42)[0] == (5 * 10 ** 41, 5 * 10 ** 41)
    print(f'enclosure checks OK (pi enclosure contains the independent Decimal pi, width <= 1e-(G+16); cos/sin enclosures <= 3e-G wide and contain the independent Decimal values, 11 k)', flush=True)

    for k in (5, 15, 16):
        aD, _ = _dcos_sin(PI / k, 110)
        a = F(aD)                                                  # a Decimal-exact rational within ~1e-110 of cos(pi/k)
        e30, e40 = F(1, 10 ** 30), F(1, 10 ** 40)
        # 1. N = 1 at the centre, r = a -/+ 1e-30 and -/+ 1e-40, -/+ 1e-60
        for e in (e30, e40, F(1, 10 ** 60)):
            assert run(k, a - e, [('0', '0')]) == 'IMPROVES', (k, e)
            v = run(k, a + e, [('0', '0')]); assert v.startswith('INVALID circle 1 crosses'), (k, e, v)
        assert run(k, a - e30, [('0', '0')], fmt_dec(a - e30, 90)) == 'VALID_NOT_BETTER'     # r == record: not better
        # 2. one circle tangent to side j (centre accurate to ~1e-100), r = r0 -/+ 1e-30 and -/+ 1e-40
        r0 = Decimal('0.1')
        for j in range(k):
            cD, sD = _dcos_sin(PI * (Decimal(3) / 2 + Decimal(2 * j) / k), 110)
            t = aD - r0
            x, y = (t * cD).quantize(Decimal('1e-100')), (t * sD).quantize(Decimal('1e-100'))
            for e in ('1e-30', '1e-40'):
                v1 = run(k, format(r0 - Decimal(e), 'f'), [(str(x), str(y))]); v2 = run(k, format(r0 + Decimal(e), 'f'), [(str(x), str(y))])
                assert v1 == 'IMPROVES', (k, j, e, v1)
                assert v2.startswith('INVALID circle 1 crosses the side with normal at'), (k, j, e, v2)
                assert f'normal at {angle_str(j, k)} deg' in v2, (k, j, v2)
    print('N = 1 at the centre (r = a -/+ 1e-30, 1e-40, 1e-60) and a circle tangent to every side (r -/+ 1e-30, 1e-40): '
          'k = 5, 15, 16 OK', flush=True)
    # 3. two circles at exactly 2r (horizontal, and a 3-4-5 direction), and 1e-25 overlaps; overlap far apart in the list
    for k in (5, 15, 16):
        assert run(k, '0.1', [('-0.1', '0'), ('0.1', '0')]) == 'IMPROVES'
        assert run(k, '0.1', [('0.06', '0.08'), ('-0.06', '-0.08')]) == 'IMPROVES'
        assert run(k, '0.1', [('-0.1', '0'), ('0.0999999999999999999999999', '0')]) == 'INVALID overlap 1 2'
        assert run(k, '0.1', [('0.06', '0.08'), ('-0.06', '-0.0799999999999999999999999')]) == 'INVALID overlap 1 2'
        assert run(k, '0.1000000000000000000000001', [('-0.1', '0'), ('0.1', '0')]) == 'INVALID overlap 1 2'
        far = [('-0.3', '-0.3'), ('0.3', '0.3'), ('0', '0.4'), ('-0.3', '-0.2999999999999999999999999')]
        assert run(k, '0.05', far) == 'INVALID overlap 1 4'
        # 4. wrong count, format guards
        assert check_values(k, '0.1', [('0', '0')], 2, '0') == 'INVALID count 1 != 2'
        assert check_values(k, '0.1', [('0', '0'), ('0.5', '0')], 1, '0') == 'INVALID count 2 != 1'
        assert run(k, '0', [('0', '0')]) == 'INVALID r <= 0'
        assert run(k, '-0.1', [('0', '0')]) == 'INVALID r <= 0'
        for junk in ('nan', 'inf', '.', '-', '1/3', '0x1', '٣', '1e', '--1', '1.2.3', ''):
            assert run(k, '0.1', [('0', junk)]).startswith('INVALID not a decimal'), junk
        assert run(k, '1.0e-1', [('-1e-1', '0'), ('0.1', '0E0')]) == 'IMPROVES'
        assert run(k, '1.5e-1', [('-1e-1', '0'), ('0.1', '0E0')]) == 'INVALID overlap 1 2'
    assert run(2, '0.1', [('0', '0')]).startswith('INVALID unsupported k')
    # 5. orientation: k = 15 has a vertex at (0, 1) and a flat bottom; k = 16 is flat at top and bottom
    assert run(15, '0.01', [('0', '0.985')]) == 'IMPROVES'              # under the top vertex (0,1): inside
    assert run(15, '0.01', [('0', '-0.985')]).startswith('INVALID')     # below the flat bottom (y = -0.97815)
    assert run(16, '0.01', [('0', '0.975')]).startswith('INVALID') and run(16, '0.01', [('0', '-0.975')]).startswith('INVALID')
    assert run(16, '0.01', [('0', '0.97')]) == 'IMPROVES' and run(16, '0.01', [('0', '-0.97')]) == 'IMPROVES'
    # 6. exact rational slacks (Niven values) are decided exactly: k = 3, a = 1/2 exactly, N = 1 at the centre
    assert run(3, '0.5', [('0', '0')]) == 'IMPROVES'                   # exact tangency to all three sides
    assert run(3, '0.5' + '0' * 300 + '1', [('0', '0')]).startswith('INVALID circle 1 crosses')
    #    below the resolution: a true slack of +1e-250 cannot be proved >= 0 at 1e-200 -> INVALID 'undecided' (conservative,
    #    never VALID by accident); -1e-250 -> INVALID as well.
    PI3 = _dpi(320)
    for k in (15, 16):
        with localcontext() as ctx:
            ctx.prec = 340; x3 = PI3 / k                              # the argument itself must carry 320+ digits
        a3, _ = _dcos_sin(x3, 320)
        for sg in (-1, 1):
            v = run(k, fmt_dec(F(a3) + sg * F(1, 10 ** 250), 300), [('0', '0')])
            assert v.startswith('INVALID circle 1'), (k, sg, v)
            if sg < 0: assert 'undecided' in v, v
        #    ... while +/- 1e-190 is resolved both ways at the last level (G = 202)
        assert run(k, fmt_dec(F(a3) - F(1, 10 ** 190), 300), [('0', '0')]) == 'IMPROVES', k
        assert run(k, fmt_dec(F(a3) + F(1, 10 ** 190), 300), [('0', '0')]).startswith('INVALID circle 1 crosses'), k
    print('pairs, count, format, orientation, exact-value and resolution tests OK (slack 1e-190 decided both ways; '
          '1e-250 -> INVALID undecided)', flush=True)
    # 7. float cross-check of the side tables for many k (float is only the TEST reference): wherever the float slack is
    #    clearly nonzero, the proved sign must agree.
    rng = random.Random(20260925); agree = 0
    for k in (3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 16, 17, 20):
        a = math.cos(math.pi / k)
        for _ in range(400):
            S = 10 ** 9; X, Y, R = rng.randrange(-S, S), rng.randrange(-S, S), rng.randrange(1, S // 2)
            fs = min(a - R / S - (math.cos(math.radians(270 + 360 * j / k)) * X / S + math.sin(math.radians(270 + 360 * j / k)) * Y / S)
                     for j in range(k))
            if abs(fs) < 1e-9: continue
            st = wall_status(k, S, R, X, Y)
            assert (st is None) == (fs > 0), (k, X, Y, R, fs, st)
            agree += 1
    print(f'float cross-check (13 k, random circles): {agree} verdicts agree, 0 disagree', flush=True)
    # 8. k = 3 and 5 against certify_poly.py (exact field arithmetic; same frame for odd k) on random rational circles
    sys.path.insert(0, ROOT)
    import certify_poly as cp
    agree = 0
    for k in (3, 5):
        for _ in range(3000):
            S = 10 ** 9; X, Y, R = rng.randrange(-S, S), rng.randrange(-S, S), rng.randrange(1, S // 2)
            assert (wall_status(k, S, R, X, Y) is None) == (cp.wall_ok(k, S, R, X, Y) is None), (k, X, Y, R)
            agree += 1
    print(f'k = 3, 5 vs certify_poly.py (random rational circles): {agree} agree, 0 disagree', flush=True)
    # 9. k = 5 against certify_poly.py on all 200 published cpt packings, at printed - 2e-12, printed, printed + 1e-11 and
    #    mirrored (y -> -y, the wrong orientation): verdict categories must be identical.
    Tc = cp.page_radii(); same = diff = 0; cats = {}
    for n in range(1, 201):
        pts = cp.load_pub(n)
        mir = [(x, y[1:] if y.startswith('-') else '-' + y) for x, y in pts]
        for off in (F(-2, 10 ** 12), F(0), F(1, 10 ** 11)):
            r = fmt_dec(dec_fraction(Tc[n]) + off, 15)
            for P in (pts, mir):
                v1 = check_values(5, r, P, n, Tc[n]); v2 = cp.check_values(5, r, P, n, Tc[n])
                c1, c2 = v1.split()[0], v2.split()[0]
                if c1 == c2 and (c1 != 'INVALID' or v1.split()[1] == v2.split()[1]): same += 1
                else: diff += 1; print('DISAGREE', n, off, v1, '|', v2)
                cats[c1] = cats.get(c1, 0) + 1
    print(f'k = 5 vs certify_poly.py on the 200 published cpt packings (3 radii x 2 orientations = 1200 checks): '
          f'{same} identical, {diff} different; verdicts {cats}', flush=True)
    # 10. every published cxd / cpd packing, via certificate files at printed - 2e-12
    with tempfile.TemporaryDirectory(prefix='kgon_selftest_') as tmp:        # removed afterwards
        res = published_check(via_files=tmp)
    bad = res[('cxd', '0.000000000002')] + res[('cpd', '0.000000000002')]
    # 11. negative control: the same packings rotated by 180/k deg (the other orientation) must be INVALID, except packings
    #     that are themselves symmetric under that rotation (a mirror for odd k); float check_geom found cpd N = 7, 91.
    still = {}
    for code, k in TABLES.items():
        T = page_radii(code); PIk = _dpi(60); cD, sD = _dcos_sin(PIk / k, 60); still[code] = []
        for n in sorted(T):
            if n == 1: continue
            P = []
            for x, y in load_pub(code, n):
                xd, yd = Decimal(x), Decimal(y)
                P.append((str((cD * xd - sD * yd).quantize(Decimal('1e-30'))), str((sD * xd + cD * yd).quantize(Decimal('1e-30')))))
            r = fmt_dec(dec_fraction(T[n]) - F(2, 10 ** 12), 15)
            if not check_values(k, r, P, n, T[n]).startswith('INVALID'): still[code].append(n)
        print(f'negative control {code}: rotated by 180/{k} deg, N = 2-{max(T)}: {max(T) - 1 - len(still[code])}/{max(T) - 1} '
              f'INVALID; still valid: {still[code] or "none"}', flush=True)
    ok = not bad and diff == 0 and not still['cxd'] and set(still['cpd']) <= {7, 91}
    print(f'SELFTEST {"OK" if ok else "FAILED"}  ({time.time() - t0:.0f} s)' +
          ('' if ok else f'  published not VALID_NOT_BETTER at printed - 2e-12: {bad}; cpt disagreements {diff}; rotated still valid {still}'))
    return ok


if __name__ == '__main__':
    if sys.argv[1:2] == ['--selftest']:
        sys.exit(0 if selftest() else 1)
    elif sys.argv[1:2] == ['--published']:
        r = published_check()
        sys.exit(0 if not (r[('cxd', '0.000000000002')] or r[('cpd', '0.000000000002')]) else 1)
    elif len(sys.argv) == 5:
        try:
            k = int(sys.argv[1]); n = int(sys.argv[3])
        except ValueError:
            print('VERDICT: INVALID bad k or N argument'); sys.exit(0)
        print('VERDICT:', check_file(k, sys.argv[2], n, sys.argv[4]))
    else:
        print(__doc__); sys.exit(2)
