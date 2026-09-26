# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Exact checker A for equal d-dimensional balls in the unit d-BALL (Packomania hsp4 / hsp5 / hsp6), any dimension d.
Authors: Moki&Julio (2026-09-25). Python 3.11, standard library only. No floating point number takes part in any decision.

FRAME: container = the unit ball centred at the origin in R^d. A ball of radius r centred at c is inside iff
       |c|^2 <= (1 - r)^2 and r <= 1;  balls i != j do not overlap iff |c_i - c_j|^2 >= (2r)^2.

FILE: line 1 "r <radius>", then exactly N lines of exactly d numbers (blank lines are ignored). Every number must be a plain
decimal, optionally with an exponent (e.g. 1.5e-17); anything else (nan, inf, fractions, hex, non-ASCII digits) is refused.
The numbers are read EXACTLY: value = m * 10^e with integers m, e, and all of them are put on ONE integer grid
value = I / S with S = 10^D, D = the largest number of decimal places needed. Every test below is on these integers.

VALIDITY (all exact, all integers):
  count  : exactly N centres, each with exactly d coordinates; 0 < r <= 1 (0 < R <= S).
  wall   : for every centre, sum_k X_k^2 <= (S - R)^2          (|c|^2 <= (1 - r)^2, times S^2; S - R >= 0).
  pairs  : for EVERY pair i < j (plain O(N^2 d) loop; N <= a few thousand is instant), sum_k (X_ik - X_jk)^2 >= (2R)^2.
           Tangency (equality) is valid in both tests.
VERDICT: IMPROVES iff valid at radius r and r > record_radius (both exact decimals); VALID_NOT_BETTER if valid and
         r <= record; else INVALID <why>. LOSS IN RADIUS: zero (the tests are exact, so VALID means the balls fit).

usage: python checkers/certify_ball.py <d> <file> <N> <record_radius>   -> prints "VERDICT: ..."
       python checkers/certify_ball.py --selftest   -> unit tests (d = 2..6) + every published hsp4 file (data/big/hsp4/coords):
            VALID_NOT_BETTER at printed - 3e-12, INVALID at printed + 1e-9, the exact deficit of each 12-decimal file at the
            printed radius, negative controls at printed - 3e-12 (closest pair pushed together by 1e-10 of its separation; outermost centre
            pushed out by 1e-10 of its norm), and the truncated hsp5 / hsp6 files refused at their own dimension.
"""
import sys, os, re, random
from fractions import Fraction as F
from math import isqrt

# ---------------------------------------------------------------- exact decimal parsing
DEC = re.compile(r'^([+-]?)([0-9]+)?(?:\.([0-9]*))?(?:[eE]([+-]?[0-9]+))?$', re.ASCII)


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
def grid(r_s, pts_s):
    """radius string + coordinate strings -> (R, [[X_ik]], S) on one exact integer grid (value = I / S)."""
    vals = [dec(r_s)] + [dec(t) for p in pts_s for t in p]
    D = max(0, max(-e for _, e in vals)); S = 10 ** D
    I = [m * 10 ** (e + D) for m, e in vals]
    d = len(pts_s[0]) if pts_s else 0
    return I[0], [I[1 + k * d: 1 + (k + 1) * d] for k in range(len(pts_s))], S


def check_values(d, r_s, pts_s, n, rec_s):
    """dimension, radius string, list of d-tuples of strings, expected count, record radius string -> verdict string."""
    if len(pts_s) != n: return f'INVALID count {len(pts_s)} != {n}'
    for k, p in enumerate(pts_s):
        if len(p) != d: return f'INVALID centre {k + 1} has {len(p)} coordinates, expected {d}'
    try:
        R, P, S = grid(r_s, pts_s)
        rec = dec_fraction(rec_s)
    except ValueError as ex:
        return f'INVALID {ex}'
    if R <= 0: return 'INVALID r <= 0'
    if R > S: return 'INVALID r > 1'
    W = (S - R) * (S - R)
    for i, x in enumerate(P):
        if sum(v * v for v in x) > W: return f'INVALID ball {i + 1} crosses the container wall'
    R4 = 4 * R * R
    for i in range(n):
        xi = P[i]
        for j in range(i + 1, n):
            xj = P[j]; s = 0
            for k in range(d):
                t = xi[k] - xj[k]; s += t * t
            if s < R4: return f'INVALID overlap {i + 1} {j + 1}'
    return 'IMPROVES' if F(R, S) > rec else 'VALID_NOT_BETTER'


def read_cert(path, d):
    L = [l.split() for l in open(path, encoding='ascii') if l.strip()]
    if not L or L[0][0] != 'r' or len(L[0]) != 2: raise ValueError('header must be "r <radius>"')
    for k, t in enumerate(L[1:]):
        if len(t) != d: raise ValueError(f'centre line {k + 1} has {len(t)} numbers, expected {d}')
    return L[0][1], [tuple(t) for t in L[1:]]


def check_file(d, path, n, rec_s):
    try:
        r_s, pts = read_cert(path, d)
    except (ValueError, OSError, UnicodeDecodeError) as ex:
        return f'INVALID {ex}'
    return check_values(d, r_s, pts, n, rec_s)


# ---------------------------------------------------------------- exact min-radius (for the published-file report)
def exact_rmin_bounds(pts_s, K=10 ** 30):
    """(lower, upper) Fractions bracketing the TRUE min-radius min(1 - |c_i|, min |c_i - c_j| / 2) of the exact
    decimals, to within 2 / K (integer square roots)."""
    _, P, S = grid('1', pts_s); n = len(P); d = len(P[0])
    q = max(sum(v * v for v in x) for x in P)                        # max |c|^2 * S^2
    nlo, nhi = F(isqrt(q * K * K // (S * S)), K), F(isqrt(-(-q * K * K // (S * S))) + 1, K)
    wall_lo, wall_hi = 1 - nhi, 1 - nlo
    best = None
    for i in range(n):
        for j in range(i + 1, n):
            s = sum((P[i][k] - P[j][k]) ** 2 for k in range(d))
            if best is None or s < best: best = s
    if best is None: return wall_lo, wall_hi
    plo = F(isqrt(best * K * K // (S * S)), 2 * K); phi = F(isqrt(-(-best * K * K // (S * S))) + 1, 2 * K)
    return min(wall_lo, plo), min(wall_hi, phi)


# ---------------------------------------------------------------- published hsp4 files
ROOT = os.path.dirname(os.path.abspath(__file__))
PUB = os.path.join(ROOT, 'data', 'big', 'hsp{d}', 'coords')
PAGE = os.path.join(ROOT, 'data', 'refs', 'packomania_hsp{d}_2026-09-25.html')


def page_radii(d):
    T = {}
    for line in open(PAGE.format(d=d), encoding='utf-8', errors='replace'):
        if not line.startswith('<tr><td'): continue
        cells = [re.sub(r'<[^>]+>', '', c).strip() for c in re.findall(r'<td[^>]*>(.*?)</td>', line)]
        if len(cells) >= 2 and cells[0].isdigit(): T[int(cells[0])] = cells[1]
    return T


def load_pub(d, n):
    """(radius string of line 1, [coordinate tuples]) from hsp<d>-<N>.txt ('idx x1 .. x4' lines, as served)."""
    L = [l.split() for l in open(os.path.join(PUB.format(d=d), f'hsp{d}-{n}.txt')) if l.strip()]
    return L[0][0], [tuple(t[1:]) for t in L[1:]]


def fmt_dec(q, places=40):
    """Fraction -> decimal string rounded toward zero to `places` digits."""
    neg = q < 0; q = abs(q); m = (q.numerator * 10 ** places) // q.denominator
    s = str(m).rjust(places + 1, '0'); s = s[:-places] + '.' + s[-places:]
    return ('-' if neg else '') + s


def published_check(verbose=True):
    """Every hsp4 file: VALID_NOT_BETTER at printed - 3e-12 (record = printed), INVALID at printed + 1e-9; the exact deficit
    printed - rmin(file) (positive = the 12-decimal file does not reach its printed radius exactly)."""
    T = page_radii(4); d = 4
    Ns = sorted(int(m.group(1)) for f in os.listdir(PUB.format(d=4)) for m in [re.fullmatch(r'hsp4-(\d+)\.txt', f)] if m)
    lo_e, hi_e = F(3, 10 ** 12), F(1, 10 ** 9); fails = []; ok = 0; line1_mismatch = []
    at_printed_invalid = []; worst = (F(-1), None); deficits = {}; sq_worst = [(F(-1), None)]
    for n in Ns:
        r1, pts = load_pub(4, n)
        if r1 != T[n]: line1_mismatch.append(n)
        rp = dec_fraction(T[n])
        v1 = check_values(d, fmt_dec(rp - lo_e), pts, n, T[n])
        v2 = check_values(d, fmt_dec(rp + hi_e), pts, n, T[n])
        v0 = check_values(d, T[n], pts, n, T[n])
        if v1 == 'VALID_NOT_BETTER' and v2.startswith('INVALID'): ok += 1
        else: fails.append((n, v1, v2))
        if v0.startswith('INVALID'): at_printed_invalid.append(n)
        lo, hi = exact_rmin_bounds(pts); deficit = rp - lo; deficits[n] = float(deficit)
        R_, P_, S_ = grid(T[n], pts)                                  # violations in SQUARED terms at the printed radius
        wv = max(sum(v * v for v in x) for x in P_) - (S_ - R_) ** 2
        pv = 4 * R_ * R_ - min(sum((P_[i][k] - P_[j][k]) ** 2 for k in range(d)) for i in range(n) for j in range(i + 1, n))
        sq_worst[0] = max(sq_worst[0], (F(max(wv, pv), S_ * S_), n))
        if deficit > worst[0]: worst = (deficit, n)
    top = sorted(deficits, key=lambda k: -deficits[k])[:3]
    if verbose:
        print(f'published hsp4 files: {len(Ns)} (N = {Ns[0]}..{Ns[-1]}); line-1 radius == page radius in {len(Ns) - len(line1_mismatch)} '
              f'(mismatch: {line1_mismatch or "none"}); VALID_NOT_BETTER at printed - 3e-12 and INVALID at printed + 1e-9: {ok}; '
              f'failures: {fails[:10] or "none"}')
        print(f'  at the printed radius itself: {len(at_printed_invalid)} files INVALID (12-decimal rounding), {len(Ns) - len(at_printed_invalid)} VALID; '
              f'largest exact deficit printed - rmin(file) = {float(worst[0]):.3e} at N = {worst[1]} '
              f'(top 3: {", ".join(f"N={k} {deficits[k]:.2e}" for k in top)}); largest violation in SQUARED terms '
              f'(|c|^2 - (1 - r)^2 or 4r^2 - |ci - cj|^2) = {float(sq_worst[0][0]):.3e} at N = {sq_worst[0][1]}')
    return Ns, fails, dict(n_files=len(Ns), ok=ok, fails=fails, invalid_at_printed=len(at_printed_invalid),
                            worst_deficit=float(worst[0]), worst_N=worst[1], top3={k: deficits[k] for k in top},
                            worst_squared_violation=float(sq_worst[0][0]), worst_squared_N=sq_worst[0][1])


# ---------------------------------------------------------------- self-tests
def selftest():
    e30 = F(1, 10 ** 30)
    def run(d, r, pts, rec='0'):
        return check_values(d, r if isinstance(r, str) else fmt_dec(r, 70), [tuple(p) for p in pts], len(pts), rec)
    def z(d): return ['0'] * d
    for d in (2, 3, 4, 5, 6):
        # 1. one ball at the centre: r = 1 valid; 1 + 1e-30 invalid; 1 - 1e-30 valid
        assert run(d, '1', [z(d)]) == 'IMPROVES'
        assert run(d, F(1) + e30, [z(d)]) == 'INVALID r > 1'
        assert run(d, F(1) - e30, [z(d)]) == 'IMPROVES'
        # 2. two balls at +-(1/2) e_1 with r = 1/2: tangent to each other and to the wall; +1e-30 -> invalid (the wall)
        p, q = z(d), z(d); p[0], q[0] = '0.5', '-0.5'
        assert run(d, '0.5', [p, q]) == 'IMPROVES', d
        assert run(d, F(1, 2) + e30, [p, q]) == 'INVALID ball 1 crosses the container wall', d
        # two balls at +-(1/4) e_k: pair tangency only; +1e-30 -> overlap
        for k in range(d):
            p, q = z(d), z(d); p[k], q[k] = '0.25', '-0.25'
            assert run(d, '0.25', [p, q]) == 'IMPROVES'
            assert run(d, F(1, 4) + e30, [p, q]) == 'INVALID overlap 1 2'
        # 3. wall tangency along an exact rational unit vector (3/5, 4/5) and (2, 3, 6) / 7, r = 1/5; 1e-30 outside -> invalid
        for u, den in (((3, 4), 5), ((2, 3, 6), 7)):
            if len(u) > d: continue
            r = F(1, 5); L = 1 - r
            for sg in (1, -1):
                p = z(d)
                for k, t in enumerate(u): p[k] = fmt_dec(sg * L * t / den, 70)
                assert run(d, '0.2', [p]) == 'IMPROVES', (d, u)
                p2 = list(p); p2[0] = fmt_dec(sg * (L * u[0] / den + e30), 70)
                assert run(d, '0.2', [p2]) == 'INVALID ball 1 crosses the container wall', (d, u)
                p3 = list(p); p3[0] = fmt_dec(sg * (L * u[0] / den - e30), 70)
                assert run(d, '0.2', [p3]) == 'IMPROVES'
        # 4. format guards, count, width, record comparison, exponents
        for junk in ('nan', 'inf', '.', '-', '1/3', '0x1', '٣', '1e', '--1', '1.2.3', '1,5', ' '):
            pp = z(d); pp[-1] = junk
            assert run(d, '0.1', [pp]).startswith('INVALID not a decimal'), junk
        assert check_values(d, '0.1', [tuple(z(d))], 2, '0') == 'INVALID count 1 != 2'
        assert check_values(d, '0.1', [tuple(z(d - 1))], 1, '0') == f'INVALID centre 1 has {d - 1} coordinates, expected {d}'
        assert check_values(d, '0.1', [tuple(z(d) + ['0'])], 1, '0') == f'INVALID centre 1 has {d + 1} coordinates, expected {d}'
        assert run(d, '0', [z(d)]) == 'INVALID r <= 0'
        assert run(d, '-0.1', [z(d)]) == 'INVALID r <= 0'
        assert run(d, '0.1', [z(d)], '0.1') == 'VALID_NOT_BETTER'
        assert run(d, '0.1', [z(d)], '0.0999999999999999999999999999999999999999') == 'IMPROVES'
        assert run(d, '0.1', [z(d)], '0.1000000000000000000000000000000000000001') == 'VALID_NOT_BETTER'
        assert run(d, '0.1', [z(d)], '1e-1') == 'VALID_NOT_BETTER'
        p, q = z(d), z(d); p[0], q[0] = '-2.5e-1', '0.25E0'
        assert run(d, '2.5e-1', [p, q]) == 'IMPROVES'
        assert run(d, '2.5000000000000000000000000000001e-1', [p, q]) == 'INVALID overlap 1 2'
    # 5. the 24-cell (d = 4): (+-1, +-1, 0, 0) / sqrt 2 scaled to 1 - r with r = 1/3 - 1e-12, written rounded INWARD
    #    (towards zero, 40 digits, so the vectors are slightly short), + one ball at the centre -> valid; one centre moved
    #    1e-25 outside -> invalid; one line with 3 numbers -> invalid; wrong count -> invalid
    r = F(1, 3) - F(1, 10 ** 12); L = 1 - r
    # L / sqrt 2 rounded toward zero at 40 digits: isqrt of (L^2 / 2) * 10^80
    a = F(isqrt(L.numerator ** 2 * 10 ** 80 // (2 * L.denominator ** 2)), 10 ** 40); A = fmt_dec(a, 40)
    V = []
    for i in range(4):
        for j in range(i + 1, 4):
            for s1 in (A, '-' + A):
                for s2 in (A, '-' + A):
                    v = ['0'] * 4; v[i] = s1; v[j] = s2; V.append(v)
    rs = fmt_dec(r, 40)
    assert run(4, rs, V + [['0'] * 4]) == 'IMPROVES'
    assert run(4, rs, V + [['0'] * 4], '0.3333333333') == 'IMPROVES'
    assert run(4, rs, V + [['0'] * 4], '0.333333333333') == 'VALID_NOT_BETTER'           # 1/3 - 1e-12 < 0.333333333333
    W2 = [list(v) for v in V]
    k = next(i for i, t in enumerate(V[5]) if t != '0'); sg = -1 if V[5][k].startswith('-') else 1
    W2[5][k] = fmt_dec(sg * (a + F(1, 10 ** 25)), 40)
    assert run(4, rs, W2 + [['0'] * 4]) == 'INVALID ball 6 crosses the container wall', run(4, rs, W2 + [['0'] * 4])
    assert check_values(4, rs, [tuple(v) for v in V[:-1]] + [tuple(V[-1][:3])], 24, '0').startswith('INVALID centre 24 has 3')
    assert check_values(4, rs, [tuple(v) for v in V], 25, '0') == 'INVALID count 24 != 25'
    # at r = 1/3 exactly with the inward-rounded vectors the centre ball overlaps (|c| < 2/3): the rounding is seen
    assert run(4, '0.3333333333333333333333333333333333333333333333', V + [['0'] * 4]).startswith('INVALID')
    # 6. random exact cross-check against Fractions (d = 2..6): the integer-grid verdict equals a Fraction recomputation
    rng = random.Random(20260925); agree = 0
    for trial in range(300):
        d = rng.randrange(2, 7); n = rng.randrange(1, 12); Sx = 10 ** rng.randrange(2, 20)
        P = [[F(rng.randrange(-Sx, Sx + 1), d * Sx) for _ in range(d)] for _ in range(n)]      # |p| <= 1 / sqrt(d)
        walls = [1 - F(isqrt(int(sum(x * x for x in pp) * 10 ** 40)) + 1, 10 ** 20) for pp in P]
        pairs = [F(isqrt(int(sum((pp[k] - qq[k]) ** 2 for k in range(d)) * 10 ** 40)), 2 * 10 ** 20)
                 for i, pp in enumerate(P) for qq in P[i + 1:]]
        dm = min(walls + pairs)
        for R in (dm - F(1, 10 ** 15), dm + F(1, 10 ** 15), dm * 2):
            if R <= 0: continue
            rs_ = fmt_dec(R, 30); Rq = dec_fraction(rs_)
            pts = [tuple(fmt_dec(x, 30) for x in p) for p in P]; Pq = [[dec_fraction(t) for t in p] for p in pts]
            ok_ = Rq <= 1 and all(sum(x * x for x in p) <= (1 - Rq) ** 2 for p in Pq) and \
                all(sum((p[k] - q[k]) ** 2 for k in range(d)) >= 4 * Rq * Rq for i, p in enumerate(Pq) for q in Pq[i + 1:])
            v = check_values(d, rs_, pts, n, '0')
            assert (v == 'IMPROVES') == ok_, (trial, v, ok_); agree += 1
    print(f'random exact cross-check (integer grid vs Fractions, d = 2..6): {agree} verdicts agree')
    print('unit self-tests OK (d = 2..6)')
    # 7. published hsp4 files
    Ns, fails, pubstats = published_check()
    T = page_radii(4); neg_bad = []; wall_bad = []
    for n in Ns:
        _, pts = load_pub(4, n); rp = dec_fraction(T[n]); r_s = fmt_dec(rp - F(3, 10 ** 12))
        R, P, S = grid('1', pts)
        if n >= 2:
            best = min(((sum((P[i][k] - P[j][k]) ** 2 for k in range(4)), i, j) for i in range(n) for j in range(i + 1, n)))
            _, i, j = best; t = F(1, 10 ** 10); moved = list(pts)
            moved[i] = tuple(fmt_dec(F(P[i][k], S) + t * F(P[j][k] - P[i][k], S), 40) for k in range(4))
            v = check_values(4, r_s, moved, n, T[n])
            if not v.startswith('INVALID overlap'): neg_bad.append((n, v))
        io = max(range(n), key=lambda i: sum(x * x for x in P[i])); moved = list(pts)
        moved[io] = tuple(fmt_dec(F(P[io][k], S) * (1 + F(1, 10 ** 10)), 40) for k in range(4))
        v = check_values(4, r_s, moved, n, T[n])
        if v != f'INVALID ball {io + 1} crosses the container wall': wall_bad.append((n, v))
    print(f'negative controls on {len(Ns)} published files (at printed - 3e-12): closest pair pushed together by 1e-10 of its '
          f'separation -> INVALID overlap in {len(Ns) - len(neg_bad)}/{len(Ns)} (not: {neg_bad[:5] or "none"}); outermost centre '
          f'pushed out by 1e-10 of its norm -> INVALID wall in {len(Ns) - len(wall_bad)}/{len(Ns)} (not: {wall_bad[:5] or "none"})')
    # 8. hsp5 / hsp6 published files are truncated to 4 columns: refused at their own dimension
    trunc = {}
    for d in (5, 6):
        T5 = page_radii(d); bad = 0; tot = 0
        for n in (10, 50, 100, 200):
            _, pts = load_pub(d, n); tot += 1
            v = check_values(d, fmt_dec(dec_fraction(T5[n]) - F(3, 10 ** 12)), pts, n, T5[n])
            bad += v.startswith('INVALID centre 1 has 4 coordinates')
        trunc[d] = (bad, tot)
    print(f'truncated hsp5 / hsp6 files refused at d = 5 / 6: {trunc}')
    ok = not fails and not neg_bad and not wall_bad and all(b == t for b, t in trunc.values())
    print('SELFTEST ' + ('OK' if ok else f'FAILED: {fails[:5]} {neg_bad[:5]} {wall_bad[:5]} {trunc}'))
    return ok, pubstats


if __name__ == '__main__':
    args = sys.argv[1:]
    if args[:1] == ['--selftest']:
        sys.exit(0 if selftest()[0] else 1)
    elif len(args) == 4:
        print('VERDICT:', check_file(int(args[0]), args[1], int(args[2]), args[3]))
    else:
        print(__doc__); sys.exit(2)
