#!/usr/bin/env python3
# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Exact equal-circle checker for a regular k-gon, circumradius 1.

Frame: centre at the origin, one side horizontal at the bottom.
Side j has outward unit normal (cos t_j, sin t_j),
t_j = 270 deg + 360 deg * j / k, and sits at distance a = cos(pi/k).

No binary floating point is used in any decision. Decisions are Fraction
comparisons on exact pair gaps and on rational enclosures of wall slacks.

Usage:
  python verify_exact_kgon.py <k> <file> <N> <record_radius>
  python verify_exact_kgon.py --selftest
"""

from __future__ import annotations

import io
import os
import sys
import tempfile
from contextlib import redirect_stdout
from decimal import Decimal, InvalidOperation
from fractions import Fraction


# Stop and refuse rather than guess once a slack enclosure is this narrow.
SLACK_LIMIT = Fraction(1, 10**200)

_PI = None  # (lo, hi) decimal bracket, lo <= pi <= hi
_TRIG = {}  # reduced angle f -> (width, cos_lo, cos_hi, sin_lo, sin_hi)


def as_fraction(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return Fraction(value)
    raise TypeError('not an exact rational')


def parse_number(token):
    """Decimal string, including scientific notation, to an exact Fraction."""
    try:
        dec = Decimal(str(token).strip())
    except InvalidOperation as exc:
        raise ValueError(f'bad number {token}') from exc
    if not dec.is_finite():
        raise ValueError(f'non-finite {token}')
    return Fraction(dec)


def floor_to_scale(x, scale):
    """Greatest multiple of 1/scale that is <= x."""
    return Fraction((x.numerator * scale) // x.denominator, scale)


def ceil_to_scale(x, scale):
    """Least multiple of 1/scale that is >= x."""
    num = x.numerator * scale
    den = x.denominator
    return Fraction(-((-num) // den), scale)


def mod2(q):
    """q - 2*floor(q/2), in [0, 2). Angle q*pi reduced by full turns."""
    cycles = q.numerator // (q.denominator * 2)
    f = q - 2 * cycles
    if f < 0 or f >= 2:
        raise RuntimeError('mod2 failed')
    return f


def atan_bounds(x, n_terms):
    """Rigorous bracket of arctan(x) for 0 <= x <= 1.

    The series is alternating and termwise decreasing, so the remainder after
    n_terms has the sign of the first omitted term and absolute value at most
    that term.
    """
    if x < 0 or x > 1:
        raise RuntimeError('atan domain')
    if n_terms < 1:
        raise RuntimeError('atan terms')
    acc = Fraction(0)
    mag = x
    x2 = x * x
    for n in range(n_terms):
        if n % 2 == 0:
            acc += mag
        else:
            acc -= mag
        mag = mag * x2 * Fraction(2 * n + 1, 2 * n + 3)
    if mag < 0:
        raise RuntimeError('atan term sign')
    if n_terms % 2 == 0:
        return acc, acc + mag
    return acc - mag, acc


def machin_pi(max_width):
    """pi = 16*arctan(1/5) - 4*arctan(1/239), width at most max_width."""
    x1 = Fraction(1, 5)
    x2 = Fraction(1, 239)
    n = 2
    while n <= 20000:
        a_lo, a_hi = atan_bounds(x1, n)
        b_lo, b_hi = atan_bounds(x2, n)
        lo = 16 * a_lo - 4 * b_hi
        hi = 16 * a_hi - 4 * b_lo
        if lo > hi:
            raise RuntimeError('pi bracket crossed')
        if hi - lo <= max_width and Fraction(3) < lo and hi < Fraction(4):
            return lo, hi
        n *= 2
    raise RuntimeError('pi precision exhausted')


def pi_bounds(max_width):
    """Decimal rational bracket of pi with width <= max_width.

    Machin produces an exact bracket; flooring and ceiling onto a short
    decimal scale keeps later powers small and still contains pi.
    """
    global _PI
    if max_width <= 0:
        raise RuntimeError('pi width')
    if _PI is not None and _PI[1] - _PI[0] <= max_width:
        return _PI
    need = 4 / max_width  # 10^places >= 4/max_width => snap adds <= max_width/2
    power = 1
    while power < need:
        power *= 10
    lo_m, hi_m = machin_pi(max_width / 2)
    lo = floor_to_scale(lo_m, power)
    hi = ceil_to_scale(hi_m, power)
    if not (lo <= lo_m <= hi_m <= hi):
        raise RuntimeError('pi snap escaped the series bracket')
    if hi - lo > max_width or not (Fraction(3) < lo and hi < Fraction(4)):
        raise RuntimeError('pi snap too wide')
    _PI = (lo, hi)
    return _PI


def pow_over_fact(x, m):
    """x**m / m! for x >= 0, as an exact Fraction."""
    if x < 0 or m < 0:
        raise RuntimeError('pow_over_fact')
    acc = Fraction(1)
    for i in range(1, m + 1):
        acc *= x
        acc /= i
    return acc


def taylor_sin_cos(t, n_terms):
    """Taylor partial sums of sin and cos through sin degree 2N+1 and cos degree 2N."""
    sin_acc = Fraction(0)
    cos_acc = Fraction(0)
    power = Fraction(1)
    fact = 1
    for n in range(n_terms + 1):
        term = power / fact
        if n % 2 == 0:
            cos_acc += term
        else:
            cos_acc -= term
        power *= t
        fact *= 2 * n + 1
        term = power / fact
        if n % 2 == 0:
            sin_acc += term
        else:
            sin_acc -= term
        if n < n_terms:
            power *= t
            fact *= 2 * n + 2
    return sin_acc, cos_acc


def enclose_sin_cos(w, req_width, depth=0):
    """Enclosures of sin(w*pi) and cos(w*pi) for w in (0, 1/4], each of width <= req_width.

    The partial sum is taken at one rational point t*. The true angle differs
    from t* by at most delta, and |sin| and |cos| are 1-Lipschitz, so the
    values move by at most delta. The Lagrange remainder is at most
    x**(m+1)/(m+1)! because every derivative of sin or cos is bounded by 1.
    """
    if req_width <= 0:
        raise RuntimeError('trig width')
    if depth > 8:
        raise RuntimeError('trig enclosure did not tighten')
    # delta <= W/4 and each remainder <= W/4 => enclosure width 2*(R+delta) <= W.
    pi_lo, pi_hi = pi_bounds((req_width / 2) / w)
    if w * pi_hi >= 1:
        pi_lo, pi_hi = pi_bounds(Fraction(1, 1000))
        if w * pi_hi >= 1:
            raise RuntimeError('reduced angle is not below 1 radian')
    t_star = w * (pi_lo + pi_hi) / 2
    delta = w * (pi_hi - pi_lo) / 2
    xbound = w * pi_hi
    if t_star < 0 or xbound < t_star:
        raise RuntimeError('taylor centre')
    # A slightly larger 3-decimal bound keeps the remainder cheap and valid.
    x_simple = ceil_to_scale(xbound, 1000)
    if x_simple < xbound or x_simple >= 1:
        raise RuntimeError('remainder bound left [xbound, 1)')
    rem_budget = req_width / 4
    n_terms = 1
    r_sin = None
    r_cos = None
    while n_terms <= 5000:
        r_cos = pow_over_fact(x_simple, 2 * n_terms + 1)
        r_sin = pow_over_fact(x_simple, 2 * n_terms + 2)
        if r_sin <= rem_budget and r_cos <= rem_budget:
            break
        n_terms += 1
    else:
        raise RuntimeError('taylor did not converge')
    p_sin, p_cos = taylor_sin_cos(t_star, n_terms)
    sin_lo = p_sin - r_sin - delta
    sin_hi = p_sin + r_sin + delta
    cos_lo = p_cos - r_cos - delta
    cos_hi = p_cos + r_cos + delta
    if (
        sin_lo > sin_hi
        or cos_lo > cos_hi
        or sin_hi - sin_lo > req_width
        or cos_hi - cos_lo > req_width
    ):
        return enclose_sin_cos(w, req_width / 2, depth + 1)
    return sin_lo, sin_hi, cos_lo, cos_hi


def quadrant(f):
    """Reduce f in [0, 2) to w in [0, 1/4] plus signs and a sin/cos swap.

    sin(f*pi) = sin_sign * (cos(w*pi) if swap else sin(w*pi))
    cos(f*pi) = cos_sign * (sin(w*pi) if swap else cos(w*pi))
    """
    half = Fraction(1, 2)
    quarter = Fraction(1, 4)
    if f < 0 or f >= 2:
        raise RuntimeError('quadrant range')
    if f <= half:
        if f <= quarter:
            return f, 1, 1, False
        return half - f, 1, 1, True
    if f <= 1:
        u = 1 - f
        if u <= quarter:
            return u, 1, -1, False
        return half - u, 1, -1, True
    if f <= Fraction(3, 2):
        u = f - 1
        if u <= quarter:
            return u, -1, -1, False
        return half - u, -1, -1, True
    u = 2 - f
    if u <= quarter:
        return u, -1, 1, False
    return half - u, -1, 1, True


def _apply_sign(lo, hi, sign):
    if sign >= 0:
        return lo, hi
    return -hi, -lo


def cos_sin_qpi(q, width):
    """Enclosures of cos(q*pi) and sin(q*pi), each of width <= `width`."""
    if width < 0:
        raise RuntimeError('width')
    f = mod2(as_fraction(q))
    cached = _TRIG.get(f)
    if cached is not None and cached[0] <= width:
        return cached[1], cached[2], cached[3], cached[4]
    w, sin_sign, cos_sign, swap = quadrant(f)
    if w == 0:
        sin_lo = sin_hi = Fraction(0)
        cos_lo = cos_hi = Fraction(1)
        achieved = Fraction(0)
    else:
        sin_lo, sin_hi, cos_lo, cos_hi = enclose_sin_cos(w, width)
        achieved = width
    if swap:
        sin_lo, sin_hi, cos_lo, cos_hi = cos_lo, cos_hi, sin_lo, sin_hi
    sin_lo, sin_hi = _apply_sign(sin_lo, sin_hi, sin_sign)
    cos_lo, cos_hi = _apply_sign(cos_lo, cos_hi, cos_sign)
    if cos_lo > cos_hi or sin_lo > sin_hi:
        raise RuntimeError('trig bracket crossed')
    if cached is None or achieved < cached[0]:
        _TRIG[f] = (achieved, cos_lo, cos_hi, sin_lo, sin_hi)
    return cos_lo, cos_hi, sin_lo, sin_hi


def inradius_bounds(k, width):
    """a = cos(pi/k), the distance from the origin to each side."""
    cos_lo, cos_hi, _sin_lo, _sin_hi = cos_sin_qpi(Fraction(1, k), width)
    return cos_lo, cos_hi


def _mul_scalar(lo, hi, scale):
    if scale >= 0:
        return lo * scale, hi * scale
    return hi * scale, lo * scale


def slack_bounds(k, j, x, y, radius, slack_width):
    """Enclosure of a - n_j·c - r of width <= slack_width."""
    scale = Fraction(1) + abs(x) + abs(y)
    trig_width = slack_width / scale
    q = Fraction(3 * k + 4 * j, 2 * k)
    cos_lo, cos_hi, sin_lo, sin_hi = cos_sin_qpi(q, trig_width)
    a_lo, a_hi = inradius_bounds(k, trig_width)
    cx_lo, cx_hi = _mul_scalar(cos_lo, cos_hi, x)
    sy_lo, sy_hi = _mul_scalar(sin_lo, sin_hi, y)
    lo = a_lo - cx_hi - sy_hi - radius
    hi = a_hi - cx_lo - sy_lo - radius
    if lo > hi or hi - lo > slack_width:
        raise RuntimeError('slack enclosure failed')
    return lo, hi


def side_status(k, j, x, y, radius):
    """'in' if slack >= 0 is proved, 'out' if slack < 0 is proved, else 'undecided'."""
    width = Fraction(1, 10**12)
    while True:
        lo, hi = slack_bounds(k, j, x, y, radius, width)
        if lo >= 0:
            return 'in'
        if hi < 0:
            return 'out'
        if width <= SLACK_LIMIT:
            return 'undecided'
        nxt = width / 10000
        width = nxt if nxt > SLACK_LIMIT else SLACK_LIMIT


def judge(k, radius, centers, n_centres, record):
    """Return the VERDICT line. Sound: an invalid packing is never accepted."""
    if isinstance(k, bool) or not isinstance(k, int) or k < 3:
        return 'VERDICT: INVALID bad k'
    if isinstance(n_centres, bool) or not isinstance(n_centres, int):
        return 'VERDICT: INVALID bad count'
    try:
        radius = as_fraction(radius)
        record = as_fraction(record)
        centres = [(as_fraction(x), as_fraction(y)) for x, y in centers]
    except (TypeError, ValueError):
        return 'VERDICT: INVALID bad number'
    if radius <= 0:
        return 'VERDICT: INVALID radius not positive'
    if len(centres) != n_centres:
        return f'VERDICT: INVALID count {len(centres)} != {n_centres}'
    four_r2 = 4 * radius * radius
    for i in range(n_centres):
        xi, yi = centres[i]
        for j in range(i + 1, n_centres):
            dx = xi - centres[j][0]
            dy = yi - centres[j][1]
            if dx * dx + dy * dy < four_r2:
                return f'VERDICT: INVALID pair {i} {j} overlap'
    for i in range(n_centres):
        xi, yi = centres[i]
        for j in range(k):
            status = side_status(k, j, xi, yi, radius)
            if status == 'out':
                return f'VERDICT: INVALID circle {i} side {j} outside'
            if status == 'undecided':
                return f'VERDICT: INVALID circle {i} side {j} slack undecided'
    if radius > record:
        return 'VERDICT: IMPROVES'
    return 'VERDICT: VALID_NOT_BETTER'


def load_packing(path):
    with open(path, encoding='utf-8') as handle:
        raw_lines = handle.read().splitlines()
    lines = [line.strip() for line in raw_lines if line.strip()]
    if not lines:
        raise ValueError('empty file')
    head = lines[0].split()
    if len(head) != 2 or head[0] != 'r':
        raise ValueError('bad header')
    radius = parse_number(head[1])
    centres = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != 2:
            raise ValueError('bad centre line')
        centres.append((parse_number(parts[0]), parse_number(parts[1])))
    return radius, centres


def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv) != 5:
        print('usage: python verify_exact_kgon.py <k> <file> <N> <record_radius>')
        return 2
    try:
        k = int(argv[1])
        n_centres = int(argv[3])
        record = parse_number(argv[4])
        radius, centres = load_packing(argv[2])
    except (ValueError, InvalidOperation, OSError) as exc:
        message = ' '.join(str(exc).split())
        print(f'VERDICT: INVALID {message}')
        return 0
    print(judge(k, radius, centres, n_centres, record))
    return 0


def _expect(label, got, want):
    print(f'{label}: {got}')
    if got != want:
        raise AssertionError(f'{label}: got {got!r} want {want!r}')


def _expect_invalid(label, got, token):
    print(f'{label}: {got}')
    if not (got.startswith('VERDICT: INVALID') and token in got):
        raise AssertionError(f'{label}: got {got!r}, missing {token!r}')


def _check_known_angles():
    width = Fraction(1, 10**30)
    cos_lo, cos_hi, sin_lo, sin_hi = cos_sin_qpi(Fraction(1, 6), width)
    if not (sin_lo <= Fraction(1, 2) <= sin_hi):
        raise AssertionError(f'sin(pi/6) {[sin_lo, sin_hi]}')
    if sin_hi - sin_lo > width or cos_hi - cos_lo > width:
        raise AssertionError('sin(pi/6) width')
    if not (Fraction(4, 5) < cos_lo and cos_hi < Fraction(9, 10)):
        raise AssertionError(f'cos(pi/6) {[cos_lo, cos_hi]}')
    cos_lo, cos_hi, sin_lo, sin_hi = cos_sin_qpi(Fraction(1, 3), width)
    if not (cos_lo <= Fraction(1, 2) <= cos_hi and sin_lo > 0):
        raise AssertionError(f'cos(pi/3) {[cos_lo, cos_hi, sin_lo, sin_hi]}')
    cos_lo, cos_hi, sin_lo, sin_hi = cos_sin_qpi(Fraction(0), width)
    if not (cos_lo == 1 and cos_hi == 1 and sin_lo == 0 and sin_hi == 0):
        raise AssertionError('angle 0')
    cos_lo, cos_hi, sin_lo, sin_hi = cos_sin_qpi(Fraction(5, 2), width)
    if not (cos_lo == 0 and cos_hi == 0 and sin_lo == 1 and sin_hi == 1):
        raise AssertionError('angle 5/2 pi')
    print('known angles: OK')


def _expect_axis(k, j, cos_value, sin_value):
    q = Fraction(3 * k + 4 * j, 2 * k)
    cos_lo, cos_hi, sin_lo, sin_hi = cos_sin_qpi(q, Fraction(1, 10**6))
    if not (cos_lo == cos_value == cos_hi and sin_lo == sin_value == sin_hi):
        raise AssertionError(
            f'axis k={k} j={j}: cos {[cos_lo, cos_hi]} sin {[sin_lo, sin_hi]}'
        )


def _check_axes(k):
    _expect_axis(k, 0, Fraction(0), Fraction(-1))
    if k % 4 == 0:
        _expect_axis(k, k // 4, Fraction(1), Fraction(0))
        _expect_axis(k, k // 2, Fraction(0), Fraction(1))
        _expect_axis(k, (3 * k) // 4, Fraction(-1), Fraction(0))
    print(f'k={k} axes: OK')


def _check_centre(k):
    width = Fraction(1, 10**45)
    eps = Fraction(1, 10**30)
    a_lo, a_hi = inradius_bounds(k, Fraction(1, 10**8))
    floor = Fraction(1, 2) if k <= 6 else Fraction(9, 10)
    if not (floor < a_lo <= a_hi < 1 and a_hi - a_lo <= Fraction(1, 10**8)):
        raise AssertionError(f'k={k} inradius {[a_lo, a_hi]}')
    a_lo, a_hi = inradius_bounds(k, width)
    a_mid = (a_lo + a_hi) / 2
    origin = [(Fraction(0), Fraction(0))]
    r_under = a_mid - eps
    r_over = a_mid + eps
    _expect(f'k={k} N=1 under', judge(k, r_under, origin, 1, Fraction(0)), 'VERDICT: IMPROVES')
    _expect(
        f'k={k} N=1 record equal',
        judge(k, r_under, origin, 1, r_under),
        'VERDICT: VALID_NOT_BETTER',
    )
    _expect(
        f'k={k} N=1 record above',
        judge(k, r_under, origin, 1, r_under + 1),
        'VERDICT: VALID_NOT_BETTER',
    )
    _expect_invalid(f'k={k} N=1 over', judge(k, r_over, origin, 1, Fraction(0)), 'outside')


def _check_touch(k):
    """Centre built from a 1e-45 enclosure of side j=3, then radius moved by 1e-30."""
    width = Fraction(1, 10**45)
    eps = Fraction(1, 10**30)
    r_nom = Fraction(1, 100)
    q = Fraction(3 * k + 4 * 3, 2 * k)
    cos_lo, cos_hi, sin_lo, sin_hi = cos_sin_qpi(q, width)
    a_lo, a_hi = inradius_bounds(k, width)
    gap = (a_lo + a_hi) / 2 - r_nom
    centre = (gap * (cos_lo + cos_hi) / 2, gap * (sin_lo + sin_hi) / 2)
    r_under = r_nom - eps
    r_over = r_nom + eps
    _expect(f'k={k} touch under', judge(k, r_under, [centre], 1, Fraction(0)), 'VERDICT: IMPROVES')
    _expect(
        f'k={k} touch record equal',
        judge(k, r_under, [centre], 1, r_under),
        'VERDICT: VALID_NOT_BETTER',
    )
    _expect_invalid(f'k={k} touch over', judge(k, r_over, [centre], 1, Fraction(0)), 'outside')


def _check_pairs(k):
    radius = Fraction(1, 10)
    eps = Fraction(1, 10**25)
    zero = Fraction(0)
    touching = [(-radius, zero), (radius, zero)]
    _expect(f'k={k} pair touch', judge(k, radius, touching, 2, Fraction(0)), 'VERDICT: IMPROVES')
    _expect(
        f'k={k} pair record equal',
        judge(k, radius, touching, 2, radius),
        'VERDICT: VALID_NOT_BETTER',
    )
    overlap = [(-radius, zero), (radius - eps, zero)]
    _expect_invalid(f'k={k} overlap', judge(k, radius, overlap, 2, Fraction(0)), 'overlap')
    _expect_invalid(f'k={k} count', judge(k, radius, touching, 3, Fraction(0)), 'count')
    _expect_invalid(
        f'k={k} zero radius',
        judge(k, Fraction(0), [(zero, zero)], 1, Fraction(0)),
        'radius',
    )


def _check_file():
    text = 'r 1e-3\n-1e-3 1.5e-17\n1e-3 0\n'
    fd, path = tempfile.mkstemp(suffix='.txt')
    os.close(fd)
    try:
        with open(path, 'w', encoding='utf-8', newline='\n') as handle:
            handle.write(text)
        cases = (
            (['verify_exact_kgon.py', '16', path, '2', '0'], 'VERDICT: IMPROVES', None),
            (['verify_exact_kgon.py', '15', path, '2', '1e-3'], 'VERDICT: VALID_NOT_BETTER', None),
            (['verify_exact_kgon.py', '16', path, '1', '0'], None, 'count'),
        )
        for argv, exact, token in cases:
            buf = io.StringIO()
            with redirect_stdout(buf):
                code = main(argv)
            got = buf.getvalue().strip()
            print(f'file {argv[3]} {argv[4]}: {got}')
            if code != 0:
                raise AssertionError(got)
            if exact is not None and got != exact:
                raise AssertionError(got)
            if token is not None and token not in got:
                raise AssertionError(got)
    finally:
        os.remove(path)


def run_selftests():
    print('SELFTEST: running')
    _check_known_angles()
    for k in (5, 15, 16):
        _check_axes(k)
        _check_centre(k)
        _check_touch(k)
        _check_pairs(k)
    _check_file()


if __name__ == '__main__':
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] == '--selftest'):
        try:
            run_selftests()
            print('SELFTEST: OK')
        except Exception as exc:
            print(f'SELFTEST: FAIL {exc}')
            sys.exit(1)
    else:
        sys.exit(main())
