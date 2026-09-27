#!/usr/bin/env python3
# check_local_blind.py
# Blind checker for LOCAL.md (N = 16 local stage). Python 3 standard library only.
# Decisions use Fraction / int only. Floats are never read into a decision.
#
# Stricter readings where LOCAL.md is ambiguous:
# 1. LOCAL.md names fields, not keys. The schema below is required; any other
#    layout FAILs. Optional derived keys, when present as rationals, must equal
#    the recomputation: r, Lambda, A, B, E, R, delta, Phi, m_lo2, m_hi2, epsilons.
#    Keys m_lo and m_hi are ignored (the roots need not be rational).
# 2. Every ||d|| in Lemma L is the skew sup-norm on the 30 coordinates (u, v).
#    The third barycentric coordinate is dependent and is not an extra coordinate.
# 3. Lemma V is instantiated at delta(m_lo). m_t <= m_lo and m_lo <= m_hi are
#    decided by comparing squares.
# 4. Each listed wall is exactly 0 at c~. Kind "s" means 1 - u - v.
# 5. Multiplier order is the 20 pairs (stored endpoint order; either order;
#    duplicate unordered pairs rejected; both ends in 0..14) then the 13 walls.
#    Indices are 0-based. Pairs do not involve the rattler. Counts are exact.
# 6. Cluster frame order is certificate frame order. No symmetry and no
#    permutation search is applied. The rattler is excluded from the ball test.
#    Exactly 3 clusters, matched by integer id.
# 7. rho_q is one positive integer grid radius for that cluster (top-level rho
#    is the default). G is one positive integer. The closed test on each frame
#    point is max(|A/G - u|, |B/G - v|) + rho/G <= R.
# 8. m_t is required and positive, from the optional third argument and/or an
#    m_t field; every supplied copy must agree. Every cluster needs m_t^2 <= m_lo^2.
# 9. A <= 0 or E >= 1 FAILs (the stated formulas divide by 2A and by 1 - E).
#    The pair remainder uses the spec constant Q <= 12 ||d||_inf^2.
# 10. The mirror claim is not re-derived. Each cluster's own rationals must carry it.
# 11. delta > R does not fail: ||d||_inf <= R still implies ||d||_inf <= delta,
#     so Lemma V stays a valid upper bound (no tighter than the ball).

import json
import re
import sys
from fractions import Fraction

N_FRAME = 15
N_COORD = 30
N_PAIRS = 20
N_WALLS = 13
N_ROWS = N_PAIRS + N_WALLS
N_CONFIG = 16
N_CONFIG_PAIRS = 120
N_Y = N_COORD * 2
Q_COEFF = 12  # sharp: |du| = |dv| = 2 ||d||_inf gives Q = 12 ||d||_inf^2

_INT = re.compile(r"[+-]?\d+")
_RAT = re.compile(r"^([+-]?\d+)(?:/([+-]?\d+)(?:\^(\d+))?)?$")


class CheckerError(Exception):
    def __init__(self, message):
        super().__init__(" ".join(str(message).split()))


def need(obj, key, label):
    if not isinstance(obj, dict):
        raise CheckerError(f"{label} is not an object")
    if key not in obj:
        raise CheckerError(f"{label} missing '{key}'")
    return obj[key]


def parse_integer(x, label):
    if isinstance(x, bool) or isinstance(x, float):
        raise CheckerError(f"{label}: refuse non-exact float/bool")
    if isinstance(x, int):
        return x
    if isinstance(x, str) and _INT.fullmatch(x.strip()):
        return int(x.strip())
    raise CheckerError(f"{label}: expected integer, got {x!r}")


def parse_rational(x, label):
    if isinstance(x, bool) or isinstance(x, float):
        raise CheckerError(f"{label}: refuse non-exact float/bool")
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, str):
        m = _RAT.fullmatch(x.strip().replace("**", "^"))
        if not m:
            raise CheckerError(f"{label}: bad rational {x!r}")
        num = int(m.group(1))
        if m.group(2) is None:
            return Fraction(num)
        base = int(m.group(2))
        if m.group(3) is None:
            den = base
        else:
            exp = int(m.group(3))
            if exp > 100000:
                raise CheckerError(f"{label}: exponent too large")
            den = base ** exp
        try:
            return Fraction(num, den)
        except ZeroDivisionError as exc:
            raise CheckerError(f"{label}: zero denominator") from exc
    if isinstance(x, list):
        if len(x) != 2:
            raise CheckerError(f"{label}: rational pair must have length 2")
        return parse_rational_pair(x[0], x[1], label)
    raise CheckerError(f"{label}: unsupported type {type(x).__name__}")


def parse_rational_pair(num, den, label):
    try:
        return Fraction(parse_integer(num, label), parse_integer(den, label))
    except ZeroDivisionError as exc:
        raise CheckerError(f"{label}: zero denominator") from exc


def parse_vector(x, n, label):
    if not isinstance(x, list) or len(x) != n:
        raise CheckerError(f"{label}: expected a list of length {n}")
    return [parse_rational(v, f"{label}[{i}]") for i, v in enumerate(x)]


def require_nonnegative(vec, label):
    for i, v in enumerate(vec):
        if v < 0:
            raise CheckerError(f"{label}[{i}] is negative")


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Fraction(0))


def l1(vec):
    return sum((abs(v) for v in vec), Fraction(0))


def Q_diff(p, q):
    du = p[0] - q[0]
    dv = p[1] - q[1]
    return du * du + du * dv + dv * dv


def in_triangle(p, label):
    u, v = p
    if u < 0 or v < 0 or u + v > 1:
        raise CheckerError(f"{label} is outside the triangle")


def parse_point(x, label):
    if not isinstance(x, list) or len(x) != 2:
        raise CheckerError(f"{label}: expected [u, v]")
    return (parse_rational(x[0], label + ".u"), parse_rational(x[1], label + ".v"))


def parse_int_point(x, label):
    if not isinstance(x, list) or len(x) != 2:
        raise CheckerError(f"{label}: expected [A, B]")
    return (parse_integer(x[0], label + ".A"), parse_integer(x[1], label + ".B"))


def pair_row(i, j, frame):
    du = frame[i][0] - frame[j][0]
    dv = frame[i][1] - frame[j][1]
    gu = 2 * du + dv
    gv = du + 2 * dv
    row = [Fraction(0)] * N_COORD
    row[2 * i] = gu
    row[2 * i + 1] = gv
    row[2 * j] = -gu
    row[2 * j + 1] = -gv
    return row


def wall_value(frame, point, kind):
    u, v = frame[point]
    if kind == "u":
        return u
    if kind == "v":
        return v
    return 1 - u - v


def wall_row(point, kind):
    row = [Fraction(0)] * N_COORD
    if kind == "u":
        row[2 * point] = Fraction(1)
    elif kind == "v":
        row[2 * point + 1] = Fraction(1)
    else:
        row[2 * point] = Fraction(-1)
        row[2 * point + 1] = Fraction(-1)
    return row


def jt_dot(J, y):
    out = [Fraction(0)] * N_COORD
    for row, yk in zip(J, y):
        if yk == 0:
            continue
        for i, jki in enumerate(row):
            if jki != 0:
                out[i] += jki * yk
    return out


def load_json(path):
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except OSError as exc:
        raise CheckerError(f"cannot read {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise CheckerError(f"invalid JSON in {path}: {exc}") from exc


def match_stored_rational(obj, key, expected, label):
    if key not in obj or isinstance(obj[key], float):
        return
    got = parse_rational(obj[key], f"{label}.{key}")
    if got != expected:
        raise CheckerError(f"{label}.{key} stored {got} != recomputed {expected}")


def match_stored_vector(obj, key, expected, label):
    if key not in obj:
        return
    if isinstance(obj[key], float):
        raise CheckerError(f"{label}.{key}: refuse float")
    got = parse_vector(obj[key], len(expected), f"{label}.{key}")
    if got != list(expected):
        raise CheckerError(f"{label}.{key} does not match recomputation")


def parse_pairs(raw, label):
    if not isinstance(raw, list) or len(raw) != N_PAIRS:
        raise CheckerError(f"{label}: expected {N_PAIRS} pairs")
    pairs = []
    seen = set()
    for n, item in enumerate(raw):
        if not isinstance(item, list) or len(item) != 2:
            raise CheckerError(f"{label} pair {n} is not a pair")
        i = parse_integer(item[0], f"{label} pair {n}")
        j = parse_integer(item[1], f"{label} pair {n}")
        if i == j or not (0 <= i < N_FRAME and 0 <= j < N_FRAME):
            raise CheckerError(f"{label} pair {n} is not two distinct frame points")
        key = (i, j) if i < j else (j, i)
        if key in seen:
            raise CheckerError(f"{label} repeats unordered pair {key}")
        seen.add(key)
        pairs.append((i, j))
    return pairs


def parse_walls(raw, label):
    if not isinstance(raw, list) or len(raw) != N_WALLS:
        raise CheckerError(f"{label}: expected {N_WALLS} walls")
    walls = []
    seen = set()
    for n, item in enumerate(raw):
        where = f"{label} wall {n}"
        point = parse_integer(need(item, "point", where), where)
        kind = need(item, "kind", where)
        if kind not in ("u", "v", "s") or not (0 <= point < N_FRAME):
            raise CheckerError(f"{where} has a bad point or kind")
        if (point, kind) in seen:
            raise CheckerError(f"{where} repeats {(point, kind)}")
        seen.add((point, kind))
        walls.append((point, kind))
    return walls


def parse_y(raw, label):
    if not isinstance(raw, list) or len(raw) != N_Y:
        raise CheckerError(f"{label}: expected {N_Y} multipliers y")
    ymap = {}
    for n, item in enumerate(raw):
        where = f"{label} y[{n}]"
        coord = parse_integer(need(item, "coord", where), where)
        sign = parse_integer(need(item, "sign", where), where)
        if not (0 <= coord < N_COORD) or sign not in (1, -1):
            raise CheckerError(f"{where} has coord/sign outside 0..29 and ±1")
        if (coord, sign) in ymap:
            raise CheckerError(f"{label} repeats y[{coord},{sign}]")
        values = parse_vector(need(item, "values", where), N_ROWS, where)
        require_nonnegative(values, where)
        ymap[(coord, sign)] = values
    for coord in range(N_COORD):
        for sign in (1, -1):
            if (coord, sign) not in ymap:
                raise CheckerError(f"{label} missing y[{coord},{sign}]")
    return ymap


def check_epsilons(obj, eps_map, label):
    if "epsilons" not in obj:
        return
    raw = obj["epsilons"]
    if not isinstance(raw, list) or len(raw) != N_Y:
        raise CheckerError(f"{label}.epsilons: expected {N_Y} entries")
    seen = set()
    for n, item in enumerate(raw):
        where = f"{label}.epsilons[{n}]"
        coord = parse_integer(need(item, "coord", where), where)
        sign = parse_integer(need(item, "sign", where), where)
        values = parse_vector(need(item, "values", where), N_COORD, where)
        key = (coord, sign)
        if key not in eps_map or values != eps_map[key]:
            raise CheckerError(f"{where} does not match recomputation")
        seen.add(key)
    if seen != set(eps_map):
        raise CheckerError(f"{label}.epsilons does not cover every (coord, sign)")


def min_pair_Q(points, label):
    best = None
    count = 0
    for a in range(N_CONFIG):
        for b in range(a + 1, N_CONFIG):
            q = Q_diff(points[a], points[b])
            if q <= 0:
                raise CheckerError(f"{label}: coincident points {a}, {b}")
            if best is None or q < best:
                best = q
            count += 1
    if count != N_CONFIG_PAIRS or best is None:
        raise CheckerError(f"{label}: internal pair count {count}")
    return best


def check_certificate_cluster(obj, cid):
    label = f"cluster {cid}"
    frame_raw = need(obj, "frame", label)
    if not isinstance(frame_raw, list) or len(frame_raw) != N_FRAME:
        raise CheckerError(f"{label}: frame must be {N_FRAME} points")
    frame = [parse_point(p, f"{label} frame[{i}]") for i, p in enumerate(frame_raw)]
    for i, p in enumerate(frame):
        in_triangle(p, f"{label} frame[{i}]")
    rattler = parse_point(need(obj, "rattler", label), f"{label} rattler")
    in_triangle(rattler, f"{label} rattler")

    pairs = parse_pairs(need(obj, "pairs", label), label)
    walls = parse_walls(need(obj, "walls", label), label)
    for point, kind in walls:
        if wall_value(frame, point, kind) != 0:
            raise CheckerError(f"{label}: wall {kind} of point {point} is not exactly 0")

    lam = parse_vector(need(obj, "lambda", label), N_ROWS, f"{label} lambda")
    require_nonnegative(lam, f"{label} lambda")
    Lambda = sum(lam[:N_PAIRS], Fraction(0))
    if Lambda <= 0:
        raise CheckerError(f"{label}: pair-stress sum Lambda is not positive")

    m_lo2 = min_pair_Q(list(frame) + [rattler], label)
    gvals = [Q_diff(frame[i], frame[j]) - m_lo2 for i, j in pairs]
    gvals.extend(wall_value(frame, point, kind) for point, kind in walls)
    phi = dot(lam, [Q_diff(frame[i], frame[j]) for i, j in pairs] + gvals[N_PAIRS:])
    if phi - Lambda * m_lo2 != dot(lam, gvals):
        raise CheckerError(f"{label}: internal Phi identity failed")

    J = [pair_row(i, j, frame) for i, j in pairs]
    J.extend(wall_row(point, kind) for point, kind in walls)
    ymap = parse_y(need(obj, "y", label), label)

    # ε = J^T y + σ e_i, the sign that makes σ d_i = ε·d - y·(J d).
    pair_mass_max = Fraction(0)
    B = None
    E = Fraction(0)
    eps_map = {}
    for coord in range(N_COORD):
        for sign in (1, -1):
            y = ymap[(coord, sign)]
            mass = sum(y[:N_PAIRS], Fraction(0))
            if mass > pair_mass_max:
                pair_mass_max = mass
            jty = jt_dot(J, y)
            extra = Fraction(sign)
            eps = jty[:]
            eps[coord] += extra
            for t, (ev, jv) in enumerate(zip(eps, jty)):
                add = extra if t == coord else Fraction(0)
                if ev - jv != add:
                    raise CheckerError(f"{label}: internal multiplier identity failed")
            e_l1 = l1(eps)
            if e_l1 > E:
                E = e_l1
            slack = dot(y, gvals)
            if B is None or slack > B:
                B = slack
            eps_map[(coord, sign)] = eps

    A = Q_COEFF * pair_mass_max
    if A <= 0:
        raise CheckerError(f"{label}: A <= 0, so R = (1 - E) / (2A) is undefined")
    if E >= 1:
        raise CheckerError(f"{label}: E >= 1, so the localisation ball is empty")
    if B < 0:
        raise CheckerError(f"{label}: B(m_lo) < 0")
    one_minus_E = 1 - E
    R = one_minus_E / (2 * A)
    delta = (2 * B) / one_minus_E
    if R <= 0:
        raise CheckerError(f"{label}: recomputed R is not positive")

    r = jt_dot(J, lam)
    m_hi2 = (phi + l1(r) * delta + Q_COEFF * Lambda * delta * delta) / Lambda
    if m_hi2 < m_lo2:
        raise CheckerError(f"{label}: recomputed m_hi^2 is below m_lo^2")

    match_stored_vector(obj, "r", r, label)
    match_stored_rational(obj, "Lambda", Lambda, label)
    match_stored_rational(obj, "A", A, label)
    match_stored_rational(obj, "B", B, label)
    match_stored_rational(obj, "E", E, label)
    match_stored_rational(obj, "R", R, label)
    match_stored_rational(obj, "delta", delta, label)
    match_stored_rational(obj, "Phi", phi, label)
    match_stored_rational(obj, "m_lo2", m_lo2, label)
    match_stored_rational(obj, "m_hi2", m_hi2, label)
    check_epsilons(obj, eps_map, label)
    return {"frame": frame, "R": R, "m_lo2": m_lo2, "m_hi2": m_hi2}


def index_clusters(raw, label):
    if not isinstance(raw, list) or len(raw) != 3:
        raise CheckerError(f"{label} must contain exactly 3 clusters")
    out = {}
    for item in raw:
        cid = parse_integer(need(item, "id", label), f"{label} id")
        if cid in out:
            raise CheckerError(f"{label} repeats cluster id {cid}")
        out[cid] = item
    return out


def collect_m_t(cli, cert, boxes):
    found = []
    if cli is not None:
        found.append(parse_rational(cli, "m_t argv"))
    for src, obj in (("cert", cert), ("clusters", boxes)):
        if isinstance(obj, dict) and "m_t" in obj:
            found.append(parse_rational(obj["m_t"], f"{src}.m_t"))
    if not found:
        raise CheckerError("m_t missing: pass it as the third argument or as m_t")
    if any(v != found[0] for v in found[1:]):
        raise CheckerError("m_t values disagree")
    if found[0] <= 0:
        raise CheckerError("m_t is not positive")
    return found[0]


def cluster_rho(obj, default_rho, cid):
    label = f"cluster {cid}"
    if "rho" in obj:
        rho = parse_integer(obj["rho"], f"{label}.rho")
    elif default_rho is not None:
        rho = default_rho
    else:
        raise CheckerError(f"{label}: rho missing")
    if rho <= 0:
        raise CheckerError(f"{label}: rho is not positive")
    return rho


def check_box(obj, certified, G, rho, cid):
    label = f"cluster {cid}"
    if "G" in obj and parse_integer(obj["G"], f"{label}.G") != G:
        raise CheckerError(f"{label}: G differs from the file G")
    raw = need(obj, "frame", label)
    if not isinstance(raw, list) or len(raw) != N_FRAME:
        raise CheckerError(f"{label}: box frame must be {N_FRAME} integer points")
    rad = Fraction(rho, G)
    R = certified["R"]
    for i, item in enumerate(raw):
        A, B = parse_int_point(item, f"{label} box[{i}]")
        u, v = certified["frame"][i]
        dist = max(abs(Fraction(A, G) - u), abs(Fraction(B, G) - v))
        reach = dist + rad
        if reach > R:
            raise CheckerError(
                f"{label} point {i}: box reaches {reach} which is past R = {R}"
            )


def run(argv):
    if len(argv) not in (3, 4):
        raise CheckerError(
            "usage: check_local_blind.py <cert.json> <clusters.json> [m_t]"
        )
    cert = load_json(argv[1])
    boxes = load_json(argv[2])
    if not isinstance(cert, dict) or not isinstance(boxes, dict):
        raise CheckerError("both JSON roots must be objects")
    m_t = collect_m_t(argv[3] if len(argv) == 4 else None, cert, boxes)
    G = parse_integer(need(boxes, "G", "clusters"), "clusters.G")
    if G <= 0:
        raise CheckerError("G is not positive")
    default_rho = None
    if "rho" in boxes:
        default_rho = parse_integer(boxes["rho"], "clusters.rho")
        if default_rho <= 0:
            raise CheckerError("clusters.rho is not positive")

    cert_map = index_clusters(need(cert, "clusters", "cert"), "cert")
    box_map = index_clusters(need(boxes, "clusters", "clusters"), "clusters")
    if set(cert_map) != set(box_map):
        raise CheckerError("certificate cluster ids and box cluster ids differ")

    m_t2 = m_t * m_t
    for cid in sorted(cert_map):
        certified = check_certificate_cluster(cert_map[cid], cid)
        if m_t2 > certified["m_lo2"]:
            raise CheckerError(
                f"cluster {cid}: m_t^2 = {m_t2} exceeds m_lo^2 = {certified['m_lo2']}"
            )
        rho = cluster_rho(box_map[cid], default_rho, cid)
        check_box(box_map[cid], certified, G, rho, cid)


def main(argv):
    try:
        run(argv)
    except CheckerError as exc:
        print(f"VERDICT: FAIL {exc}")
        return 1
    print("VERDICT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
