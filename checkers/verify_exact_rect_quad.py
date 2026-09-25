# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
# verify_exact2.py — equal circles in rect:<h> or quad. No float in any decision.
import sys
from decimal import Decimal
from fractions import Fraction

def die():
    print("VERDICT: FAIL")
    raise SystemExit(1)

def Q(s):
    # Finite decimal -> exact rational. Decimal stores every base-10 digit.
    try:
        return Fraction(Decimal(s))
    except Exception:
        die()

def main():
    if len(sys.argv) != 5:
        die()
    spec, path, n_s, rec_s = sys.argv[1:5]
    try:
        N = int(n_s)
    except ValueError:
        die()
    rec = Q(rec_s)
    try:
        rows = [ln.split() for ln in open(path) if ln.strip()]
    except OSError:
        die()
    if not rows or len(rows[0]) != 2 or rows[0][0] != "r":
        die()
    r, pts = Q(rows[0][1]), []
    for row in rows[1:]:
        if len(row) != 2:
            die()
        pts.append((Q(row[0]), Q(row[1])))
    if spec == "quad":
        gap = 1 - r  # gap >= 0 before the square: disk of radius 1 - r
        def inside(x, y):
            return x >= r and y >= r and gap >= 0 and x * x + y * y <= gap * gap
    elif spec.startswith("rect:"):
        h = Q(spec[5:])
        x0, x1 = -Fraction(1, 2) + r, Fraction(1, 2) - r
        y0, y1 = -h / 2 + r, h / 2 - r
        def inside(x, y):
            return x0 <= x <= x1 and y0 <= y <= y1
    else:
        die()
    fit = [] if r < 0 else [i for i, p in enumerate(pts) if not inside(*p)]
    need = 4 * r * r  # >= (2r)^2; equal distance is contact
    hits = sum(
        (pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2 < need
        for i in range(len(pts)) for j in range(i + 1, len(pts))
    )
    margin = r - rec
    rel = "IMPROVES" if margin > 0 else "TIES" if margin == 0 else "WORSE"
    print("count OK" if len(pts) == N else f"count FAIL got {len(pts)} expected {N}")
    print("fit FAIL r<0" if r < 0 else ("fit OK" if not fit else "fit FAIL " + " ".join(map(str, fit[:12]))))
    print("pairs OK" if hits == 0 else f"pairs FAIL {hits}")
    print(f"{rel} margin={margin}")
    passed = len(pts) == N and r >= 0 and not fit and hits == 0
    print(f"VERDICT: {rel}" if passed else "VERDICT: FAIL")
    return 0 if passed else 1

if __name__ == "__main__":
    sys.exit(main())
