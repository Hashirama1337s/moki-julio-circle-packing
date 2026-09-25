# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
# verify_exact.py — equal circles, unit isosceles right triangle, Packomania crt (circle form).
# CLI: verify_exact.py file N record_radius_decimal
# File: first line "r <decimal>", then N lines "x y". Decisions use Fraction only.
import sys
from decimal import Decimal
from fractions import Fraction

def frac(s):
    # Decimal(s) keeps every digit; Fraction of that is exact. float() is never called.
    return Fraction(Decimal(s))

def main():
    path, n_s, rec_s = sys.argv[1:]
    N, record = int(n_s), frac(rec_s)
    rows = [ln.split() for ln in open(path, encoding="utf-8") if ln.strip()]
    good_r = bool(rows) and rows[0][0] == "r" and len(rows[0]) == 2
    r = frac(rows[0][1]) if good_r else None
    pts = [(frac(a), frac(b)) for a, b in rows[1:]] if all(len(t) == 2 for t in rows[1:]) else []
    ok_n = r is not None and r > 0 and len(rows) == N + 1 and len(pts) == N
    print("a) count:", "PASS" if ok_n else "FAIL")
    ok_b = ok_n and all(x >= r and y >= r for x, y in pts)
    print("b) legs:", "PASS" if ok_b else "FAIL")
    # Hypotenuse: need L = 1-x-y >= r*sqrt(2). sqrt(2) is irrational.
    # L < 0 is already outside. For L >= 0 and r > 0, squaring is monotone on [0, inf):
    # L*L >= 2*r*r  is exactly L >= r*sqrt(2). Equality (touching) passes.
    ok_c = ok_n
    if ok_n:
        for x, y in pts:
            L = 1 - x - y
            if L < 0 or L * L < 2 * r * r:
                ok_c = False
                break
    print("c) hypotenuse:", "PASS" if ok_c else "FAIL")
    ok_d = ok_n
    if ok_n:
        for i in range(N):
            for xj, yj in pts[i + 1:]:
                dx, dy = pts[i][0] - xj, pts[i][1] - yj
                if dx * dx + dy * dy < 4 * r * r:
                    ok_d = False
                    break
            if not ok_d:
                break
    print("d) separation:", "PASS" if ok_d else "FAIL")
    if r is None:
        rel, margin = "FAIL", None
    else:
        margin = r - record
        rel = "IMPROVES" if r > record else "TIES" if r == record else "WORSE"
    print(f"e) record: {rel} margin={margin}")
    print("VERDICT:", rel if (ok_n and ok_b and ok_c and ok_d) else "FAIL")

if __name__ == "__main__":
    main()
