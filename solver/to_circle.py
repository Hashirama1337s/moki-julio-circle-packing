# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Point-form candidate (rational points in the unit triangle) -> circle-form file for the independent circle-form checker.
Circles: r = d/(2+(2+sqrt2)d), centre = r + L*p with L = 1-(2+sqrt2)r (computed at 80 digits), centres rounded to 50 digits,
then r_c = the largest radius the ROUNDED centres admit (walls, hypotenuse, half pair distances) rounded DOWN at 45 digits.
usage: py -3.11 to_circle.py cert/cand_N.txt  ->  cert/circ_N.txt"""
import sys, mpmath as mp
from fractions import Fraction as F
mp.mp.dps = 80
def conv(path, out=None):
    L = [l.split() for l in open(path) if l.strip()]; n = int(L[0][1])
    P = [(F(a), F(b)) for a, b in L[1:]]
    d2 = min((p[0]-q[0])**2 + (p[1]-q[1])**2 for i, p in enumerate(P) for q in P[i+1:])
    d = mp.sqrt(mp.mpf(d2.numerator) / d2.denominator); s2 = mp.sqrt(2)
    r = d / (2 + (2 + s2) * d); Lg = 1 - (2 + s2) * r
    C = [(mp.nstr(r + Lg * mp.mpf(x.numerator) / x.denominator, 50, strip_zeros=False),
          mp.nstr(r + Lg * mp.mpf(y.numerator) / y.denominator, 50, strip_zeros=False)) for x, y in P]
    Cm = [(mp.mpf(a), mp.mpf(b)) for a, b in C]
    rmax = min(min(x, y, (1 - x - y) / s2) for x, y in Cm)
    rmax = min(rmax, min(mp.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2) for i, p in enumerate(Cm) for q in Cm[i+1:]) / 2)
    rc = mp.floor(rmax * mp.mpf(10) ** 45) / mp.mpf(10) ** 45
    import os
    out = out or os.path.join(os.path.dirname(path), "circ_" + os.path.basename(path))
    assert os.path.abspath(out) != os.path.abspath(path), "refusing to overwrite the input certificate"
    with open(out, "w") as f:
        f.write(f"r {mp.nstr(rc, 45, strip_zeros=False)}\n")
        for a, b in C: f.write(f"{a} {b}\n")
    return out, n
if __name__ == "__main__":
    for p in sys.argv[1:]: print(conv(p))
