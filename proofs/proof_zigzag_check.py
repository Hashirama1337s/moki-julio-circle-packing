"""Standalone check of PROOFS.md (standard library only): for every row of PROOF_ZIGZAG.csv (table, h, N, Packomania's published
radius) compute the zig-zag optimum r_zz exactly enough (80-digit decimals), confirm it lies in the theorem's window
h(2 - sqrt3) <= r_zz <= h/2, and confirm Packomania's radius equals it to 28 digits.   usage: python proof_zigzag_check.py
"""
import csv, os
from decimal import Decimal as D, getcontext
getcontext().prec = 80
HERE = os.path.dirname(os.path.abspath(__file__))
ok = bad = 0
for row in csv.DictReader(open(os.path.join(HERE, "PROOF_ZIGZAG.csv"))):
    h, n, pub = D(row["h"]), int(row["N"]), D(row["packomania_radius"])
    k = D((n - 1) ** 2)
    rz = ((1 + k * h) - ((1 + k * h) ** 2 - 1 - k * h * h).sqrt()) / 2
    window = h * (2 - D(3).sqrt()) <= rz <= h / 2
    eq = abs((n - 1) * (4 * h * rz - h * h).sqrt() + 2 * rz - 1) < D(10) ** -60
    match = abs(rz - pub) < D(10) ** -28
    if window and eq and match: ok += 1
    else: bad += 1; print("FAIL", row["table"], n, window, eq, match)
print(f"{ok} / {ok + bad} entries: published radius == proven zig-zag optimum (28 digits), theorem window holds")
