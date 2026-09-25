"""Flagship csq, stage 2: hex row-lattice families (vision/recon2_families.py: alt / shift / sq in the unit square) for a count
N' >= N, minus N' - N least-contacted circles -> certificate -> certify_big.py. Replaces cand_big/csq/csq_<N>.txt only if the
certified radius is larger than the file there. Targets: out/hats_csq_families.json (N, gain_vs_best_so_far, gain_vs_table, src).
usage: py -3.11 hats_csq_fam.py   -> out/hats_csq_fam.jsonl
"""
import os, sys, json, io, contextlib
import numpy as np
from decimal import Decimal, ROUND_FLOOR, getcontext
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "vision"))
with contextlib.redirect_stdout(io.StringIO()): import recon2_families as F
import certify_big, hats_flagship as HF
getcontext().prec = 60
done = {json.loads(l)["N"] for l in open(os.path.join(HERE, "out", "hats_csq_fam.jsonl"))} if os.path.exists(os.path.join(HERE, "out", "hats_csq_fam.jsonl")) else set()
T = HF.table("csq"); out = open(os.path.join(HERE, "out", "hats_csq_fam.jsonl"), "a"); ok = 0
for n, g_best, g_tab, src in json.load(open(os.path.join(HERE, "out", "hats_csq_families.json"))):
    if n in done: continue
    m, desc = src[0], src[1]; kind = desc.split()[0]; k = int(desc.split("k=")[1].split()[0]); nn = int(desc.split("n=")[1])
    N2, r = F.family(1.0, 1.0, k, nn, kind); assert N2 == m
    c = np.array(F.build(1.0, 1.0, k, nn, kind, r)) - 0.5
    r_m = HF.clearance("csq", c)
    from scipy.spatial import cKDTree
    near = cKDTree(c).query_pairs(2 * r_m * (1 + 1e-9), output_type="ndarray"); deg = np.bincount(near.ravel(), minlength=len(c))
    deg = deg + (0.5 - np.abs(c) <= r_m * (1 + 1e-9)).sum(1)
    keep = np.sort(np.argsort(deg, kind="stable")[m - n:]); ck = c[keep]; r_f = HF.clearance("csq", ck)
    r_claim = (Decimal(repr(float(r_f))) * (1 - Decimal(10) ** -12)).quantize(Decimal(10) ** -25, rounding=ROUND_FLOOR)
    f = os.path.join(HERE, "cand_big", "csq", f"csq_{n}.txt"); tmp = f + ".fam"
    open(tmp, "w", newline="\n").write(f"r {r_claim}\n" + "".join(f"{format(Decimal(repr(float(x))), 'f')} {format(Decimal(repr(float(y))), 'f')}\n" for x, y in ck))
    v = certify_big.check("square", tmp, n, T[n])
    old = Decimal(open(f).readline().split()[1]) if os.path.exists(f) else Decimal(0)
    keepit = v == "IMPROVES" and r_claim > old
    if keepit: os.replace(tmp, f); ok += 1
    else: os.remove(tmp)
    row = {"table": "csq", "N": n, "from": f"{desc} ({m}) minus {m - n}", "r_claim": str(r_claim), "r_pub": T[n], "checker_a": v,
           "gain_vs_packomania": float(r_claim / Decimal(T[n]) - 1), "replaced": bool(keepit), "was": str(old)}
    out.write(json.dumps(row) + "\n")
print("family certificates kept:", ok)
