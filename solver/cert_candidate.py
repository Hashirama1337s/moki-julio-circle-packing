"""Certify a candidate packing (float centres .npy) through the full chain: hpslp.converge -> refine_circ Newton -> 45-digit
certificate (admissible r rounded DOWN) -> checker A + checker B (both exact) vs the PACKOMANIA radius -> lopt.
Keeps it only if it beats our current best (v1.1 converged file, else v1.0 submission).  -> cand_hp/<shelf>/<shelf>_<N>.txt + .json
usage: py -3.11 cert_candidate.py <shelf> <N> <candidate.npy> [tag]
"""
import os, sys, json, numpy as np, mpmath as mp
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hpslp, refine_circ, certify_circ, lopt, finalize_circ
from lopt_fix import checker_b_verdict

def ours(shelf, n):
    for p in (os.path.join(HERE, "cand_hp", shelf, f"{shelf}_{n}.txt"), os.path.join(HERE, "lopt_hp", shelf, f"{shelf}_{n}.txt"),
              os.path.join(HERE, "submission", shelf, f"{shelf}_{n}.txt")):
        if os.path.exists(p): return F(open(p).readline().split()[1]), p
    return None, None

def run(shelf, n, npy, tag="cand"):
    cont, cs, pat, rp = finalize_circ.info(shelf)
    rec = {int(l.split()[0]): l.split()[1] for l in open(rp) if l.strip()}[n]
    c = np.load(npy); assert c.shape == (n, 2), c.shape
    mp.mp.dps = 80
    C1, r1, info = hpslp.converge([(mp.mpf(repr(float(x))), mp.mpf(repr(float(y)))) for x, y in c], cs, floor_rel=1e-47, t_cap=600)
    mp.mp.dps = 60
    Cr, rr, res, shape = refine_circ.refine(np.array([[float(x), float(y)] for x, y in C1]), cont, cs)
    d = os.path.join(HERE, "cand_hp", shelf); os.makedirs(d, exist_ok=True)
    tmp = os.path.join(d, f"{shelf}_{n}.{tag}.tmp")
    best = None
    for name, CC in ([("newton", Cr)] if res < mp.mpf(10) ** -50 else []) + [("slp", C1)]:
        refine_circ.write_hp(CC, cs, tmp); rn = F(open(tmp).readline().split()[1])
        if best is None or rn > best[1]: best = (name, rn); os.replace(tmp, tmp + ".best")
    os.replace(tmp + ".best", tmp); how, r_new = best
    r_ours, src = ours(shelf, n)
    row = {"shelf": shelf, "N": n, "tag": tag, "how": how, "newton_residual": float(res), "r_new": str(r_new), "r_ours_before": str(r_ours),
           "ours_src": os.path.relpath(src, HERE) if src else None, "gain_vs_ours": float(r_new / r_ours - 1) if r_ours else None,
           "gain_vs_packomania": float(r_new / F(rec) - 1)}
    row["checker_a"] = certify_circ.check(cs, tmp, rec, verbose=False)[0]
    row["checker_b"] = checker_b_verdict(shelf, cs, tmp, n, rec)
    lo = lopt.certify(cs, tmp); row["lopt"] = lo["verdict"]
    row["lopt_detail"] = {k: v for k, v in lo.items() if k in ("rho", "Delta_rel", "backbone", "rattlers", "redundancy")}
    keep = r_ours is None or r_new > r_ours
    row["kept"] = keep
    if keep:
        os.replace(tmp, os.path.join(d, f"{shelf}_{n}.txt"))
        json.dump(row, open(os.path.join(d, f"{shelf}_{n}.json"), "w"), indent=1)
    else: os.remove(tmp)
    print(json.dumps(row, indent=1)); return row

if __name__ == "__main__":
    run(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "cand")
