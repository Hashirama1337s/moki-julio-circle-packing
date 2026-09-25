# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""H1 (vision/MERGED-HATS.md): FLAGSHIP monotone harvest on Packomania's csq (unit square) and cci (unit disc).
For every listed N whose radius is beaten by some larger listed M's radius: take Packomania's own M packing (data/big/, downloaded
2026-09-24), delete M - N circles (fewest near contacts first), write a certificate (radius = the exact min clearance of what
remains, rounded down) and check it with certify_big.py. Nothing is polished here (that is the next stage).
Reference radii: the LIVE table pages saved 2026-09-24 (data/refs/packomania_<s>_2026-09-24.html).
usage: py -3.11 hats_flagship.py [--sample=10]   -> cand_big/<s>/<s>_<N>.txt + out/hats_flagship.jsonl
"""
import os, sys, re, json, math, random
import numpy as np
from scipy.spatial import cKDTree
from decimal import Decimal, ROUND_FLOOR, getcontext
HERE = os.path.dirname(os.path.abspath(__file__)); P = lambda *a: os.path.join(HERE, *a); sys.path.insert(0, HERE)
import certify_big
getcontext().prec = 60
CONT = {"csq": "square", "cci": "circle"}

def table(s):
    t = open(P("data", "refs", f"packomania_{s}_2026-09-24.html"), encoding="utf-8", errors="ignore").read()
    return {int(n): v for n, v in re.findall(r'name="[a-z]+(\d+)">\s*<strong>\s*\d+</strong></a></td>\s*<td>(?:<strong>)?\s*([0-9.]+)', t)}

def coords(s, m):
    return [l.split()[1:3] for l in open(P("data", "big", s, f"{s}{m}.txt")) if l.strip() and not l.lstrip().startswith("#")]

def clearance(s, c):
    d, _ = cKDTree(c).query(c, k=2); half = d[:, 1].min() / 2
    wall = (0.5 - np.abs(c)).min() if s == "csq" else (1 - np.linalg.norm(c, axis=1)).min()
    return min(half, wall)

def build(s, n, m, T):
    rows = coords(s, m); c = np.array([[float(x), float(y)] for x, y in rows]); r_m = clearance(s, c)
    near = cKDTree(c).query_pairs(2 * r_m * (1 + 1e-9), output_type="ndarray")
    deg = np.bincount(near.ravel(), minlength=len(c))
    if s == "csq": deg = deg + (0.5 - np.abs(c) <= r_m * (1 + 1e-9)).sum(1)
    else: deg = deg + (1 - np.linalg.norm(c, axis=1) <= r_m * (1 + 1e-9))
    keep = np.sort(np.argsort(deg, kind="stable")[m - n:])               # drop the m - n least-contacted circles
    ck = c[keep]; r_f = clearance(s, ck)
    r_claim = (Decimal(repr(float(r_f))) * (1 - Decimal(10) ** -12)).quantize(Decimal(10) ** -25, rounding=ROUND_FLOOR)
    d = P("cand_big", s); os.makedirs(d, exist_ok=True); f = os.path.join(d, f"{s}_{n}.txt")
    open(f, "w", newline="\n").write(f"r {r_claim}\n" + "".join(f"{rows[i][0]} {rows[i][1]}\n" for i in keep))
    v = certify_big.check(CONT[s], f, n, T[n])
    return {"table": s, "N": n, "from_M": m, "deleted": m - n, "r_pub": T[n], "r_claim": str(r_claim),
            "gain_vs_packomania": float(Decimal(str(r_claim)) / Decimal(T[n]) - 1), "checker_a": v, "file": os.path.relpath(f, HERE)}

def targets(s):
    T = table(s); ns = sorted(T); out = []; run = (0.0, None)
    for n in reversed(ns):
        if run[1] is not None and run[0] > float(T[n]) * (1 + 1e-10): out.append((n, run[1]))
        if float(T[n]) > run[0]: run = (float(T[n]), n)
    return T, sorted(out)

if __name__ == "__main__":
    k = next((int(a.split("=")[1]) for a in sys.argv if a.startswith("--sample=")), None)
    with open(P("out", "hats_flagship.jsonl"), "a") as log:
        for s in ("csq", "cci"):
            T, tg = targets(s)
            if k: random.seed(20260924); tg = sorted(random.sample(tg, k))
            print(s, "targets", len(tg), flush=True)
            for n, m in tg:
                row = build(s, n, m, T); log.write(json.dumps(row) + "\n"); print(json.dumps(row), flush=True)
