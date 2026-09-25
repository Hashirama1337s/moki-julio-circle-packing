"""csq cells N <= 2000 whose best packing so far is only a DELETION from a larger published N (not ours to claim): polish it with
the KD-tree SLP (slp_big) and keep it only if it rises ABOVE Packomania's monotone envelope and Amore (lai_compare.csq_gate) — then
the packing is our work. Certificate: float centres as exact decimals, radius = computed clearance * (1 - 1e-12), rounded down;
BOTH O(N) exact checkers. -> cand_big_hp/csq/csq_<N>.txt + .json (the build's third source).
usage: py -3.11 hats_csq_polish.py [--workers=2] [--cap=150]   -> out/hats_csq_polish.jsonl
"""
import os, sys
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import json, time, subprocess
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ARG = lambda k, d: next((a.split("=", 1)[1] for a in sys.argv if a.startswith(f"--{k}=")), d)

def job(n):
    import numpy as np, mpmath as mp, slp_big, certify_big, lai_compare, hats_flagship as HF
    from decimal import Decimal, ROUND_FLOOR, getcontext
    getcontext().prec = 60; t0 = time.time()
    L = [l.split() for l in open(os.path.join(HERE, "cand_big", "csq", f"csq_{n}.txt")) if l.strip()]
    c = np.array([[float(a), float(b)] for a, b in L[1:]]); r_seed = float(L[0][1])
    c2, r2 = slp_big.polish(c, ("rect", 1.0), t_cap=float(ARG("cap", 150)))
    row = {"N": n, "r_seed": r_seed, "r_polish": r2, "rel": r2 / r_seed - 1, "secs": round(time.time() - t0, 1)}
    r_claim = (Decimal(repr(float(slp_big.rmin(c2, ("rect", 1.0))))) * (1 - Decimal(10) ** -12)).quantize(Decimal(10) ** -25, rounding=ROUND_FLOOR)
    ok, why = lai_compare.csq_gate(n, mp.mpf(str(r_claim))); row["gate"] = why or "claimable"
    if not ok: return row
    d = os.path.join(HERE, "cand_big_hp", "csq"); os.makedirs(d, exist_ok=True); f = os.path.join(d, f"csq_{n}.txt")
    if os.path.exists(f) and Decimal(open(f).readline().split()[1]) >= r_claim: row["kept"] = False; return row
    tmp = f + ".polish.tmp"
    open(tmp, "w", newline="\n").write(f"r {r_claim}\n" + "".join(f"{format(Decimal(repr(float(x))), 'f')} {format(Decimal(repr(float(y))), 'f')}\n" for x, y in c2))
    rec = HF.table("csq")[n]
    a = certify_big.check("square", tmp, n, rec)
    out = subprocess.run([sys.executable, os.path.join(HERE, "grok", "verify_exact_big.py"), "square", tmp, str(n), rec], capture_output=True, text=True).stdout
    b = [l for l in out.splitlines() if l.startswith("VERDICT")]; b = b[0].split(":")[1].strip() if b else "ERROR"
    row.update({"claude": a, "grok": b, "r_claim": str(r_claim)})
    if a == "IMPROVES" and b == "IMPROVES":
        os.replace(tmp, f); row["kept"] = True
        json.dump({"shelf": "csq", "N": n, "tag": "csq_polish", "how": "deletion seed from Packomania's larger N, polished above the envelope (slp_big)",
                   "r_new": str(r_claim), "claude": a, "grok": b, "lopt": "NOT_ATTEMPTED"}, open(f[:-4] + ".json", "w"), indent=1)
    else: os.remove(tmp); row["kept"] = False
    return row

if __name__ == "__main__":
    import glob
    fam = {int(os.path.basename(p)[4:-4]) for p in glob.glob(os.path.join(HERE, "cand_big_hp", "csq", "csq_*.txt"))}
    T = sorted(n for n in (int(os.path.basename(p)[4:-4]) for p in glob.glob(os.path.join(HERE, "cand_big", "csq", "csq_*.txt")))
               if n <= 2000 and n not in fam)
    print(len(T), "csq deletion cells N <= 2000 to polish", flush=True)
    with Pool(int(ARG("workers", 2))) as p, open(os.path.join(HERE, "out", "hats_csq_polish.jsonl"), "a") as f:
        for row in p.imap_unordered(job, T):
            f.write(json.dumps(row) + "\n"); f.flush()
            print(f"N={row['N']} rel {row['rel']:+.2e} {row['gate'][:60]} {'KEPT' if row.get('kept') else ''} {row['secs']}s", flush=True)
    print("DONE", flush=True)
