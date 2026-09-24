"""Four hats, zero-cost check results turned into certified candidates (vision/CLAUDE-SEALED-HATS.md, out/hats_zero.json):
 - TAO monotone: crc_700 N in {373, 367, 401, 249, 214} are beaten by our own N+1 packing minus one circle -> transplant_sweep.job
   (delete / insert seeds from ours and Packomania's neighbours, float polish, full exact chain if it beats ours).
 - EINSTEIN mirror / TAO stacking: crc_800 N=282 = our crc_400 N=141 record + its mirror image across the long side, polished,
   then the full exact chain (cert_candidate.run: both checkers, lopt; kept only if it beats ours).
usage: py -3.11 hats_seeds.py   -> out/hats_seeds.jsonl
"""
import os, sys
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import json, time, io, contextlib, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "out", "hats_seeds.jsonl")

def stack_job(_):
    import slp_circ, finalize_circ, transplant_probe as tp, cert_candidate
    c4, r4, src = tp.packing("crc_400", 141)
    q = np.vstack([c4 + [0.0, 0.2], c4 * [1.0, -1.0] - [0.0, 0.2]])      # 1 x 0.4 at y in [0, 0.4] and its mirror in [-0.4, 0]
    cont = finalize_circ.info("crc_800")[0]
    r_seed = slp_circ.rmin(q, cont)
    c, r = slp_circ.polish(q, cont, t_cap=120.0)
    if r < r_seed: c, r = q, r_seed
    npy = os.path.join(HERE, "out", "hats_crc_800_282.npy"); np.save(npy, c)
    with contextlib.redirect_stdout(io.StringIO()): row = cert_candidate.run("crc_800", 282, npy, "hats_stack")
    return {"shelf": "crc_800", "N": 282, "seed": f"crc_400(141, {src}) x2 mirror", "r_seed": r_seed, "r_float": r, **{k: row.get(k) for k in
            ("r_new", "r_ours_before", "gain_vs_ours", "gain_vs_packomania", "claude", "grok", "lopt", "kept")}}

def mono_job(n):
    import transplant_sweep
    with contextlib.redirect_stdout(io.StringIO()): return {"seed": "monotone N+1 delete / N-1 insert", **transplant_sweep.job(("crc_700", n, time.time() + 3600))}

def run(a):
    t = time.time()
    try: row = stack_job(a) if a == "stack" else mono_job(a)
    except Exception as e: row = {"job": a, "ERROR": repr(e)}
    row["secs"] = round(time.time() - t, 1); return row

if __name__ == "__main__":
    jobs = [373, 367, "stack", 401, 249, 214]
    with Pool(2) as p:
        for row in p.imap_unordered(run, jobs):
            open(OUT, "a").write(json.dumps(row) + "\n"); print(json.dumps(row), flush=True)
    print("DONE", flush=True)
