"""Sealed probe (vision/CLAUDE-SEALED-TRANSPLANT.md): seed N from our N+1 (delete a circle) and N-1 (insert into a hole).
-> out/transplant_probe.jsonl (one line per seed), candidate packings out/transplant/<shelf>_<N>_<tag>.npy
usage: py -3.11 transplant_probe.py [--workers=3]
"""
import os, sys
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import json, glob, time, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SHELVES = ["crt", "ccq", "crc_500", "crc_800"]
MODE2 = "--mode=2" in sys.argv                      # sealed replication: vision/CLAUDE-SEALED-TRANSPLANT2.md
OUT = os.path.join(HERE, "out", "transplant2_probe.jsonl" if MODE2 else "transplant_probe.jsonl")

def load_txt(p):
    L = [l.split() for l in open(p) if l.strip()]
    return np.array([[float(a), float(b)] for a, b in L[1:]]), float(L[0][1])

def packing(shelf, n):
    """Our best circle-form packing at n: v1.1 converged file, else v1.0 submission, else Packomania coords."""
    import finalize_circ
    for p in (os.path.join(HERE, "cand_hp", shelf, f"{shelf}_{n}.txt"), os.path.join(HERE, "lopt_hp", shelf, f"{shelf}_{n}.txt"),
              os.path.join(HERE, "submission", shelf, f"{shelf}_{n}.txt")):
        if os.path.exists(p): return load_txt(p) + ("ours",)
    pat = finalize_circ.info(shelf)[2]
    if not os.path.exists(pat.format(n)): return None
    c = np.array([[float(l.split()[1]), float(l.split()[2])] for l in open(pat.format(n)) if l.strip()])
    if len(c) != n: return None          # 09-24: csc sizes whose coordinates Packomania has not published (empty files)
    return c, None, "packomania"

def targets(shelves=SHELVES, offset=False, exclude=()):
    out = []
    for s in shelves:
        ns = sorted(int(os.path.basename(p)[:-4].split("_")[-1]) for p in glob.glob(os.path.join(HERE, "submission", s, f"{s}_*.txt")))
        ok = [n for n in ns if packing(s, n - 1) is not None and packing(s, n + 1) is not None]
        k = max(1, len(ok) // 8); pick = ok[k // 2::k] if offset else ok[::k]
        out += [(s, n) for n in pick if (s, n) not in exclude][:8]
    return out

def near_counts(c, r, cont):
    import slp_circ
    g, _ = slp_circ.walls(c, cont); D = np.sqrt(((c[:, None] - c[None]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf)
    return (D <= 2 * r * (1 + 1e-7)).sum(1) + (g <= r * (1 + 1e-7)).sum(1)

def holes(c, cont, k=3, seed=0):
    import slp_circ
    rng = np.random.default_rng(seed); lo, hi = c.min(0) - 0.05, c.max(0) + 0.05
    P = slp_circ.repair(rng.uniform(lo, hi, (20000, 2)), cont)
    g, _ = slp_circ.walls(P, cont); D = np.sqrt(((P[:, None] - c[None]) ** 2).sum(-1)).min(1)
    clr = np.minimum(g.min(1), D / 2); pick = []
    for i in np.argsort(-clr):
        if all(np.hypot(*(P[i] - P[j])) > 2 * clr[i] for j in pick): pick.append(i)
        if len(pick) == k: break
    return P[pick]

def job(a):
    shelf, n = a
    import slp_circ, finalize_circ
    cont = finalize_circ.info(shelf)[0]; c_n, r_best, _ = packing(shelf, n)
    if r_best is None: r_best = slp_circ.rmin(c_n, cont)
    res = []
    up = packing(shelf, n + 1); c_up = up[0]; r_up = up[1] or slp_circ.rmin(c_up, cont)
    nc = near_counts(c_up, r_up, cont); order = [int(i) for i in np.argsort(nc, kind="stable")][:6]
    for i in order:
        q = np.delete(c_up, i, 0); t = time.process_time(); c, r = slp_circ.polish(q, cont)
        res.append({"shelf": shelf, "n": n, "seed": f"del{i}(nc={int(nc[i])})", "src": up[2], "r": float(r), "r_best": float(r_best),
                    "rel": float(r / r_best - 1), "cpu": time.process_time() - t})
        if r > r_best * (1 + 1e-10): np.save(os.path.join(HERE, "out", "transplant", f"{shelf}_{n}_del{i}.npy"), c)
    dn = packing(shelf, n - 1); c_dn = dn[0]
    for j, h in enumerate(holes(c_dn, cont)):
        q = np.vstack([c_dn, h]); t = time.process_time(); c, r = slp_circ.polish(q, cont)
        res.append({"shelf": shelf, "n": n, "seed": f"ins{j}", "src": dn[2], "r": float(r), "r_best": float(r_best),
                    "rel": float(r / r_best - 1), "cpu": time.process_time() - t})
        if r > r_best * (1 + 1e-10): np.save(os.path.join(HERE, "out", "transplant", f"{shelf}_{n}_ins{j}.npy"), c)
    return res

if __name__ == "__main__":
    W = int(next((x.split("=")[1] for x in sys.argv if x.startswith("--workers=")), 3))
    os.makedirs(os.path.join(HERE, "out", "transplant"), exist_ok=True)
    if MODE2:
        prev = {tuple(t) for t in json.load(open(os.path.join(HERE, "out", "transplant_targets.json")))}
        T = targets([f"crc_{k}" for k in range(100, 900, 100)], offset=True, exclude=prev)
        json.dump(T, open(os.path.join(HERE, "out", "transplant2_targets.json"), "w"))
    else:
        T = targets(); json.dump(T, open(os.path.join(HERE, "out", "transplant_targets.json"), "w"))
    print(f"{len(T)} targets: {T}", flush=True); t0 = time.time()
    with Pool(W) as pool, open(OUT, "a") as f:
        for res in pool.imap_unordered(job, T):
            for o in res: f.write(json.dumps(o) + "\n")
            f.flush(); b = max(res, key=lambda o: o["rel"])
            print(f"[{time.time() - t0:.0f}s] {b['shelf']} N={b['n']}: best seed {b['seed']} rel {b['rel']:+.2e}"
                  + ("  HIT" if b["rel"] > 1e-10 else ""), flush=True)
    print("DONE", flush=True)
