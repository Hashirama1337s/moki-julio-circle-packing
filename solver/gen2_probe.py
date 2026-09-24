"""Second-generation transplants (Claude stage-2 #3 + Grok #3 force-ranked relocation). For sizes where the one-step sweep found
NOTHING (rel_float <= 0 in the stage-2 logs), try:
  A 'two-step': N-2 + two circles into the two largest holes;  N+2 minus the two fewest-contact circles.
  B 'relocate1': move the weakest circle of OUR packing at N into the largest hole (weakest = smallest total contact force from the
    local-optimality certificate when one exists, else fewest near contacts).
  C 'relocate2': the same with the two weakest circles into the two largest holes.
Each seed -> slp_circ.polish; the best strict gain -> full chain (cert_candidate.run).  Targets are fixed by the sealed rule
(vision/CLAUDE-SEALED-GEN2.md) and written to out/gen2_targets.json before any polish.   -> out/gen2_probe.jsonl
usage: py -3.11 gen2_probe.py [--workers=12] [--per-shelf=6]
"""
import os, sys
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import json, glob, time, io, contextlib, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SHELVES = ["crt", "ccq"] + [f"crc_{k}" for k in range(100, 900, 100)]

def weakest(shelf, n, c, r, cont, k):
    """Indices of the k weakest circles: least total stress from lopt_certs (if certified), else fewest near contacts.
    NOTE (found 2026-09-24): lopt_certs belong to the v1.1 submission files; where our packing at n changed after v1.1 the
    stress indices do not match these circles, so the pick there is effectively arbitrary. Results are unaffected (every
    candidate is certified from scratch); night_probe.order() checks that the certificate belongs to the packing."""
    import transplant_probe as tp
    lc = os.path.join(HERE, "lopt_certs", shelf, f"{shelf}_{n}.json")
    if os.path.exists(lc):
        d = json.load(open(lc)); load = np.zeros(n)
        for s, lamh in zip(d["S"], d["lambda"]):
            lam = float.fromhex(lamh)
            for i in (s[1:3] if s[0] == "p" else s[1:2]): load[i] += lam
        return [int(i) for i in np.argsort(load, kind="stable")[:k]], "force"
    nc = tp.near_counts(c, r, cont)
    return [int(i) for i in np.argsort(nc, kind="stable")[:k]], "contacts"

def job(a):
    shelf, n = a
    import slp_circ, finalize_circ, transplant_probe as tp, cert_candidate
    cont = finalize_circ.info(shelf)[0]
    cur = tp.packing(shelf, n); r_best = cur[1] or slp_circ.rmin(cur[0], cont)
    best = (r_best, None, None); tried = []
    def trial(q, tag):
        nonlocal best
        c, r = slp_circ.polish(q, cont, t_cap=60.0); tried.append((tag, float(r / r_best - 1)))
        if r > best[0]: best = (r, c, tag)
    dn2, up2 = tp.packing(shelf, n - 2), tp.packing(shelf, n + 2)
    if dn2 is not None:
        trial(np.vstack([dn2[0], tp.holes(dn2[0], cont, k=2)]), "A_ins2")
    if up2 is not None:
        c2 = up2[0]; w, _ = weakest(shelf, n + 2, c2, up2[1] or slp_circ.rmin(c2, cont), cont, 2)
        trial(np.delete(c2, w, 0), "A_del2")
    c0 = cur[0]; r0 = cur[1] or slp_circ.rmin(c0, cont)
    for k in (1, 2):
        w, how = weakest(shelf, n, c0, r0, cont, k)
        base = np.delete(c0, w, 0); trial(np.vstack([base, tp.holes(base, cont, k=k)]), f"{'B' if k == 1 else 'C'}_reloc{k}_{how}")
    row = {"shelf": shelf, "N": n, "tried": tried, "seed": best[2], "rel_float": float(best[0] / r_best - 1)}
    if best[1] is None or best[0] <= r_best * (1 + 1e-10): return row
    npy = os.path.join(HERE, "out", "transplant", f"gen2_{shelf}_{n}_{best[2]}.npy"); np.save(npy, best[1])
    try:
        with contextlib.redirect_stdout(io.StringIO()): cr = cert_candidate.run(shelf, n, npy, tag="gen2_" + best[2])
        row.update({k: cr.get(k) for k in ("gain_vs_ours", "gain_vs_packomania", "claude", "grok", "lopt", "kept")})
    except Exception as e: row["error"] = repr(e)
    return row

if __name__ == "__main__":
    W = int(next((x.split("=")[1] for x in sys.argv if x.startswith("--workers=")), 12))
    per = int(next((x.split("=")[1] for x in sys.argv if x.startswith("--per-shelf=")), 6))
    only = next((x.split("=")[1].split(",") for x in sys.argv if x.startswith("--shelves=")), None)   # run the probe in halves
    tag = next((x.split("=")[1] for x in sys.argv if x.startswith("--tag=")), "")                      # (no two writers per shelf)
    logs = [os.path.join(HERE, "out", f) for f in ("transplant_sweep.jsonl", "transplant_sweep_tq.jsonl", "sweep2_A.jsonl", "sweep2_B.jsonl")]
    miss = {}
    for p in logs:
        if os.path.exists(p):
            for l in open(p):
                r = json.loads(l)
                if "rel_float" in r and r["rel_float"] <= 0 and r.get("base") == "ours": miss[(r["shelf"], r["N"])] = 1
    T = []
    for s in (only or SHELVES):
        ns = sorted(n for (sh, n) in miss if sh == s)
        if "--all" in sys.argv: T += [(s, n) for n in ns]                       # full run after the sealed PASS
        else: k = max(1, len(ns) // per); T += [(s, n) for n in ns[k // 2::k][:per]]
    if "--all" in sys.argv:                                                     # skip sizes the probe already tried
        tried = set()
        for f in glob.glob(os.path.join(HERE, "out", "gen2_probe*.jsonl")):
            tried |= {(json.loads(l)["shelf"], json.loads(l)["N"]) for l in open(f)}
        T = [t for t in T if t not in tried]
    json.dump(T, open(os.path.join(HERE, "out", f"gen2_targets{tag}.json"), "w"))
    print(f"{len(T)} targets (one-step misses on our own packings), {W} workers", flush=True); t0 = time.time(); hits = 0
    with Pool(W) as pool, open(os.path.join(HERE, "out", f"gen2_probe{tag}.jsonl"), "a") as f:
        for o in pool.imap_unordered(job, T):
            f.write(json.dumps(o) + "\n"); f.flush()
            kept = o.get("kept") and o.get("claude") == "IMPROVES" and o.get("grok") == "IMPROVES"; hits += bool(kept)
            if kept or o.get("error"):
                print(f"[{time.time() - t0:.0f}s] {o['shelf']} N={o['N']} {o['seed']}: "
                      + (f"KEPT {o['gain_vs_packomania']:+.2e} vs Packomania" if kept else f"ERROR {o['error']}"), flush=True)
    print(f"DONE hits {hits}/{len(T)} -> {'PASS' if hits >= 3 else 'NULL'}", flush=True)
