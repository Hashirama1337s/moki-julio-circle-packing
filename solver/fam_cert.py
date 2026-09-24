"""Certify the closed-form lattice-family and monotone-deletion candidates found by the proofs scout (vision/families_now.out,
re-scanned against our CURRENT records). Per candidate: rebuild the exact family packing (scout's build(), frame [0,1]x[0,h] ->
shifted to our centred frame), or take our packing at the larger N (monotone); delete the surplus circles with the fewest
contacts; SLP polish; if the float radius beats our current one by > 1e-10 -> full chain (cert_candidate.run) into cand_hp/.
Single writer: sizes still queued in the running full sweep (out/night_full_targets.json, not yet in its jsonl) are DEFERRED.
usage: py -3.11 fam_cert.py [--workers=3] [--deferred]   -> out/fam_cert.jsonl
"""
import os, sys
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import re, json, io, contextlib, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SRC = open(os.path.join(HERE, "lattice_families.py"), encoding="utf-8").read()
NS = {}; exec(compile(SRC[:SRC.index("for hh in range(1, 9):")], "recon2_families_defs", "exec"), NS)

def parse():
    out = []
    for l in open(os.path.join(HERE, "vision", "families_now.out")):
        m = re.match(r"(crc_\d+) N=(\d+)\s+FAMILY (rows|cols) k=(\d+) n=(\d+) (\w+) \((\d+) circles\)", l)
        if m: out.append({"shelf": m[1], "N": int(m[2]), "src": "family", "orient": m[3], "k": int(m[4]), "n": int(m[5]), "kind": m[6], "M": int(m[7])})
        m = re.match(r"(\w+) N=(\d+)\s+MONOTONE: delete (\d+) circle\(s\) from N=(\d+)", l)
        if m: out.append({"shelf": m[1], "N": int(m[2]), "src": "monotone", "from": int(m[4])})
    return out

def job(c):
    import slp_circ, finalize_circ, transplant_probe as tp, cert_candidate
    s, N = c["shelf"], c["N"]; cont = finalize_circ.info(s)[0]; h = cont[1] if cont[0] == "rect" else None
    cur = tp.packing(s, N); r_cur = slp_circ.rmin(cur[0], cont)
    if c["src"] == "family":
        W, Hh = (1.0, h) if c["orient"] == "rows" else (h, 1.0)
        N2, r = NS["family"](W, Hh, c["k"], c["n"], c["kind"]); P = NS["build"](W, Hh, c["k"], c["n"], c["kind"], r)
        if c["orient"] == "cols": P = [(y, x) for x, y in P]
        P = np.array(P) - [0.5, h / 2]
    else:
        big = tp.packing(s, c["from"]); P = big[0]
    r0 = slp_circ.rmin(P, cont); row = {**c, "r_family": float(r0), "r_ours": float(r_cur)}
    while len(P) > N:                                             # delete the surplus: fewest near contacts first
        P = np.delete(P, int(np.argmin(tp.near_counts(P, slp_circ.rmin(P, cont), cont))), 0)
    q, rq = slp_circ.polish(P, cont, t_cap=60.0)
    row.update({"r_polished": float(rq), "rel_vs_ours": float(rq / r_cur - 1)})
    if rq > r_cur * (1 + 1e-10):
        npy = os.path.join(HERE, "out", "transplant", f"fam_{s}_{N}_{c['src']}.npy"); np.save(npy, q)
        try:
            with contextlib.redirect_stdout(io.StringIO()): cr = cert_candidate.run(s, N, npy, tag="family_" + c["src"])
            row.update({k: cr.get(k) for k in ("r_new", "gain_vs_ours", "gain_vs_packomania", "claude", "grok", "lopt", "kept")})
        except Exception as e: row["error"] = repr(e)
    return row

if __name__ == "__main__":
    W = int(next((a.split("=")[1] for a in sys.argv if a.startswith("--workers=")), 3))
    C = parse()
    queued = {tuple(t) for t in json.load(open(os.path.join(HERE, "out", "night_full_targets.json")))}
    fp = os.path.join(HERE, "out", "night_probe_full.jsonl")
    done = {(json.loads(l)["shelf"], json.loads(l)["N"]) for l in open(fp)} if os.path.exists(fp) else set()
    prev = os.path.join(HERE, "out", "fam_cert.jsonl")
    tried = {(json.loads(l)["shelf"], json.loads(l)["N"]) for l in open(prev)} if os.path.exists(prev) else set()
    if "--force" in sys.argv: queued = set()               # only after the RP sweep has EXITED (no other writer left)
    run = [c for c in C if (c["shelf"], c["N"]) not in tried and ((c["shelf"], c["N"]) not in queued or (c["shelf"], c["N"]) in done)]
    defer = [(c["shelf"], c["N"]) for c in C if (c["shelf"], c["N"]) in queued and (c["shelf"], c["N"]) not in done]
    print(f"{len(C)} candidates: running {len(run)}, deferred (still queued in the full sweep) {len(defer)} {defer}", flush=True)
    with Pool(W) as pool, open(prev, "a") as f:
        for o in pool.imap_unordered(job, run):
            f.write(json.dumps(o) + "\n"); f.flush()
            ok = o.get("kept") and o.get("claude") == "IMPROVES" and o.get("grok") == "IMPROVES"
            print(f"{o['shelf']} N={o['N']} {o['src']}: polished {o['rel_vs_ours']:+.2e} vs ours"
                  + (f" -> RECORD {o['gain_vs_packomania']:+.2e} vs Packomania (lopt {o['lopt']})" if ok else f" {o.get('claude', '')} {o.get('grok', '')} {o.get('error', '')}"), flush=True)
    print("DONE", flush=True)
