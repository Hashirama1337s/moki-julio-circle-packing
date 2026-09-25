"""Lattice-family seed sweep over every rectangle size (the proofs scout's closed-form row families: alt / shift / sq, rows or
columns; recon2_families.py). Per size N: the best families with N <= M <= N+6 circles (surplus deleted, fewest contacts first)
and M = N-1, N-2 (circles added in the largest holes), keeping only families within 5e-3 of our current radius; SLP polish
each; the best strict float gain over our packing -> full chain (cert_candidate.run) into cand_hp/.
Single writer: sizes still queued in the running full sweep are skipped (re-run with --all-queued after it ends).
usage: py -3.11 fam_sweep.py [--workers=6] [--stop=HH:MM] [--shelves=crc_700,...] [--include-queued] [--exclude-pending]
       [--tol=5e-3] [--mmax=6] [--out=fam_sweep.jsonl]
"""
import os, sys
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import json, time, io, contextlib, datetime, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SRC = open(os.path.join(HERE, "lattice_families.py"), encoding="utf-8").read()
NS = {}; exec(compile(SRC[:SRC.index("for hh in range(1, 9):")], "recon2_families_defs", "exec"), NS)
ARG = lambda k, d: next((x.split("=", 1)[1] for x in sys.argv if x.startswith(f"--{k}=")), d)
STOP = datetime.datetime.fromisoformat(os.environ["FAM_STOP"]) if os.environ.get("FAM_STOP") else None

def families(shelf):
    """All (M, r, spec) for one rectangle table, M <= max table N + 8."""
    h = int(shelf.split("_")[1]) / 1000; nmax = max(int(l.split()[0]) for l in open(os.path.join(HERE, "shelves", shelf, "radius.txt")) if l.strip())
    out = []
    for orient, (W, H) in (("rows", (1.0, h)), ("cols", (h, 1.0))):
        for k in range(1, 60):
            for n in range(1, 480):
                if k * n > 2 * (nmax + 8) + 60: break
                for kind in ("alt", "shift", "sq"):
                    res = NS["family"](W, H, k, n, kind)
                    if res and res[0] <= nmax + 8: out.append((res[0], res[1], (orient, k, n, kind)))
    return shelf, out

def job(a):
    shelf, N, cands = a
    if STOP and datetime.datetime.now() >= STOP: return {"shelf": shelf, "N": N, "skipped": "stop"}
    import slp_circ, finalize_circ, transplant_probe as tp, cert_candidate
    cont = finalize_circ.info(shelf)[0]; h = cont[1]
    cur = tp.packing(shelf, N); r_cur = slp_circ.rmin(cur[0], cont); best = (r_cur, None, None); tried = []
    for M, r, (orient, k, n, kind) in cands:
        W, Hh = (1.0, h) if orient == "rows" else (h, 1.0)
        if kind == "alt" and n < 2: continue                        # 02:08: build() divides by n - 1 (k = 1 families)
        try: P = NS["build"](W, Hh, k, n, kind, r)
        except (ZeroDivisionError, ValueError) as e: tried.append((f"{orient}_{kind}_k{k}_n{n}", "build error " + repr(e))); continue
        if orient == "cols": P = [(y, x) for x, y in P]
        P = np.array(P, dtype=float) - [0.5, h / 2]
        while len(P) > N: P = np.delete(P, int(np.argmin(tp.near_counts(P, slp_circ.rmin(P, cont), cont))), 0)
        if len(P) < N: P = np.vstack([P, tp.holes(P, cont, k=N - len(P))])
        q, rq = slp_circ.polish(P, cont, t_cap=60.0); tag = f"{orient}_{kind}_k{k}_n{n}_M{M}"
        tried.append((tag, round(float(rq / r_cur - 1), 9)))
        if rq > best[0]: best = (rq, q, tag)
    row = {"shelf": shelf, "N": N, "rel_float": float(best[0] / r_cur - 1), "seed": best[2], "tried": tried}
    if best[1] is not None and best[0] > r_cur * (1 + 1e-10):
        npy = os.path.join(HERE, "out", "transplant", f"famsw_{shelf}_{N}.npy"); np.save(npy, best[1])
        try:
            with contextlib.redirect_stdout(io.StringIO()): cr = cert_candidate.run(shelf, N, npy, tag="famsweep_" + best[2])
            row.update({k: cr.get(k) for k in ("r_new", "gain_vs_ours", "gain_vs_packomania", "checker_a", "checker_b", "lopt", "kept")})
        except Exception as e: row["error"] = repr(e)
    return row

def safe_job(a):
    """Any failure in one size becomes an error row; it never takes the pool (and in-flight certifications) down."""
    try: return job(a)
    except Exception as e: return {"shelf": a[0], "N": a[1], "error": repr(e)}

if __name__ == "__main__":
    W = int(ARG("workers", 6)); shelves = ARG("shelves", ",".join(f"crc_{k}" for k in (700, 800, 600, 400, 300, 200, 500, 100))).split(",")
    if ARG("stop", None):
        now = datetime.datetime.now(); hh, mm = map(int, ARG("stop", None).split(":"))
        s = now.replace(hour=hh, minute=mm, second=0, microsecond=0); s += datetime.timedelta(days=1) if s <= now else datetime.timedelta(0)
        os.environ["FAM_STOP"] = s.isoformat()
    import transplant_probe as tp, slp_circ, finalize_circ
    queued = {tuple(t) for t in json.load(open(os.path.join(HERE, "out", "night_full_targets.json")))}
    fp = os.path.join(HERE, "out", "night_probe_full.jsonl")
    done_full = {(json.loads(l)["shelf"], json.loads(l)["N"]) for l in open(fp)} if os.path.exists(fp) else set()
    OUTF = ARG("out", "fam_sweep.jsonl"); TOL = float(ARG("tol", "5e-3")); MX = int(ARG("mmax", "6"))
    NLO, NHI = map(int, ARG("nrange", "1,100000").split(","))
    prev = {(json.loads(l)["shelf"], json.loads(l)["N"]) for p in ("fam_cert.jsonl", OUTF) if os.path.exists(os.path.join(HERE, "out", p))
            for l in open(os.path.join(HERE, "out", p))}
    pending = set()                                         # --exclude-pending: sizes another running sweep has queued but not done
    if "--exclude-pending" in sys.argv:
        for tl, jl in (("night_full_targets.json", "night_probe_full.jsonl"), ("night_P_targets.json", "night_probe_Ponly.jsonl")):
            q = {tuple(x) for x in json.load(open(os.path.join(HERE, "out", tl)))}
            d = {(json.loads(l)["shelf"], json.loads(l)["N"]) for l in open(os.path.join(HERE, "out", jl))} if os.path.exists(os.path.join(HERE, "out", jl)) else set()
            pending |= q - d
    t0 = time.time()
    with Pool(W) as pool: fams = dict(pool.map(families, shelves))
    print(f"families built in {time.time() - t0:.0f}s: " + ", ".join(f"{s} {len(v)}" for s, v in fams.items()), flush=True)
    T = []
    for s in shelves:
        cont = finalize_circ.info(s)[0]
        for N in sorted({int(l.split()[0]) for l in open(os.path.join(HERE, "shelves", s, "radius.txt")) if l.strip()}, reverse=True):
            if not (NLO <= N <= NHI) or (s, N) in prev or (s, N) in pending or ((s, N) in queued and (s, N) not in done_full and "--include-queued" not in sys.argv): continue
            cur = tp.packing(s, N); r_cur = slp_circ.rmin(cur[0], cont)
            near = sorted([f for f in fams[s] if N - 2 <= f[0] <= N + MX and f[1] >= r_cur * (1 - TOL)], key=lambda f: -f[1])
            pick, seen = [], set()
            for f in near:
                key = (f[0], round(f[1], 12))
                if key not in seen: seen.add(key); pick.append(f)
                if len(pick) == 4: break
            if pick: T.append((s, N, pick))
    print(f"{len(T)} sizes with a family within {TOL:g} of our radius (M <= N+{MX}); {W} workers; pending elsewhere {len(pending)}", flush=True); hits = 0
    with Pool(W) as pool, open(os.path.join(HERE, "out", OUTF), "a") as f:
        for o in pool.imap_unordered(safe_job, T, chunksize=1):
            f.write(json.dumps(o) + "\n"); f.flush()
            ok = o.get("kept") and o.get("checker_a") == "IMPROVES" and o.get("checker_b") == "IMPROVES"; hits += bool(ok)
            if ok or o.get("error"):
                go = o.get("gain_vs_ours")                        # None at a size we had no record for (02:05 crash, fixed)
                print(f"[{time.time() - t0:.0f}s] {o['shelf']} N={o['N']} {o.get('seed')}: " + (f"RECORD {o['gain_vs_packomania']:+.2e} vs Packomania, "
                      + (f"{go:+.2e} vs ours" if go is not None else "NEW SIZE") + f" (lopt {o['lopt']})" if ok else f"ERROR {o['error']}"), flush=True)
    print(f"DONE {hits} records / {len(T)} sizes, {time.time() - t0:.0f}s", flush=True)
