"""Converge every record to its EXACT local optimum, re-certify it, and run the local-optimality certificate (lopt.py).
Per packing:  submission cert -> hpslp.converge (mixed-precision SLP, ~1e-17) -> refine_circ.refine (Newton on the identified
contacts, residual < 1e-50) -> write_hp (45 digits, admissible r rounded DOWN) -> keep ONLY if r_new >= r_old ->
Claude's exact checker + Grok's exact checker (both must say IMPROVES vs the Packomania record) -> lopt.certify.
Nothing in submission/ is touched: results go to lopt_hp/<shelf>/<shelf>_<N>.txt + .json (atomic; existing .json = skip).
usage: py -3.11 lopt_fix.py [--workers=12] [--only=shelf:N,shelf:N] [shelf ...]
"""
import os, sys
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import glob, json, time, shutil, subprocess, traceback
from multiprocessing import Pool
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SHELVES = ["crt", "ccq"] + [f"crc_{k}" for k in range(100, 900, 100)]
ESCAPE_ROUNDS = int(os.environ.get("ESCAPE_ROUNDS", "0"))

def grok_verdict(shelf, cs, path, n, rec):
    g = os.path.join(HERE, "grok", "verify_exact.py" if shelf == "crt" else "verify_exact2.py")
    argv = [sys.executable, g, path, str(n), rec] if shelf == "crt" else [sys.executable, g, cs, path, str(n), rec]
    out = subprocess.run(argv, capture_output=True, text=True).stdout
    vl = [l for l in out.splitlines() if l.startswith("VERDICT")]
    return vl[0].split(":")[1].strip() if vl else "ERROR"

def job(args):
    shelf, n = args
    import numpy as np, mpmath as mp
    import hpslp, refine_circ, certify_circ, lopt, finalize_circ
    d = os.path.join(HERE, "lopt_hp", shelf); os.makedirs(d, exist_ok=True)
    out_txt, out_json = os.path.join(d, f"{shelf}_{n}.txt"), os.path.join(d, f"{shelf}_{n}.json")
    if os.path.exists(out_json): return shelf, n, "skip", None, 0.0
    t0 = time.time(); row = {"shelf": shelf, "N": n}
    try:
        cont, cs, pat, rp = finalize_circ.info(shelf)
        rec = {int(l.split()[0]): l.split()[1] for l in open(rp) if l.strip()}[n]
        src = os.path.join(HERE, "submission", shelf, f"{shelf}_{n}.txt")
        L = [l.split() for l in open(src) if l.strip()]; r_old = F(L[0][1])
        tmp = out_txt + ".tmp"

        def finish(C1, path, r_min):
            """Newton on the identified contacts, then 45-digit certificate; returns (how, info) or ('none', ...) if r < r_min."""
            mp.mp.dps = 60                                           # refine_circ's working precision (as in finalize_circ)
            Cr, rr, res, shape = refine_circ.refine(np.array([[float(x), float(y)] for x, y in C1]), cont, cs)
            for name, CC in ([("newton", Cr)] if res < mp.mpf(10) ** -50 else []) + [("slp", C1)]:
                refine_circ.write_hp(CC, cs, path)
                if F(open(path).readline().split()[1]) >= r_min: return name, {"newton_residual": float(res), "contacts": shape}
            return "none", {"newton_residual": float(res), "contacts": shape}

        mp.mp.dps = 80
        C1, r1, info = hpslp.converge([(mp.mpf(a), mp.mpf(b)) for a, b in L[1:]], cs, floor_rel=1e-47)
        row["slp"] = {"iters": info["iters"], "accepted": info["accepted"], "gain_rel": float(info["gain"] / r1)}
        how, extra = finish(C1, tmp, r_old); row.update(extra)
        if how == "none": shutil.copyfile(src, tmp)                  # never lose: keep the published certificate
        lo = lopt.certify(cs, tmp)
        # FLEX ESCAPE: in straight-walled containers a flex that must carry stress fails the second-order condition (the pair
        # terms of the Lagrangian Hessian are positive), so a better packing lies next to it. Nudge +-s along the flex, converge,
        # keep only a strict gain; repeat while NOT_RIGID (max 2 rounds; s = 1e-3, 1e-5 r). Applied to the quadrant too (empirical there).
        row["escape"] = []
        for rnd in range(ESCAPE_ROUNDS):   # 18:25 restart: 0 (escapes moved to the dedicated flex walk)
            if lo["verdict"] != "NOT_RIGID" or "vec" not in (lo.get("flex") or {}): break
            Lc = [l.split() for l in open(tmp) if l.strip()]; r_cur = F(Lc[0][1]); mp.mp.dps = 80
            Cc = [[mp.mpf(a), mp.mpf(b)] for a, b in Lc[1:]]; best = (mp.mpf(Lc[0][1]), None)
            for srel in (1e-3, 1e-5):
                for sg in (1, -1):
                    Cp = [list(c) for c in Cc]
                    for i, xy, w in lo["flex"]["vec"]: Cp[i][xy] += sg * srel * mp.mpf(Lc[0][1]) * mp.mpf(w)
                    C2, r2, _ = hpslp.converge([tuple(c) for c in Cp], cs, floor_rel=1e-47)
                    row["escape"].append([rnd, srel, sg, float(r2 / mp.mpf(Lc[0][1]) - 1)])
                    if r2 > best[0]: best = (r2, C2)
            if best[1] is None or best[0] <= mp.mpf(Lc[0][1]) * (1 + mp.mpf(10) ** -15): break
            tmp2 = tmp + "2"; how2, extra2 = finish(best[1], tmp2, r_cur)
            if how2 == "none" or F(open(tmp2).readline().split()[1]) <= r_cur: break
            os.replace(tmp2, tmp); how = "escape+" + how2; row.update(extra2); lo = lopt.certify(cs, tmp)
        r_new = F(open(tmp).readline().split()[1])
        row.update({"how": how, "r_old": str(r_old), "r_new": str(r_new), "gain_vs_old": float(r_new / r_old - 1)})
        row["claude"] = certify_circ.check(cs, tmp, rec, verbose=False)[0]
        row["grok"] = grok_verdict(shelf, cs, tmp, n, rec)
        row["lopt"] = lo["verdict"]
        row["lopt_detail"] = {k: v for k, v in lo.items() if k not in ("S", "lambda", "r", "file", "container")}
        with open(out_json + ".cert.tmp", "w") as f: json.dump({"S": lo.get("S"), "lambda": lo.get("lambda")}, f)
        os.replace(out_json + ".cert.tmp", os.path.join(d, f"{shelf}_{n}.lopt.json"))
        os.replace(tmp, out_txt)
    except Exception as e:
        row["error"] = repr(e); row["trace"] = traceback.format_exc(); row["lopt"] = "ERROR"
    row["seconds"] = round(time.time() - t0, 1)
    with open(out_json + ".tmp", "w") as f: json.dump(row, f)
    os.replace(out_json + ".tmp", out_json)
    return shelf, n, row["lopt"], row.get("gain_vs_old"), row["seconds"]

if __name__ == "__main__":
    W = int(next((a.split("=")[1] for a in sys.argv if a.startswith("--workers=")), 12))
    only = next((a.split("=")[1] for a in sys.argv if a.startswith("--only=")), None)
    if only: jobs = [(s.split(":")[0], int(s.split(":")[1])) for s in only.split(",")]
    else:
        shelves = [a for a in sys.argv[1:] if not a.startswith("--")] or SHELVES
        jobs = [(s, int(os.path.basename(p)[:-4].split("_")[-1])) for s in shelves
                for p in glob.glob(os.path.join(HERE, "submission", s, f"{s}_*.txt"))]
    jobs.sort(key=lambda j: -j[1])
    print(f"{len(jobs)} packings, {W} workers", flush=True); t0 = time.time(); k = 0
    with Pool(W) as pool:
        for shelf, n, v, g, sec in pool.imap_unordered(job, jobs):
            k += 1
            print(f"[{k}/{len(jobs)} {time.time() - t0:.0f}s] {shelf} N={n}: {v} gain {g} ({sec}s)", flush=True)
    print(f"DONE {len(jobs)} in {time.time() - t0:.0f}s", flush=True)
