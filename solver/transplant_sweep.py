# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""FULL transplant sweep (sealed plan: vision/SEALED-TRANSPLANT2.md PASSED 10/64). For every size N of each shelf's
Packomania table: seed from our current best at N+1 (delete one of the 4 circles with fewest near contacts) and N-1 (insert into
the 3 largest holes), float polish; if a seed beats our current best at N by > 1e-10, run the FULL chain in the same worker
(cert_candidate.run: 80-digit converge, Newton, 45-digit certificate, BOTH exact checkers, lopt) -> cand_hp/ if strictly better.
Cascade: sizes next to a kept hit are re-queued for the next pass (max 3 passes). Hard stop: --stop=HH:MM (no new jobs after).
usage: py -3.11 transplant_sweep.py --workers=3 --stop=21:30 shelf [shelf ...]   -> out/transplant_sweep.jsonl
"""
import os, sys
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import json, time, datetime, io, contextlib, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "out", next((a.split("=")[1] for a in sys.argv if a.startswith("--out=")), "transplant_sweep.jsonl"))  # one file per process

def job(a):
    shelf, n, stop_ts = a
    if time.time() > stop_ts: return {"shelf": shelf, "N": n, "skipped": "hard stop"}
    import slp_circ, finalize_circ, transplant_probe as tp, cert_candidate
    cont = finalize_circ.info(shelf)[0]
    cur = tp.packing(shelf, n); up = tp.packing(shelf, n + 1); dn = tp.packing(shelf, n - 1)
    if cur is None or up is None or dn is None: return {"shelf": shelf, "N": n, "skipped": "no neighbours"}
    r_best = cur[1] or slp_circ.rmin(cur[0], cont)
    best = (r_best, None, None)
    # ping-pong guard (stage-2 review): if our N+1 packing was itself made by INSERTING a circle (appended last), deleting that circle
    # just restores the old N — skip it.
    jup = os.path.join(HERE, "cand_hp", shelf, f"{shelf}_{n + 1}.json")
    no_del = {n} if os.path.exists(jup) and "_ins" in json.load(open(jup)).get("tag", "") else set()
    c_up = up[0]; nc = tp.near_counts(c_up, up[1] or slp_circ.rmin(c_up, cont), cont)
    for i in [int(i) for i in np.argsort(nc, kind="stable") if int(i) not in no_del][:4]:
        c, r = slp_circ.polish(np.delete(c_up, i, 0), cont, t_cap=60.0)
        if r > best[0]: best = (r, c, f"del{i}")
    for j, h in enumerate(tp.holes(dn[0], cont)):
        c, r = slp_circ.polish(np.vstack([dn[0], h]), cont, t_cap=60.0)
        if r > best[0]: best = (r, c, f"ins{j}")
    # SECOND SEED (stage-2 review): Packomania's own N-1 / N+1 packings are a different basin family from ours — one insertion and
    # one deletion from them as well, when ours differ from theirs.
    pat = finalize_circ.info(shelf)[2]
    for m, kind in ((n - 1, "ins"), (n + 1, "del")):
        mine = up if m == n + 1 else dn
        if mine[2] != "ours" or not os.path.exists(pat.format(m)): continue
        cp = np.array([[float(l.split()[1]), float(l.split()[2])] for l in open(pat.format(m)) if l.strip()])
        if kind == "ins": q = np.vstack([cp, tp.holes(cp, cont, k=1)[0]])
        else: q = np.delete(cp, int(np.argsort(tp.near_counts(cp, slp_circ.rmin(cp, cont), cont), kind="stable")[0]), 0)
        c, r = slp_circ.polish(q, cont, t_cap=60.0)
        if r > best[0]: best = (r, c, f"pk{kind}")
    row = {"shelf": shelf, "N": n, "rel_float": best[0] / r_best - 1, "seed": best[2], "base": cur[2]}
    if best[1] is None or best[0] <= r_best * (1 + 1e-10): return row
    npy = os.path.join(HERE, "out", "transplant", f"sweep_{shelf}_{n}_{best[2]}.npy"); np.save(npy, best[1])
    try:
        with contextlib.redirect_stdout(io.StringIO()): cr = cert_candidate.run(shelf, n, npy, tag="sweep_" + best[2])
        row.update({k: cr.get(k) for k in ("gain_vs_ours", "gain_vs_packomania", "checker_a", "checker_b", "lopt", "kept")})
    except Exception as e: row["error"] = repr(e)
    return row

if __name__ == "__main__":
    W = int(next((a.split("=")[1] for a in sys.argv if a.startswith("--workers=")), 3))
    hh, mm = map(int, next(a.split("=")[1] for a in sys.argv if a.startswith("--stop=")).split(":"))
    _now = datetime.datetime.now(); _s = _now.replace(hour=hh, minute=mm, second=0, microsecond=0)
    if _s <= _now: _s += datetime.timedelta(days=1)          # 21:36 bug: "--stop=06:00" meant this morning -> ended at once
    stop = _s.timestamp(); print(f"hard stop {_s.isoformat(timespec='minutes')}", flush=True)
    shelves = [a for a in sys.argv[1:] if not a.startswith("--")]
    import finalize_circ
    queue = []
    for s in shelves:
        rp = finalize_circ.info(s)[3]; ns = sorted(int(l.split()[0]) for l in open(rp) if l.strip())
        queue += [(s, n) for n in ns if n >= 3]
    # --skip-done=a.jsonl,b.jsonl : drop sizes already swept (pass 0) in earlier runs
    # --seed-cascade=a.jsonl,...  : put the neighbours (N-1, N+1) of every KEPT hit in those runs back at the FRONT of the queue
    def _rows(spec):
        out = []
        for f in (spec.split(",") if spec else []):
            p = os.path.join(HERE, "out", f)
            if os.path.exists(p): out += [json.loads(l) for l in open(p) if l.strip()]
        return out
    skip = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--skip-done=")), "")
    seed = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--seed-cascade=")), "")
    done = {(r["shelf"], r["N"]) for r in _rows(skip) if "rel_float" in r or r.get("skipped") == "no neighbours"}
    tables = {s: {int(l.split()[0]) for l in open(finalize_circ.info(s)[3]) if l.strip()} for s in shelves}
    front = []
    for r in _rows(seed):
        if r.get("kept") and r.get("checker_a") == "IMPROVES" and r.get("checker_b") == "IMPROVES" and r["shelf"] in tables:
            front += [(r["shelf"], m) for m in (r["N"] - 1, r["N"] + 1) if m >= 3 and m in tables[r["shelf"]]]
    front = list(dict.fromkeys(front))
    queue = front + [q for q in queue if q not in done and q not in set(front)]
    print(f"queue: {len(front)} cascade seeds first, then {len(queue) - len(front)} unswept sizes (skipped {len(done)} done)", flush=True)
    t0 = time.time()
    for p in range(3):
        if not queue or time.time() > stop: break
        print(f"pass {p}: {len(queue)} sizes", flush=True); nxt = set()
        with Pool(W) as pool, open(OUT, "a") as f:
            for row in pool.imap_unordered(job, [(s, n, stop) for s, n in queue], chunksize=1):
                row["pass"] = p; f.write(json.dumps(row) + "\n"); f.flush()
                if row.get("kept") and row.get("checker_a") == "IMPROVES" and row.get("checker_b") == "IMPROVES":
                    go = row.get("gain_vs_ours")                # None = a size where we had NO record before (new record)
                    print(f"[{time.time() - t0:.0f}s] KEPT {row['shelf']} N={row['N']} {row['seed']} "
                          + (f"+{go:.3e} vs ours, " if go is not None else "NEW SIZE (no previous record of ours), ")
                          + f"{row['gain_vs_packomania']:+.3e} vs Packomania, lopt {row['lopt']}", flush=True)
                    nxt |= {(row["shelf"], row["N"] - 1), (row["shelf"], row["N"] + 1)}
                elif row.get("error"): print(f"[{time.time() - t0:.0f}s] ERROR {row['shelf']} N={row['N']}: {row['error']}", flush=True)
        queue = sorted(nxt)
    print(f"DONE {time.time() - t0:.0f}s", flush=True)
