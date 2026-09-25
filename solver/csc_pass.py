"""First pass on a NEW table: Packomania csc (equal circles in the unit semicircle, y >= 0), N = 1..250.
Per size, one worker: polish Packomania's own packing; Packomania's N-1 packing + one circle in each of its 3 largest holes;
Packomania's N+1 packing minus each of its 3 fewest-contact circles; then monotonic basin hopping (campaign_circ.perturb kicks,
big kick after a stall) from the best of those for the rest of the time budget. A float gain > 1e-10 over the published radius
goes through the full chain (cert_candidate.run: 80-digit SLP -> Newton -> 45-digit certificate, radius rounded down -> checker A
+ checker B (independent verify_exact3.py) -> local optimality) into cand_hp/csc/.
Not attempted: the 10 sizes whose coordinates Packomania has not published (Hogan, 09-Sep-2026) — the table-radius gate cannot
run there. Claims for N = 151..200 are additionally gated by Lai et al. 2025 (Table 8) at build time.
usage: py -3.11 csc_pass.py [--workers=6] [--sec=0.4,20] [--stop=HH:MM] [--from=1] [--to=250] [--out=csc_pass.jsonl]
"""
import os, sys
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import json, time, io, contextlib, datetime, zlib, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ARG = lambda k, d: next((x.split("=", 1)[1] for x in sys.argv if x.startswith(f"--{k}=")), d)
SEC = tuple(map(float, ARG("sec", "0.4,20").split(",")))
STOP = datetime.datetime.fromisoformat(os.environ["CSC_STOP"]) if os.environ.get("CSC_STOP") else None
SH = os.path.join(HERE, "shelves", "csc")

def pk(n):
    p = os.path.join(SH, "coords", f"csc{n}.txt")
    if not os.path.exists(p): return None
    L = [l.split() for l in open(p) if l.strip()]
    return np.array([[float(a), float(b)] for _, a, b in L]) if len(L) == n else None

def job(n):
    if STOP and datetime.datetime.now() >= STOP: return {"shelf": "csc", "N": n, "skipped": "stop"}
    import slp_circ, campaign_circ, transplant_probe as tp, cert_candidate
    cont = ("semi",); rec = {int(l.split()[0]): l.split()[1] for l in open(os.path.join(SH, "radius.txt")) if l.strip()}[n]
    P = pk(n)
    if P is None: return {"shelf": "csc", "N": n, "skipped": "no published coordinates"}
    r_rec = float(rec); t0 = time.time(); budget = SEC[0] * n + SEC[1]; tried = []
    rng = np.random.default_rng(zlib.crc32(f"csc_{n}".encode()))
    best = (0.0, None, None)
    def trial(q, tag):
        nonlocal best
        c, r = slp_circ.polish(q, cont, t_cap=60.0); tried.append((tag, float(r / r_rec - 1)))
        if r > best[0]: best = (r, c, tag)
    trial(P, "pk")
    dn, up = (pk(n - 1) if n > 1 else None), pk(n + 1)
    if dn is not None:
        for k, h in enumerate(tp.holes(dn, cont, k=3)): trial(np.vstack([dn, h]), f"pkins{k}")
    if up is not None:
        ru = slp_circ.rmin(up, cont)
        for k, i in enumerate(np.argsort(tp.near_counts(up, ru, cont), kind="stable")[:3]): trial(np.delete(up, i, 0), f"pkdel{k}")
    c, r = best[1], best[0]; steps = 0; bad = 0; seed_tag = best[2]
    while time.time() - t0 < budget:                                   # MBH from the best seed
        q, rq = slp_circ.polish(campaign_circ.perturb(rng, c, r, cont), cont, t_cap=60.0); steps += 1
        if rq > r + 1e-15: c, r, bad = q, rq, 0
        else: bad += 1
        if r > best[0]: best = (r, c, seed_tag + "+mbh")
        if bad >= 40:
            p = best[1][rng.integers(n)]; m = ((best[1] - p) ** 2).sum(1) < (10 * best[0]) ** 2
            k = best[1].copy(); k[m] += rng.uniform(-1.4 * best[0], 1.4 * best[0], (m.sum(), 2))
            c, r = slp_circ.polish(slp_circ.repair(k, cont), cont, t_cap=60.0); bad = 0
    row = {"shelf": "csc", "N": n, "rel_float": float(best[0] / r_rec - 1), "seed": best[2], "tried": tried, "mbh_steps": steps,
           "secs": round(time.time() - t0, 1)}
    if best[0] > r_rec * (1 + 1e-10):
        npy = os.path.join(HERE, "out", "transplant", f"csc_{n}_{best[2].replace('+', '_')}.npy"); np.save(npy, best[1])
        try:
            with contextlib.redirect_stdout(io.StringIO()): cr = cert_candidate.run("csc", n, npy, tag="csc_" + best[2].replace("+", "_"))
            row.update({k: cr.get(k) for k in ("r_new", "gain_vs_packomania", "checker_a", "checker_b", "lopt", "kept")})
        except Exception as e: row["error"] = repr(e)
    return row

if __name__ == "__main__":
    W = int(ARG("workers", 6)); out = ARG("out", "csc_pass.jsonl"); n1, n2 = int(ARG("from", 1)), int(ARG("to", 250))
    if ARG("stop", None):
        now = datetime.datetime.now(); h, m = map(int, ARG("stop", None).split(":"))
        s = now.replace(hour=h, minute=m, second=0, microsecond=0); s += datetime.timedelta(days=1) if s <= now else datetime.timedelta(0)
        os.environ["CSC_STOP"] = s.isoformat(); print(f"stop {s:%m-%d %H:%M}", flush=True)
    T = sorted(range(n1, n2 + 1), key=lambda n: -n)                     # big N first (longest jobs start early)
    print(f"csc first pass: {len(T)} sizes, {W} workers, budget {SEC[0]}*N+{SEC[1]} s", flush=True); t0 = time.time(); hits = 0
    with Pool(W) as pool, open(os.path.join(HERE, "out", out), "a") as f:
        for o in pool.imap_unordered(job, T, chunksize=1):
            f.write(json.dumps(o) + "\n"); f.flush()
            ok = o.get("kept") and o.get("checker_a") == "IMPROVES" and o.get("checker_b") == "IMPROVES"; hits += bool(ok)
            if ok or o.get("error") or (o.get("rel_float", 0) > 1e-10):
                print(f"[{time.time() - t0:.0f}s] csc N={o['N']} {o.get('seed')}: rel_float {o.get('rel_float', 0):+.2e} "
                      + (f"RECORD {o['gain_vs_packomania']:+.2e} (lopt {o['lopt']})" if ok else f"checker_a {o.get('checker_a')} checker_b {o.get('checker_b')} {o.get('error', '')}"), flush=True)
    print(f"DONE {hits} certified records / {len(T)} sizes, {time.time() - t0:.0f}s", flush=True)
