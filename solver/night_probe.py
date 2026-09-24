"""Night attack #2 (vision/CLAUDE-SEALED-NIGHT.md): deeper graph-edit moves, two arms, the same wall-clock budget per size.
  R 'rmbh'  : iterated relocation = basin hopping whose kick is a graph edit: remove k in {1,2,3} circles drawn from the 6 weakest
              (least contact force from the v1.1 certificate when it belongs to this exact packing, else fewest near contacts),
              put them into k holes drawn from the 5 largest, polish; accept if r drops by < 1e-6 relative; up to 20 kicks.
  T 'chain' : deep transplant chains: grow N-3 -> N-2 -> N-1 -> N (insert into a hole drawn from the 3 largest) or shrink
              N+3 -> N (delete one of the 3 fewest-contact circles), polishing at every step; up to 4 chains (2 up, 2 down).
Each arm's best strict float gain over OUR packing at N goes through the full chain (cert_candidate.run); an arm scores a hit when
its certified radius beats our radius at the START of the job and both exact checkers say IMPROVES (vs Packomania).
Targets fixed by the sealed rule and written to out/night_targets{tag}.json before any polish.   -> out/night_probe{tag}.jsonl
usage: py -3.11 night_probe.py [--workers=12] [--budget=240] [--tag=] [--targets=file.json] [--arms=RT] [--stop=HH:MM]
"""
import os, sys
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import json, glob, time, io, contextlib, datetime, zlib, numpy as np
from fractions import Fraction as F
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SHELVES = ["crt", "ccq"] + [f"crc_{k}" for k in range(100, 900, 100)]
ARG = lambda k, d: next((x.split("=", 1)[1] for x in sys.argv if x.startswith(f"--{k}=")), d)
BUDGET = float(ARG("budget", 240)); ARMS = ARG("arms", "RT")
STOP = datetime.datetime.fromisoformat(os.environ["NIGHT_STOP"]) if os.environ.get("NIGHT_STOP") else None  # set by the parent

def src_path(shelf, n):
    for p in (os.path.join(HERE, "cand_hp", shelf, f"{shelf}_{n}.txt"), os.path.join(HERE, "lopt_hp", shelf, f"{shelf}_{n}.txt"),
              os.path.join(HERE, "submission", shelf, f"{shelf}_{n}.txt")):
        if os.path.exists(p): return p
    return None

def order(shelf, n, c, r, cont):
    """Circle indices, weakest first. Force only when the v1.1 certificate belongs to THIS packing (same file as submission/)."""
    import transplant_probe as tp
    lc = os.path.join(HERE, "lopt_certs", shelf, f"{shelf}_{n}.json"); sub = os.path.join(HERE, "submission", shelf, f"{shelf}_{n}.txt")
    sp = src_path(shelf, n)
    if os.path.exists(lc) and sp and os.path.exists(sub) and open(sp).readline() == open(sub).readline():
        d = json.load(open(lc)); load = np.zeros(n)
        for s, lamh in zip(d["S"], d["lambda"]):
            lam = float.fromhex(lamh)
            for i in (s[1:3] if s[0] == "p" else s[1:2]): load[i] += lam
        return list(np.argsort(load, kind="stable")), "force"
    return list(np.argsort(tp.near_counts(c, r, cont), kind="stable")), "contacts"

def arm_R(shelf, n, c0, r0, cont, rng, t_end):
    import slp_circ, transplant_probe as tp
    cur_c, cur_r = c0, r0; best = (r0, None); o, how = order(shelf, n, c0, r0, cont); log = []
    for t in range(20):
        if time.time() > t_end: break
        k = int(rng.choice([1, 2, 3], p=[.5, .3, .2])); rm = rng.choice(o[:6], size=k, replace=False)
        base = np.delete(cur_c, rm, 0); H = tp.holes(base, cont, k=5, seed=int(rng.integers(1 << 30)))
        q = np.vstack([base, H[rng.choice(len(H), size=k, replace=False)]])
        c, r = slp_circ.polish(q, cont, t_cap=60.0); log.append(round(r / r0 - 1, 12))
        if r > best[0]: best = (r, c)
        if r >= cur_r * (1 - 1e-6):
            cur_c, cur_r = c, r; o = list(np.argsort(tp.near_counts(c, r, cont), kind="stable"))
    return best, {"weak": how, "kicks": len(log), "trace": log}

def arm_T(shelf, n, r0, cont, rng, t_end):
    import slp_circ, transplant_probe as tp
    best = (r0, None); log = []
    for ch in range(4):
        up = ch < 2; start = tp.packing(shelf, n - 3 if up else n + 3)
        if start is None or time.time() > t_end: continue
        c = start[0]; r = start[1] or slp_circ.rmin(c, cont)
        for step in range(3):
            if up: H = tp.holes(c, cont, k=3, seed=int(rng.integers(1 << 30))); q = np.vstack([c, H[rng.integers(len(H))]])
            else: w = np.argsort(tp.near_counts(c, r, cont), kind="stable")[:3]; q = np.delete(c, rng.choice(w), 0)
            c, r = slp_circ.polish(q, cont, t_cap=60.0)
        log.append(("up" if up else "down", round(r / r0 - 1, 12)))
        if r > best[0]: best = (r, c)
    return best, {"chains": log}

def lines(c, r):
    """Lattice lines: cluster centres by the coordinate ACROSS the line direction; pick the direction whose consecutive line
    spacings look hexagonal (1.55..1.95 r). -> (score, ax = coordinate the lines run along, list of index lists)."""
    best = None
    for ax in (0, 1):
        v = c[:, 1 - ax]; o = np.argsort(v); cl = [[o[0]]]
        for a, b in zip(o[:-1], o[1:]):
            if v[b] - v[a] < 0.5 * r: cl[-1].append(b)
            else: cl.append([b])
        sp = np.diff([v[x].mean() for x in cl]) / r; score = float(np.mean((sp > 1.55) & (sp < 1.95))) if len(sp) else 0.0
        if best is None or score > best[0]: best = (score, ax, cl)
    return best

def arm_P(shelf, n, c0, r0, cont, rng, t_end):
    """Grok's row-phase flip (MERGED-NIGHT): collective moves of whole lattice lines, worst interfaces first."""
    import slp_circ
    if cont[0] != "rect": return (r0, None), {"skip": "not a rectangle"}
    score, ax, L = lines(c0, r0)
    if score < 0.5: return (r0, None), {"skip": f"no lattice lines ({score:.2f})"}
    info = []
    for i in range(len(L) - 1):
        A, B = c0[L[i], ax], c0[L[i + 1], ax]; d = B[:, None] - A[None, :]; s = d[np.arange(len(B)), np.abs(d).argmin(1)]
        info.append((float(np.std(s) / r0 + abs(abs(s.mean()) - r0) / r0), i, float(s.mean())))
    info.sort(reverse=True); best = (r0, None); log = []
    for defect, i, sm in info:
        one, above = np.array(L[i + 1]), np.concatenate(L[i + 1:])
        for idx, delta in ((one, r0), (one, -r0), (above, -2 * sm), (above, r0), (above, -r0)):
            if time.time() > t_end: break
            if abs(delta) < 0.05 * r0: continue
            q = c0.copy(); q[idx, ax] += delta; lim = (0.5 if ax == 0 else cont[1] / 2) - r0
            for line in (L[j] for j in range(len(L)) if set(L[j]) <= set(idx.tolist())):   # wrap: a circle pushed through the
                ln = np.array(line); out = ln[np.abs(q[ln, ax]) > lim]                     # wall re-enters at the free end
                for k in out:
                    rest = ln[np.abs(q[ln, ax]) <= lim]
                    if len(rest) == 0: break                                                # short line: clip below
                    q[k, ax] = (q[rest, ax].min() - 2 * r0) if q[k, ax] > 0 else (q[rest, ax].max() + 2 * r0)
                q[ln, ax] = np.clip(q[ln, ax], -lim, lim)
            c, r = slp_circ.polish(q, cont, t_cap=60.0); log.append((i, round(r / r0 - 1, 9)))
            if r > best[0]: best = (r, c)
        if time.time() > t_end: break
    return best, {"axis": "xy"[ax], "lattice": round(score, 2), "lines": len(L), "tries": len(log), "trace": log[:40]}

def job(a):
    shelf, n = a
    if STOP and datetime.datetime.now() >= STOP: return {"shelf": shelf, "N": n, "skipped": "stop"}
    import slp_circ, finalize_circ, transplant_probe as tp, cert_candidate
    cont = finalize_circ.info(shelf)[0]; cur = tp.packing(shelf, n)
    r_start_exact, _ = cert_candidate.ours(shelf, n)
    c0 = cur[0]; r0 = slp_circ.rmin(c0, cont)          # float radius of our packing as loaded (the arms' baseline)
    rng = np.random.default_rng(zlib.crc32(f"{shelf}_{n}".encode())); row = {"shelf": shelf, "N": n, "r0": r0}
    for arm in ARMS:
        t_end = time.time() + BUDGET; t0 = time.time()
        (rb, cb), info = (arm_R(shelf, n, c0, r0, cont, rng, t_end) if arm == "R" else arm_T(shelf, n, r0, cont, rng, t_end)
                          if arm == "T" else arm_P(shelf, n, c0, r0, cont, rng, t_end))
        res = {"rel_float": float(rb / r0 - 1), "secs": round(time.time() - t0, 1), **info}
        if cb is not None and rb > r0 * (1 + 1e-10):
            npy = os.path.join(HERE, "out", "transplant", f"night_{shelf}_{n}_{arm}.npy"); np.save(npy, cb)
            try:
                with contextlib.redirect_stdout(io.StringIO()): cr = cert_candidate.run(shelf, n, npy, tag="night_" + arm)
                res.update({k: cr.get(k) for k in ("r_new", "gain_vs_ours", "gain_vs_packomania", "claude", "grok", "lopt", "kept")})
                res["beats_start"] = r_start_exact is None or F(cr["r_new"]) / r_start_exact - 1 > F(1, 10 ** 10)
                res["hit"] = bool(res["beats_start"] and cr["claude"] == "IMPROVES" and cr["grok"] == "IMPROVES")
            except Exception as e: res["error"] = repr(e)
        row[arm] = res
    return row

def sealed_targets():
    """36 sizes: 12 gen2 hits, 12 gen2 misses, 12 cascade-wave (sweep3) misses; evenly spaced by N inside each group;
    only sizes with our packings at N-3..N+3; groups disjoint."""
    rows = lambda p: [json.loads(l) for l in open(os.path.join(HERE, "out", p))] if os.path.exists(os.path.join(HERE, "out", p)) else []
    g2 = rows("gen2_probe_B.jsonl") + rows("gen2_probe_Bfull.jsonl") + rows("gen2_probe_Afull.jsonl")
    hit = lambda r: r.get("kept") and r.get("claude") == "IMPROVES" and r.get("grok") == "IMPROVES"
    import transplant_probe as tp
    ok = lambda s, n: all(tp.packing(s, m) is not None and tp.packing(s, m)[2] == "ours" for m in range(n - 3, n + 4))
    groups = [sorted({(r["shelf"], r["N"]) for r in g2 if hit(r)}, key=lambda t: (t[1], t[0])),
              sorted({(r["shelf"], r["N"]) for r in g2 if not hit(r) and "rel_float" in r}, key=lambda t: (t[1], t[0])),
              sorted({(r["shelf"], r["N"]) for r in rows("sweep3.jsonl") if "rel_float" in r and r["rel_float"] <= 0},
                     key=lambda t: (t[1], t[0]))]
    T, used = [], set()
    for G in groups:
        G = [t for t in G if t not in used and ok(*t)]; k = max(1, len(G) // 12); pick = G[k // 2::k][:12]
        T += pick; used |= set(pick)
    return T

def grok_targets():
    """Grok's rule (MERGED-NIGHT): crc_600/700/800, zero loose circles (every circle >= 3 near contacts) in the packing we hold,
    N within 15 of 188 / 219 / 286; 10 per band: every unbeaten one first, then evenly spaced beaten ones (by N, table)."""
    import transplant_probe as tp, finalize_circ, slp_circ
    R = json.load(open(os.path.join(HERE, "out", "final_records.json"))); have = {(r["shelf"], r["n"]) for r in R if r["claimable"]}
    for p in glob.glob(os.path.join(HERE, "cand_hp", "*", "*.json")):
        j = json.load(open(p))
        if j.get("claude") == "IMPROVES" and j.get("grok") == "IMPROVES": have.add((j["shelf"], j["N"]))
    T = []
    for b in (188, 219, 286):
        U, B = [], []                                     # only ~11 unbeaten zero-loose sizes exist in the bands, so each band
        for n in range(b - 15, b + 16):                   # is topped up to 10 with zero-loose sizes we already hold (MERGED-NIGHT)
            for s in ("crc_600", "crc_700", "crc_800"):
                P = tp.packing(s, n)
                if P is None or (s, n) in T: continue
                cont = finalize_circ.info(s)[0]; c = P[0]; r = P[1] or slp_circ.rmin(c, cont)
                if (tp.near_counts(c, r, cont) >= 3).all(): (B if (s, n) in have else U).append((s, n))
        U = U[:10]; need = 10 - len(U); k = max(1, len(B) // max(need, 1)); T += U + (B[k // 2::k][:need] if need else [])
    return T

if __name__ == "__main__":
    W = int(ARG("workers", 12)); tag = ARG("tag", ""); tf = ARG("targets", None)
    if ARG("stop", None):                                     # next occurrence of HH:MM; children read it from the environment
        now = datetime.datetime.now(); h, m = map(int, ARG("stop", None).split(":"))
        STOP = now.replace(hour=h, minute=m, second=0, microsecond=0)
        if STOP <= now: STOP += datetime.timedelta(days=1)
        os.environ["NIGHT_STOP"] = STOP.isoformat()
    T = [tuple(t) for t in json.load(open(tf))] if tf else sealed_targets() + [t for t in grok_targets() if t not in set(sealed_targets())]
    json.dump(T, open(os.path.join(HERE, "out", f"night_targets{tag}.json"), "w"))
    print(f"{len(T)} targets, arms {ARMS}, budget {BUDGET:.0f} s/arm, {W} workers" + (f", stop {STOP:%m-%d %H:%M}" if STOP else ""), flush=True)
    t0 = time.time(); hits = {a: 0 for a in ARMS}
    with Pool(W) as pool, open(os.path.join(HERE, "out", f"night_probe{tag}.jsonl"), "a") as f:
        for o in pool.imap_unordered(job, T):
            f.write(json.dumps(o) + "\n"); f.flush()
            for a in ARMS:
                x = o.get(a, {}); hits[a] += bool(x.get("hit"))
                if x.get("hit") or x.get("error"):
                    print(f"[{time.time() - t0:.0f}s] {o['shelf']} N={o['N']} {a}: " + (f"HIT {x['gain_vs_packomania']:+.2e} vs Packomania, "
                          f"lopt {x['lopt']}, kept {x['kept']}" if x.get("hit") else f"ERROR {x['error']}"), flush=True)
    print(f"DONE hits {hits} / {len(T)} -> " + ", ".join(f"{a} {'PASS' if hits[a] >= 3 else 'NULL'}" for a in ARMS), flush=True)
