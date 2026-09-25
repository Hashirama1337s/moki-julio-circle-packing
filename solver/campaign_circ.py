# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Record campaign, circle form, any container (crc rect:h, ccq quad, crt tri). Same MBH scheme as campaign.py, seeded with the
published record; perturbation scales in units of the circle diameter 2r. Candidates d > record are written as exact circle-form
certificates (r_c = largest radius the rounded centres admit, rounded DOWN at 45 digits, computed at 80 digits) and checked with
certify_circ.py immediately.
usage: py -3.11 campaign_circ.py <shelf> N1 N2 [--sec-per-n a,b] [--hard-stop HH:MM] [--workers 12] [--step 1]
  shelf: crc_100 .. crc_800 | ccq
"""
import sys, os, time, json, argparse, datetime, numpy as np, mpmath as mp
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SH = os.path.join(HERE, "shelves")

def shelf_info(shelf):
    """(container tuple for the solver, container string for the exact checker, coords filename pattern, radius.txt path)."""
    if shelf == "ccq":
        return ("quad",), "quad", os.path.join(SH, "ccq", "coords", "ccq{}.txt"), os.path.join(SH, "ccq", "radius.txt")
    if shelf == "csc":
        return ("semi",), "semi", os.path.join(SH, "csc", "coords", "csc{}.txt"), os.path.join(SH, "csc", "radius.txt")
    if shelf == "csq":                               # 09-24 (four hats): Packomania's unit square, centred = rect with h = 1
        return ("rect", 1.0), "rect:1", os.path.join(SH, "csq", "coords", "csq{}.txt"), os.path.join(SH, "csq", "radius.txt")
    k = int(shelf.split("_")[1])                     # crc_300 -> k = 300 -> height 0.3
    h_str = "0." + f"{k:03d}".rstrip("0")            # "0.3" (exact decimal for the checker)
    pat = os.path.join(SH, shelf, "coords", "crc{}_0." + f"{k:03d}" + "000000000.txt")
    return ("rect", float(h_str)), "rect:" + h_str, pat, os.path.join(SH, shelf, "radius.txt")

def perturb(rng, c, r, cont):
    import slp_circ
    q = c.copy(); n = len(c); d = 2 * r
    if rng.random() < 0.5:
        s = rng.choice([0.01, 0.02, 0.05, 0.1, 0.2]) * d; q += rng.uniform(-s, s, q.shape)
    else:
        p = c[rng.integers(n)]; rho = rng.choice([1.5, 2.5, 4.0]) * d; s = rng.choice([0.2, 0.4, 0.7]) * d
        m = ((c - p) ** 2).sum(1) < rho * rho; q[m] += rng.uniform(-s, s, (m.sum(), 2))
    return slp_circ.repair(q, cont)

def chain(args):
    n, cont, seconds, seed, seed_c, stall = args
    import slp_circ
    rng = np.random.default_rng(seed); t0 = time.time()
    c, r = slp_circ.polish(seed_c, cont); best_c, best_r = c, r; bad = 0; steps = 0
    while time.time() - t0 < seconds:
        q, rq = slp_circ.polish(perturb(rng, c, r, cont), cont); steps += 1
        if rq > r + 1e-15: c, r, bad = q, rq, 0
        else: bad += 1
        if r > best_r: best_c, best_r = c, r
        if bad >= stall:
            p = best_c[rng.integers(n)]; rho = 10 * best_r; m = ((best_c - p) ** 2).sum(1) < rho * rho
            k = best_c.copy(); k[m] += rng.uniform(-1.4 * best_r, 1.4 * best_r, (m.sum(), 2))
            c, r = slp_circ.polish(slp_circ.repair(k, cont), cont); bad = 0
    return best_r, best_c, steps

def write_cert(cont_s, c, path):
    """Centres as the floats' decimal reprs; r_c = largest admissible radius of those exact decimals, rounded DOWN at 45 digits.
    Only NEAR pairs are evaluated at 80 digits: a pair whose float distance exceeds 2 r_float (1 + 1e-6) cannot be the minimum
    (float distances are accurate to ~1e-15 relative), so the minimum over near pairs IS the minimum over all pairs.
    finalize_circ.py re-checks every pair exactly anyway."""
    mp.mp.dps = 80
    c = np.asarray(c, dtype=np.float64); n = len(c); iu = np.triu_indices(n, 1)
    df = np.sqrt(((c[:, None] - c[None]) ** 2).sum(-1))[iu]
    C = [(mp.mpf(repr(float(x))), mp.mpf(repr(float(y)))) for x, y in c]
    if cont_s.startswith("rect:"):
        h = mp.mpf(cont_s.split(":")[1]); w = min(min(x + mp.mpf(1) / 2, mp.mpf(1) / 2 - x, y + h / 2, h / 2 - y) for x, y in C)
    elif cont_s == "quad":
        w = min(min(x, y, 1 - mp.sqrt(x * x + y * y)) for x, y in C)
    elif cont_s == "semi":
        w = min(min(y, 1 - mp.sqrt(x * x + y * y)) for x, y in C)
    else:
        w = min(min(x, y, (1 - x - y) / mp.sqrt(2)) for x, y in C)
    near = np.nonzero(df <= df.min() * (1 + 1e-6))[0]
    pr = min(mp.sqrt((C[i][0] - C[j][0]) ** 2 + (C[i][1] - C[j][1]) ** 2) for i, j in zip(iu[0][near], iu[1][near])) / 2
    rc = mp.floor(min(w, pr) * mp.mpf(10) ** 45) / mp.mpf(10) ** 45
    with open(path, "w") as f:
        f.write(f"r {mp.nstr(rc, 45, strip_zeros=False, min_fixed=-mp.inf, max_fixed=mp.inf)}\n")
        for x, y in c: f.write(f"{repr(float(x))} {repr(float(y))}\n")

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("shelf"); ap.add_argument("n1", type=int); ap.add_argument("n2", type=int)
    ap.add_argument("--sec-per-n", default="0.08,4"); ap.add_argument("--hard-stop", default=None)
    ap.add_argument("--workers", type=int, default=12); ap.add_argument("--stall", type=int, default=60)
    ap.add_argument("--step", type=int, default=1); ap.add_argument("--tag", default="run")
    a = ap.parse_args()
    import certify_circ
    cont, cont_s, pat, rpath = shelf_info(a.shelf)
    R = {int(l.split()[0]): l.split()[1] for l in open(rpath) if l.strip()}
    os.makedirs(os.path.join(HERE, "cert_" + a.shelf), exist_ok=True)
    stop = None
    if a.hard_stop:
        hh, mm = map(int, a.hard_stop.split(":")); stop = datetime.datetime.now().replace(hour=hh, minute=mm, second=0, microsecond=0)
    ka, kb = map(float, a.sec_per_n.split(","))
    log = open(os.path.join(HERE, "out", f"campaign_{a.shelf}.jsonl"), "a")
    for n in range(a.n1, a.n2 + 1, a.step):
        if n not in R: continue
        if stop and datetime.datetime.now() >= stop: print("hard stop reached", flush=True); break
        sec = ka * n + kb
        if stop: sec = min(sec, max(5.0, (stop - datetime.datetime.now()).total_seconds() - 30))
        c0 = np.array([[float(l.split()[1]), float(l.split()[2])] for l in open(pat.format(n)) if l.strip()])
        assert len(c0) == n, (n, len(c0))
        rr = float(R[n]); t0 = time.time()
        with Pool(a.workers) as pool:
            res = pool.map(chain, [(n, cont, sec, 7000 * n + w, c0, a.stall) for w in range(a.workers)])
        res.sort(key=lambda z: -z[0]); r, c = res[0][0], res[0][1]; rel = r / rr - 1
        row = {"shelf": a.shelf, "n": n, "r": r, "r_rec": R[n], "rel": rel, "sec": round(time.time() - t0, 1),
               "steps": sum(z[2] for z in res), "tag": a.tag, "time": datetime.datetime.now().isoformat(timespec="seconds")}
        v = ""
        if rel > 1e-13:
            path = os.path.join(HERE, "cert_" + a.shelf, f"circ_{n}.txt")
            prev = certify_circ.check(cont_s, path, R[n], verbose=False) if os.path.exists(path) else None
            if prev is None or prev[2] < rel:
                write_cert(cont_s, c, path); np.save(path.replace(".txt", ".npy"), c)
            vv = certify_circ.check(cont_s, path, R[n], verbose=False); v = vv[0]; row["cert"] = vv[0]; row["cert_rel"] = vv[2]
        log.write(json.dumps(row) + "\n"); log.flush()
        print(f"{a.shelf} N={n:3d} rel {rel:+.3e} ({row['steps']} steps, {row['sec']}s) {v}", flush=True)

if __name__ == "__main__":
    main()
