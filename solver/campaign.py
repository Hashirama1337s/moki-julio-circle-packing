# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Record campaign for Packomania crt (points in the unit isosceles right triangle, maximise min distance).
For each N: 12 CPU workers run monotonic basin hopping (SLP local solver) SEEDED with the published record (and, for small N,
random starts too). Perturbations: global U(+-s d) with small s, or regional (points within rho of a random spot, U(+-s d)).
On stall: big regional kick from the worker's incumbent. Every result is appended to out/campaign.jsonl; any d > record is
written as an exact-certificate candidate (cert/cand_N.txt) and checked with certify.py immediately.
usage: py -3.11 campaign.py N1 N2 [--step k] [--sec S | --sec-per-n a,b] [--hard-stop HH:MM] [--random-frac f]
"""
import sys, os, time, json, argparse, datetime, numpy as np
from multiprocessing import Pool
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)

def clipT(q):
    q = np.clip(q, 0, 1); o = q.sum(1) - 1; q[o > 0] -= o[o > 0, None] / 2; return q

def rand_T(rng, n):
    u = rng.random((n, 2)); f = u.sum(1) > 1; u[f] = 1 - u[f]; return u

def perturb(rng, p, d):
    q = p.copy(); n = len(p)
    if rng.random() < 0.5:
        s = rng.choice([0.01, 0.02, 0.05, 0.1, 0.2]) * d
        q += rng.uniform(-s, s, q.shape)
    else:
        c = p[rng.integers(n)]; rho = rng.choice([1.5, 2.5, 4.0]) * d; s = rng.choice([0.2, 0.4, 0.7]) * d
        m = ((p - c) ** 2).sum(1) < rho * rho; q[m] += rng.uniform(-s, s, (m.sum(), 2))
    return clipT(q)

def chain(args):
    n, seconds, seed, seed_pts, random_start, stall = args
    import slp as search
    rng = np.random.default_rng(seed); t0 = time.time()
    start = rand_T(rng, n) if random_start else seed_pts
    p, d = search.polish(start); best_p, best_d = p, d; bad = steps = kicks = 0
    while time.time() - t0 < seconds:
        q, dq = search.polish(perturb(rng, p, d)); steps += 1
        if dq > d + 1e-15: p, d, bad = q, dq, 0
        else: bad += 1
        if d > best_d: best_p, best_d = p, d
        if bad >= stall:                       # kick from the incumbent (non-monotonic), keep global best
            c = best_p[rng.integers(n)]; rho = 5 * best_d; m = ((best_p - c) ** 2).sum(1) < rho * rho
            k = best_p.copy(); k[m] += rng.uniform(-0.7 * best_d, 0.7 * best_d, (m.sum(), 2))
            p, d = search.polish(clipT(k)); bad = 0; kicks += 1
    return best_d, best_p, steps, kicks

def write_candidate(n, p, path):
    """Exact-rational-safe writer: decimal repr of each float, then exact repair so every point is inside T as a rational."""
    with open(path, "w") as f:
        f.write(f"N {n}\n")
        for x, y in p:
            X, Y = F(repr(float(x))), F(repr(float(y)))
            X = max(X, F(0)); Y = max(Y, F(0))
            if X + Y > 1: Y = 1 - X
            f.write(f"{float(X)!r} {float(Y)!r}\n" if (F(repr(float(X))) == X and F(repr(float(Y))) == Y) else f"{X.numerator}/{X.denominator} {Y.numerator}/{Y.denominator}\n")

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("n1", type=int); ap.add_argument("n2", type=int)
    ap.add_argument("--step", type=int, default=1); ap.add_argument("--sec", type=float, default=None)
    ap.add_argument("--sec-per-n", default="0.1,0"); ap.add_argument("--hard-stop", default=None)
    ap.add_argument("--random-frac", type=float, default=0.0); ap.add_argument("--stall", type=int, default=60)
    ap.add_argument("--workers", type=int, default=12); ap.add_argument("--tag", default="run")
    a = ap.parse_args()
    from records import table, coords, circle_radius, to_points
    import certify
    Dt = table("distance.txt"); os.makedirs(os.path.join(HERE, "out"), exist_ok=True); os.makedirs(os.path.join(HERE, "cert"), exist_ok=True)
    stop = None
    if a.hard_stop:
        hh, mm = map(int, a.hard_stop.split(":")); now = datetime.datetime.now()
        stop = now.replace(hour=hh, minute=mm, second=0, microsecond=0)
    ka, kb = map(float, a.sec_per_n.split(","))
    log = open(os.path.join(HERE, "out", "campaign.jsonl"), "a")
    for n in range(a.n1, a.n2 + 1, a.step):
        if stop and datetime.datetime.now() >= stop: print("hard stop reached", flush=True); break
        sec = a.sec if a.sec else ka * n + kb
        if stop: sec = min(sec, max(5.0, (stop - datetime.datetime.now()).total_seconds() - 30))
        dr = float(Dt[n]); pts = coords(n); r = circle_radius(pts)
        R = np.array([[float(x), float(y)] for x, y in to_points(pts, r)])
        nrand = int(round(a.random_frac * a.workers))
        jobs = [(n, sec, 1000 * n + w, R, w < nrand, a.stall) for w in range(a.workers)]
        t0 = time.time()
        with Pool(a.workers) as pool: res = pool.map(chain, jobs)
        res.sort(key=lambda z: -z[0]); d, p = res[0][0], res[0][1]
        rel = d / dr - 1; steps = sum(z[2] for z in res)
        rec = {"n": n, "d": d, "d_rec": dr, "rel": rel, "sec": round(time.time() - t0, 1), "steps": steps, "tag": a.tag,
               "time": datetime.datetime.now().isoformat(timespec="seconds"), "workers_d": [z[0] for z in res]}
        verdict = ""
        if rel > 1e-13:
            path = os.path.join(HERE, "cert", f"cand_{n}.txt")
            prev = None
            if os.path.exists(path):
                prev = certify.check(path, verbose=False)
            if prev is None or prev[2] < rel:
                write_candidate(n, p, path); np.save(os.path.join(HERE, "cert", f"cand_{n}.npy"), p)
            v = certify.check(path, verbose=False); verdict = v[0]; rec["cert"] = v[0]; rec["cert_rel"] = v[2]
        log.write(json.dumps(rec) + "\n"); log.flush()
        print(f"N={n:3d} rel {rel:+.3e} ({steps} steps, {rec['sec']}s) {verdict}", flush=True)

if __name__ == "__main__":
    main()
