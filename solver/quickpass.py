"""Quick pass over a whole shelf: ONE MBH chain per N (seeded with the published record), 12 N in parallel (no per-N pool
start-up). Catches records that are merely unpolished or one move away. Writes improvements as exact circle-form certificates
(cert_<shelf>/circ_N.txt, checked by certify_circ.py) and appends every result to out/quick_<shelf>.jsonl.
usage: py -3.11 quickpass.py <shelf> [--n1 a] [--n2 b] [--sec-per-n a,b] [--hard-stop HH:MM] [--workers 12]
  shelf: crc_100..crc_800 | ccq | crt (crt uses the circle form of the triangle)
"""
import sys, os, time, json, argparse, datetime, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import campaign_circ as cc

def info(shelf):
    if shelf == "crt":
        return ("tri",), "tri", os.path.join(HERE, "data", "coords", "crt{}.txt"), os.path.join(HERE, "data", "radius.txt")
    return cc.shelf_info(shelf)

def job(args):
    shelf, n, sec, seed, stop_ts = args
    cont, cs, pat, rp = info(shelf)
    c0 = np.array([[float(l.split()[1]), float(l.split()[2])] for l in open(pat.format(n)) if l.strip()])
    if stop_ts: sec = min(sec, max(1.0, stop_ts - time.time()))
    r, c, steps = cc.chain((n, cont, sec, seed, c0, 60))
    return n, r, c, steps

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("shelf"); ap.add_argument("--n1", type=int, default=1); ap.add_argument("--n2", type=int, default=10 ** 6)
    ap.add_argument("--sec-per-n", default="0.05,5"); ap.add_argument("--hard-stop", default=None); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--tag", default="quick"); ap.add_argument("--only", default=None, help="comma list of N (overrides n1/n2)")
    a = ap.parse_args()
    import certify_circ
    cont, cs, pat, rp = info(a.shelf)
    R = {int(l.split()[0]): l.split()[1] for l in open(rp) if l.strip()}
    ns = [n for n in sorted(R) if a.n1 <= n <= a.n2 and n >= 2]
    if a.only: ns = [n for n in (int(v) for v in a.only.split(",")) if n in R and n >= 2]
    stop_ts = None
    if a.hard_stop:
        hh, mm = map(int, a.hard_stop.split(":"))
        stop_ts = datetime.datetime.now().replace(hour=hh, minute=mm, second=0, microsecond=0).timestamp()
    ka, kb = map(float, a.sec_per_n.split(","))
    os.makedirs(os.path.join(HERE, f"cert_{a.shelf}"), exist_ok=True)
    log = open(os.path.join(HERE, "out", f"quick_{a.shelf}.jsonl"), "a"); T0 = time.time(); found = 0
    jobs = [(a.shelf, n, ka * n + kb, 99000 + n, stop_ts) for n in sorted(ns, key=lambda n: -n)]   # big first (balance)
    with Pool(a.workers) as pool:
        for n, r, c, steps in pool.imap_unordered(job, jobs):
            rel = r / float(R[n]) - 1; v = ""
            if rel > 1e-13:                       # candidate only; exact checks happen once, in finalize_circ.py (09-23 speed fix)
                path = os.path.join(HERE, f"cert_{a.shelf}", f"circ_{n}.txt"); npy = path.replace(".txt", ".npy")
                import slp_circ
                prev = slp_circ.rmin(np.load(npy), cont) if os.path.exists(npy) else -1.0
                if r > prev:
                    cc.write_cert(cs, c, path); np.save(npy, c)
                v = "CAND"; found += 1
            log.write(json.dumps({"shelf": a.shelf, "n": n, "r": r, "r_rec": R[n], "rel": rel, "steps": steps, "cert": v, "tag": a.tag,
                                  "time": datetime.datetime.now().isoformat(timespec="seconds")}) + "\n"); log.flush()
            if v: print(f"{a.shelf} N={n:3d} rel {rel:+.3e} ({steps} steps) {v}", flush=True)
    print(f"{a.shelf} QUICKPASS_DONE: {len(ns)} N, {found} candidates, {time.time() - T0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
