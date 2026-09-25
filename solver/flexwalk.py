# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""FLEX WALK (attack #3 of plan A = attack #1 of plan B, converged 09-23): a NOT_RIGID certificate has a load-bearing first-order flex
v (G_S v = 0, dr = 0). Along v the stressed contacts open at second order, so r can grow once the structure re-relaxes.
For each kernel direction and sign: s_max = the distance along v at which the first OTHER constraint (any pair or wall, linearised)
becomes tight (capped at 5% of r); screen s = s_max x (0.1, 0.3, 0.6, 1.0) with the float polish; the best strict gain goes through
the full chain (cert_candidate.run: 80-digit converge -> Newton -> 45-digit certificate -> BOTH exact checkers -> lopt).
Repeat while NOT_RIGID (max 3 rounds).
usage: py -3.11 flexwalk.py <shelf> <N> [certificate.txt]
"""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import lopt, slp_circ, finalize_circ, cert_candidate

def kernel_dirs(o, n):
    S, lam, bcols, G, ker = o["_so"]
    dirs = []
    for k in range(ker.shape[1]):
        v = np.zeros(2 * n + 1)
        for q, col in enumerate(bcols): v[int(col)] = ker[q, k]
        v /= np.abs(v[:-1]).max(); dirs.append(v)
    return dirs

def s_max(c, r, v, cont):
    """Largest s such that no pair / wall constraint (linearised in s) drops below the current radius r, capped 0.05 r."""
    n = len(c); d = v[:-1].reshape(n, 2); iu = np.triu_indices(n, 1)
    diff = c[:, None] - c[None]; D = np.sqrt((diff ** 2).sum(-1)); dd = d[:, None] - d[None]
    rate = (diff * dd).sum(-1) / np.maximum(D, 1e-300)          # d/ds |ci - cj|
    sl = D / 2 - r; rt = rate / 2
    g, Gw = slp_circ.walls(c, cont); wsl = g - r; wrt = (Gw * d[:, None, :]).sum(-1)
    cand = [0.05 * r]
    for s_, r_ in ((sl[iu], rt[iu]), (wsl.ravel(), wrt.ravel())):
        m = (r_ < -1e-12) & (s_ > 1e-12 * r)
        if m.any(): cand.append(float((s_[m] / -r_[m]).min()))
    return min(cand)

def walk(shelf, n, path=None, rounds=3, log=print, max_dirs=4, fracs=(0.3, 1.0)):
    cont, cs, pat, rp = finalize_circ.info(shelf)
    if path is None: path = cert_candidate.ours(shelf, n)[1]
    hist = []
    for rnd in range(rounds):
        o = lopt.certify(cs, path, want_flex="second")
        if o["verdict"] != "NOT_RIGID" or "_so" not in o:
            log(f"{shelf} N={n} round {rnd}: {o['verdict']} -> stop"); break
        r0, C = lopt.load(path); c0 = np.array([[float(x), float(y)] for x, y in C]); rf = float(r0)
        best = (rf, None, None)
        for k, v in enumerate(kernel_dirs(o, n)[:max_dirs]):
            for sg in (1, -1):
                sm = s_max(c0, rf, sg * v, cont)
                for frac in fracs:
                    q = c0 + frac * sm * sg * v[:-1].reshape(n, 2)
                    cq, rq = slp_circ.polish(q, cont, max_iter=300, t_cap=60.0)
                    if rq > best[0]: best = (rq, cq, (k, sg, frac, sm))
        gain = best[0] / rf - 1
        log(f"{shelf} N={n} round {rnd}: dims {len(kernel_dirs(o, n))}, best float gain {gain:+.3e} {best[2]}")
        hist.append({"round": rnd, "float_gain": gain, "move": best[2]})
        if best[1] is None or gain <= 1e-12: break
        os.makedirs(os.path.join(HERE, "out", "flexwalk"), exist_ok=True)
        npy = os.path.join(HERE, "out", "flexwalk", f"{shelf}_{n}_r{rnd}.npy"); np.save(npy, best[1])
        row = cert_candidate.run(shelf, n, npy, tag=f"flexwalk{rnd}")
        hist[-1].update({k: row[k] for k in ("gain_vs_ours", "gain_vs_packomania", "checker_a", "checker_b", "lopt", "kept")})
        if not row["kept"]: break
        path = os.path.join(HERE, "cand_hp", shelf, f"{shelf}_{n}.txt")
    return hist

if __name__ == "__main__":
    h = walk(sys.argv[1], int(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else None)
    print(json.dumps(h, indent=1, default=str))
