"""v1.7: GPU basin search for BIG N (square csq up to 10,000 circles). Same idea as gpu_mbh.py (B perturbed copies relaxed together at
a target radius that shrinks, then regrows past our record), but pair forces come from a Verlet neighbour list (KD-tree on the CPU,
rebuilt every 40 steps; cut-off 2.6 r) instead of an N x N matrix, so memory is O(B N). The best copy is polished on the CPU
(slp_big), gated (lai_compare.csq_gate: Packomania envelope + Amore) and certified by BOTH O(N) exact checkers; kept in
cand_big_hp/csq only if it beats the file there.
usage: py -3.11 gpu_big.py --targets=out/x.json [--batch=32] [--rounds=2] [--tag=]   -> out/gpu_big<tag>.jsonl
"""
import os, sys, json, time, glob, subprocess
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import numpy as np
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ARG = lambda k, d: next((a.split("=", 1)[1] for a in sys.argv if a.startswith(f"--{k}=")), d)

def load_best(n):
    for p in (os.path.join(HERE, "cand_big_hp", "csq", f"csq_{n}.txt"), os.path.join(HERE, "cand_hp", "csq", f"csq_{n}.txt"),
              os.path.join(HERE, "cand_big", "csq", f"csq_{n}.txt")):
        if os.path.exists(p):
            L = [l.split() for l in open(p) if l.strip()]
            return np.array([[float(a), float(b)] for a, b in L[1:]]), float(L[0][1]), p
    # 09-25: open cells ON Packomania's envelope have no file of ours -> start from Packomania's own packing (its printed radius).
    p = os.path.join(HERE, "data", "big", "csq", f"csq{n}.txt")
    if os.path.exists(p):
        import hats_flagship as HF
        c = np.array([[float(v) for v in l.split()[1:3]] for l in open(p) if l.strip() and not l.lstrip().startswith("#")])
        if len(c) == n: return c, float(HF.table("csq")[n]), p
    return None

def certify(n, c, tag):
    import mpmath as mp, slp_big, certify_big, lai_compare, hats_flagship as HF
    from decimal import Decimal, ROUND_FLOOR, getcontext
    getcontext().prec = 60
    c2, r2 = slp_big.polish(c, ("rect", 1.0), t_cap=float(ARG("cap", 150)))
    r_claim = (Decimal(repr(float(slp_big.rmin(c2, ("rect", 1.0))))) * (1 - Decimal(10) ** -12)).quantize(Decimal(10) ** -25, rounding=ROUND_FLOOR)
    ok, why = lai_compare.csq_gate(n, mp.mpf(str(r_claim)))
    row = {"N": n, "r_claim": str(r_claim), "gate": why or "claimable"}
    f = os.path.join(HERE, "cand_big_hp", "csq", f"csq_{n}.txt"); old = Decimal(open(f).readline().split()[1]) if os.path.exists(f) else Decimal(0)
    row["vs_file"] = float(r_claim / old - 1) if old else None
    if not ok or r_claim <= old: row["kept"] = False; return row
    tmp = f + f".{tag}.tmp"; rec = HF.table("csq")[n]
    open(tmp, "w", newline="\n").write(f"r {r_claim}\n" + "".join(f"{format(Decimal(repr(float(x))), 'f')} {format(Decimal(repr(float(y))), 'f')}\n" for x, y in c2))
    a = certify_big.check("square", tmp, n, rec)
    o = subprocess.run([sys.executable, os.path.join(HERE, "grok", "verify_exact_big.py"), "square", tmp, str(n), rec], capture_output=True, text=True).stdout
    b = [l for l in o.splitlines() if l.startswith("VERDICT")]; b = b[0].split(":")[1].strip() if b else "ERROR"
    row.update(claude=a, grok=b)
    if a == "IMPROVES" and b == "IMPROVES":
        os.replace(tmp, f); row["kept"] = True
        json.dump({"shelf": "csq", "N": n, "tag": "csq_gpu_big", "how": "GPU basin search (gpu_big.py) from our best packing, slp_big polish",
                   "r_new": str(r_claim), "claude": a, "grok": b, "lopt": "NOT_ATTEMPTED"}, open(f[:-4] + ".json", "w"), indent=1)
    else: os.remove(tmp); row["kept"] = False
    return row

def main():
    import torch
    from concurrent.futures import ProcessPoolExecutor
    dev = torch.device("cuda"); DT = torch.float64
    torch.cuda.set_per_process_memory_fraction(0.70)
    T = [int(t[1]) if isinstance(t, list) else int(t) for t in json.load(open(ARG("targets", None)))]
    B0 = int(ARG("batch", 32)); ROUNDS = int(ARG("rounds", 2)); tag = ARG("tag", "")
    out = open(os.path.join(HERE, "out", f"gpu_big{tag}.jsonl"), "a"); pool = ProcessPoolExecutor(int(ARG("cpu", 2))); futs = []; t_start = time.time()

    def walls(P): x, y = P[..., 0], P[..., 1]; return torch.stack([0.5 + x, 0.5 - x, 0.5 + y, 0.5 - y], -1)

    def vlist(Pc, cut):
        B, n = Pc.shape[:2]; I, J = [], []
        for b in range(B):
            pr = cKDTree(Pc[b]).query_pairs(cut, output_type="ndarray"); I.append(pr[:, 0] + b * n); J.append(pr[:, 1] + b * n)
        return torch.tensor(np.concatenate(I), device=dev), torch.tensor(np.concatenate(J), device=dev)

    for n in T:
        got = load_best(n)
        if got is None: continue
        c0, r_file, src = got; t0 = time.time(); r0 = float(min(cKDTree(c0).query(c0, k=2)[0][:, 1].min() / 2, (0.5 - np.abs(c0)).min()))
        rng = np.random.default_rng(n); best = (0.0, None)
        B = max(8, min(B0, int(4.0e5 // n)))
        for rnd in range(ROUNDS):
            P0 = np.repeat(c0[None], B, 0).copy()
            for b in range(B):                                   # several local shakes per copy + a small global jitter
                for _ in range(1 + n // 800):
                    cen = c0[rng.integers(n)]; rho = rng.choice([3.0, 5.0, 8.0]) * 2 * r0
                    m = ((c0 - cen) ** 2).sum(1) < rho * rho
                    P0[b, m] += rng.normal(0, rng.choice([0.1, 0.3, 0.6]) * r0, (m.sum(), 2))
                P0[b] += rng.normal(0, 0.02 * r0, (n, 2))
            P0 = np.clip(P0, -0.5, 0.5)
            P = torch.tensor(P0, device=dev, dtype=DT).requires_grad_(True)
            shrink = torch.tensor([0.99, 0.995, 0.998], device=dev, dtype=DT)[torch.arange(B, device=dev) % 3]
            opt = torch.optim.Adam([P], lr=0.02 * r0); T1, T2 = 300, 400; I = J = None; bidx = None
            for t in range(T1 + T2):
                if t % 40 == 0:
                    I, J = vlist(P.detach().cpu().numpy(), 2.6 * r0); bidx = I // n
                s = shrink + (1 + 2e-5 - shrink) * min(1.0, t / T1); rt = r0 * s
                Pf = P.reshape(-1, 2); d = torch.sqrt(((Pf[I] - Pf[J]) ** 2).sum(-1) + 1e-300)
                E = (torch.relu(2 * rt[bidx] - d) ** 2).sum() + (torch.relu(rt[:, None, None] - walls(P)) ** 2).sum()
                opt.zero_grad(); E.backward(); opt.step()
                with torch.no_grad(): P.clamp_(-0.5, 0.5)
                if t == T1: opt = torch.optim.Adam([P], lr=0.004 * r0)
            with torch.no_grad():
                Pc = P.detach().cpu().numpy()
                for b in range(B):
                    rb = min(cKDTree(Pc[b]).query(Pc[b], k=2)[0][:, 1].min() / 2, float((0.5 - np.abs(Pc[b])).min()))
                    if rb > best[0]: best = (rb, Pc[b].copy())
            del P, opt; torch.cuda.empty_cache()
        row = {"N": n, "src": os.path.relpath(src, HERE), "r0": r0, "gpu_best_rel": best[0] / r0 - 1, "B": B, "secs": round(time.time() - t0, 1)}
        out.write(json.dumps(row) + "\n"); out.flush()
        print(f"[{time.time() - t_start:.0f}s] csq N={n} gpu best {row['gpu_best_rel']:+.2e} (B={B}, {row['secs']}s)", flush=True)
        if best[1] is not None and best[0] > r0 * (1 - float(ARG("send", 3e-3))):   # 09-25: was 2e-4 -> almost nothing polished (raw GPU output sits ~3e-4 below until polish)
            futs.append(pool.submit(certify, n, best[1], "gpubig" + tag))
        for f in [f for f in futs if f.done()]:
            futs.remove(f); r = f.result(); out.write(json.dumps({"cpu": r}) + "\n"); out.flush()
            print(f"    CPU N={r['N']} vs file {r['vs_file']} {r['gate'][:50]}" + (" HIT" if r.get("kept") else ""), flush=True)
    for f in futs:
        r = f.result(); out.write(json.dumps({"cpu": r}) + "\n"); out.flush()
        print(f"    CPU N={r['N']} vs file {r['vs_file']} {r['gate'][:50]}" + (" HIT" if r.get("kept") else ""), flush=True)
    rows = [json.loads(l) for l in open(os.path.join(HERE, "out", f"gpu_big{tag}.jsonl"))]
    print(f"DONE {sum(1 for r in rows if 'cpu' in r and r['cpu'].get('kept'))} kept, {time.time() - t_start:.0f}s", flush=True)

if __name__ == "__main__":
    main()
