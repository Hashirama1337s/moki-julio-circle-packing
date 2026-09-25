"""GPU basin search (2026-09-24). For each target (shelf, N): B perturbed copies of OUR best packing relax
together on the GPU (float64) under an overlap energy at a target radius that first SHRINKS (room to rearrange) and then GROWS a
little PAST our record; the true radius of every copy is measured at the end. The best few distinct copies are polished on the CPU
(slp_circ.polish) and, if they beat our best by > 1e-10, sent through the full exact chain (cert_candidate.run: both checkers,
lopt; kept only if better than the file on disk).
Perturbations (thirds of the batch): relocate 1-3 low-contact circles to random points; shake a random disc of circles; global jitter.
usage: py -3.11 gpu_mbh.py --targets=out/x.json [--batch=96] [--rounds=2] [--cpu=2] [--tag=] [--stop=HH:MM]  -> out/gpu_mbh<tag>.jsonl
"""
import os, sys, json, time, datetime, math
for v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"): os.environ[v] = "1"
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ARG = lambda k, d: next((a.split("=", 1)[1] for a in sys.argv if a.startswith(f"--{k}=")), d)

def cpu_job(a):
    # ONE job per size (probe bug 17:4x: two candidates of the same size certified in parallel shared a temp file): polish every
    # candidate, then certify only the best one.
    shelf, n, npys, r_best, tag = a
    import io, contextlib, slp_circ, finalize_circ, cert_candidate
    cont = finalize_circ.info(shelf)[0]; r2, c2, npy = -1.0, None, npys[0]
    for p in npys:
        cc, rr = slp_circ.polish(np.load(p), cont, t_cap=90.0)
        if rr > r2: r2, c2, npy = rr, cc, p
    row = {"shelf": shelf, "N": n, "r_polish": r2, "rel_polish": r2 / r_best - 1, "n_cand": len(npys)}
    if r2 > r_best * (1 + 1e-10):
        np.save(npy, c2)
        try:
            with contextlib.redirect_stdout(io.StringIO()): cr = cert_candidate.run(shelf, n, npy, tag=tag)
            row.update({k: cr.get(k) for k in ("gain_vs_ours", "gain_vs_packomania", "checker_a", "checker_b", "lopt", "kept")})
        except Exception as e: row["error"] = repr(e)
    return row

def gpu_main():
    import torch
    from concurrent.futures import ProcessPoolExecutor
    dev = torch.device("cuda"); DT = torch.float64
    import finalize_circ, transplant_probe as tp, slp_circ
    T = [tuple(t) for t in json.load(open(ARG("targets", None)))]
    B0 = int(ARG("batch", 96)); ROUNDS = int(ARG("rounds", 2)); tag = ARG("tag", ""); W = int(ARG("cpu", 2))
    stop = ARG("stop", None)
    if stop:
        now = datetime.datetime.now(); hh, mm = map(int, stop.split(":")); stop = now.replace(hour=hh, minute=mm, second=0)
        if stop <= now: stop += datetime.timedelta(days=1)
    out = open(os.path.join(HERE, "out", f"gpu_mbh{tag}.jsonl"), "a"); cdir = os.path.join(HERE, "out", "gpu_cand"); os.makedirs(cdir, exist_ok=True)
    S2 = math.sqrt(2.0)

    def walls(P, cont):
        x, y = P[..., 0], P[..., 1]
        if cont[0] == "rect": h = cont[1]; return torch.stack([0.5 + x, 0.5 - x, h / 2 + y, h / 2 - y], -1)
        if cont[0] == "tri": return torch.stack([x, y, (1 - x - y) / S2], -1)
        rr = torch.sqrt(x * x + y * y + 1e-300)
        if cont[0] == "quad": return torch.stack([x, y, 1 - rr], -1)
        if cont[0] == "semi": return torch.stack([y, 1 - rr], -1)
        return (1 - rr)[..., None]

    def pair_d(P):
        dx = P[:, :, None, :] - P[:, None, :, :]; return torch.sqrt((dx * dx).sum(-1) + 1e-300)

    def radius(P, cont):
        n = P.shape[1]; D = pair_d(P) + torch.eye(n, device=dev, dtype=DT) * 9
        return torch.minimum(D.amin((1, 2)) / 2, walls(P, cont).amin((1, 2)))

    def repair(P, cont):
        if cont[0] == "rect":
            return torch.stack([P[..., 0].clamp(-0.5, 0.5), P[..., 1].clamp(-cont[1] / 2, cont[1] / 2)], -1)
        if cont[0] == "tri":
            P = P.clamp(0, 1); o = (P.sum(-1) - 1).clamp(min=0); return P - o[..., None] / 2
        rr = torch.sqrt((P * P).sum(-1) + 1e-300)[..., None]; P = torch.where(rr > 1, P / rr, P)
        if cont[0] == "quad": return P.clamp(min=0)
        if cont[0] == "semi": return torch.stack([P[..., 0], P[..., 1].clamp(min=0)], -1)
        return P

    torch.cuda.set_per_process_memory_fraction(0.70)          # 18:1x: the card filled up (16.0/16.3 GB) and spilled to system RAM (20x slower)
    pool = ProcessPoolExecutor(W); futs = []; t_start = time.time(); nt = 0
    for shelf, n in T:
        if stop and datetime.datetime.now() >= stop: break
        cont = finalize_circ.info(shelf)[0]; cur = tp.packing(shelf, n)
        if cur is None: continue
        c0 = np.asarray(cur[0], dtype=np.float64); r0 = float(cur[1] or slp_circ.rmin(c0, cont)); t0 = time.time()
        g = torch.Generator(device=dev); g.manual_seed(hash((shelf, n)) % 2 ** 31)
        C0 = torch.tensor(c0, device=dev, dtype=DT)
        # contact degree (low-contact circles are relocated first)
        D0 = pair_d(C0[None])[0] + torch.eye(n, device=dev, dtype=DT) * 9
        deg = (D0 < 2 * r0 * (1 + 1e-6)).sum(1) + (walls(C0[None], cont)[0] < r0 * (1 + 1e-6)).sum(1)
        low = torch.argsort(deg.to(DT) + torch.rand(n, device=dev, dtype=DT, generator=g) * 0.5)[: max(3, n // 10)]
        best_rows = []; Bn = max(8, min(B0, int(3.0e7 // (n * n))))      # batch scaled to N^2 (memory)
        for rnd in range(ROUNDS):
          B = Bn
          while True:
           try:
            P = C0.expand(B, n, 2).clone(); k = B // 3
            # 1) relocation of 1-3 low-contact circles to random points (repair() pulls them inside)
            for b in range(k):
                m = 1 + b % 3; idx = low[torch.randint(len(low), (m,), device=dev, generator=g)]
                P[b, idx] = (torch.rand(m, 2, device=dev, dtype=DT, generator=g) - 0.5) * 1.2 + C0.mean(0)
            # 2) shake a random disc of radius rho * 2r0 by sigma * r0
            cen = C0[torch.randint(n, (B - k,), device=dev, generator=g)]
            rho = torch.tensor([2.0, 3.0, 5.0], device=dev, dtype=DT)[torch.randint(3, (B - k,), device=dev, generator=g)] * 2 * r0
            sig = torch.tensor([0.1, 0.3, 0.6], device=dev, dtype=DT)[torch.randint(3, (B - k,), device=dev, generator=g)] * r0
            inside = ((C0[None] - cen[:, None]) ** 2).sum(-1) < rho[:, None] ** 2
            half = (B - k) // 2
            inside[half:] = True; sig[half:] = sig[half:] * 0.15            # 3) the last part: global jitter
            P[k:] = P[k:] + inside[..., None] * sig[:, None, None] * torch.randn(B - k, n, 2, device=dev, dtype=DT, generator=g)
            P = repair(P, cont).requires_grad_(True)
            shrink = torch.tensor([0.97, 0.985, 0.995], device=dev, dtype=DT)[torch.arange(B, device=dev) % 3]
            opt = torch.optim.Adam([P], lr=0.02 * r0); T1, T2 = 250, 350
            for t in range(T1 + T2):
                s = shrink + (1 + 2e-5 - shrink) * min(1.0, t / T1)
                rt = r0 * s
                D = pair_d(P); ov = torch.relu(2 * rt[:, None, None] - D).triu(1)
                E = (ov * ov).sum((1, 2)) + (torch.relu(rt[:, None, None] - walls(P, cont)) ** 2).sum((1, 2))
                opt.zero_grad(); E.sum().backward(); opt.step()
                with torch.no_grad(): P.copy_(repair(P, cont))
                if t == T1: opt = torch.optim.Adam([P], lr=0.004 * r0)
            with torch.no_grad():
                R = radius(P, cont); mv = (P - C0[None]).norm(dim=-1).amax(1) / r0
                for b in torch.argsort(-R)[:3].tolist():
                    if mv[b] > 0.05: best_rows.append((float(R[b]), P[b].detach().cpu().numpy()))
            break
           except torch.cuda.OutOfMemoryError:
            P = None; torch.cuda.empty_cache(); B = max(4, B // 2)
            if B == 4: break
        torch.cuda.empty_cache()
        best_rows.sort(key=lambda x: -x[0]); sent = 0; npys = []
        for rg, cb in best_rows[:2]:
            if rg < r0 * (1 - 3e-3): continue
            npy = os.path.join(cdir, f"{shelf}_{n}_{sent}.npy"); np.save(npy, cb); sent += 1; npys.append(npy)
        if npys: futs.append(pool.submit(cpu_job, (shelf, n, npys, r0, "gpu" + tag)))
        nt += 1
        row = {"shelf": shelf, "N": n, "r0": r0, "gpu_best_rel": (best_rows[0][0] / r0 - 1) if best_rows else None, "sent": sent,
               "secs": round(time.time() - t0, 1)}
        out.write(json.dumps(row) + "\n"); out.flush()
        print(f"[{time.time() - t_start:.0f}s] {shelf} N={n} gpu best {row['gpu_best_rel']} sent {sent} ({row['secs']}s)", flush=True)
        done = [f for f in futs if f.done()]
        for f in done:
            futs.remove(f); r = f.result(); out.write(json.dumps({"cpu": r}) + "\n"); out.flush()
            ok = r.get("kept") and r.get("checker_a") == "IMPROVES" and r.get("checker_b") == "IMPROVES"
            print(f"    CPU {r['shelf']} N={r['N']} polish {r['rel_polish']:+.2e}" + (f" HIT {r['gain_vs_ours'] if r['gain_vs_ours'] is None else format(r['gain_vs_ours'], '+.2e')} vs ours (None = new size), {r['gain_vs_packomania']:+.2e} vs Packomania" if ok else ""), flush=True)
    for f in futs:
        r = f.result(); out.write(json.dumps({"cpu": r}) + "\n"); out.flush()
        ok = r.get("kept") and r.get("checker_a") == "IMPROVES" and r.get("checker_b") == "IMPROVES"
        print(f"    CPU {r['shelf']} N={r['N']} polish {r['rel_polish']:+.2e}" + (f" HIT {r['gain_vs_ours']} vs ours" if ok else ""), flush=True)
    rows = [json.loads(l) for l in open(os.path.join(HERE, "out", f"gpu_mbh{tag}.jsonl"))]
    hits = sum(1 for r in rows if "cpu" in r and r["cpu"].get("kept") and r["cpu"].get("checker_a") == "IMPROVES" and r["cpu"].get("checker_b") == "IMPROVES")
    print(f"DONE {nt} targets, {hits} certified kept hits, {time.time() - t_start:.0f}s", flush=True)

if __name__ == "__main__":
    gpu_main()
