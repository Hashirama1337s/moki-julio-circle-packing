# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""cpt probe (2026-09-25): can a local optimiser beat Packomania cpt on Specht-program-only cells, N 101-200?

  population = N in 101..200 whose page reference is exactly [1] (E. Specht, program cpt, 2023)
  sample     = random.Random(20260925).sample(sorted(population), 16)
  per N      = start from the published packing -> SLP polish -> basin hopping {relocate 1-3 fewest-contact circles to the
               largest holes | shake a random disc of circles} -> polish, keep the best (sideways moves within 1e-12 accepted),
               until a CPU budget (process_time) of BUDGET s or a wall cap.
  BEATEN     = best > published * (1 + 1e-9)   (float margin)
  SEALED BAR = >= 8 of 16 BEATEN -> PASS
  positive control = first 2 sampled N: coords * 0.999 + U(+-1e-3 r) -> polish must recover >= published * (1 - 1e-9).
usage: py -3.11 probe.py [budget_cpu_s=300] [workers=4] [wall_cap_s=480]
"""
import os, sys, time, json, random
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geom
from scipy.spatial import cKDTree

OUT = os.path.join(geom.HERE, 'out'); SEED = 20260925


def _lower_priority():
    try:
        import psutil; psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass


def population():
    T = geom.page_table()
    return T, sorted(n for n in range(101, 201) if T[n][1] == ('1',))


def sample_pts(rng, m):
    lo = geom.VERT.min(0); hi = geom.VERT.max(0); out = []
    while sum(len(x) for x in out) < m:
        p = lo + (hi - lo) * rng.random((2 * m, 2)); out.append(p[geom.inside(p)])
    return np.concatenate(out)[:m]


def clearance(p, cs, r):
    g = (geom.A - p @ geom.NRM.T).min(1)
    if len(cs) == 0: return g
    d, _ = cKDTree(cs).query(p, k=1); return np.minimum(g, d - r)


def find_holes(cs, r, k, rng, m=6000):
    """k greedy largest holes (max clearance points) for the remaining centres cs; refine by local pattern search."""
    p = sample_pts(rng, m); rho = clearance(p, cs, r); holes = []
    for _ in range(k):
        top = np.argsort(-rho)[:8]; bestp, bestv = None, -np.inf
        ref = np.vstack([cs] + ([np.array(holes)] if holes else []))
        for t in top:
            x, v, s = p[t].copy(), rho[t], 0.5 * r
            while s > 1e-4 * r:
                cand = x + s * np.array([[1, 0], [-1, 0], [0, 1], [0, -1], [.7, .7], [-.7, .7], [.7, -.7], [-.7, -.7]])
                cand = cand[geom.inside(cand)]
                if len(cand) == 0: s /= 2; continue
                cv = clearance(cand, ref, r)
                j = int(np.argmax(cv))
                if cv[j] > v: x, v = cand[j], cv[j]
                else: s /= 2
            if v > bestv: bestp, bestv = x, v
        holes.append(bestp)
        rho = np.minimum(rho, np.sqrt(((p - bestp) ** 2).sum(1)) - r)
    return np.array(holes)


def move_relocate(c, r, rng):
    n = len(c); k = int(rng.choice([1, 2, 3], p=[0.5, 0.3, 0.2]))
    cnt = geom.contacts(c, r); key = cnt + rng.random(n) * 0.5          # random tie-break among equal contact counts
    pool = np.argsort(key)[:max(k, 2 * k)]; pick = rng.choice(pool, size=k, replace=False)
    keep = np.setdiff1d(np.arange(n), pick)
    h = find_holes(c[keep], r, k, rng)
    q = c.copy(); q[pick] = h + rng.uniform(-1e-3 * r, 1e-3 * r, h.shape); return geom.repair(q), f'reloc{k}'


def move_shake(c, r, rng):
    n = len(c); i0 = int(rng.integers(n)); R = r * rng.uniform(2.5, 6.0); s = r * rng.uniform(0.1, 0.6)
    idx = np.nonzero(((c - c[i0]) ** 2).sum(1) <= R * R)[0]
    ang = rng.uniform(0, 2 * np.pi, len(idx)); rad = s * np.sqrt(rng.random(len(idx)))
    q = c.copy(); q[idx] += np.stack([rad * np.cos(ang), rad * np.sin(ang)], 1); return geom.repair(q), 'shake'


def run_one(args):
    n, rpub_s, budget, wall_cap, seed = args
    _lower_priority()
    rng = np.random.default_rng(seed); rpub = float(rpub_s)
    t_cpu0 = time.process_time(); t_w0 = time.time()
    c0 = geom.load_coords(n); r_coords = geom.rmin(c0)
    best_c, best_r, _ = geom.polish(c0); r_pol0 = best_r
    hops = impr = side = 0; kinds = {}; hist = [(0.0, best_r)]
    while time.process_time() - t_cpu0 < budget and time.time() - t_w0 < wall_cap:
        if rng.random() < 0.5: q, kind = move_relocate(best_c, best_r, rng)
        else: q, kind = move_shake(best_c, best_r, rng)
        q, rq, _ = geom.polish(q, t_cap=60.0); hops += 1
        kinds.setdefault(kind, [0, 0]); kinds[kind][0] += 1
        if rq > best_r * (1 + 1e-13):
            best_c, best_r = q, rq; impr += 1; kinds[kind][1] += 1; hist.append((time.process_time() - t_cpu0, best_r))
        elif rq >= best_r * (1 - 1e-12):
            best_c = q; side += 1
    fc, fr, _ = geom.polish(best_c, max_iter=800)                    # final tight polish
    if fr >= best_r: best_c, best_r = fc, fr
    g, _ = geom.walls(best_c); iu = np.triu_indices(n, 1)
    d = np.sqrt(((best_c[:, None] - best_c[None]) ** 2).sum(-1))[iu]
    feas = float(max((best_r - g).max(), (2 * best_r - d).max()))   # float self-check of the claimed radius
    np.save(os.path.join(OUT, f'best_{n}.npy'), best_c)
    rel = best_r / rpub - 1
    return dict(n=n, r_pub=rpub_s, r_coords=r_coords, r_polish=r_pol0, rel_polish=r_pol0 / rpub - 1, r_best=best_r, rel=rel,
                beaten=bool(best_r > rpub * (1 + 1e-9)), hops=hops, improvements=impr, sideways=side, kinds=kinds,
                feas_violation=feas, cpu=time.process_time() - t_cpu0, wall=time.time() - t_w0, hist=hist[-12:])


def control(ns, T, seeds=(1, 2, 3)):
    rows = []
    for n in ns:
        rpub = float(T[n][0]); c = geom.load_coords(n); _, rp0, _ = geom.polish(c)
        for sd in seeds:
            rng = np.random.default_rng(1000 * n + sd)
            s = geom.repair(0.999 * c + rng.uniform(-1e-3 * rpub, 1e-3 * rpub, c.shape)); r_start = geom.rmin(s)
            t0 = time.process_time(); _, rr, it = geom.polish(s); cpu = time.process_time() - t0
            rows.append(dict(n=n, seed=sd, r_pub=T[n][0], r_start=r_start, r_rec=rr, rel_vs_pub=rr / rpub - 1,
                             rel_vs_polished_pub=rr / rp0 - 1, iters=it, cpu=cpu, ok=bool(rr >= rpub * (1 - 1e-9))))
    return rows


if __name__ == '__main__':
    _lower_priority()
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    wall_cap = float(sys.argv[3]) if len(sys.argv) > 3 else 480.0
    os.makedirs(OUT, exist_ok=True); T0 = time.time(); C0 = time.process_time()
    T, pop = population()
    samp = random.Random(SEED).sample(pop, 16)
    print(f'population (Specht-program-only, N 101-200): {len(pop)}\n{pop}', flush=True)
    print(f'sample (random.Random({SEED}).sample(pop, 16)): {samp}', flush=True)
    ctl = control(samp[:2], T)
    for x in ctl:
        print(f"CONTROL N={x['n']} seed {x['seed']}: start {x['r_start']:.9f} -> {x['r_rec']:.15f}  vs pub {x['rel_vs_pub']:+.3e}"
              f"  vs polished pub {x['rel_vs_polished_pub']:+.3e}  iters {x['iters']}  {'OK' if x['ok'] else 'FAIL'}", flush=True)
    rows = []
    with Pool(W) as pool:
        for x in pool.imap_unordered(run_one, [(n, T[n][0], budget, wall_cap, SEED + n) for n in samp]):
            rows.append(x)
            print(f"N={x['n']:3d} pub {x['r_pub']}  polish {x['r_polish']:.15f} ({x['rel_polish']:+.3e})  best {x['r_best']:.15f} "
                  f"({x['rel']:+.3e}) {'BEATEN' if x['beaten'] else '-'}  hops {x['hops']} impr {x['improvements']} "
                  f"side {x['sideways']} {x['kinds']} feas {x['feas_violation']:.1e} cpu {x['cpu']:.0f}s wall {x['wall']:.0f}s", flush=True)
            json.dump(dict(population=pop, sample=samp, control=ctl, rows=rows), open(os.path.join(OUT, 'probe_results.json'), 'w'), indent=1)
    nb = sum(r['beaten'] for r in rows)
    print(f'BEATEN {nb}/16 -> {"PASS" if nb >= 8 else "FAIL"}  (bar >= 8)   total wall {time.time() - T0:.0f}s  '
          f'worker cpu {sum(r["cpu"] for r in rows):.0f}s', flush=True)
