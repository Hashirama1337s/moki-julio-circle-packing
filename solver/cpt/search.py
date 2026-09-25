# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""cpt sweep (2026-09-25): polish + basin-hop every Packomania cpt cell N in 49..200 credited ONLY to ref [1]
(E. Specht, program cpt, 2023), starting from the published packing. Float search only; exact certification is
checkers/certify_poly.py + checkers/verify_exact_poly.py.

Frame (see geom.py, verified there on all 200 published packings): regular pentagon, circumradius 1, centred at
the origin, vertex at (0, 1), flat bottom side, every side at distance a = cos 36 deg.

Per N:  published packing -> SLP polish (geom.polish) -> basin hopping with the probe's two moves
        (relocate 1-3 fewest-contact circles into the largest holes | shake a random disc of circles), each followed by
        a polish. The walk keeps a CURRENT state (sideways moves within 1e-12 of the best are accepted) and a separate
        BEST state (updated only on a gain > 1e-13 relative), so the reported best radius is always the true float
        min-radius of the saved coordinates. Stops at the CPU budget (process time), a per-N wall cap, or the global
        deadline, whichever comes first; then a final tight polish; then a dense O(N^2) float recheck.
Order:  the 108 cells NOT sampled by the probe first (descending N: the stale block is at large N), then the 16 probe
        cells (they already have a probe float best in cpt_probe/out/). If the global deadline cuts the sweep, the
        cells not started are reported as skipped.
Outputs: cpt/out/best_<N>.npy, cpt/out/sweep_rows.jsonl (one JSON row per finished N; --resume skips those N).

usage: python search.py [--budget 240] [--workers 4] [--wall-cap 600] [--deadline-h 2.5] [--only 150,151] [--no-resume]
Launch long runs at BelowNormal (Start-Process + PriorityClass); every worker also lowers itself to BelowNormal.
"""
import os, sys, time, json, argparse
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)                          # geom.py, probe.py sit next to this file
import geom                                     # noqa: E402  (geom.py: frame, loaders, SLP polish)
import probe                                    # noqa: E402  (probe.py: move_relocate, move_shake)

OUT = os.path.join(HERE, 'out')
ROWS = os.path.join(OUT, 'sweep_rows.jsonl')
SEED = 20260925
PROBE_SAMPLE = [151, 122, 175, 121, 164, 182, 194, 188, 139, 172, 159, 115, 176, 142, 102, 146]


def lower_priority():
    try:
        import psutil; psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass


def population():
    """N in 49..200 whose Packomania cpt reference list is exactly [1] (Specht's program)."""
    T = geom.page_table()
    return T, [n for n in range(49, 201) if T[n][1] == ('1',)]


def order(pop):
    first = sorted((n for n in pop if n not in PROBE_SAMPLE), reverse=True)
    last = sorted((n for n in pop if n in PROBE_SAMPLE), reverse=True)
    return first + last


def dense_r(c):
    """Float min-radius from ALL pairs (O(N^2)) and the five wall slacks: an independent recheck of geom.rmin."""
    g, _ = geom.walls(c); n = len(c)
    if n == 1: return float(g.min())
    iu = np.triu_indices(n, 1)
    d = np.sqrt(((c[:, None, :] - c[None, :, :]) ** 2).sum(-1))[iu]
    return float(min(g.min(), d.min() / 2))


def run_one(args):
    n, rpub_s, budget, wall_cap, deadline, seed, out = args
    lower_priority()
    if time.time() >= deadline:
        return dict(n=n, r_pub=rpub_s, skipped=True)
    rng = np.random.default_rng(seed); rpub = float(rpub_s)
    t_cpu0 = time.process_time(); t_w0 = time.time()
    c0 = geom.load_coords(n); r_coords = dense_r(c0)
    cur_c, cur_r, _ = geom.polish(c0); r_pol0 = cur_r
    best_c, best_r = cur_c.copy(), cur_r
    hops = impr = side = 0; kinds = {}; hist = [(0.0, best_r)]; stop = 'budget'
    while True:
        if time.process_time() - t_cpu0 >= budget: stop = 'budget'; break
        if time.time() - t_w0 >= wall_cap: stop = 'wall_cap'; break
        if time.time() >= deadline: stop = 'deadline'; break
        if rng.random() < 0.5: q, kind = probe.move_relocate(cur_c, cur_r, rng)
        else: q, kind = probe.move_shake(cur_c, cur_r, rng)
        q, rq, _ = geom.polish(q, t_cap=60.0); hops += 1
        kinds.setdefault(kind, [0, 0]); kinds[kind][0] += 1
        if rq > best_r * (1 + 1e-13):
            best_c, best_r = q.copy(), rq; cur_c, cur_r = q, rq; impr += 1; kinds[kind][1] += 1
            hist.append((round(time.process_time() - t_cpu0, 1), best_r))
        elif rq >= best_r * (1 - 1e-12):
            cur_c, cur_r = q, rq; side += 1                         # sideways: move the walk, keep the best separate
    fc, fr, _ = geom.polish(best_c, max_iter=800)                   # final tight polish
    if fr > best_r: best_c, best_r = fc, fr
    r_dense = dense_r(best_c)
    np.save(os.path.join(out, f'best_{n}.npy'), best_c)
    return dict(n=n, r_pub=rpub_s, skipped=False, r_coords=r_coords, r_polish=r_pol0, rel_polish=r_pol0 / rpub - 1,
                r_best=best_r, r_dense=r_dense, rel=r_dense / rpub - 1, beaten_float=bool(r_dense > rpub * (1 + 1e-9)),
                hops=hops, improvements=impr, sideways=side, kinds=kinds, stop=stop,
                cpu=round(time.process_time() - t_cpu0, 1), wall=round(time.time() - t_w0, 1), hist=hist[-12:])


def done_set(rows):
    s = set()
    if os.path.exists(rows):
        for line in open(rows):
            try:
                x = json.loads(line)
                if not x.get('skipped'): s.add(x['n'])
            except Exception:
                pass
    return s


if __name__ == '__main__':
    lower_priority()
    ap = argparse.ArgumentParser()
    ap.add_argument('--budget', type=float, default=240.0)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--wall-cap', type=float, default=600.0)
    ap.add_argument('--deadline-h', type=float, default=2.5)
    ap.add_argument('--only', type=str, default='')
    ap.add_argument('--no-resume', action='store_true')
    ap.add_argument('--out', type=str, default=OUT)
    a = ap.parse_args()
    W = min(a.workers, 4)                                          # house rule: at most 4 worker processes
    out = os.path.abspath(a.out); rows_path = os.path.join(out, 'sweep_rows.jsonl')
    os.makedirs(out, exist_ok=True)
    T0 = time.time(); deadline = T0 + 3600.0 * a.deadline_h
    T, pop = population()
    todo = order(pop)
    if a.only: todo = [n for n in todo if n in {int(x) for x in a.only.split(',')}]
    if not a.no_resume:
        d = done_set(rows_path); todo = [n for n in todo if n not in d]
    print(f'population (Specht-program-only, N 49-200): {len(pop)}  (49-100: {sum(n <= 100 for n in pop)}, '
          f'101-200: {sum(n > 100 for n in pop)})', flush=True)
    print(f'to run now: {len(todo)}  budget {a.budget:.0f}s CPU/N  wall cap {a.wall_cap:.0f}s/N  workers {W}  '
          f'deadline {a.deadline_h} h  order {todo}', flush=True)
    jobs = [(n, T[n][0], a.budget, a.wall_cap, deadline, SEED + n, out) for n in todo]
    nb = nd = ns = 0
    with Pool(W) as pool:
        for x in pool.imap_unordered(run_one, jobs):
            with open(rows_path, 'a') as f: f.write(json.dumps(x) + '\n')
            if x.get('skipped'):
                ns += 1; print(f"N={x['n']:3d} SKIPPED (deadline)", flush=True); continue
            nd += 1; nb += x['beaten_float']
            print(f"N={x['n']:3d} pub {x['r_pub']}  polish {x['r_polish']:.15f} ({x['rel_polish']:+.3e})  best {x['r_dense']:.15f} "
                  f"({x['rel']:+.3e}) {'BEATEN' if x['beaten_float'] else '-'}  hops {x['hops']} impr {x['improvements']} "
                  f"side {x['sideways']} stop {x['stop']} cpu {x['cpu']:.0f}s wall {x['wall']:.0f}s   "
                  f"[{nd} done, {nb} beaten, {ns} skipped, elapsed {time.time() - T0:.0f}s]", flush=True)
    print(f'SWEEP END: {nd} run, {nb} beaten in float (> 1e-9 rel), {ns} skipped; wall {time.time() - T0:.0f}s', flush=True)
