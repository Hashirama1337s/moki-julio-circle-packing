# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""chx float search (equal circles in a regular HEXAGON), 2026-09-26.  Moki&Julio.
Driver around the k-gon machinery (not edited): ../kgon/geom.py Frame(6) + SLP polish, ../kgon/search.py moves (relocate the
fewest-contact circles into the largest holes | shake a disc), imported by file path.

Per cell N:
  starts : Packomania's own packing (data/big/chx/chx<N>.txt) -> SLP polish;  Amore 2023's packing mapped into our frame
           (amore/chx_<N>.npy, the Amore mapping step) -> SLP polish;  the better one starts the walk.
  bar    : bar_eff = max(data/refs/chx/chx_bar.json bar, POLISHED Amore radius) -- a polish of HIS packing is not ours to claim.
  walk   : basin hopping exactly as ../kgon/search.run_one (sideways within 1e-12, best updated on > 1e-13 relative), CPU budget
           (process time) + per-cell wall cap + global deadline; final tight polish; dense O(N^2) float recheck.
  beaten : claim rule in float:  r > bar_eff (1 + 1e-10)  AND  r >= printed + 2e-12.
Every worker lowers itself to BelowNormal (chx_prio.lower) and records the class read back from the OS.
Output: <out>/best_<N>.npy, <out>/rows.jsonl (resumable: finished cells are skipped).
usage: python search_chx.py --cells 7,33,... --out out/probe [--budget 240] [--workers 4] [--deadline 14:40]
"""
import os, sys, time, json, argparse, datetime
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import chx_prio                   # noqa: E402
import chx_common as C            # noqa: E402
ks = C._load('kgon_search', os.path.join(ROOT, 'kgon', 'search.py'))

SEED = 20260926


def bars():
    return {int(n): g for n, g in json.load(open(os.path.join(C.REFS, 'chx_bar.json'))).items()}


def beats(r, bar_eff, printed):
    return bool(r > bar_eff * (1 + 1e-10) and r >= printed + 2e-12)


def run_one(args):
    n, budget, wall_cap, deadline, seed, out = args
    prio = chx_prio.lower()
    if time.time() >= deadline:
        return dict(n=n, skipped=True)
    F = C.frame(); G = bars()[n]
    rng = np.random.default_rng(seed); rpub = float(G['printed']); bar = float(G['bar'])
    t_cpu0 = time.process_time(); t_w0 = time.time()
    c0 = C.load_coords(n); r_coords = F.dense_r(c0)
    cp, rp, _ = F.polish(c0, max_iter=800)
    start, cur_c, cur_r = 'published', cp, rp
    ra_raw = ra_pol = None
    apath = os.path.join(HERE, 'amore', f'chx_{n}.npy')
    if os.path.exists(apath):
        ca = np.load(apath); ra_raw = F.dense_r(ca)
        ca, ra_pol, _ = F.polish(F.repair(ca), max_iter=800)
        if ra_pol > cur_r: start, cur_c, cur_r = 'amore', ca, ra_pol
    bar_eff = max(bar, ra_pol or 0.0)
    bar_eff_src = 'amore_polished' if (ra_pol or 0.0) > bar else G['bar_source']
    best_c, best_r = cur_c.copy(), cur_r
    hops = impr = side = 0; kinds = {}; hist = [(round(time.process_time() - t_cpu0, 1), best_r)]; stop = 'budget'
    while True:
        if time.process_time() - t_cpu0 >= budget: stop = 'budget'; break
        if time.time() - t_w0 >= wall_cap: stop = 'wall_cap'; break
        if time.time() >= deadline: stop = 'deadline'; break
        if rng.random() < 0.5: q, kind = ks.move_relocate(F, cur_c, cur_r, rng)
        else: q, kind = ks.move_shake(F, cur_c, cur_r, rng)
        q, rq, _ = F.polish(q, t_cap=60.0); hops += 1
        kinds.setdefault(kind, [0, 0]); kinds[kind][0] += 1
        if rq > best_r * (1 + 1e-13):
            best_c, best_r = q.copy(), rq; cur_c, cur_r = q, rq; impr += 1; kinds[kind][1] += 1
            hist.append((round(time.process_time() - t_cpu0, 1), best_r))
        elif rq >= best_r * (1 - 1e-12):
            cur_c, cur_r = q, rq; side += 1
    fc, fr, _ = F.polish(best_c, max_iter=800)
    if fr > best_r: best_c, best_r = fc, fr
    r_dense = F.dense_r(best_c)
    np.save(os.path.join(out, f'best_{n}.npy'), best_c)
    return dict(n=n, skipped=False, priority_class=prio, r_pub=G['printed'], page_refs=G['page_refs'], bar=G['bar'],
                bar_src=G['bar_source'], bar_eff=bar_eff, bar_eff_src=bar_eff_src, r_coords=r_coords, r_polish_pub=rp,
                rel_polish_pub=rp / rpub - 1, r_amore_raw=ra_raw, r_amore_polished=ra_pol, start=start, r_best=best_r,
                r_dense=r_dense, rel_vs_printed=r_dense / rpub - 1, rel_vs_bar=r_dense / bar_eff - 1,
                beaten_float=beats(r_dense, bar_eff, rpub), hops=hops, improvements=impr, sideways=side, kinds=kinds,
                stop=stop, cpu=round(time.process_time() - t_cpu0, 1), wall=round(time.time() - t_w0, 1), hist=hist[-12:])


def done_set(out):
    s = set(); p = os.path.join(out, 'rows.jsonl')
    if os.path.exists(p):
        for line in open(p):
            try:
                x = json.loads(line)
                if not x.get('skipped'): s.add(x['n'])
            except Exception:
                pass
    return s


if __name__ == '__main__':
    print('main priority class:', chx_prio.lower(), flush=True)
    ap = argparse.ArgumentParser()
    ap.add_argument('--cells', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--budget', type=float, default=240.0)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--wall-cap', type=float, default=0.0)
    ap.add_argument('--deadline', type=str, default='14:40', help='local HH:MM today')
    ap.add_argument('--seed-offset', type=int, default=0)
    a = ap.parse_args()
    W = min(a.workers, 4)                                          # house rule for this task: at most 4 workers
    out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    wall_cap = a.wall_cap or 2 * a.budget + 120
    hh, mm = map(int, a.deadline.split(':'))
    deadline = datetime.datetime.now().replace(hour=hh, minute=mm, second=0, microsecond=0).timestamp()
    cells = [int(x) for x in a.cells.split(',') if x]
    d = done_set(out); todo = [n for n in cells if n not in d]
    T0 = time.time()
    print(f'to run: {len(todo)} cells {todo}; budget {a.budget:.0f}s CPU/cell, wall cap {wall_cap:.0f}s, workers {W}, '
          f'deadline {a.deadline}', flush=True)
    jobs = [(n, a.budget, wall_cap, deadline, SEED + n + a.seed_offset, out) for n in todo]
    nb = nd = ns = 0
    with Pool(W) as pool:
        for x in pool.imap_unordered(run_one, jobs):
            with open(os.path.join(out, 'rows.jsonl'), 'a') as f: f.write(json.dumps(x) + '\n')
            if x.get('skipped'):
                ns += 1; print(f"N={x['n']:3d} SKIPPED (deadline)", flush=True); continue
            nd += 1; nb += x['beaten_float']
            print(f"N={x['n']:3d} [{','.join(x['page_refs'])}] pub {x['r_pub']} bar {x['bar_eff']:.15f} ({x['bar_eff_src']}) "
                  f"start {x['start']}  best {x['r_dense']:.15f} (vs printed {x['rel_vs_printed']:+.3e}, vs bar {x['rel_vs_bar']:+.3e}) "
                  f"{'BEATEN' if x['beaten_float'] else '-'}  hops {x['hops']} impr {x['improvements']} stop {x['stop']} "
                  f"cpu {x['cpu']:.0f}s wall {x['wall']:.0f}s prio {x['priority_class']}  [{nd} done, {nb} beaten, {ns} skipped, "
                  f"elapsed {time.time() - T0:.0f}s]", flush=True)
    print(f'END: {nd} run, {nb} beaten in float (claim rule vs bar_eff), {ns} skipped; wall {time.time() - T0:.0f}s', flush=True)
