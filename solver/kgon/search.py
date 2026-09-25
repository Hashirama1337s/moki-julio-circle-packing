# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""k-gon sweep (cxd = regular 16-gon, cpd = regular 15-gon), 2026-09-25.  Moki&Julio.
k/table-parametrised version of ../cpt/search.py.  Float search only; exact certification is
checkers/certify_kgon.py + checkers/verify_exact_kgon.py.

Frame: geom.py (regular k-gon, circumradius 1, flat bottom side, verified on every published packing).

Per cell (table, N):
  starts : the published packing -> SLP polish;  Amore's 2023 packing mapped into our frame (amore/<table>_<N>.npy,
           made by the gate script from Zenodo 7574070) -> SLP polish;  the better one is the start of the walk.
  bar    : bar_eff = max(gate bar (out/gate_<table>.json: printed, Amore recomputed, Lai 2025 upper), POLISHED Amore
           radius).  The polished Amore radius enters the bar because his centres carry ~1e-10 noise: a polish of HIS
           packing is not ours to claim.
  walk   : basin hopping with the cpt moves (relocate 1-3 fewest-contact circles into the largest holes | shake a random
           disc of circles), each followed by a polish.  CURRENT state accepts sideways moves within 1e-12; BEST is updated
           only on a gain > 1e-13 relative, so the reported radius is always the true float min-radius of the saved centres.
           Stops at the CPU budget (process time), a per-cell wall cap, or the global deadline; then a final tight polish;
           then a dense O(N^2) float recheck.
  beaten : claim rule in float:  r > bar_eff (1 + 1e-10)  AND  r >= printed + 2e-12.
Cells:
  --mode probe : the 12 probe cells (6 per table, evenly spread over the Specht-program-only cells N >= 49) -> out/probe/
  --mode sweep : every Specht-program-only cell: cxd N 49-150, then cpd N 49-122, then cxd N 151-200 (Lai-gated), descending
                 N within each block; the probe cells go last.  -> out/<table>/best_<N>.npy, out/<table>/sweep_rows.jsonl
                 (--resume is the default: cells with a finished row are skipped).
usage: python search.py --mode probe|sweep [--budget 300] [--workers 5] [--deadline-h 3.5] [--only cxd:150,cpd:99]
Launch long runs at BelowNormal (Start-Process + PriorityClass); every worker also lowers itself to BelowNormal.
"""
import os, sys, time, json, argparse
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import geom                                     # noqa: E402

OUT = os.path.join(HERE, 'out')
SEED = 20260925
TOFF = {'cxd': 0, 'cpd': 100000}


def lower_priority():
    try:
        import psutil; psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass


# ---------------------------------------------------------------------- moves (../cpt/probe.py, generalised to a Frame)
def sample_pts(F, rng, m):
    lo = F.VERT.min(0); hi = F.VERT.max(0); out = []
    while sum(len(x) for x in out) < m:
        p = lo + (hi - lo) * rng.random((2 * m, 2)); out.append(p[F.inside(p)])
    return np.concatenate(out)[:m]


def clearance(F, p, cs, r):
    g = (F.A - p @ F.NRM.T).min(1)
    if len(cs) == 0: return g
    d, _ = cKDTree(cs).query(p, k=1); return np.minimum(g, d - r)


def find_holes(F, cs, r, k, rng, m=6000):
    """k greedy largest holes (max clearance points) for the remaining centres cs; refined by a local pattern search."""
    p = sample_pts(F, rng, m); rho = clearance(F, p, cs, r); holes = []
    for _ in range(k):
        top = np.argsort(-rho)[:8]; bestp, bestv = None, -np.inf
        ref = np.vstack([cs] + ([np.array(holes)] if holes else []))
        for t in top:
            x, v, s = p[t].copy(), rho[t], 0.5 * r
            while s > 1e-4 * r:
                cand = x + s * np.array([[1, 0], [-1, 0], [0, 1], [0, -1], [.7, .7], [-.7, .7], [.7, -.7], [-.7, -.7]])
                cand = cand[F.inside(cand)]
                if len(cand) == 0: s /= 2; continue
                cv = clearance(F, cand, ref, r)
                j = int(np.argmax(cv))
                if cv[j] > v: x, v = cand[j], cv[j]
                else: s /= 2
            if v > bestv: bestp, bestv = x, v
        holes.append(bestp)
        rho = np.minimum(rho, np.sqrt(((p - bestp) ** 2).sum(1)) - r)
    return np.array(holes)


def move_relocate(F, c, r, rng):
    n = len(c); k = int(rng.choice([1, 2, 3], p=[0.5, 0.3, 0.2]))
    cnt = F.contacts(c, r); key = cnt + rng.random(n) * 0.5          # random tie-break among equal contact counts
    pool = np.argsort(key)[:max(k, 2 * k)]; pick = rng.choice(pool, size=k, replace=False)
    keep = np.setdiff1d(np.arange(n), pick)
    h = find_holes(F, c[keep], r, k, rng)
    q = c.copy(); q[pick] = h + rng.uniform(-1e-3 * r, 1e-3 * r, h.shape); return F.repair(q), f'reloc{k}'


def move_shake(F, c, r, rng):
    n = len(c); i0 = int(rng.integers(n)); R = r * rng.uniform(2.5, 6.0); s = r * rng.uniform(0.1, 0.6)
    idx = np.nonzero(((c - c[i0]) ** 2).sum(1) <= R * R)[0]
    ang = rng.uniform(0, 2 * np.pi, len(idx)); rad = s * np.sqrt(rng.random(len(idx)))
    q = c.copy(); q[idx] += np.stack([rad * np.cos(ang), rad * np.sin(ang)], 1); return F.repair(q), 'shake'


# ---------------------------------------------------------------------- cells
def gate(table):
    return {int(n): g for n, g in json.load(open(os.path.join(OUT, f'gate_{table}.json'))).items()}


def populations():
    """{table: sorted Specht-program-only N >= 49}."""
    return {t: geom.program_only(t, 49) for t in geom.TABLES}


def probe_cells(pop=None):
    pop = pop or populations(); out = []
    for t in ('cxd', 'cpd'):
        P = pop[t]; L = len(P)
        out += [(t, P[min(L - 1, int((i + 0.5) * L / 6))]) for i in range(6)]
    return out


def sweep_order(pop=None):
    pop = pop or populations(); pr = set(probe_cells(pop))
    A = [('cxd', n) for n in sorted(pop['cxd'], reverse=True) if n <= 150]
    B = [('cpd', n) for n in sorted(pop['cpd'], reverse=True)]
    C = [('cxd', n) for n in sorted(pop['cxd'], reverse=True) if n > 150]
    allc = A + B + C
    return [c for c in allc if c not in pr] + [c for c in allc if c in pr]


def beats(r, bar_eff, printed):
    return bool(r > bar_eff * (1 + 1e-10) and r >= printed + 2e-12)


def run_one(args):
    table, n, budget, wall_cap, deadline, seed, out = args
    lower_priority()
    if time.time() >= deadline:
        return dict(table=table, n=n, skipped=True)
    k = geom.TABLES[table]; F = geom.frame(k); G = gate(table)[n]
    rng = np.random.default_rng(seed); rpub = float(G['printed']); bar = float(G['bar'])
    t_cpu0 = time.process_time(); t_w0 = time.time()
    c0 = geom.load_coords(table, n); r_coords = F.dense_r(c0)
    cp, rp, _ = F.polish(c0, max_iter=800)
    start, cur_c, cur_r = 'published', cp, rp
    ra_raw = ra_pol = None
    apath = os.path.join(HERE, 'amore', f'{table}_{n}.npy')
    if os.path.exists(apath):
        ca = np.load(apath); ra_raw = F.dense_r(ca)
        ca, ra_pol, _ = F.polish(F.repair(ca), max_iter=800)
        if ra_pol > cur_r: start, cur_c, cur_r = 'amore', ca, ra_pol
    bar_eff = max(bar, ra_pol or 0.0)
    bar_eff_src = 'amore_polished' if (ra_pol or 0.0) > bar else G['bar_src']
    best_c, best_r = cur_c.copy(), cur_r
    hops = impr = side = 0; kinds = {}; hist = [(round(time.process_time() - t_cpu0, 1), best_r)]; stop = 'budget'
    while True:
        if time.process_time() - t_cpu0 >= budget: stop = 'budget'; break
        if time.time() - t_w0 >= wall_cap: stop = 'wall_cap'; break
        if time.time() >= deadline: stop = 'deadline'; break
        if rng.random() < 0.5: q, kind = move_relocate(F, cur_c, cur_r, rng)
        else: q, kind = move_shake(F, cur_c, cur_r, rng)
        q, rq, _ = F.polish(q, t_cap=60.0); hops += 1
        kinds.setdefault(kind, [0, 0]); kinds[kind][0] += 1
        if rq > best_r * (1 + 1e-13):
            best_c, best_r = q.copy(), rq; cur_c, cur_r = q, rq; impr += 1; kinds[kind][1] += 1
            hist.append((round(time.process_time() - t_cpu0, 1), best_r))
        elif rq >= best_r * (1 - 1e-12):
            cur_c, cur_r = q, rq; side += 1                         # sideways: move the walk, keep the best separate
    fc, fr, _ = F.polish(best_c, max_iter=800)                      # final tight polish
    if fr > best_r: best_c, best_r = fc, fr
    r_dense = F.dense_r(best_c)
    os.makedirs(os.path.join(out, table), exist_ok=True)
    np.save(os.path.join(out, table, f'best_{n}.npy'), best_c)
    return dict(table=table, n=n, k=k, skipped=False, r_pub=G['printed'], bar=G['bar'], bar_src=G['bar_src'],
                bar_eff=bar_eff, bar_eff_src=bar_eff_src, r_coords=r_coords, r_polish_pub=rp, rel_polish_pub=rp / rpub - 1,
                r_amore_raw=ra_raw, r_amore_polished=ra_pol, start=start, r_best=best_r, r_dense=r_dense,
                rel_vs_printed=r_dense / rpub - 1, rel_vs_bar=r_dense / bar_eff - 1, beaten_float=beats(r_dense, bar_eff, rpub),
                hops=hops, improvements=impr, sideways=side, kinds=kinds, stop=stop,
                cpu=round(time.process_time() - t_cpu0, 1), wall=round(time.time() - t_w0, 1), hist=hist[-12:])


def done_set(out):
    s = set()
    for t in geom.TABLES:
        p = os.path.join(out, t, 'sweep_rows.jsonl')
        if os.path.exists(p):
            for line in open(p):
                try:
                    x = json.loads(line)
                    if not x.get('skipped'): s.add((x['table'], x['n']))
                except Exception:
                    pass
    return s


if __name__ == '__main__':
    lower_priority()
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['probe', 'sweep'], required=True)
    ap.add_argument('--budget', type=float, default=300.0)
    ap.add_argument('--workers', type=int, default=5)
    ap.add_argument('--wall-cap', type=float, default=0.0, help='per-cell wall cap in s (default 2 x budget + 120)')
    ap.add_argument('--deadline-h', type=float, default=3.5)
    ap.add_argument('--only', type=str, default='', help='e.g. cxd:150,cpd:99')
    ap.add_argument('--no-resume', action='store_true')
    ap.add_argument('--out', type=str, default='')
    a = ap.parse_args()
    W = min(a.workers, 5)                                          # house rule: at most 5 worker processes
    out = os.path.abspath(a.out or (os.path.join(OUT, 'probe') if a.mode == 'probe' else OUT))
    wall_cap = a.wall_cap or 2 * a.budget + 120
    T0 = time.time(); deadline = T0 + 3600.0 * a.deadline_h
    pop = populations()
    todo = probe_cells(pop) if a.mode == 'probe' else sweep_order(pop)
    if a.only:
        want = {(x.split(':')[0], int(x.split(':')[1])) for x in a.only.split(',')}; todo = [c for c in todo if c in want]
    if not a.no_resume:
        d = done_set(out); todo = [c for c in todo if c not in d]
    print(f"mode {a.mode}; populations (Specht-program-only, N >= 49): cxd {len(pop['cxd'])} "
          f"(49-150: {sum(n <= 150 for n in pop['cxd'])}, 151-200: {sum(n > 150 for n in pop['cxd'])}), cpd {len(pop['cpd'])}", flush=True)
    print(f'to run now: {len(todo)} cells; budget {a.budget:.0f}s CPU/cell, wall cap {wall_cap:.0f}s, workers {W}, '
          f'deadline {a.deadline_h} h; order {todo}', flush=True)
    soff = 0 if a.mode == 'probe' else 5000                        # sweep walks are independent of the probe walks
    jobs = [(t, n, a.budget, wall_cap, deadline, SEED + TOFF[t] + n + soff, out) for t, n in todo]
    nb = nd = ns = 0
    for t in geom.TABLES: os.makedirs(os.path.join(out, t), exist_ok=True)
    with Pool(W) as pool:
        for x in pool.imap_unordered(run_one, jobs):
            with open(os.path.join(out, x['table'], 'sweep_rows.jsonl'), 'a') as f: f.write(json.dumps(x) + '\n')
            if x.get('skipped'):
                ns += 1; print(f"{x['table']} N={x['n']:3d} SKIPPED (deadline)", flush=True); continue
            nd += 1; nb += x['beaten_float']
            print(f"{x['table']} N={x['n']:3d} pub {x['r_pub']} bar {x['bar_eff']:.15f} ({x['bar_eff_src']})  start {x['start']}  "
                  f"best {x['r_dense']:.15f} (vs printed {x['rel_vs_printed']:+.3e}, vs bar {x['rel_vs_bar']:+.3e}) "
                  f"{'BEATEN' if x['beaten_float'] else '-'}  hops {x['hops']} impr {x['improvements']} side {x['sideways']} "
                  f"stop {x['stop']} cpu {x['cpu']:.0f}s wall {x['wall']:.0f}s   [{nd} done, {nb} beaten, {ns} skipped, "
                  f"elapsed {time.time() - T0:.0f}s]", flush=True)
    print(f'{a.mode.upper()} END: {nd} run, {nb} beaten in float (claim rule vs bar_eff), {ns} skipped; wall {time.time() - T0:.0f}s',
          flush=True)
