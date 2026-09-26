# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""hsp PROBE (2026-09-25; Moki&Julio): equal balls in the unit d-ball, Packomania hsp4 / hsp5 / hsp6. Float search only;
exact certification is the certification step (checker A = checkers/certify_ball.py, checker B = checkers/verify_exact_ball.py).

The 20 cells of recon 7 ("Proposed 20-cell probe", the recon of the hsp tables), each with ~3 min CPU:
   1-7  hsp4 87 119 154 242 259 285 299   warm: Packomania's coordinates -> polish -> basin hopping  (rigidity deficit)
   8-11 hsp4 70 149 175 213               warm, the same                                            (rigid: proxy control)
   12   hsp4 200                          cold: code shell of N - k + k inner (k near the page's core count)
   13   hsp4 30                           code 29 + centre -> polish     POSITIVE control: must hold >= +6.05e-6 over the page
   14   hsp4 25                           24-cell + centre -> polish     NEGATIVE control: must show no gain
   15   hsp5 42                           code 41 + centre -> polish     (does polish add beyond the code seed?)
   16-17 hsp5 160 250                     cold: code shell + inner cluster
   18   hsp6 200                          code 199 + centre -> polish
   19   hsp6 155                          cold
   20   hsp6 66                           E6 72-point kissing configuration minus 6 -> polish (+ the Cohn 66-point code)
Per cell: every seed of its plan is polished (slpd.polish, pretest), the best one gets the remaining CPU as basin hopping
(relocate the 1-3 least-contacted balls into the largest holes | shake a region | jiggle all; each + penalty relax + polish,
kept if the TRUE float min-radius grows). Warm cells also get two neighbour seeds from the table's own N - 1 (insert one
ball into the largest hole) and N + 1 (delete the fewest-contact ball), and our own packings at N +- 1..4 when we hold
them. Float pre-screen = the claim rule: r > bar (1 + 1e-10) and r >= printed + 2e-12 (bar: out/bars.json).
Outputs: out/probe.jsonl (one row per cell), out/best/hsp<d>_<N>.npy (best float packing, never overwritten by a
worse one; a run with --tag X other than 'probe' writes out/X.jsonl and out/best_X/). Every process runs BelowNormal (launch python.exe directly; workers lower themselves); OMP threads = 1.

usage: python search.py [--workers 4] [--cpu 180] [--only 13,14] [--tag probe]
"""
import os, sys, time, json, argparse
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import geomd, slpd, seeds, bars                                     # noqa: E402

BEST = os.path.join(geomd.OUT, 'best')
SEED = 20260925
PROBE = [(1, 4, 87, 'warm'), (2, 4, 119, 'warm'), (3, 4, 154, 'warm'), (4, 4, 242, 'warm'), (5, 4, 259, 'warm'),
         (6, 4, 285, 'warm'), (7, 4, 299, 'warm'), (8, 4, 70, 'warm'), (9, 4, 149, 'warm'), (10, 4, 175, 'warm'),
         (11, 4, 213, 'warm'), (12, 4, 200, 'cold'), (13, 4, 30, 'code_centre'), (14, 4, 25, 'cell24'),
         (15, 5, 42, 'code_centre'), (16, 5, 160, 'cold'), (17, 5, 250, 'cold'), (18, 6, 200, 'code_centre'),
         (19, 6, 155, 'cold'), (20, 6, 66, 'e6')]


def beats(r, B):
    return bool(r > float(B['bar']) * (1 + 1e-10) and r >= float(B['printed']) + 2e-12)


def set_tag(tag):
    global BEST
    BEST = os.path.join(geomd.OUT, 'best' if tag == 'probe' else f'best_{tag}')


def best_path(d, n):
    return os.path.join(BEST, f'hsp{d}_{n}.npy')


def save_best(d, n, c):
    os.makedirs(BEST, exist_ok=True); p = best_path(d, n); r = geomd.rmin(c)
    if os.path.exists(p):
        old = np.load(p)
        if old.shape == c.shape and geomd.rmin(old) >= r: return False
    tmp = p + '.tmp.npy'; np.save(tmp, c); os.replace(tmp, p); return True


def ours(d, n):
    p = best_path(d, n)
    if os.path.exists(p):
        c = np.load(p)
        if c.shape == (n, d): return c
    return None


# ---------------------------------------------------------------- basin-hopping moves (d-D)
def move_relocate(c, r, rng):
    n = len(c); k = int(rng.choice([1, 2, 3], p=[0.5, 0.3, 0.2]))
    key = geomd.contacts(c, r) + rng.random(n) * 0.5
    pick = rng.choice(np.argsort(key)[:2 * k], size=k, replace=False); keep = np.setdiff1d(np.arange(n), pick)
    h = seeds.find_holes(c[keep], r, k, rng)
    q = c.copy(); q[pick] = h + rng.uniform(-1e-3 * r, 1e-3 * r, h.shape)
    return seeds.relax(q, r * 1.0005, 300)[0], f'reloc{k}'


def move_shake(c, r, rng):
    n, d = c.shape; i0 = int(rng.integers(n)); R = r * rng.uniform(2.5, 5.0); s = r * rng.uniform(0.1, 0.5)
    idx = np.nonzero(geomd.norms(c - c[i0]) <= R)[0]
    v = geomd.random_sphere(len(idx), d, rng)
    q = c.copy(); q[idx] += v * (s * rng.random(len(idx)) ** (1 / d))[:, None]
    return seeds.relax(q, r * 1.0005, 300)[0], 'shake'


def move_jiggle(c, r, rng):
    a = r * 10 ** rng.uniform(-5, -3)
    return c + rng.uniform(-a, a, c.shape), 'jiggle'


# ---------------------------------------------------------------- one cell
def job(args):
    pid, d, n, plan, cpu, seed, tag = args
    geomd.lower_priority(); set_tag(tag)
    rng = np.random.default_rng(seed); T0 = time.process_time(); W0 = time.time()
    B = bars.bar(d, n); pr, bar = float(B['printed']), float(B['bar']); r_ref = max(pr, bar)
    left = lambda: cpu - (time.process_time() - T0)
    tried = []; best = [None, -1.0, None]

    def polish(lab, c, cap, pretest=True):
        if c is None or left() <= 1: return None
        t = time.process_time(); r0 = geomd.rmin(c); st = {}
        q, r = slpd.polish(c, t_cap=max(1.0, min(cap, left())), stats=st, pretest=pretest, seed=int(rng.integers(1 << 30)))
        tried.append(dict(seed=lab, r_seed=r0, r=r, gain_vs_printed=r - pr, rel_vs_bar=r / bar - 1,
                          cpu=round(time.process_time() - t, 1), stop=st.get('stop'), lp=st.get('lp'), soc=st.get('soc')))
        if r > best[1]: best[:] = [q, r, lab]
        return q

    # ---- seeds by plan
    if plan == 'warm':
        polish('packomania', seeds.packomania(d, n), 0.35 * cpu)
        donors = []
        if n - 1 >= 2:
            try: donors.append((f'insert1_into_page{n - 1}', seeds.insert(seeds.packomania(d, n - 1), float(bars.bar(d, n - 1)['printed']), 1, rng)))
            except Exception: pass
        if n + 1 <= geomd.NMAX[d]:
            try:
                cd = seeds.packomania(d, n + 1); donors.append((f'delete1_from_page{n + 1}', seeds.delete(cd, geomd.rmin(cd), 1, rng)))
            except Exception: pass
        for k in range(1, 5):                                     # our own packings at N +- k
            for m, how in ((n - k, 'insert'), (n + k, 'delete')):
                co = ours(d, m) if 2 <= m <= geomd.NMAX[d] else None
                if co is not None and geomd.rmin(co) > best[1]:
                    rc = geomd.rmin(co)
                    donors.append((f'{how}{k}_ours{m}', seeds.insert(co, rc, k, rng) if how == 'insert' else seeds.delete(co, rc, k, rng)))
        for lab, c in donors: polish(lab, c, 0.1 * cpu)
    elif plan == 'cold':
        T = geomd.table(d); core = T[n]['core']
        ks = sorted({max(0, core + dk) for dk in (0, -2, 2)})
        cands = [(f'shell{n - k}+core{k}:code', k, 'code') for k in ks] + [(f'shell{n - core}+core{core}:random', core, 'random')]
        for lab, k, inner in cands:
            if left() < 0.5 * cpu: break
            polish(lab, seeds.shell_core(d, n, k, r_ref, rng, inner), 0.06 * cpu)
        if left() > 0.5 * cpu: polish('cold_random', seeds.cold(d, n, r_ref, rng), 0.06 * cpu)
        if best[0] is not None: polish(best[2] + '+long', best[0], 0.2 * cpu, pretest=False)
    elif plan == 'code_centre':
        cs = seeds.code_seed(d, n - 1, centre=True)
        polish(f'code{n - 1}+centre', cs[0] if cs else None, 0.3 * cpu)
        ca = seeds.code_seed(d, n, centre=False)
        if ca: polish(f'code{n}', ca[0], 0.15 * cpu)
    elif plan == 'cell24':
        polish('24cell+centre', seeds.place_code(seeds.cell24(), True)[0], 0.3 * cpu)
    elif plan == 'e6':
        E = seeds.e6_roots(); U0 = E
        for t in range(4):
            if left() < 0.5 * cpu: break
            keep = np.sort(rng.permutation(72)[:n]) if t else np.sort(np.argsort(E @ E[0])[6:])   # t = 0: a structured cut
            cE = seeds.place_code(E[keep])[0]
            if t: cE = move_jiggle(cE, 1.0 / 3.0, rng)[0]         # random cuts: + a 1e-5..1e-3 r jiggle (the cut is stationary)
            polish(f'E6-{72 - n}' + ('' if t else ':cap') + (f':rand{t}+jiggle' if t else ''), cE, 0.1 * cpu)
        ca = seeds.code_seed(d, n)
        if ca: polish(f'code{n}', ca[0], 0.1 * cpu)
    # ---- basin hopping on the best
    hops = impr = side = 0; kinds = {}
    if best[0] is not None:
        cur_c, cur_r = best[0].copy(), best[1]; bh_r0 = best[1]
        while left() > 8:
            u = rng.random()
            if u < 0.4: q, kind = move_relocate(cur_c, cur_r, rng)
            elif u < 0.8: q, kind = move_shake(cur_c, cur_r, rng)
            else: q, kind = move_jiggle(cur_c, cur_r, rng)
            q, rq = slpd.polish(q, t_cap=max(1.0, min(15.0, left() - 4)), lp_time=10.0, seed=int(rng.integers(1 << 30)))
            hops += 1; kinds.setdefault(kind, [0, 0]); kinds[kind][0] += 1
            if rq > best[1] * (1 + 1e-13):
                best[:] = [q.copy(), rq, best[2].split('+bh')[0] + '+bh']; cur_c, cur_r = q, rq; impr += 1; kinds[kind][1] += 1
            elif rq >= best[1] * (1 - 1e-12):
                cur_c, cur_r = q, rq; side += 1
        tried.append(dict(seed='basin_hopping', r_seed=bh_r0, r=best[1], hops=hops, improvements=impr, sideways=side, kinds=kinds))
    if best[0] is None:
        return dict(id=pid, d=d, N=n, plan=plan, error='no seed', tried=tried)
    c = best[0]; r = geomd.float_r(c); saved = save_best(d, n, c)
    return dict(id=pid, d=d, N=n, plan=plan, tag=tag, printed=B['printed'], bar=B['bar'], bar_source=B['bar_source'],
                ref=B['ref'], best_seed=best[2], r_best=r, gain_abs_vs_printed=r - pr, gain_rel_vs_printed=r / pr - 1,
                gain_rel_vs_bar=r / bar - 1, beats_bar_float=beats(r, B), saved=saved, cpu=round(time.process_time() - T0, 1),
                wall=round(time.time() - W0, 1), tried=tried)


def main():
    geomd.lower_priority()
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4); ap.add_argument('--cpu', type=float, default=180.0)
    ap.add_argument('--only', default=''); ap.add_argument('--tag', default='probe')
    a = ap.parse_args(); W = min(a.workers, 4); set_tag(a.tag)
    only = {int(x) for x in a.only.split(',')} if a.only else None
    cells = [x for x in PROBE if only is None or x[0] in only]
    cells.sort(key=lambda x: -(x[2] * x[1]))                       # big cells first (better packing of the pool)
    os.makedirs(geomd.OUT, exist_ok=True); rows = os.path.join(geomd.OUT, f'{a.tag}.jsonl')
    print(f'hsp probe start {time.strftime("%Y-%m-%d %H:%M:%S")}  workers {W}  cpu/cell {a.cpu:.0f}s  cells {[x[0] for x in cells]}', flush=True)
    T0 = time.time(); nb = 0
    with Pool(W) as pool:
        for x in pool.imap_unordered(job, [(i, d, n, plan, a.cpu, SEED + 1000 * d + n, a.tag) for i, d, n, plan in cells]):
            with open(rows, 'a') as f: f.write(json.dumps(x) + '\n')
            if 'error' in x:
                print(f"#{x['id']} hsp{x['d']} N={x['N']} ERROR {x['error']}", flush=True); continue
            nb += x['beats_bar_float']
            print(f"#{x['id']:2d} hsp{x['d']} N={x['N']:3d} {x['plan']:11s} best {x['best_seed']:34s} r {x['r_best']:.15f}  vs printed "
                  f"{x['gain_abs_vs_printed']:+.3e} abs  vs bar({x['bar_source']}) {x['gain_rel_vs_bar']:+.3e} rel  "
                  f"{'BEATS' if x['beats_bar_float'] else '-'}  cpu {x['cpu']}s wall {x['wall']}s  [{time.time() - T0:.0f}s]", flush=True)
            for t in x['tried']:
                print('      ' + json.dumps(t), flush=True)
    print(f'PROBE END {time.strftime("%H:%M:%S")}: {nb}/{len(cells)} beat the bar in float; wall {time.time() - T0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
