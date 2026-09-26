# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""ssp PROBE (2026-09-26; Moki&Julio): equal spheres in the unit sphere, Packomania ssp. Float search only; exact
certification is the certification step (checker A = checkers/certify_ball.py d = 3, checker B = checkers/verify_exact_ball.py d = 3).
The bar and the verdict rule are stated in the probe plan BEFORE the run.

Recipe = ../hsp/sweep4.py Cell / job_warm, unchanged in substance (../hsp/slpd.polish, ../hsp/search.move_relocate / move_shake /
move_jiggle, imported from ../hsp/), with the ssp seed:
  warm  : Packomania's packing = the file's N - 1 centres + the reconstructed N-th centre (missing.py; every ssp file
          has row N blank) -> slpd.polish (pretest, cap 60 s) -> basin hopping from the best (relocate 1-3 least-contacted
          spheres into the largest holes | shake a region | jiggle all; each + penalty relax + polish (cap 15 s), kept if
          the TRUE float min-radius grows), --cpu seconds per cell, EXTENDED by --ext up to twice when the walk improved in
          the last 45 s of its budget.
  del   : (positive control N = 918; the page's N = 919 radius is ABOVE the printed 918 radius) Packomania's packing of the
          larger M whose radius is the bar (else M = N + 1) minus M - N spheres, four choices (the fewest-contact ones, random
          most-contacted ones, the innermost ones, random wall-contact ones), each polished (cap 20 s), + Packomania's own N;
          then basin hopping from the best.
Float pre-screen = the claim rule r > bar (1 + 1e-10) and r >= printed + 2e-12 (bar: data/refs/ssp/ssp_bar.json).
Outputs: out/<tag>.jsonl (one row per cell), out/best_<tag>/ssp_<N>.npy (best float packing, never overwritten by a
worse one). Every process is BelowNormal (python.exe launched directly; every worker calls prio.lower() and reports the
class it read back); OMP / BLAS threads = 1.

usage: python probe.py [--workers 4] [--cpu 240] [--ext 60] [--only 918,12] [--tag probe]
       [--cells N1,N2,...  (ad-hoc warm cells instead of the probe list)] [--del-cells N,...  (deletion from the larger N
       whose printed radius is the bar, or from N + 1)] [--deadline-at <unix time>]
"""
import os, sys, time, json, argparse
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sspgeo                                               # noqa: E402  (puts ../hsp/ on sys.path)
import prio                                                 # noqa: E402
prio.patch_geomd(); prio.lower()
import numpy as np
from multiprocessing import Pool
import geomd, slpd, seeds, search                           # noqa: E402  (hsp modules, unchanged)

BARS = os.path.join(sspgeo.REFS, 'ssp_bar.json')
SEED = 20260926
# (id, N, plan, group)
PROBE = [(1, 918, 'del', 'pos_control'), (2, 12, 'warm', 'neg_control'),
         (3, 150, 'warm', 'mid'), (4, 200, 'warm', 'mid'), (5, 250, 'warm', 'mid'), (6, 305, 'warm', 'mid'), (7, 350, 'warm', 'mid'),
         (8, 420, 'warm', 'tail'), (9, 457, 'warm', 'tail'), (10, 500, 'warm', 'tail'), (11, 543, 'warm', 'tail'),
         (12, 600, 'warm', 'tail'), (13, 650, 'warm', 'tail'), (14, 700, 'warm', 'tail'), (15, 750, 'warm', 'tail'),
         (16, 800, 'warm', 'tail'), (17, 850, 'warm', 'tail'), (18, 900, 'warm', 'tail'), (19, 950, 'warm', 'tail'),
         (20, 1000, 'warm', 'tail')]
_B = None


def bars():
    global _B
    if _B is None: _B = json.load(open(BARS))['cells']
    return _B


def beats(r, B):
    return bool(r > float(B['bar']) * (1 + 1e-10) and r >= float(B['printed']) + 2e-12)


def best_dir(tag):
    return os.path.join(sspgeo.OUT, f'best_{tag}')


def save_best(tag, n, c):
    d = best_dir(tag); os.makedirs(d, exist_ok=True); p = os.path.join(d, f'ssp_{n}.npy'); r = geomd.rmin(c)
    if os.path.exists(p):
        old = np.load(p)
        if old.shape == c.shape and geomd.rmin(old) >= r: return False
    tmp = p + '.tmp.npy'; np.save(tmp, c); os.replace(tmp, p); return True


class Cell:
    """Budgeted search state for one cell (../hsp/sweep4.Cell with the ssp bar)."""
    def __init__(self, n, cpu, ext, deadline, seed):
        self.n, self.B = n, bars()[str(n)]
        self.pr, self.bar = float(self.B['printed']), float(self.B['bar'])
        self.rng = np.random.default_rng(seed); self.T0 = time.process_time(); self.W0 = time.time()
        self.cpu, self.ext, self.deadline = cpu, ext, deadline; self.next = 0
        self.tried = []; self.best = [None, -1.0, None]

    def left(self):
        return min(self.cpu - (time.process_time() - self.T0), self.deadline - time.time())

    def polish(self, lab, c, cap, pretest=True, lp_time=30.0):
        if c is None or self.left() <= 1: return None
        t = time.process_time(); r0 = geomd.rmin(c); st = {}
        q, r = slpd.polish(c, t_cap=max(1.0, min(cap, self.left())), stats=st, pretest=pretest, lp_time=lp_time,
                           seed=int(self.rng.integers(1 << 30)))
        self.tried.append(dict(seed=lab, r_seed=r0, r=r, gain_vs_printed=r - self.pr, rel_vs_bar=r / self.bar - 1,
                               cpu=round(time.process_time() - t, 1), stop=st.get('stop'), lp=st.get('lp'), kept=st.get('kept')))
        if r > self.best[1]: self.best[:] = [q, r, lab]
        return q

    def bh(self, extend=True, hop_cap=15.0):
        if self.best[0] is None: return
        rng = self.rng; cur_c, cur_r = self.best[0].copy(), self.best[1]; r0 = cur_r
        hops = impr = side = 0; kinds = {}; last_impr = -1e9
        while True:
            if self.left() <= 6:
                used = time.process_time() - self.T0
                if (extend and self.next < 2 and used - last_impr <= 45.0 and self.deadline - time.time() > self.ext + 10):
                    self.cpu += self.ext; self.next += 1; continue
                break
            u = rng.random()
            if u < 0.4: q, kind = search.move_relocate(cur_c, cur_r, rng)
            elif u < 0.8: q, kind = search.move_shake(cur_c, cur_r, rng)
            else: q, kind = search.move_jiggle(cur_c, cur_r, rng)
            q, rq = slpd.polish(q, t_cap=max(1.0, min(hop_cap, self.left() - 3)), lp_time=10.0, seed=int(rng.integers(1 << 30)))
            hops += 1; kinds.setdefault(kind, [0, 0]); kinds[kind][0] += 1
            if rq > self.best[1] * (1 + 1e-13):
                self.best[:] = [q.copy(), rq, self.best[2].split('+bh')[0] + '+bh']; cur_c, cur_r = q, rq; impr += 1
                kinds[kind][1] += 1; last_impr = time.process_time() - self.T0
            elif rq >= self.best[1] * (1 - 1e-12):
                cur_c, cur_r = q, rq; side += 1
        self.tried.append(dict(seed='basin_hopping', r_seed=r0, r=self.best[1], hops=hops, improvements=impr, sideways=side,
                               kinds=kinds, extensions=self.next))

    def row(self, tag, extra=None):
        base = dict(N=self.n, printed=self.B['printed'], bar=self.B['bar'], bar_source=self.B['bar_source'], ref=self.B['ref'])
        if self.best[0] is None: return dict(base, error='no seed', tried=self.tried)
        c = self.best[0]; r = geomd.float_r(c); saved = save_best(tag, self.n, c)
        base.update(best_seed=self.best[2], r_best=r, gain_abs_vs_printed=r - self.pr, gain_rel_vs_bar=r / self.bar - 1,
                    beats_bar_float=beats(r, self.B), saved=saved, cpu=round(time.process_time() - self.T0, 1),
                    wall=round(time.time() - self.W0, 1), tried=self.tried)
        if extra: base.update(extra)
        return base


def del_seeds(c, r, rng, m=1):
    """Four m-sphere deletions of a packing c (radius r): the m fewest-contact spheres, m random most-contacted ones, the m
    innermost, m random wall-contact ones."""
    k = geomd.contacts(c, r); nr = geomd.norms(c); n = len(c); out = []
    wall = np.nonzero((1.0 - nr) <= r * (1 + 1e-7))[0]
    top = np.argsort(-(k + rng.random(n) * 0.5))[:max(m, 8)]
    picks = [('fewest', np.argsort(k + rng.random(n) * 0.5)[:m]),
             ('mostcontacted', rng.choice(top, size=m, replace=False)),
             ('innermost', np.argsort(nr)[:m]),
             ('wall', rng.choice(wall, size=m, replace=False) if len(wall) >= m else np.argsort(-nr)[:m])]
    for lab, I in picks:
        I = [int(i) for i in I]
        out.append((f"{lab}#{'+'.join(map(str, I))}(k={'+'.join(str(int(k[i])) for i in I)})", np.delete(c, I, axis=0)))
    return out


def job(args):
    pid, n, plan, group, cpu, ext, deadline, seed, tag = args
    pc = prio.lower()
    if time.time() >= deadline - 20: return dict(id=pid, N=n, skipped=True)
    C = Cell(n, cpu, ext, deadline, seed)
    if plan == 'del':
        M = C.B['bar_M'] if C.B['bar_source'] == 'larger_N' and C.B['bar_M'] else n + 1
        big = sspgeo.load_full(M)
        if big is not None:
            for lab, c in del_seeds(big, geomd.rmin(big), C.rng, M - n):
                C.polish(f'page{M}-del{M - n}:{lab}', c, 20.0)
    own = sspgeo.load_full(n)
    if own is None: return dict(id=pid, N=n, error='no reconstructed packing (missing.py)')
    C.polish('packomania', own, 60.0)
    C.bh()
    return C.row(tag, dict(id=pid, plan=plan, group=group, prio=pc, r_packomania_full=geomd.rmin(own)))


def main():
    pc = prio.lower()
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4); ap.add_argument('--cpu', type=float, default=240.0)
    ap.add_argument('--ext', type=float, default=60.0); ap.add_argument('--only', default='')
    ap.add_argument('--cells', default=''); ap.add_argument('--tag', default='probe'); ap.add_argument('--del-cells', default='')
    ap.add_argument('--deadline-at', type=float, default=0.0); ap.add_argument('--hours', type=float, default=2.0)
    a = ap.parse_args(); W = min(a.workers, 4)
    T0 = time.time(); deadline = a.deadline_at if a.deadline_at > 0 else T0 + 3600 * a.hours
    if a.cells:
        cells = [(i + 1, int(x), 'warm', 'adhoc') for i, x in enumerate(a.cells.split(','))]
        if a.del_cells: cells += [(100 + i, int(x), 'del', 'nonmono') for i, x in enumerate(a.del_cells.split(','))]
    elif a.del_cells:
        cells = [(100 + i, int(x), 'del', 'nonmono') for i, x in enumerate(a.del_cells.split(','))]
    else:
        only = {int(x) for x in a.only.split(',')} if a.only else None
        cells = [x for x in PROBE if only is None or x[1] in only]
    cells.sort(key=lambda x: -x[1])
    os.makedirs(sspgeo.OUT, exist_ok=True); rows = os.path.join(sspgeo.OUT, f'{a.tag}.jsonl')
    print(f'ssp probe start {time.strftime("%Y-%m-%d %H:%M:%S")}  deadline {time.strftime("%H:%M:%S", time.localtime(deadline))}  '
          f'workers {W}  cpu/cell {a.cpu:.0f}+2x{a.ext:.0f}s  main prio {pc}  cells {[x[1] for x in cells]}', flush=True)
    nb = 0
    with Pool(W, maxtasksperchild=1) as pool:
        jobs = [(i, n, plan, g, a.cpu, a.ext, deadline, SEED + n, a.tag) for i, n, plan, g in cells]
        for x in pool.imap_unordered(job, jobs):
            with open(rows, 'a') as f: f.write(json.dumps(x) + '\n')
            if x.get('skipped') or 'error' in x:
                print(f"#{x['id']} N={x['N']} {'SKIPPED' if x.get('skipped') else 'ERROR ' + x['error']}", flush=True); continue
            nb += x['beats_bar_float']; bh = next((t for t in x['tried'] if t['seed'] == 'basin_hopping'), {})
            print(f"#{x['id']:2d} N={x['N']:4d} {x['group']:11s} best {x['best_seed'][:34]:34s} r {x['r_best']:.15f} vs printed "
                  f"{x['gain_abs_vs_printed']:+.3e} abs  vs bar({x['bar_source']}) {x['gain_rel_vs_bar']:+.3e} rel  "
                  f"{'BEATS' if x['beats_bar_float'] else '-'}  cpu {x['cpu']}s hops {bh.get('hops', 0)} impr {bh.get('improvements', 0)} "
                  f"ext {bh.get('extensions', 0)} prio {x.get('prio')} [{time.time() - T0:.0f}s]", flush=True)
            for t in x['tried']:
                print('      ' + json.dumps(t), flush=True)
    print(f'PROBE END {time.strftime("%H:%M:%S")}: {nb}/{len(cells)} beat the bar in float; wall {time.time() - T0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
