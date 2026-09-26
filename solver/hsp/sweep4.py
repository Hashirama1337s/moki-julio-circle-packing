# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""hsp4 SWEEP N 56-300 + the hsp5 / hsp6 code-seed arm (2026-09-25; Moki&Julio). Float search only; exact certification
is the certification step (checker A = checkers/certify_ball.py, checker B = checkers/verify_exact_ball.py).

Triggered by the probe verdict (out/verdict_probe.json: controls OK, hsp4 9/11 wins -> sweep N 56-300).
Phases (one pool of <= 4 BelowNormal workers, one global wall deadline, --resume from the jsonl rows):
  P1 warm  : every hsp4 N in 56..300 except the 9 probe wins (87 119 149 175 213 242 259 285 299; their certified
             packings are kept; N = 119 is certified from the better dry-run packing, no re-run) = 236 cells, in
             DESCENDING N (the probe's largest gains and rigidity deficits are at large N, so a deadline cut loses the
             least). Seeds: Packomania's own coordinates -> slpd.polish (pretest, cap 60 s); a code seed where the bar is
             code-derived (none has a local code file in 56..300: 120 / 121 are the 600-cell cells, warm from the page);
             our own earlier packing for N if we hold one. Then basin hopping from the best (search.move_relocate /
             move_shake / move_jiggle, each + polish, kept if the TRUE float min-radius grows). CPU budget per cell
             --cpu (150 s), EXTENDED by --ext (60 s) up to twice when the walk improved in the last 45 s of its budget.
             (The probe's page-neighbour seeds -- insert into page N-1 / delete from page N+1 -- never came within 8e-6
             of the page in 11 cells, so they are dropped here.)
  P2 nbr   : hsp4 cells not beaten in float after P1, ASCENDING N: neighbour seeds from OUR packings at N +- 1..4
             (nearest donor below: insert N - M balls into its largest holes; nearest donor above with a radius above
             our best for N: delete the M - N fewest-contact balls), each polished (cap 25 s), then 30 s of basin hopping
             from the best if a seed improved it. Wall sub-deadline --p2-min minutes after P2 starts.
  P3 codes : hsp5 / hsp6 code-derived cells (data/big/hsp_codes/code_derived_beats_2026-09-25.json, 173 cells),
             smallest N first: the code seed as recorded (M points, + centre, minus the fewest-contact deleted points)
             -> polish -> basin hopping, --cpu3 (90 s) per cell. (This is what won hsp6 N = 66 in the probe.)
Rows: out/sweep4/rows.jsonl (env HSP_SWEEP_OUT overrides the directory, for tests) (one per (phase, d, N)); best float packings out/sweep4/best/hsp<d>_<N>.npy
(never overwritten by a worse one). Float pre-screen = the claim rule r > bar (1 + 1e-10) and r >= printed + 2e-12.
Priority: launch python.exe directly and set BelowNormal on it; every worker lowers itself with prio.lower() (the
geomd / ctypes version is a no-op on 64-bit Python) and reports the class it read back (row field 'prio').

usage: python sweep4.py [--workers 4] [--hours 3.5] [--cpu 150] [--ext 60] [--cpu3 90]
       [--phases 1,2,3] [--p2-min 20] [--only 200,201] [--no-resume] [--deadline-at <unix time>]
"""
import os, sys, time, json, argparse
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prio                                                          # noqa: E402
prio.patch_geomd()
import geomd, slpd, seeds, bars, search                              # noqa: E402

OUT = os.environ.get('HSP_SWEEP_OUT', os.path.join(geomd.OUT, 'sweep4')); BEST = os.path.join(OUT, 'best'); ROWS = os.path.join(OUT, 'rows.jsonl')
PROBE_WINS = (87, 119, 149, 175, 213, 242, 259, 285, 299)
SEED = 20260926
BEST_DIRS = (BEST, os.path.join(geomd.OUT, 'best'), os.path.join(geomd.OUT, 'best_dry'))


def ours(d, n):
    """Best float packing we hold for (d, N) over the sweep, the probe and the dry run -> (c, r, dir) or None."""
    got = []
    for dd in BEST_DIRS:
        p = os.path.join(dd, f'hsp{d}_{n}.npy')
        if os.path.exists(p):
            c = np.load(p)
            if c.shape == (n, d): got.append((c, geomd.rmin(c), os.path.basename(dd)))
    return max(got, key=lambda t: t[1]) if got else None


def save_best(d, n, c):
    os.makedirs(BEST, exist_ok=True); p = os.path.join(BEST, f'hsp{d}_{n}.npy'); r = geomd.rmin(c)
    if os.path.exists(p):
        old = np.load(p)
        if old.shape == c.shape and geomd.rmin(old) >= r: return False
    tmp = p + '.tmp.npy'; np.save(tmp, c); os.replace(tmp, p); return True


class Cell:
    """Budgeted search state for one cell: polish seeds, then basin hopping with the extension rule."""
    def __init__(self, d, n, cpu, ext, deadline, seed):
        self.d, self.n, self.B = d, n, bars.bar(d, n)
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
                               cpu=round(time.process_time() - t, 1), stop=st.get('stop'), lp=st.get('lp')))
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

    def row(self, phase, extra=None):
        base = dict(phase=phase, d=self.d, N=self.n, printed=self.B['printed'], bar=self.B['bar'], bar_source=self.B['bar_source'],
                    ref=self.B['ref'])
        if self.best[0] is None: return dict(base, error='no seed', tried=self.tried)
        c = self.best[0]; r = geomd.float_r(c); saved = save_best(self.d, self.n, c)
        base.update(best_seed=self.best[2], r_best=r, gain_abs_vs_printed=r - self.pr, gain_rel_vs_bar=r / self.bar - 1,
                    beats_bar_float=search.beats(r, self.B), saved=saved, cpu=round(time.process_time() - self.T0, 1),
                    wall=round(time.time() - self.W0, 1), tried=self.tried)
        if extra: base.update(extra)
        return base


def job_warm(args):
    phase, d, n, cpu, ext, deadline, seed = args
    pc = prio.lower()
    if time.time() >= deadline - 20: return dict(phase=phase, d=d, N=n, skipped=True)
    C = Cell(d, n, cpu, ext, deadline, seed)
    C.polish('packomania', seeds.packomania(d, n), 60.0)
    if C.B['bar_source'] in ('code', 'code+centre') and C.B['bar_M']:
        cs = seeds.code_seed(d, C.B['bar_M'], centre=C.B['bar_source'] == 'code+centre')
        if cs is not None:
            c = cs[0]
            if len(c) > n: c = seeds.delete(c, cs[1], len(c) - n, C.rng)
            if len(c) == n: C.polish(f"code{C.B['bar_M']}" + ('+centre' if C.B['bar_source'] == 'code+centre' else ''), c, 20.0)
    o = ours(d, n)
    if o is not None and o[1] > C.best[1]: C.polish(f'ours:{o[2]}', o[0], 20.0)
    C.bh()
    return C.row(phase, dict(prio=pc))


def job_nbr(args):
    phase, d, n, cpu, deadline, seed = args
    pc = prio.lower()
    if time.time() >= deadline - 20: return dict(phase=phase, d=d, N=n, skipped=True)
    C = Cell(d, n, cpu, 0, deadline, seed); o = ours(d, n)
    if o is None: return dict(phase=phase, d=d, N=n, error='no packing for N')
    C.best[:] = [o[0].copy(), o[1], f'ours:{o[2]}']; r_in = o[1]
    lo = hi = None
    for k in range(1, 5):
        if lo is None and n - k >= 2:
            x = ours(d, n - k)
            if x is not None and x[1] > r_in: lo = (k, n - k, x)
        if hi is None and n + k <= geomd.NMAX[d]:
            x = ours(d, n + k)
            if x is not None and x[1] > r_in: hi = (k, n + k, x)
    if lo is None and hi is None: return dict(phase=phase, d=d, N=n, no_donor=True, r_start=r_in)
    if lo:
        k, m, (c, r, _) = lo; C.polish(f'insert{k}_ours{m}', seeds.insert(c, r, k, C.rng), 25.0, pretest=False)
    if hi:
        k, m, (c, r, _) = hi; C.polish(f'delete{k}_ours{m}', seeds.delete(c, r, k, C.rng), 25.0, pretest=False)
    if C.best[1] > r_in: C.bh(extend=False)
    return C.row(phase, dict(prio=pc, r_start=r_in, donors=[x[:2] for x in (lo, hi) if x]))


def job_code(args):
    phase, d, n, M, centre, deleted, cpu, deadline, seed = args
    pc = prio.lower()
    if time.time() >= deadline - 20: return dict(phase=phase, d=d, N=n, skipped=True)
    C = Cell(d, n, cpu, 0, deadline, seed)
    cs = seeds.code_seed(d, M, centre=centre)
    if cs is None: return dict(phase=phase, d=d, N=n, error=f'no local code c{d}_{M}')
    c = cs[0]
    if len(c) > n: c = seeds.delete(c, cs[1], len(c) - n, C.rng)
    lab = f'code{M}' + ('+centre' if centre else '') + (f'-{deleted}' if deleted else '')
    C.polish(lab, c, 20.0)
    o = ours(d, n)
    if o is not None and o[1] > C.best[1]: C.polish(f'ours:{o[2]}', o[0], 10.0)
    C.bh(extend=False)
    return C.row(phase, dict(prio=pc))


def done_pairs():
    s = set()
    if os.path.exists(ROWS):
        for line in open(ROWS):
            try: x = json.loads(line)
            except Exception: continue
            if not x.get('skipped') and 'error' not in x: s.add((x['phase'], x['d'], x['N']))
    return s


def run_pool(pool, fn, jobs, T0, tally):
    for x in pool.imap_unordered(fn, jobs):
        with open(ROWS, 'a') as f: f.write(json.dumps(x) + '\n')
        if x.get('skipped'): tally['skipped'] += 1; continue
        if x.get('no_donor'): tally['no_donor'] += 1; continue
        if 'error' in x: tally['error'] += 1; print(f"{x['phase']} hsp{x['d']} N={x['N']} ERROR {x['error']}", flush=True); continue
        tally['done'] += 1; tally['beats'] += x['beats_bar_float']
        bh = next((t for t in x['tried'] if t['seed'] == 'basin_hopping'), {})
        print(f"{x['phase']} hsp{x['d']} N={x['N']:3d} best {x['best_seed']:26s} vs printed {x['gain_abs_vs_printed']:+.3e} abs  vs bar({x['bar_source'][:6]}) "
              f"{x['gain_rel_vs_bar']:+.3e} rel {'BEATS' if x['beats_bar_float'] else '-    '} cpu {x['cpu']:.0f}s hops {bh.get('hops', 0)} impr {bh.get('improvements', 0)} "
              f"ext {bh.get('extensions', 0)} prio {x.get('prio')}  [{tally['done']} done, {tally['beats']} beat, {time.time() - T0:.0f}s]", flush=True)


def main():
    pc = prio.lower()
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4); ap.add_argument('--hours', type=float, default=3.5)
    ap.add_argument('--cpu', type=float, default=150.0); ap.add_argument('--ext', type=float, default=60.0)
    ap.add_argument('--cpu3', type=float, default=90.0); ap.add_argument('--phases', default='1,2,3')
    ap.add_argument('--p2-min', type=float, default=20.0); ap.add_argument('--only', default='')
    ap.add_argument('--no-resume', action='store_true'); ap.add_argument('--deadline-at', type=float, default=0.0)
    a = ap.parse_args(); W = min(a.workers, 4); os.makedirs(BEST, exist_ok=True)
    T0 = time.time(); deadline = a.deadline_at if a.deadline_at > 0 else T0 + 3600 * a.hours
    phases = {int(p) for p in a.phases.split(',')}; only = {int(x) for x in a.only.split(',')} if a.only else None
    done = set() if a.no_resume else done_pairs()
    sel = lambda ph, d, L: [n for n in L if (ph, d, n) not in done and (only is None or n in only)]
    p1 = sel('P1', 4, [n for n in range(300, 55, -1) if n not in PROBE_WINS])
    print(f'hsp sweep start {time.strftime("%Y-%m-%d %H:%M:%S")}  deadline {time.strftime("%H:%M:%S", time.localtime(deadline))}  '
          f'workers {W}  cpu {a.cpu:.0f}+2x{a.ext:.0f}s  main prio {pc}  P1 cells {len(p1)}  resume: {len(done)} done', flush=True)
    tally = dict(done=0, beats=0, skipped=0, error=0, no_donor=0)
    with Pool(W) as pool:
        if 1 in phases:
            run_pool(pool, job_warm, [('P1', 4, n, a.cpu, a.ext, deadline, SEED + n) for n in p1], T0, tally)
        if 2 in phases and time.time() < deadline - 60:
            d2 = min(deadline, time.time() + 60 * a.p2_min); Bs = bars.load(); L = []
            for n in range(56, 301):
                if n in PROBE_WINS or ('P2', 4, n) in done or (only and n not in only): continue
                o = ours(4, n)
                if o is None or not search.beats(o[1], Bs['d4'][str(n)]): L.append(n)
            print(f'P2 neighbour seeds (ascending N, sub-deadline {time.strftime("%H:%M:%S", time.localtime(d2))}): {len(L)} hsp4 cells not beaten', flush=True)
            run_pool(pool, job_nbr, [('P2', 4, n, 80.0, d2, SEED + 7 * n) for n in L], T0, tally)
        if 3 in phases and time.time() < deadline - 60:
            J = json.load(open(os.path.join(geomd.CODES, 'code_derived_beats_2026-09-25.json'), encoding='utf-8'))
            J = sorted([x for x in J if x['table'] in ('hsp5', 'hsp6')], key=lambda x: (x['N'], x['table']))
            jobs = [('P3', int(x['table'][3]), x['N'], x['code_points'], x['centre'], x['deleted'], a.cpu3, deadline, SEED + 13 * x['N'] + int(x['table'][3]))
                    for x in J if ('P3', int(x['table'][3]), x['N']) not in done and (only is None or x['N'] in only)]
            print(f'P3 hsp5/hsp6 code seed + basin hopping (smallest N first): {len(jobs)} cells', flush=True)
            run_pool(pool, job_code, jobs, T0, tally)
    print(f'SWEEP END {time.strftime("%H:%M:%S")}: {tally}; wall {time.time() - T0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
