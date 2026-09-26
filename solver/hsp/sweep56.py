# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""hsp5 / hsp6 SWEEP from Packomania's RECOVERED packings (2026-09-26; Moki&Julio). Float search only; exact certification is
the certification step (checker A = checkers/certify_ball.py, checker B = checkers/verify_exact_ball.py), as for hsp4.

Packomania's hsp5 / hsp6 files are truncated to 4 columns; recon56.py completes them (the missing coordinates are pinned by
the wall and the contacts) and the completed packing, polished, reaches the printed radius (out/recon56/test.jsonl). Per cell
(the hsp4 recipe of sweep4.py P1, whose Cell class, basin hopping and row format are reused unchanged): seeds = the recovered
packing (recon56.complete -> slpd.polish, cap 60 s), the code construction when the bar is code-derived, our own earlier packing
for the cell if any (out/sweep4/best: the P3 code arm); then basin hopping (relocate / shake / jiggle, each + polish) with the
--cpu / --ext budget. Cells in DESCENDING N (the largest gains sat at large N in hsp4). Float pre-screen = the claim rule
r > bar (1 + 1e-10) and r >= printed + 2e-12 (bar: out/bars.json, which includes Cohn's codes and larger N).
Rows: out/sweep56/rows.jsonl; best float packings out/sweep56/best/hsp<d>_<N>.npy (never overwritten by a worse one).
Priority: launch python.exe directly; every worker lowers itself (prio.lower()).

usage: python sweep56.py [--tables 5,6] [--nmin 20] [--workers 4] [--cpu 150] [--ext 60]
       [--deadline-at <unix time>] [--only 5:100,6:200]
"""
import os, sys, time, json, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('HSP_SWEEP_OUT', os.path.join(HERE, 'out', 'sweep56'))   # 6-D runs use HSP_SWEEP_OUT=out/sweep6
import numpy as np                                                   # noqa: E402
import prio                                                          # noqa: E402
import geomd, seeds, sweep4, recon56                                 # noqa: E402
from multiprocessing import Pool                                     # noqa: E402
assert sweep4.OUT == os.environ['HSP_SWEEP_OUT'], "sweep4 must read HSP_SWEEP_OUT at import"
SW4_BEST = os.path.join(geomd.OUT, 'sweep4', 'best')


def job_rec(args):
    phase, d, n, cpu, ext, deadline, seed = args
    pc = prio.lower()
    if time.time() >= deadline - 20: return dict(phase=phase, d=d, N=n, skipped=True)
    C = sweep4.Cell(d, n, cpu, ext, deadline, seed)
    r_pr, A = recon56.trunc(d, n)
    rp = os.path.join(geomd.OUT, 'recon56', f'hsp{d}_{n}.npy')          # 6-D: the recovered + polished packing (recover6.py)
    if d == 6 and os.path.exists(rp) and np.load(rp).shape == (n, d): comp, viol, log = np.load(rp), 0.0, [dict(method='recover6')]
    else: comp, viol, log = recon56.complete(d, A, r_pr, C.rng, starts=8)
    C.polish('packomania_recovered', comp, 60.0)
    if C.B['bar_source'] in ('code', 'code+centre') and C.B['bar_M']:
        cs = seeds.code_seed(d, C.B['bar_M'], centre=C.B['bar_source'] == 'code+centre')
        if cs is not None:
            c = cs[0]
            if len(c) > n: c = seeds.delete(c, cs[1], len(c) - n, C.rng)
            if len(c) == n: C.polish(f"code{C.B['bar_M']}" + ('+centre' if C.B['bar_source'] == 'code+centre' else ''), c, 20.0)
    for lab, p in (('ours:sweep56', os.path.join(sweep4.BEST, f'hsp{d}_{n}.npy')), ('ours:sweep4', os.path.join(SW4_BEST, f'hsp{d}_{n}.npy')),
                   ('ours:first', os.path.join(geomd.OUT, 'sweep56', 'best', f'hsp{d}_{n}.npy')),     # a second pass also starts from the first
                   ('ours:first6', os.path.join(geomd.OUT, 'sweep6', 'best', f'hsp{d}_{n}.npy')),
                   ('ours:second5', os.path.join(geomd.OUT, 'sweep5b', 'best', f'hsp{d}_{n}.npy'))):
        if os.path.exists(p):
            c = np.load(p)
            if c.shape == (n, d) and geomd.rmin(c) > C.best[1]: C.polish(lab, c, 20.0)
    C.bh()
    return C.row(phase, dict(prio=pc, recon_viol=viol, recon_starts=len(log)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tables', default='5,6'); ap.add_argument('--nmin', type=int, default=20); ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--cpu', type=float, default=150.0); ap.add_argument('--ext', type=float, default=60.0)
    ap.add_argument('--hours', type=float, default=3.5); ap.add_argument('--deadline-at', type=float, default=None); ap.add_argument('--only', default='')
    ap.add_argument('--seed0', type=int, default=20260926)                  # a second pass: a new seed stream
    a = ap.parse_args(); prio.lower()
    deadline = a.deadline_at or time.time() + 3600 * a.hours
    os.makedirs(sweep4.OUT, exist_ok=True)
    if a.only: cells = [tuple(map(int, x.split(':'))) for x in a.only.split(',')]
    else:
        cells = [(d, n) for d in map(int, a.tables.split(',')) for n in range(a.nmin, geomd.NMAX[d] + 1)
                 if os.path.exists(geomd.pub_path(d, n))]
        cells.sort(key=lambda t: (-t[1], t[0]))                      # descending N, hsp5 before hsp6 at equal N
    done = sweep4.done_pairs(); todo = [(d, n) for d, n in cells if ('R', d, n) not in done]
    T0 = time.time(); tally = dict(done=0, beats=0, skipped=0, no_donor=0, error=0)
    print(f"sweep56 start {time.strftime('%Y-%m-%d %H:%M:%S')}  deadline {time.strftime('%H:%M:%S', time.localtime(deadline))}  workers {a.workers}  "
          f"cells {len(cells)} (todo {len(todo)})", flush=True)
    with Pool(a.workers) as pool:
        sweep4.run_pool(pool, job_rec, [('R', d, n, a.cpu, a.ext, deadline, a.seed0 + 1000 * d + n) for d, n in todo], T0, tally)
    print(f"SWEEP END {time.strftime('%H:%M:%S')} {tally}", flush=True)


if __name__ == '__main__':
    main()
