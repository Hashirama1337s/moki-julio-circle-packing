# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""ssp SWEEP pass (2026-09-26; Moki&Julio): the probe recipe at a short budget over a band of N, with the exact
certification done in the worker right after the search. Runs only after the probe verdict (out/verdict_probe.json)
says SWEEP for that band.

Per cell (probe.Cell, unchanged): Packomania's packing (file N - 1 centres + the reconstructed N-th, missing.py)
-> ../hsp/slpd.polish (pretest, cap --polish s) -> basin hopping (../hsp/search moves, hop cap 10 s) until --cpu s of CPU
(no extensions). If the float radius clears the claim rule, the certification step writes cert/ssp_<N>.txt and runs
checker A (checkers/certify_ball.py) and checker B (checkers/verify_exact_ball.py), each vs the bar and vs the printed radius;
WIN = the claim rule in exact decimals + all four IMPROVES.
Cells already WON in out/cand_ssp.json are skipped. One global wall deadline (--deadline-at); <= 4 BelowNormal workers.
Rows: out/<tag>.jsonl; best float packings: out/best_<tag>/; wins merged into out/cand_ssp.json at the end
(and after every 20 rows).

usage: python sweep.py --lo 401 --hi 1000 --workers 4 --cpu 40 --polish 15 --deadline-at <t> [--tag sweep_tail]
       [--order desc|asc|interleave] [--skip N,...]
CENSUS mode (--cpu 8 --polish 6: the polish only, no basin hopping; 1 worker, before the probe): which of Packomania's
packings are first-order stationary (slpd stop 'stationary*') and which have an improving direction (rigidity deficit).
"""
import os, sys, time, json, argparse
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sspgeo                                               # noqa: E402
import prio                                                 # noqa: E402
prio.patch_geomd(); prio.lower()
import numpy as np
from decimal import Decimal
from multiprocessing import Pool
import geomd                                                # noqa: E402
import probe, certify                                       # noqa: E402  (probe.py, the certification step)

SEED = 20260927


def job(args):
    n, cpu, pol, deadline, seed, tag = args
    pc = prio.lower()
    if time.time() >= deadline - 30: return dict(N=n, skipped=True)
    C = probe.Cell(n, cpu, 0, deadline - 15, seed)
    own = sspgeo.load_full(n)
    if own is None: return dict(N=n, error='no reconstructed packing')
    C.polish('packomania', own, pol)
    C.bh(extend=False, hop_cap=10.0)
    row = C.row(tag, dict(prio=pc, r_packomania_full=geomd.rmin(own)))
    if row.get('beats_bar_float'):
        c = C.best[0]; B = C.B
        try:
            row.update(certify.certify(n, c, B)); row['status'] = 'WIN' if row['win'] else 'certified_not_win'
        except Exception as e:
            row['status'] = 'cert_error'; row['cert_error'] = str(e)[:300]
    else:
        row['status'] = 'below_float_bar'
    return row


def merge(rows, tag):
    old = {c['N']: c for c in json.load(open(certify.CAND))['cells']} if os.path.exists(certify.CAND) else {}
    Bs = json.load(open(certify.BARS))['cells']
    for x in rows:
        if 'r_new' not in x: continue
        n = x['N']; B = Bs[str(n)]; cv = certify.code_value(B)
        r = dict(N=n, d=3, printed=B['printed'], bar=B['bar'], bar_source=B['bar_source'], ref=B['ref'], r_float=x['r_best'],
                 file=tag, route=f"{tag}:{x.get('best_seed')}", group='sweep', code_value=cv,
                 gain_float_rel_vs_bar=x['r_best'] / float(B['bar']) - 1, status=x['status'])
        for k in ('cert', 'r_new', 'checker_a_vs_bar', 'checker_a_vs_printed', 'checker_b_vs_bar', 'checker_b_vs_printed', 'claim_rule',
                  'win', 'gain_abs_vs_printed', 'gain_rel_vs_printed', 'gain_rel_vs_bar', 'secs_a', 'secs_b'):
            r[k] = x.get(k)
        r['pure_code'] = False; r['gain_rel_vs_code'] = (float(Decimal(x['r_new'])) / cv - 1) if cv else None
        o = old.get(n)
        if o and o.get('status') == 'WIN' and (r['status'] != 'WIN' or Decimal(o['r_new']) >= Decimal(r['r_new'])): continue
        old[n] = r
    json.dump(dict(table='ssp', built=time.strftime('%Y-%m-%d %H:%M:%S'),
                   rule='win = r_new > bar (1 + 1e-10) and r_new >= printed + 2e-12, checkers A and B IMPROVES vs bar and vs printed, '
                        'not a pure code construction', cells=[old[k] for k in sorted(old)]), open(certify.CAND, 'w'), indent=1)


def main():
    pc = prio.lower()
    ap = argparse.ArgumentParser()
    ap.add_argument('--lo', type=int, default=401); ap.add_argument('--hi', type=int, default=1000)
    ap.add_argument('--workers', type=int, default=4); ap.add_argument('--cpu', type=float, default=40.0)
    ap.add_argument('--polish', type=float, default=15.0); ap.add_argument('--deadline-at', type=float, required=True)
    ap.add_argument('--tag', default='sweep_tail'); ap.add_argument('--order', default='desc')
    ap.add_argument('--skip', default='')
    ap.add_argument('--seed0', type=int, default=SEED)                        # a second pass: a new seed stream
    a = ap.parse_args(); W = min(a.workers, 4)
    won = set()
    if os.path.exists(certify.CAND):
        won = {c['N'] for c in json.load(open(certify.CAND))['cells'] if c.get('status') == 'WIN'}
    rows_p = os.path.join(sspgeo.OUT, f'{a.tag}.jsonl'); done = set()
    if os.path.exists(rows_p):
        for l in open(rows_p):
            x = json.loads(l)
            if not x.get('skipped') and 'error' not in x: done.add(x['N'])
    L = [n for n in range(a.lo, a.hi + 1) if n not in won and n not in done]
    if a.skip: L = [n for n in L if n not in {int(x) for x in a.skip.split(',')}]
    if a.order == 'desc': L = L[::-1]
    elif a.order == 'interleave': L = sorted(L, key=lambda n: (n % 20, n))      # every band sampled evenly under a deadline
    print(f'ssp sweep {a.tag} start {time.strftime("%H:%M:%S")} deadline {time.strftime("%H:%M:%S", time.localtime(a.deadline_at))} '
          f'workers {W} cpu {a.cpu:.0f}s polish {a.polish:.0f}s prio {pc}: {len(L)} cells (skipping {len(won)} won, {len(done)} done)', flush=True)
    T0 = time.time(); acc = []; tally = dict(done=0, float_beats=0, wins=0, skipped=0, error=0)
    with Pool(W, maxtasksperchild=10) as pool:
        for x in pool.imap_unordered(job, [(n, a.cpu, a.polish, a.deadline_at, a.seed0 + n, a.tag) for n in L]):
            with open(rows_p, 'a') as f: f.write(json.dumps(x) + '\n')
            if x.get('skipped'): tally['skipped'] += 1; continue
            if 'error' in x: tally['error'] += 1; print(f"N={x['N']} ERROR {x['error']}", flush=True); continue
            tally['done'] += 1; tally['float_beats'] += bool(x.get('beats_bar_float')); tally['wins'] += x.get('status') == 'WIN'
            acc.append(x)
            print(f"N={x['N']:4d} {x['status']:17s} gain {x['gain_abs_vs_printed']:+.3e} abs {x['gain_rel_vs_bar']:+.3e} rel(bar) "
                  f"{x.get('checker_a_vs_bar', '')}/{x.get('checker_a_vs_printed', '')} {x.get('checker_b_vs_bar', '')}/{x.get('checker_b_vs_printed', '')} "
                  f"cpu {x['cpu']:.0f}s  [{tally['done']} done, {tally['wins']} wins, {time.time() - T0:.0f}s]", flush=True)
            if len(acc) % 20 == 0: merge(acc, a.tag)
    merge(acc, a.tag)
    print(f'SWEEP END {time.strftime("%H:%M:%S")}: {tally}; wall {time.time() - T0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
