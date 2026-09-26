# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""scu sweep (2026-09-25; Moki&Julio): equal spheres in the unit cube (Packomania scu). Float search only; exact
certification is checkers/certify_cube.py (checker A) + checkers/verify_exact_cube.py (checker B), exact.

Bars: data/refs/scu_bar.json (per N: the largest of Packomania's printed radius, Lai et al. 2023, roc-climate N = 203
  and the printed radius of any larger N). Float pre-screen used here and in the rows:
  beats_bar = r > bar (1 + 1e-10) and r >= printed + 2e-12   (the claim rule; the exact decision is the checkers').

Phases (one pool of <= 4 BelowNormal workers, global deadline, --resume from the jsonl rows):
  P1 polish : every NON-derivative cell N 201..1008 (405 free + 27 roots = 432): Specht-program free cells first, then the
              Cantrell / Pack'n'tile free cells, then the roots; the 20 cells of the probe last. Start = Packomania's own
              coordinates (507 and 619: the stale file, which reproduces an older radius). Repaired SLP polish
              (slp3d.polish, pretest=True: a first-order stationary start is detected in one small LP), cap 300 s.
  P2 polish : N 94..200 (107 cells) vs max(HTML, Lai 2023); start from the better (float min-radius) of Packomania's
              and Lai's coordinates (Lai rescaled by 1/L). Same polish.
  P3 bh     : basin hopping on the non-derivative cells whose polish did NOT beat the bar, largest stale blocks first
              (a block = a run of consecutive N with the same credit string). Moves: relocate 1-3 fewest-contact
              spheres into the largest holes | shake a random ball of spheres | jiggle every sphere by 1e-5..1e-3 r
              (escapes first-order stationary points that are not second-order optimal); each + polish (60 s cap).
              The walk keeps a CURRENT state (sideways moves within 1e-12 of the best are accepted) and a separate BEST
              state (updated only on a gain > 1e-13 relative); final polish; dense O(N^2) float recheck.
  PD delete : (phase id 5, runs after P1/P2 and before P3) for N = 1007 down to 94: if a packing we hold for some M in
              N+1..N+8 has a larger float radius than our best for N, delete the M - N fewest-contact spheres and polish
              (120 s cap). Covers derivative cells too.
  PI seeds  : (phase id 6, added 17:40 on the first run after P3 showed ~6 hops per 10 min at N ~ 760 and no gain):
              for every cell (P1, P2, P4 lists) we have not beaten yet, ascending N: take the best packing we hold for
              M in N-8..N-1 and INSERT N - M spheres into its largest holes, and the best for M in N+1..N+8 and DELETE
              M - N fewest-contact spheres (donor only if its radius exceeds our best for N); polish each (150 s cap).
              Rows are phase 'PN' (a first 17:31 attempt, rows 'PI', took the LARGEST-radius donor = the farthest one,
              8 inserts; hopeless, so the NEAREST qualifying donor on each side is used).
              Phase 6 = the non-derivative cells (before P3); phase 7 = the derivative cells (after P3, before P4).
  P4 polish : the derivative cells N 201..1008 (376), if time remains (Packomania's files for those N).
Outputs: out/sweep/best_<N>.npy (best float packing over all phases), out/sweep/rows.jsonl (one row per
(phase, N); --resume skips pairs already there), out/sweep.log.

usage: python search.py [--workers 4] [--deadline-h 4.5] [--cap 300] [--bh-wall 600]
       [--bh-cpu 300] [--phases 1,2,5,3,4] [--only 201,228] [--no-resume] [--deadline-at <unix time>]  (5 = PD, 6 = PI non-derivative, 7 = PI derivative)
Launch at BelowNormal (Start-Process + PriorityClass on python.exe itself); every worker also lowers itself.
"""
import os, sys, time, json, argparse
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common, slp3d                                       # noqa: E402

OUT = os.environ.get('SCU_SWEEP_OUT', os.path.join(HERE, 'out', 'sweep'))    # workers inherit --out via the env
ROWS = os.path.join(OUT, 'rows.jsonl')
PROBE_OUT = os.path.join(HERE, 'out')
SEED = 20260925
PROBE20 = [201, 228, 262, 291, 321, 347, 382, 426, 464, 483, 521, 571, 605, 628, 678, 738, 757, 799, 866, 885]


def bars():
    return json.load(open(os.path.join(common.REFS, 'scu_bar.json')))['cells']


def beats(r, e):
    return bool(r > float(e['bar']) * (1 + 1e-10) and r >= float(e['printed']) + 2e-12)


# ---------------------------------------------------------------- populations
def phase_lists(T=None):
    T = T or common.table(); der, roots, free = common.classes(T)
    spe = [n for n in range(201, 1009) if n in free and T[n][1] == '[31]' and n not in PROBE20]
    oth = [n for n in range(201, 1009) if n in free and T[n][1] != '[31]' and n not in PROBE20]
    rts = [n for n in range(201, 1009) if n in roots and n not in PROBE20]
    p1 = spe + oth + rts + [n for n in PROBE20 if n not in der]
    p2 = list(range(94, 201))
    p4 = [n for n in range(201, 1009) if n in der]
    return p1, p2, p4


def stale_blocks(T=None, lo=94, hi=1008):
    """{N: length of the run of consecutive N sharing N's credit string} over lo..hi."""
    T = T or common.table(); L = {}; n = lo
    while n <= hi:
        m = n
        while m + 1 <= hi and T[m + 1][1] == T[n][1]: m += 1
        for k in range(n, m + 1): L[k] = m - n + 1
        n = m + 1
    return L


# ---------------------------------------------------------------- start coordinates
def start_coords(n):
    """(coordinates, label): Packomania's own file; for N in 11..200 the better of it and Lai 2023 (float rmin)."""
    c = common.load_pub(n); lab = 'packomania_stale' if n in common.STALE else 'packomania'
    if 94 <= n <= 200 and os.path.exists(os.path.join(common.LAI, f'{n}_CubeSol.txt')):
        cl, _, _ = common.load_lai(n)
        if slp3d.rmin(cl) > slp3d.rmin(c): return cl, 'lai2023'
    return c, lab


def best_known(n):
    """Best float packing we already hold for N (sweep file, probe v2 file, probe v1 file) -> (c, r, source) or None."""
    got = []
    for p, src in ((os.path.join(OUT, f'best_{n}.npy'), 'sweep'), (os.path.join(PROBE_OUT, f'best_v3_{n}.npy'), 'probe_v3'),
                   (os.path.join(PROBE_OUT, f'best_v2_{n}.npy'), 'probe_v2'), (os.path.join(PROBE_OUT, f'best_{n}.npy'), 'probe_v1')):
        if os.path.exists(p):
            c = np.load(p)
            if c.shape == (n, 3): got.append((c, slp3d.rmin(c), src))
    return max(got, key=lambda t: t[1]) if got else None


def save_best(n, c):
    """Keep best_<N>.npy = the best float packing over all phases (never overwrite with a worse one)."""
    p = os.path.join(OUT, f'best_{n}.npy'); r = slp3d.rmin(c)
    if os.path.exists(p):
        old = np.load(p)
        if old.shape == c.shape and slp3d.rmin(old) >= r: return False
    tmp = p + '.tmp.npy'; np.save(tmp, c); os.replace(tmp, p); return True


# ---------------------------------------------------------------- basin-hopping moves (3-D)
def inside(p, r):
    return np.all(np.abs(p) <= 0.5 - r, axis=1)


def clearance(p, cs, r):
    """Largest radius a sphere centred at p could have given the centres cs (radius r) and the walls."""
    d, _ = cKDTree(cs).query(p, k=1)
    return np.minimum(d - r, (0.5 - np.abs(p)).min(1))


def find_holes(cs, r, k, rng, m=20000):
    """k greedy largest holes for the remaining centres cs; refine by a 3-D pattern search."""
    p = rng.uniform(-0.5 + r, 0.5 - r, (m, 3)); rho = clearance(p, cs, r); holes = []
    D = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]] +
                 [[a, b, c] for a in (-1, 1) for b in (-1, 1) for c in (-1, 1)], float)
    D /= np.linalg.norm(D, axis=1)[:, None]
    for _ in range(k):
        ref = np.vstack([cs] + ([np.array(holes)] if holes else []))
        top = np.argsort(-rho)[:6]; bestp, bestv = None, -np.inf
        for t in top:
            x, v, s = p[t].copy(), rho[t], 0.5 * r
            while s > 1e-4 * r:
                cand = x + s * D; cand = cand[inside(cand, 0.0)]
                if len(cand) == 0: s /= 2; continue
                cv = clearance(cand, ref, r); j = int(np.argmax(cv))
                if cv[j] > v: x, v = cand[j], cv[j]
                else: s /= 2
            if v > bestv: bestp, bestv = x, v
        holes.append(bestp)
        rho = np.minimum(rho, np.sqrt(((p - bestp) ** 2).sum(1)) - r)
    return np.array(holes)


def contacts(c, r, rel=1e-7):
    n = len(c); cnt = np.zeros(n, int)
    P = cKDTree(c).query_pairs(2 * r * (1 + rel), output_type='ndarray')
    if len(P): np.add.at(cnt, P[:, 0], 1); np.add.at(cnt, P[:, 1], 1)
    g, _ = slp3d.walls(c); cnt += (g <= r * (1 + rel)).sum(1)
    return cnt


def move_relocate(c, r, rng):
    n = len(c); k = int(rng.choice([1, 2, 3], p=[0.5, 0.3, 0.2]))
    key = contacts(c, r) + rng.random(n) * 0.5
    pool = np.argsort(key)[:2 * k]; pick = rng.choice(pool, size=k, replace=False)
    keep = np.setdiff1d(np.arange(n), pick)
    h = find_holes(c[keep], r, k, rng)
    q = c.copy(); q[pick] = h + rng.uniform(-1e-3 * r, 1e-3 * r, h.shape); return slp3d.repair(q), f'reloc{k}'


def move_shake(c, r, rng):
    n = len(c); i0 = int(rng.integers(n)); R = r * rng.uniform(2.5, 6.0); s = r * rng.uniform(0.1, 0.6)
    idx = np.nonzero(((c - c[i0]) ** 2).sum(1) <= R * R)[0]
    v = rng.normal(size=(len(idx), 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    q = c.copy(); q[idx] += v * (s * rng.random(len(idx)) ** (1 / 3))[:, None]; return slp3d.repair(q), 'shake'


def move_jiggle(c, r, rng):
    a = r * 10 ** rng.uniform(-5, -3)
    return slp3d.repair(c + rng.uniform(-a, a, c.shape)), 'jiggle'


# ---------------------------------------------------------------- jobs
def lower():
    common.lower_priority()


def job_polish(args):
    phase, n, cap, deadline, seed = args
    lower()
    B = bars()[str(n)]; now = time.time()
    if now >= deadline: return dict(phase=phase, N=n, skipped=True)
    t0 = time.time()
    try:
        c0, lab = start_coords(n)
    except Exception as ex:
        return dict(phase=phase, N=n, error=f'load: {ex}')
    r0 = slp3d.rmin(c0); st = {}
    c, r = slp3d.polish(c0, t_cap=max(5.0, min(cap, deadline - time.time())), stats=st, pretest=True, seed=seed)
    rd = common.float_r(c); saved = save_best(n, c)
    return dict(phase=phase, N=n, cls=B['cls'], ref=B['ref'], printed=B['printed'], bar=B['bar'], bar_source=B['bar_source'],
                start=lab, r_start=r0, r_best=rd, rel_printed=rd / float(B['printed']) - 1, rel_bar=rd / float(B['bar']) - 1,
                beats_bar=beats(rd, B), saved=saved, secs=round(time.time() - t0, 1),
                stats={k: (round(v, 1) if k == 'lp_secs' else v) for k, v in st.items()})


def job_bh(args):
    phase, n, cpu_budget, wall_cap, deadline, seed = args
    lower()
    B = bars()[str(n)]
    if time.time() >= deadline: return dict(phase=phase, N=n, skipped=True)
    rng = np.random.default_rng(seed); t_cpu0 = time.process_time(); t_w0 = time.time()
    bk = best_known(n)
    if bk is None:
        c0, _ = start_coords(n); c0, _ = slp3d.polish(c0, t_cap=120, pretest=True); src0 = 'fresh'
    else:
        c0, _, src0 = bk
    cur_c = c0.copy(); cur_r = slp3d.rmin(cur_c); best_c, best_r = cur_c.copy(), cur_r; r_in = best_r
    hops = impr = side = 0; kinds = {}; hist = [(0.0, best_r)]; stop = 'budget'
    while True:
        if time.process_time() - t_cpu0 >= cpu_budget: stop = 'budget'; break
        if time.time() - t_w0 >= wall_cap: stop = 'wall_cap'; break
        if time.time() >= deadline: stop = 'deadline'; break
        u = rng.random()
        if u < 0.4: q, kind = move_relocate(cur_c, cur_r, rng)
        elif u < 0.8: q, kind = move_shake(cur_c, cur_r, rng)
        else: q, kind = move_jiggle(cur_c, cur_r, rng)
        q, rq = slp3d.polish(q, t_cap=max(5.0, min(60.0, wall_cap / 3, deadline - time.time())), lp_time=20.0,
                             seed=int(rng.integers(1 << 30)))
        hops += 1; kinds.setdefault(kind, [0, 0]); kinds[kind][0] += 1
        if rq > best_r * (1 + 1e-13):
            best_c, best_r = q.copy(), rq; cur_c, cur_r = q, rq; impr += 1; kinds[kind][1] += 1
            hist.append((round(time.process_time() - t_cpu0, 1), best_r))
        elif rq >= best_r * (1 - 1e-12):
            cur_c, cur_r = q, rq; side += 1
    if time.time() < deadline:
        fc, fr = slp3d.polish(best_c, t_cap=max(5.0, min(120.0, wall_cap / 3, deadline - time.time())), max_iter=800,
                              pretest=True, lp_time=20.0)
        if fr > best_r: best_c, best_r = fc, fr
    rd = common.float_r(best_c); saved = save_best(n, best_c)
    return dict(phase=phase, N=n, cls=B['cls'], ref=B['ref'], printed=B['printed'], bar=B['bar'], bar_source=B['bar_source'],
                start=src0, r_start=r_in, r_best=rd, rel_printed=rd / float(B['printed']) - 1, rel_bar=rd / float(B['bar']) - 1,
                beats_bar=beats(rd, B), saved=saved, hops=hops, improvements=impr, sideways=side, kinds=kinds, stop=stop,
                cpu=round(time.process_time() - t_cpu0, 1), secs=round(time.time() - t_w0, 1), hist=hist[-8:])


def set_out(d):
    global OUT, ROWS
    OUT = os.path.abspath(d); ROWS = os.path.join(OUT, 'rows.jsonl'); os.environ['SCU_SWEEP_OUT'] = OUT


def job_delete(args):
    """PD: deletion seed. If a packing we hold for some M in N+1..N+8 has a larger float radius than our best for N
    (or than N's own published packing), delete M - N fewest-contact spheres from it and polish."""
    phase, n, cap, deadline, seed = args
    lower()
    B = bars()[str(n)]
    if time.time() >= deadline: return dict(phase=phase, N=n, skipped=True)
    rng = np.random.default_rng(seed); t0 = time.time()
    bk = best_known(n)
    if bk is not None: r_cur = bk[1]
    else:
        try: r_cur = slp3d.rmin(start_coords(n)[0])
        except Exception: r_cur = float(B['printed'])
    donors = []
    for m in range(n + 1, min(n + 8, 1008) + 1):
        d = best_known(m)
        if d is not None and d[1] > r_cur * (1 + 1e-12): donors.append((d[1], m, d[0]))
    if not donors:
        return dict(phase=phase, N=n, cls=B['cls'], no_donor=True, r_start=r_cur, secs=round(time.time() - t0, 1))
    rd_, m, c = max(donors, key=lambda t: t[0])
    k = m - n; key = contacts(c, rd_) + rng.random(m) * 0.5
    keep = np.sort(np.argsort(key)[k:])                       # drop the k fewest-contact spheres
    q = c[keep]; st = {}
    q, rq = slp3d.polish(q, t_cap=max(5.0, min(cap, deadline - time.time())), stats=st, seed=seed)
    rd = common.float_r(q); saved = save_best(n, q)
    return dict(phase=phase, N=n, cls=B['cls'], ref=B['ref'], printed=B['printed'], bar=B['bar'], bar_source=B['bar_source'],
                start=f'delete{k}_from_{m}', r_start=r_cur, r_donor=rd_, r_best=rd, rel_printed=rd / float(B['printed']) - 1,
                rel_bar=rd / float(B['bar']) - 1, beats_bar=beats(rd, B), saved=saved, secs=round(time.time() - t0, 1),
                stats={kk: (round(v, 1) if kk == 'lp_secs' else v) for kk, v in st.items()})


def job_seed(args):
    """PI: neighbour seeds for a cell polish did not beat. Donors = packings we hold for M in N-8..N-1 (INSERT the
    N - M missing spheres into the largest holes) and M in N+1..N+8 (DELETE the M - N fewest-contact spheres) whose
    float radius exceeds our current best for N; the best donor on each side is polished (cap each); keep the best."""
    phase, n, cap, deadline, seed = args
    lower()
    B = bars()[str(n)]
    if time.time() >= deadline: return dict(phase=phase, N=n, skipped=True)
    rng = np.random.default_rng(seed); t0 = time.time()
    bk = best_known(n)
    r_cur = bk[1] if bk is not None else slp3d.rmin(start_coords(n)[0])
    lo, hi = [], []
    for m in range(max(2, n - 8), min(n + 8, 1008) + 1):
        if m == n: continue
        d = best_known(m)
        if d is not None and d[1] > r_cur * (1 + 1e-12): (lo if m < n else hi).append((d[1], m, d[0]))
    if not lo and not hi:
        return dict(phase=phase, N=n, cls=B['cls'], no_donor=True, r_start=r_cur, secs=round(time.time() - t0, 1))
    tries = []
    if lo:
        rd_, m, c = max(lo, key=lambda t: t[1]); k = n - m          # the NEAREST qualifying donor below (fewest inserts)
        h = find_holes(c, rd_, k, rng); q = np.vstack([c, h + rng.uniform(-1e-3 * rd_, 1e-3 * rd_, h.shape)])
        tries.append((f'insert{k}_into_{m}', rd_, slp3d.repair(q)))
    if hi:
        rd_, m, c = min(hi, key=lambda t: t[1]); k = m - n          # the NEAREST qualifying donor above (fewest deletes)
        key = contacts(c, rd_) + rng.random(m) * 0.5; keep = np.sort(np.argsort(key)[k:])
        tries.append((f'delete{k}_from_{m}', rd_, c[keep]))
    best = None; st_all = {}
    for lab, rd_, q in tries:
        if time.time() >= deadline: break
        st = {}
        q, rq = slp3d.polish(q, t_cap=max(5.0, min(cap, deadline - time.time())), stats=st, seed=seed)
        st_all[lab] = dict(r=rq, stop=st.get('stop'), lp=st.get('lp'))
        if best is None or rq > best[1]: best = (lab, rq, q, rd_)
    if best is None: return dict(phase=phase, N=n, skipped=True)
    lab, _, q, rd_ = best
    rd = common.float_r(q); saved = save_best(n, q)
    return dict(phase=phase, N=n, cls=B['cls'], ref=B['ref'], printed=B['printed'], bar=B['bar'], bar_source=B['bar_source'],
                start=lab, r_start=r_cur, r_donor=rd_, r_best=rd, rel_printed=rd / float(B['printed']) - 1,
                rel_bar=rd / float(B['bar']) - 1, beats_bar=beats(rd, B), saved=saved, secs=round(time.time() - t0, 1),
                stats=dict(stop=st_all.get(lab, {}).get('stop'), lp=st_all.get(lab, {}).get('lp'), lp_fail=None, tries=st_all))


def done_pairs():
    s = set(); rows = []
    if os.path.exists(ROWS):
        for line in open(ROWS):
            try:
                x = json.loads(line)
            except Exception:
                continue
            rows.append(x)
            if not x.get('skipped') and 'error' not in x: s.add((x['phase'], x['N']))
    return s, rows


def run_pool(pool, fn, jobs, T0, tally):
    for x in pool.imap_unordered(fn, jobs):
        with open(ROWS, 'a') as f: f.write(json.dumps(x) + '\n')
        if x.get('skipped'):
            tally['skipped'] += 1; continue
        if x.get('no_donor'):
            tally['no_donor'] = tally.get('no_donor', 0) + 1; continue
        if 'error' in x:
            tally['error'] += 1; print(f"{x['phase']} N={x['N']} ERROR {x['error']}", flush=True); continue
        tally['done'] += 1; tally['beats'] += x['beats_bar']
        extra = (f"stop {x['stats'].get('stop')} lp {x['stats'].get('lp')} fail {x['stats'].get('lp_fail')}" if 'stats' in x
                 else f"hops {x['hops']} impr {x['improvements']} side {x['sideways']} {x['kinds']} stop {x['stop']}")
        print(f"{x['phase']} N={x['N']:4d} {x['cls'][:4]} {x['ref']:5s} start {x['start']:16s} vs printed {x['rel_printed']:+.3e} "
              f"vs bar({x['bar_source'][:14]}) {x['rel_bar']:+.3e} {'BEATS' if x['beats_bar'] else '-'} {x['secs']:.0f}s {extra}"
              f"   [{tally['done']} done, {tally['beats']} beat, {tally['skipped']} skipped, {time.time() - T0:.0f}s]", flush=True)


def main():
    lower()
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4); ap.add_argument('--deadline-h', type=float, default=4.5)
    ap.add_argument('--cap', type=float, default=300.0); ap.add_argument('--bh-wall', type=float, default=600.0)
    ap.add_argument('--bh-cpu', type=float, default=300.0); ap.add_argument('--phases', default='1,2,5,3,4')
    ap.add_argument('--only', default=''); ap.add_argument('--no-resume', action='store_true')
    ap.add_argument('--out', default=OUT)
    ap.add_argument('--deadline-at', type=float, default=0.0, help='absolute deadline (unix time); overrides --deadline-h')
    a = ap.parse_args(); W = min(a.workers, 4)
    set_out(a.out); os.makedirs(OUT, exist_ok=True)
    T0 = time.time(); deadline = a.deadline_at if a.deadline_at > 0 else T0 + 3600 * a.deadline_h; phases = {int(p) for p in a.phases.split(',')}
    T = common.table(); p1, p2, p4 = phase_lists(T)
    only = {int(x) for x in a.only.split(',')} if a.only else None
    done, _ = (set(), []) if a.no_resume else done_pairs()
    sel = lambda ph, L: [n for n in L if (ph, n) not in done and (only is None or n in only)]
    print(f'scu sweep start {time.strftime("%Y-%m-%d %H:%M:%S")}  deadline {time.strftime("%H:%M:%S", time.localtime(deadline))}  workers {W}  cap {a.cap:.0f}s  '
          f'P1 {len(p1)}  P2 {len(p2)}  P4 {len(p4)}  resume: {len(done)} (phase, N) done', flush=True)
    tally = dict(done=0, beats=0, skipped=0, error=0)
    with Pool(W) as pool:
        jobs = []
        if 1 in phases: jobs += [('P1', n, a.cap, deadline, SEED + n) for n in sel('P1', p1)]
        if 2 in phases: jobs += [('P2', n, a.cap, deadline, SEED + n) for n in sel('P2', p2)]
        print(f'P1+P2 jobs now: {len(jobs)}', flush=True)
        run_pool(pool, job_polish, jobs, T0, tally)
        if 5 in phases and time.time() < deadline:
            jobs = [('PD', n, 120.0, deadline, SEED + 3 * n) for n in sel('PD', list(range(1007, 93, -1)))]
            print(f'PD deletion seeds: {len(jobs)} cells (N 1007 -> 94; donors N+1..N+8)', flush=True)
            run_pool(pool, job_delete, jobs, T0, tally)
        def run_pi(L, what):
            B = bars(); jobs = []
            for n in sel('PN', L):
                bk = best_known(n)
                if bk is None or not beats(bk[1], B[str(n)]): jobs.append(n)
            jobs.sort()
            print(f'PI neighbour seeds, {what} (insert from N-8..N-1 / delete from N+1..N+8): {len(jobs)} cells not beaten', flush=True)
            run_pool(pool, job_seed, [('PN', n, 150.0, deadline, SEED + 11 * n) for n in jobs], T0, tally)
        if 6 in phases and time.time() < deadline: run_pi(p1 + p2, 'non-derivative cells')
        if 3 in phases and time.time() < deadline:
            _, rows = done_pairs(); B = bars(); L = stale_blocks(T)
            best = {}
            for x in rows:
                if x.get('skipped') or 'error' in x or x.get('no_donor'): continue
                best[x['N']] = max(best.get(x['N'], -1.0), x['r_best'])
            beaten = {n for n, r in best.items() if beats(r, B[str(n)])}
            for n in PROBE20:                                    # the probe v2 files count too
                bk = best_known(n)
                if bk is not None and beats(bk[1], B[str(n)]): beaten.add(n)
            for n in list(best):
                bk = best_known(n)
                if bk is not None and beats(bk[1], B[str(n)]): beaten.add(n)
            cand = [n for n in p1 + p2 if n in best and n not in beaten]
            cand.sort(key=lambda n: (-L[n], -n))
            cand = sel('P3', cand)
            print(f'P3 basin hopping: {len(cand)} cells not beaten after polish (largest stale blocks first): {cand[:30]} ...', flush=True)
            run_pool(pool, job_bh, [('P3', n, a.bh_cpu, a.bh_wall, deadline, SEED + 7 * n) for n in cand], T0, tally)
        if 7 in phases and time.time() < deadline: run_pi(p4, 'derivative cells')
        if 4 in phases and time.time() < deadline:
            jobs = [('P4', n, a.cap, deadline, SEED + n) for n in sel('P4', p4)]
            print(f'P4 derivative cells: {len(jobs)}', flush=True)
            run_pool(pool, job_polish, jobs, T0, tally)
    print(f'SWEEP END {time.strftime("%H:%M:%S")}: {tally}; wall {time.time() - T0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
