# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""chx float wins -> exact certificates, both checkers.  Moki&Julio 2026-09-26.

For every cell searched (out/probe/rows.jsonl, out/sweep/rows.jsonl) take the best float packing (best_<N>.npy), recompute
its float clearance r_f (dense O(N^2) pairs + all wall slacks, kgon Frame(6)) and the effective bar
  bar_eff = max(data/refs/chx/chx_bar.json bar, every POLISHED Amore radius a run recorded).
If the float claim rule holds (r_f > bar_eff (1 + 1e-10) and r_f >= printed + 2e-12):
  write cert/chx_<N>.txt : "r <R>", then N lines "x y" = the EXACT decimal expansions of the float coordinates (frame of
        certify_kgon k = 6, which is Packomania's chx frame: no rotation); R = a rigorous LOWER bound of the true min clearance of
        those exact centres (pairs: integer square roots of exact squared distances; walls: certify_kgon's interval enclosures at
        10^-82), ROUNDED DOWN to 25 significant digits (same construction as the k-gon certification step, copied here, not imported);
  run   checker A  certify_kgon.py 6 <file> <N> <rec>         for rec = bar_eff and rec = printed
        checker B  checkers/verify_exact_kgon.py 6 <file> <N> <rec> (black box) for rec = bar_eff and rec = printed
  claim rule, exact decimals: R > bar_eff (1 + 1e-10) AND R >= printed + 2e-12.
  WIN = all four verdicts IMPROVES + claim rule.
  Controls (checker A): the same centres at R (1 + 1e-15) and the centres rotated by 90 deg must both be INVALID.
Writes out/cand_chx.json (fields as the ball tables' candidate files).
usage: python certify_chx.py
"""
import os, sys, json, math, time, subprocess
from decimal import Decimal, getcontext
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import chx_prio                   # noqa: E402
import chx_common as C            # noqa: E402
ck = C.ck
getcontext().prec = 200
CERT = os.path.join(HERE, 'cert'); OUT = os.path.join(HERE, 'out')
A_PATH = os.path.join(os.path.dirname(ROOT), 'checkers', 'certify_kgon.py'); B_PATH = os.path.join(os.path.dirname(ROOT), 'checkers', 'verify_exact_kgon.py')
DIRS = tuple((nm, os.path.join(OUT, nm)) for nm in ('probe', 'sweep', 'sweep2', 'sweep3')      # sweep2/3: 09-26 afternoon (second reference's cells; seed + 1)
             if os.path.exists(os.path.join(OUT, nm, 'rows.jsonl')))


def exact_dec(v):
    s = format(Decimal(float(v)), 'f')
    return '0' if s in ('-0', '0') else s


def min_radius_lower_bound(k, pts, r_f, G=82):
    vals = [ck.dec(t) for p in pts for t in p]
    D = max(0, max(-e for _, e in vals)); S = 10 ** D
    I = [m * 10 ** (e + D) for m, e in vals]; X = I[0::2]; Y = I[1::2]; n = len(X)
    (Alo, Ahi), sides = ck.table(k, G); wl = None
    for i in range(n):
        for xl, xh, yl, yh in sides:
            a1, a2, b1, b2 = X[i] * xl, X[i] * xh, Y[i] * yl, Y[i] * yh
            lo = Alo * S - max(a1, a2) - max(b1, b2)
            if wl is None or lo < wl: wl = lo
    best = Fr(wl, S * 10 ** G)
    c = np.array([[float(x), float(y)] for x, y in pts])
    if n > 1:
        from scipy.spatial import cKDTree
        P = cKDTree(c).query_pairs(2 * r_f * (1 + 1e-6), output_type='ndarray')
        dd, jj = cKDTree(c).query(c, k=2); i0 = int(np.argmin(dd[:, 1]))
        cand = set(map(tuple, P.tolist())) | {tuple(sorted((i0, int(jj[i0, 1]))))}
        Q = 10 ** 60
        for i, j in cand:
            d2 = (X[i] - X[j]) ** 2 + (Y[i] - Y[j]) ** 2
            lb = Fr(math.isqrt(d2 * Q * Q), 2 * S * Q)
            if lb < best: best = lb
    return best


def floor_sig(q, sig=25):
    e = math.floor(math.log10(q.numerator) - math.log10(q.denominator))
    for ee in (e - 1, e, e + 1):
        if Fr(10) ** ee <= q < Fr(10) ** (ee + 1): e = ee; break
    places = sig - 1 - e
    m = (q.numerator * 10 ** places) // q.denominator if places >= 0 else q.numerator // (q.denominator * 10 ** -places)
    return format(Decimal(m).scaleb(-places), 'f')


def verdict(script, path, n, rec):
    p = subprocess.run([sys.executable, script, '6', path, str(n), rec], capture_output=True, text=True, cwd=ROOT)
    v = [l for l in p.stdout.splitlines() if l.startswith('VERDICT')]
    return v[0].split(':', 1)[1].strip() if v else 'ERROR ' + (p.stderr[-300:] or p.stdout[-300:])


def rows(d):
    out = {}; p = os.path.join(d, 'rows.jsonl')
    if os.path.exists(p):
        for line in open(p):
            x = json.loads(line)
            if not x.get('skipped'): out[x['n']] = x
    return out


def main():
    print('priority class:', chx_prio.lower(), flush=True)
    os.makedirs(CERT, exist_ok=True); os.makedirs(OUT, exist_ok=True)
    F = C.frame(); BAR = {int(n): g for n, g in json.load(open(os.path.join(C.REFS, 'chx_bar.json'))).items()}
    R_ = {name: rows(d) for name, d in DIRS}
    Ns = sorted(set().union(*[set(r) for r in R_.values()])); cells = []; t0 = time.time()
    for n in Ns:
        g = BAR[n]; rp_s = g['printed']; rp = float(rp_s); srcs = []
        for name, d in DIRS:
            f = os.path.join(d, f'best_{n}.npy')
            if n in R_[name] and os.path.exists(f):
                c = np.load(f); srcs.append((name, c, F.dense_r(c), R_[name][n]))
        bar_d = Decimal(g['bar']); bar_src = g['bar_source']
        for _, _, _, meta in srcs:
            if meta.get('r_amore_polished') is not None and Decimal(repr(meta['r_amore_polished'])) > bar_d:
                bar_d = Decimal(repr(meta['r_amore_polished'])); bar_src = 'amore2023_polished'
        bar_s = format(bar_d, 'f'); bar_f = float(bar_d)
        src, c, r_f, meta = max(srcs, key=lambda t: t[2])
        route = f"{src}:start_{meta.get('start')}:{'slp_polish' if meta.get('improvements', 0) == 0 else 'slp_polish+basin_hop'}"
        row = dict(N=n, printed=rp_s, page_refs=g['page_refs'], bar=bar_s, bar_source=bar_src, r_float=r_f, file=src, route=route,
                   gain_float_rel_vs_bar=r_f / bar_f - 1, gain_float_rel_vs_printed=r_f / rp - 1)
        if not (r_f > bar_f * (1 + 1e-10) and r_f >= rp + 2e-12):
            row['status'] = 'below_float_bar'; cells.append(row); continue
        pts = [(exact_dec(x), exact_dec(y)) for x, y in c]
        R = floor_sig(min_radius_lower_bound(6, pts, r_f), 25)
        path = os.path.join(CERT, f'chx_{n}.txt'); ck.write_cert_file(path, R, pts)
        va_bar = verdict(A_PATH, path, n, bar_s); va_pub = verdict(A_PATH, path, n, rp_s)
        vb_bar = verdict(B_PATH, path, n, bar_s); vb_pub = verdict(B_PATH, path, n, rp_s)
        Rd = Decimal(R)
        claim = bool(Rd > bar_d * (1 + Decimal('1e-10')) and Rd >= Decimal(rp_s) + Decimal('2e-12'))
        # controls (checker A, in memory): radius raised by 1e-15 relative; centres rotated by 90 degrees
        R_up = ck.fmt_dec(ck.dec_fraction(R) * Fr(10 ** 15 + 1, 10 ** 15), 40)
        ctl_up = ck.check_values(6, R_up, pts, n, '0')
        rot = [(('-' + y) if not y.startswith('-') else y[1:], x) for x, y in pts]
        ctl_rot = ck.check_values(6, R, rot, n, '0')
        win = all(v == 'IMPROVES' for v in (va_bar, va_pub, vb_bar, vb_pub)) and claim
        row.update(status='WIN' if win else 'rejected', cert=os.path.relpath(path, ROOT).replace('\\', '/'), r_new=R,
                   checker_a_vs_bar=va_bar, checker_a_vs_printed=va_pub, checker_b_vs_bar=vb_bar, checker_b_vs_printed=vb_pub,
                   claim_rule=claim, win=win, control_r_up_1e15=ctl_up, control_rotated_90=ctl_rot,
                   gain_abs_vs_printed=float(Rd - Decimal(rp_s)), gain_rel_vs_printed=float(Rd / Decimal(rp_s) - 1),
                   gain_rel_vs_bar=float(Rd / bar_d - 1))
        cells.append(row)
        print(f"N={n:3d} [{','.join(g['page_refs'])}] {route}: cert r {R} (vs bar {row['gain_rel_vs_bar']:+.3e} [{bar_src}], "
              f"vs printed {row['gain_rel_vs_printed']:+.3e})  A bar/pub {va_bar}/{va_pub}  B bar/pub {vb_bar}/{vb_pub}  "
              f"claim {claim}  controls {ctl_up[:7]}/{ctl_rot[:7]}  -> {'WIN' if win else 'rejected'}", flush=True)
    doc = dict(table='chx', k=6, built=time.strftime('%Y-%m-%d %H:%M:%S'),
               frame='certify_kgon k = 6: circumradius 1, centred at the origin, flat bottom side (= Packomania chx frame, no rotation)',
               rule='win = r_new > bar (1 + 1e-10) and r_new >= printed + 2e-12, checkers A (certify_kgon.py) and B '
                    '(checkers/verify_exact_kgon.py) IMPROVES vs bar and vs printed',
               cells=cells)
    json.dump(doc, open(os.path.join(OUT, 'cand_chx.json'), 'w'), indent=1)
    st = {}
    for x in cells: st[x['status']] = st.get(x['status'], 0) + 1
    print(f'{len(cells)} cells; status {st}; wins: {[x["N"] for x in cells if x.get("win")]}; {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
