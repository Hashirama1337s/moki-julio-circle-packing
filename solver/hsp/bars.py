# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Bars for equal balls in the unit d-ball, d = 4, 5, 6 (2026-09-25; Moki&Julio) -> out/bars.json.

bar(d, N) = max of
  page       Packomania's printed radius (HTML page of 2026-09-25; 12 decimals);
  code       the best code-derived radius for N balls from Cohn's table of spherical codes (plain-text table
             data/refs/lit/cohn_codes/spherical_codes_plain_2026-09-25.txt, d <= 32, M <= 1024):
               code alone  : M >= N points on the shell (deleting M - N points only helps): r = s / (1 + s),
               code+centre : M >= N - 1 points on the shell + one ball at the origin: r = min(s / (1 + s), 1/3),
             s = sin(theta / 2) = sqrt((1 - maxcos) / 2). The table prints maxcos ROUNDED UP at 12 decimals (e.g. its
             '3,9,0.333333333334' for the exact 1/3), so the bar uses maxcos - 1e-12 (the largest angle it allows) for a
             12-decimal value and the value itself for a shorter (exact) one: the bar is never below the true code radius;
  larger_N   the largest printed radius at any M > N in the same table (a packing of M balls minus M - N balls).
Every value is computed in Decimal (50 digits) and the bar is rounded UP to 20 significant digits.
Cross-checks: the local coordinate files data/big/hsp_codes/c<d>_<M>.txt (float max cosine, <= table maxcos + 2e-12) and
the recon's data/big/hsp_codes/code_derived_beats_2026-09-25.json (177 cells: our code value >= its float value - 1e-15).
Claim rule (used by the certification step): r_new > bar (1 + 1e-10) and r_new >= printed + 2e-12, in exact decimals.
Credit: a cell whose best packing is a code construction (no gain from our polish over the code-derived radius) is NOT
ours -- the claim rule already requires beating the code value.

usage: py -3.11 bars.py      -> out/bars.json + a summary
"""
import os, sys, json
from decimal import Decimal, getcontext, ROUND_CEILING
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import geomd                                           # noqa: E402
getcontext().prec = 50
ONE, THIRD = Decimal(1), Decimal(1) / Decimal(3)
BARS = os.path.join(geomd.OUT, 'bars.json')


def code_r(cos_s):
    c = Decimal(cos_s)
    if len(cos_s.split('.')[-1]) >= 12: c -= Decimal('1e-12')
    s = ((ONE - c) / 2).sqrt()
    return s / (ONE + s)


def up20(x):
    x = Decimal(x)
    return format(x.quantize(Decimal(1).scaleb(x.adjusted() - 19), rounding=ROUND_CEILING), 'f')


def build():
    cohn = geomd.cohn_table(); out = {}
    recon = {(int(x['table'][3:]), x['N']): x for x in json.load(open(os.path.join(geomd.CODES, 'code_derived_beats_2026-09-25.json'), encoding='utf-8'))}
    checks = dict(local_files=0, local_bad=[], recon=0, recon_bad=[])
    for d in (4, 5, 6):
        T = geomd.table(d); nmax = geomd.NMAX[d]
        R = {m: code_r(cs) for (dd, m), cs in cohn.items() if dd == d}
        Ms = sorted(R)
        # suffix maxima: best code with >= M points
        best_ge = {}; b = (Decimal(0), None)
        for m in reversed(Ms):
            if R[m] > b[0]: b = (R[m], m)
            best_ge[m] = b
        def ge(n):
            k = [m for m in Ms if m >= n]
            return best_ge[k[0]] if k else (Decimal(0), None)
        page = {n: Decimal(T[n]['r']) for n in T}
        cells = {}
        for n in range(1, nmax + 1):
            ca = ge(n); cc = ge(n - 1) if n >= 2 else (Decimal(0), None)
            cc = (min(cc[0], THIRD), cc[1])
            lg = max(((page[m], m) for m in page if m > n), default=(Decimal(0), None))
            cands = [(page[n], 'page', None), (ca[0], 'code', ca[1]), (cc[0], 'code+centre', cc[1]), (lg[0], 'larger_N', lg[1])]
            top = max(cands, key=lambda t: t[0])
            cells[str(n)] = dict(N=n, printed=T[n]['r'], ref=T[n]['ref'], contacts=T[n]['contacts'], loose=T[n]['loose'],
                                 boundary=T[n]['boundary'], core=T[n]['core'],
                                 code_alone=dict(M=ca[1], r=up20(ca[0]) if ca[1] else None),
                                 code_centre=dict(M=cc[1], r=up20(cc[0]) if cc[1] else None),
                                 larger_N=dict(M=lg[1], r=str(lg[0]) if lg[1] else None),
                                 bar=up20(top[0]) if top[1] != 'page' else T[n]['r'], bar_source=top[1], bar_M=top[2])
            rc = recon.get((d, n))
            if rc:
                checks['recon'] += 1
                mine = max(Decimal(cells[str(n)]['code_alone']['r'] or 0), Decimal(cells[str(n)]['code_centre']['r'] or 0))
                if mine < Decimal(repr(rc['r_code_float'])) - Decimal('1e-15'): checks['recon_bad'].append((d, n, str(mine), rc['r_code_float']))
        for m in Ms:
            U = geomd.load_code(d, m)
            if U is None: continue
            checks['local_files'] += 1; mc = geomd.max_cos(U)
            if mc > float(cohn[(d, m)]) + 2e-12: checks['local_bad'].append((d, m, mc, cohn[(d, m)]))
        out[f'd{d}'] = cells
    os.makedirs(geomd.OUT, exist_ok=True)
    json.dump(dict(built='2026-09-25', rule='bar = max(page, code alone M>=N, code M>=N-1 + centre (<= 1/3), larger-N page); '
                   'claim: r_new > bar (1 + 1e-10) and r_new >= printed + 2e-12', checks=checks, cells=out),
              open(BARS, 'w'), indent=1)
    return out, checks


def load():
    return json.load(open(BARS))['cells']


def bar(d, n, B=None):
    return (B or load())[f'd{d}'][str(n)]


if __name__ == '__main__':
    out, ck = build()
    print(f"checks: local code files {ck['local_files']} (max cos above the table + 2e-12: {ck['local_bad'] or 'none'}); "
          f"recon json cells {ck['recon']} (ours below the recon's float value: {ck['recon_bad'] or 'none'})")
    for d in (4, 5, 6):
        C = out[f'd{d}']; src = {}
        for x in C.values(): src[x['bar_source']] = src.get(x['bar_source'], 0) + 1
        print(f'hsp{d}: {len(C)} cells, bar sources {src}')
    for d, n in ((4, 25), (4, 30), (4, 87), (4, 119), (4, 200), (5, 42), (5, 160), (5, 250), (6, 66), (6, 155), (6, 200)):
        x = out[f'd{d}'][str(n)]
        print(f"hsp{d} N={n}: printed {x['printed']}  code {x['code_alone']}  code+centre {x['code_centre']}  bar {x['bar']} ({x['bar_source']} M={x['bar_M']})")
