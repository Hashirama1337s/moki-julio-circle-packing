# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Bars for equal spheres in the unit sphere, Packomania ssp, N = 1..1000 (2026-09-26; Moki&Julio)
-> data/refs/ssp/ssp_bar.json (+ a summary on stdout).

bar(N) = max of
  page       Packomania's printed radius (saved page of 2026-09-26, "Last update 13-Jul-2026"; 12 decimals), and the value
             of the page's own txt/radius.txt (identical for every N <= 1000, checked);
  zhou2023   Zhou, Ren, He, Liu, Li, arXiv:2305.10023v1, Tables 2-5 (N <= 400): r_up = 1 / (R_best - 5e-11), the largest
             radius their 10-decimal container ratio R allows (data/refs/ssp/lit_zhou2023_arxiv.json);
  code       Cohn's table of spherical codes, d = 3 (data/refs/lit/cohn_codes/spherical_codes_plain_2026-09-25.txt, M = 7..1024,
             + the exact small codes M = 2..6: maxcos -1, -1/2, -1/3, 0, 0), as in ../hsp/bars.py:
               code alone  : M >= N points on the shell |c| = 1 - r: r = s / (1 + s), s = sqrt((1 - maxcos) / 2);
               code+centre : M >= N - 1 on the shell + one sphere at the origin: r = min(s / (1 + s), 1/3);
             a 12-decimal maxcos (printed rounded UP) is lowered by 1e-12 first (the bar is never below the code radius);
  larger_N   the largest printed radius at any M > N on the page (rows up to N = 3 million; page and radius.txt).
No other data-level source exists (the recon of the ssp table): Lai et al. [13] is a private communication to Specht (no paper, no data);
the C&OR 2024 version of Zhou et al. and their Dec 2023 / May 2024 sets are on the page only (refs [12]).
Every value is computed in Decimal (50 digits); a non-page bar is rounded UP to 20 significant digits.
Claim rule (the certification step): r_new > bar (1 + 1e-10) and r_new >= printed + 2e-12, exact decimals.
usage: python bars.py
"""
import os, sys, json
from decimal import Decimal, getcontext, ROUND_CEILING
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sspgeo                                               # noqa: E402
import prio                                                 # noqa: E402  (../hsp/prio.py)
prio.lower()
import geomd                                                # noqa: E402  (cohn_table)
getcontext().prec = 50
ONE, THIRD = Decimal(1), Decimal(1) / Decimal(3)
BARS = os.path.join(sspgeo.REFS, 'ssp_bar.json')
SMALL = {2: '-1', 3: '-0.5', 4: str(Decimal(-1) / 3), 5: '0', 6: '0'}      # exact (code_r does not lower them: < 12 decimals or exact)


def code_r(cos_s, exact=False):
    c = Decimal(cos_s)
    if not exact and len(cos_s.split('.')[-1]) >= 12: c -= Decimal('1e-12')
    s = ((ONE - c) / 2).sqrt()
    return s / (ONE + s)


def up20(x):
    x = Decimal(x)
    return format(x.quantize(Decimal(1).scaleb(x.adjusted() - 19), rounding=ROUND_CEILING), 'f')


def radius_txt():
    T = {}
    for l in open(os.path.join(sspgeo.BIG, 'ssp_radius_2026-09-26.txt')):
        t = l.split()
        if len(t) == 2 and t[0].isdigit(): T[int(t[0])] = t[1]
    return T


def build():
    T = sspgeo.table(); TA = sspgeo.table(all_rows=True); RT = radius_txt()
    Z = json.load(open(os.path.join(sspgeo.REFS, 'lit_zhou2023_arxiv.json')))['cells']
    R = {m: code_r(cs) for (d, m), cs in geomd.cohn_table().items() if d == 3}
    for m, cs in SMALL.items(): R[m] = code_r(cs, exact=True)
    Ms = sorted(R); best_ge = {}; b = (Decimal(0), None)
    for m in reversed(Ms):
        if R[m] > b[0]: b = (R[m], m)
        best_ge[m] = b

    def ge(n):
        k = [m for m in Ms if m >= n]
        return best_ge[k[0]] if k else (Decimal(0), None)
    allr = {}
    for m in set(TA) | set(RT):
        v = [Decimal(x) for x in (TA.get(m, {}).get('r'), RT.get(m)) if x]
        allr[m] = max(v)
    mism = [(m, TA[m]['r'], RT.get(m)) for m in TA if RT.get(m) and RT[m] != TA[m]['r']]
    cells = {}; src = {}
    for n in range(1, 1001):
        page = Decimal(T[n]['r'])
        ca = ge(n); cc = ge(n - 1) if n >= 2 else (Decimal(0), None); cc = (min(cc[0], THIRD), cc[1])
        lg = max(((allr[m], m) for m in allr if m > n), default=(Decimal(0), None))
        z = Z.get(str(n)); zr = Decimal(z['r_up']) if z and z.get('r_up') else Decimal(0)
        cands = [(page, 'page', None), (zr, 'zhou2023', None), (ca[0], 'code', ca[1]), (cc[0], 'code+centre', cc[1]),
                 (lg[0], 'larger_N', lg[1])]
        top = max(cands, key=lambda t: t[0])
        # the page wins ties; a literature value counts only when it is above the page
        if top[0] <= page: top = (page, 'page', None)
        cells[str(n)] = dict(N=n, printed=T[n]['r'], radius_txt=RT.get(n), ref=T[n]['ref'], contacts=T[n]['contacts'],
                             loose=T[n]['loose'], boundary=T[n]['boundary'], core=T[n]['core'],
                             zhou2023_r_best=z.get('r_best') if z else None, zhou2023_r_up=z.get('r_up') if z else None,
                             code_alone=dict(M=ca[1], r=up20(ca[0]) if ca[1] else None),
                             code_centre=dict(M=cc[1], r=up20(cc[0]) if cc[1] else None),
                             larger_N=dict(M=lg[1], r=str(lg[0]) if lg[1] else None),
                             bar=T[n]['r'] if top[1] == 'page' else up20(top[0]), bar_source=top[1], bar_M=top[2],
                             bar_minus_printed=float(top[0] - page))
        src[top[1]] = src.get(top[1], 0) + 1
    json.dump(dict(built='2026-09-26', rule='bar = max(page (= radius.txt), Zhou 2023 arXiv r_up (N <= 400), code alone M >= N, '
                   'code M >= N-1 + centre (<= 1/3), larger-N page); claim: r_new > bar (1 + 1e-10) and r_new >= printed + 2e-12',
                   sources=src, page_vs_radius_txt_mismatch_all_rows=mism, cells=cells), open(BARS, 'w'), indent=1)
    return cells, src, mism


def load():
    return json.load(open(BARS))['cells']


def bar(n, B=None):
    return (B or load())[str(n)]


if __name__ == '__main__':
    cells, src, mism = build()
    print('bar sources:', src)
    print('page vs radius.txt mismatches (all rows; N > 1000 only matter via larger_N):', mism)
    nz = [(c['N'], c['bar_source'], c['bar_minus_printed']) for c in cells.values() if c['bar_source'] != 'page']
    print(f'{len(nz)} cells with a non-page bar:', nz[:60])
    nonmono = [(n, cells[str(n)]['printed'], cells[str(n)]['larger_N']) for n in range(1, 1000)
               if Decimal(cells[str(n)]['larger_N']['r'] or 0) > Decimal(cells[str(n)]['printed'])]
    print('page non-monotone (a larger N has a larger printed radius):', nonmono)
    for n in (2, 4, 6, 12, 13, 14, 20, 57, 100, 200, 400, 401, 600, 1000):
        x = cells[str(n)]
        print(f"N={n}: printed {x['printed']} zhou {x['zhou2023_r_up']} code {x['code_alone']} code+centre {x['code_centre']} "
              f"bar {x['bar']} ({x['bar_source']})")
