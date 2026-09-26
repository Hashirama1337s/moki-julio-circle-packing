# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""The claim bar per N for chx (equal circles in a regular hexagon).  Moki&Julio 2026-09-26.

  bar(N) = max over every N' >= N of every data-level radius at N' :
             Packomania page radius (data/refs/chx/packomania_chx_2026-09-26.html, 12 decimals),
             Packomania txt/radius.txt (data/refs/chx/chx_radius_2026-09-26.txt; differs from the page at N = 4, 84, 89, 95, 114),
             Amore 2023 (Zenodo 7574070, Sides_6, N 2-400) recomputed from his centres (data/refs/chx/lit_amore2023.json).
  (N' > N enters because deleting circles from an N'-packing gives an N-packing of the same radius.)
  claim rule: r_new > bar (1 + 1e-10)  AND  r_new >= printed + 2e-12  (printed = the page radius), decided exactly.
Other sources checked at data level and found to hold no hexagon values: data/refs/chx/lit_search_2026-09-26.json (arXiv / Zenodo /
GitHub sweep), Amore-Carrizalez-Zarate 2023 'Echoes of the hexagon' Zenodo 7507985 (folders Sides_0, 12..60 only:
data/refs/chx/echoes_zenodo7507985_names.json).
Writes data/refs/chx/chx_bar.json.
usage: python bars.py
"""
import os, sys, json
from decimal import Decimal
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import chx_common as C            # noqa: E402

if __name__ == '__main__':
    T = C.page_table(); RT = C.radius_txt()
    AM = json.load(open(os.path.join(C.REFS, 'lit_amore2023.json')))['radii']
    vals = {}                                                    # N -> [(Decimal, source)]
    for n, v in T.items(): vals.setdefault(n, []).append((Decimal(v['radius']), 'packomania_page'))
    for n, r in RT.items(): vals.setdefault(n, []).append((Decimal(r), 'packomania_radius_txt'))
    for n, v in AM.items(): vals.setdefault(int(n), []).append((Decimal(v['r']), 'amore2023'))
    Ns = sorted(vals); out = {}; best = (Decimal(0), None, None)
    for n in reversed(Ns):                                       # running max over N' >= N
        for r, src in vals[n]:
            if r > best[0]: best = (r, src, n)
        if n in T:
            here = {src: str(r) for r, src in vals[n]}
            bar, src, nfrom = best
            out[n] = dict(printed=T[n]['radius'], page_refs=T[n]['refs'], radius_txt=RT.get(n), amore2023=here.get('amore2023'),
                          bar=str(bar), bar_source=src + ('' if nfrom == n else f'@N={nfrom}'), bar_from_N=nfrom,
                          bar_vs_printed_rel=float(bar / Decimal(T[n]['radius']) - 1))
    json.dump({str(n): out[n] for n in sorted(out)}, open(os.path.join(C.REFS, 'chx_bar.json'), 'w'), indent=1)
    raised = {n: o for n, o in out.items() if o['bar_source'] != 'packomania_page'}
    print(f'{len(out)} page cells; bar above the page radius at {len(raised)}:')
    for n, o in sorted(raised.items()):
        print(f"  N={n}: printed {o['printed']}  bar {o['bar']}  ({o['bar_source']}, {o['bar_vs_printed_rel']:+.2e})")
