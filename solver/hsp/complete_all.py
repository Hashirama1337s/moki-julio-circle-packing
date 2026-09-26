# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Save the recovered (completed, NOT polished) packing of every truncated hsp5 file (recon56.py, exact MILP) or hsp6 file
(recover6.py, rigid circle-intersection growth) as out/recon56/complete/hsp<d>_<N>.npy + complete_hsp5.json (per N: r_factor, violation, interior, completed radius). These are
Specht's packings, recovered: the build's table-radius audit and classify read them (2026-09-26; Moki&Julio).
usage: python complete_all.py [--d 5]"""
import os, sys, json, argparse, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prio, geomd, recon56                                        # noqa: E402
OUT = os.path.join(geomd.OUT, 'recon56', 'complete')
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--d', type=int, default=5); a = ap.parse_args(); prio.lower()
    os.makedirs(OUT, exist_ok=True); rows = {}
    for n in range(2, geomd.NMAX[a.d] + 1):
        if not os.path.exists(geomd.pub_path(a.d, n)): continue
        t = time.process_time(); r, A = recon56.trunc(a.d, n)
        if a.d == 6:                                                   # 6-D: the rigid circle-intersection growth (recover6.py)
            import recover6
            C, v, info = recover6.recover6(A, r, np.random.default_rng(n)); log = [dict(r_factor=0.0, interior=None, method=info['method'])]
        else:
            C, v, log = recon56.complete(a.d, A, r, np.random.default_rng(n), starts=8)
        np.save(os.path.join(OUT, f'hsp{a.d}_{n}.npy'), C)
        rows[n] = dict(printed=r, r_factor=log[0].get('r_factor') if log else None, viol=v, interior=log[0].get('interior') if log else None,
                       r_completed=geomd.rmin(C), cpu=round(time.process_time() - t, 1))
        print(n, rows[n], flush=True)
    json.dump(rows, open(os.path.join(OUT, f'complete_hsp{a.d}.json'), 'w'), indent=1)
    print('done', len(rows), 'ok', sum(1 for x in rows.values() if x['viol'] < 1e-9 and x['r_factor'] is not None and x['r_factor'] <= 1e-6))
