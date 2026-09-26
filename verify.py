# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Verify every record in this repository with BOTH independent exact checkers.

    python verify.py            (Python 3.9+, standard library only; ~5-15 minutes on 12 threads)

For each row of MANIFEST.csv: unzip the certificate from certificates/<shelf>.zip, then
  1. checkers/certify_circ.py      (checker A)  - exact rational arithmetic; IMPROVES needs radius > published * (1 + 1e-10)
  2. checkers/verify_exact_*.py    (checker B, written independently) - exact rational arithmetic
The square: checkers/certify_big.py (A) + checkers/verify_exact_big.py (B). The regular pentagon (cpt): checkers/certify_poly.py (A)
+ checkers/verify_exact_poly.py (B), both exact in Q(sqrt5, sqrt(10 +- 2 sqrt5)); they say IMPROVES for any gain over the listed
radius (every pentagon record was selected only if it beats that radius by more than 1e-10 relative).
The regular 16-gon (cxd) and 15-gon (cpd): checkers/certify_kgon.py (A) + checkers/verify_exact_kgon.py (B), both exact (the
irrational walls through rigorous rational enclosures); they too say IMPROVES for any gain (every record beats its radius by > 1e-10).
The cube (scu: equal SPHERES in the cube of side 1 centred at the origin, 3-D; certificate lines "x y z"): checkers/certify_cube.py (A)
+ checkers/verify_exact_cube.py (B), both exact integer tests of every wall and every pair; they too say IMPROVES for any gain
(every cube record beats its bar, the largest published radius we know for that N, by more than 1e-10 relative).
The balls (hsp4 / hsp5 / hsp6: equal balls in the unit ball of dimension 4 / 5 / 6 centred at the origin; certificate lines of d
numbers): checkers/certify_ball.py (A) + checkers/verify_exact_ball.py (B), both exact integer tests of every wall and every pair;
they too say IMPROVES for any gain (every record beats its bar -- the published radius, the best spherical-code construction from
Cohn's table and any larger published size -- by more than 1e-10 relative).
A record passes only if BOTH say IMPROVES against the published Packomania radius listed in MANIFEST.csv.
"""
import csv, os, sys, zipfile, tempfile, subprocess
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "checkers"))

def container(shelf):
    if shelf == "crt": return "tri"
    if shelf == "ccq": return "quad"
    if shelf == "csc": return "semi"
    if shelf == "csq": return "square"
    if shelf == "cpt": return "poly:5"
    if shelf in ("cxd", "cpd"): return "poly:16" if shelf == "cxd" else "poly:15"
    if shelf == "scu": return "cube"
    if shelf in ("hsp4", "hsp5", "hsp6"): return "ball:" + shelf[3]
    return "rect:0." + f"{int(shelf.split('_')[1]):03d}".rstrip("0")          # crc_300 -> rect:0.3

def check(args):
    shelf, n, rec, path = args
    if shelf in ("hsp4", "hsp5", "hsp6"):   # equal balls in the unit ball of dimension d = 4, 5, 6: exact integer tests, wall and every pair
        import certify_ball
        d = int(shelf[3]); a = certify_ball.check_file(d, path, n, rec)
        argv = [sys.executable, os.path.join(HERE, "checkers", "verify_exact_ball.py"), str(d), path, str(n), rec]
        out = subprocess.run(argv, capture_output=True, text=True).stdout
        b = [l for l in out.splitlines() if l.startswith("VERDICT")]; b = b[0].split(":", 1)[1].strip() if b else "ERROR"
        return shelf, n, a, b
    if shelf == "scu":   # equal spheres in the cube of side 1 centred at the origin (3-D): exact integer tests, walls and every pair
        import certify_cube
        a = certify_cube.check_file(path, n, rec)
        argv = [sys.executable, os.path.join(HERE, "checkers", "verify_exact_cube.py"), path, str(n), rec]
        out = subprocess.run(argv, capture_output=True, text=True).stdout
        b = [l for l in out.splitlines() if l.startswith("VERDICT")]; b = b[0].split(":", 1)[1].strip() if b else "ERROR"
        return shelf, n, a, b
    if shelf in ("cxd", "cpd"):   # the regular 16-gon / 15-gon (circumradius 1, bottom side horizontal): rigorous exact wall enclosures
        import certify_kgon
        k = 16 if shelf == "cxd" else 15
        a = certify_kgon.check_file(k, path, n, rec)
        argv = [sys.executable, os.path.join(HERE, "checkers", "verify_exact_kgon.py"), str(k), path, str(n), rec]
        out = subprocess.run(argv, capture_output=True, text=True).stdout
        b = [l for l in out.splitlines() if l.startswith("VERDICT")]; b = b[0].split(":", 1)[1].strip() if b else "ERROR"
        return shelf, n, a, b
    if shelf == "cpt":   # the regular pentagon (circumradius 1, vertex (0, 1)): both checkers decide the irrational walls exactly
        import certify_poly
        a = certify_poly.check_file(5, path, n, rec)
        argv = [sys.executable, os.path.join(HERE, "checkers", "verify_exact_poly.py"), "5", path, str(n), rec]
        out = subprocess.run(argv, capture_output=True, text=True).stdout
        b = [l for l in out.splitlines() if l.startswith("VERDICT")]; b = b[0].split(":", 1)[1].strip() if b else "ERROR"
        return shelf, n, a, b
    if shelf == "csq":   # the square (up to 10,000 circles): the two O(N) exact checkers
        import certify_big
        a = certify_big.check("square", path, n, rec)
        argv = [sys.executable, os.path.join(HERE, "checkers", "verify_exact_big.py"), "square", path, str(n), rec]
        out = subprocess.run(argv, capture_output=True, text=True).stdout
        b = [l for l in out.splitlines() if l.startswith("VERDICT")]; b = b[0].split(":")[1].strip() if b else "ERROR"
        return shelf, n, a, b
    import certify_circ
    a = certify_circ.check(container(shelf), path, rec, verbose=False)[0]
    if shelf == "crt": argv = [sys.executable, os.path.join(HERE, "checkers", "verify_exact_crt.py"), path, str(n), rec]
    elif shelf == "csc": argv = [sys.executable, os.path.join(HERE, "checkers", "verify_exact_semi.py"), "semi", path, str(n), rec]
    else: argv = [sys.executable, os.path.join(HERE, "checkers", "verify_exact_rect_quad.py"), container(shelf), path, str(n), rec]
    out = subprocess.run(argv, capture_output=True, text=True).stdout
    b = [l for l in out.splitlines() if l.startswith("VERDICT")]; b = b[0].split(":")[1].strip() if b else "ERROR"
    return shelf, n, a, b

if __name__ == "__main__":
    rows = list(csv.DictReader(open(os.path.join(HERE, "MANIFEST.csv"))))
    tmp = tempfile.mkdtemp(prefix="packing_verify_")
    for z in sorted(os.listdir(os.path.join(HERE, "certificates"))):
        zipfile.ZipFile(os.path.join(HERE, "certificates", z)).extractall(os.path.join(tmp, z[:-4].split("_part")[0]))   # csq_part1.zip -> csq/
    jobs = [(r["shelf"], int(r["N"]), r["published_radius"], os.path.join(tmp, r["shelf"], f"{r['shelf']}_{r['N']}.txt")) for r in rows]
    with Pool() as pool: res = pool.map(check, jobs)
    ok = [r for r in res if r[2] == "IMPROVES" and r[3] == "IMPROVES"]
    per = {}
    for s, n, a, b in res: per.setdefault(s, [0, 0]); per[s][0] += 1; per[s][1] += (a == "IMPROVES" and b == "IMPROVES")
    print("per table (records, both checkers IMPROVES):", per)
    print(f"TOTAL: {len(ok)} / {len(res)} verified by both checkers")
    bad = [r for r in res if r not in ok]
    print("FAILURES:", bad if bad else "none")
    sys.exit(0 if not bad else 1)
