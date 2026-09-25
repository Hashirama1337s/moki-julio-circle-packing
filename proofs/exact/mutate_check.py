#!/usr/bin/env python3
# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Mutation test for replay_exact.py (Moki&Julio): each corruption of a certificate or tree must be REJECTED, and the
unmodified files must PASS.

For every case and every mutation below, it loads the shipped <case>_cert.json and <case>_tree.json.gz from this
directory, corrupts ONE thing, writes the pair to a fresh temporary directory (the files beside this script are never
modified), and runs replay_exact.run(case) on it. Most mutations re-bind the tree to the edited certificate (sha256), so
that the corruption itself must be caught and not just the broken binding. One mutation deliberately skips the re-binding.

Standard library only; runs from any working directory.
usage: py -3.11 mutate_check.py [substring ...]     (optional: run only the mutations whose name contains a substring)
exit code 0 iff every corruption is rejected and both unmodified cases pass.
"""
import sys, os, json, gzip, copy, shutil, hashlib, tempfile, io, contextlib, importlib.util
from fractions import Fraction as Fr

sys.dont_write_bytecode = True                     # leave no __pycache__ beside the certificates
HERE = os.path.dirname(os.path.abspath(__file__))
CASES = ("crc_800_4", "crc_600_6")
_spec = importlib.util.spec_from_file_location("replay_exact", os.path.join(HERE, "replay_exact.py"))
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)


def load(case):
    with open(os.path.join(HERE, f"{case}_cert.json"), "rb") as f:
        cert = json.loads(f.read())
    with gzip.open(os.path.join(HERE, f"{case}_tree.json.gz"), "rt", encoding="utf-8") as f:
        tree = json.load(f)
    return cert, tree


def write(tmp, case, cert, tree, rehash=True):
    blob = json.dumps(cert, indent=1).encode()
    with open(os.path.join(tmp, f"{case}_cert.json"), "wb") as f:
        f.write(blob)
    if rehash:
        tree["cert_sha256"] = hashlib.sha256(blob).hexdigest()
    with gzip.open(os.path.join(tmp, f"{case}_tree.json.gz"), "wt", encoding="utf-8") as f:
        json.dump(tree, f)


def first(tree, pred):
    for rec in tree["nodes"][1:]:
        if pred(rec):
            return rec
    raise RuntimeError("no node available for this mutation")


MUTS = []


def mut(name):
    def deco(f):
        MUTS.append((name, f))
        return f
    return deco


IDENTITY = "unmodified files (must PASS)"


@mut(IDENTITY)
def m_identity(case, cert, tree):
    return None


# ---- certificate (Lemma 1, section 4)
@mut("lambda entry + 1/10^60")
def m_lambda(case, cert, tree):
    a, b = cert["lambda"][0]; cert["lambda"][0] = [str(Fr(a) + Fr(1, 10 ** 60)), b]


@mut("lambda_lo above lambda")
def m_lambda_lo(case, cert, tree):
    cert["lambda_lo"][0] = str(Fr(cert["lambda_lo"][0]) * 2)


@mut("Lambda_hi below Lambda")
def m_Lambda(case, cert, tree):
    cert["Lambda_hi"] = str(Fr(cert["Lambda_hi"]) * Fr(99, 100))


@mut("rho larger than its formula")
def m_rho(case, cert, tree):
    cert["rho"] = str(Fr(cert["rho"]) * Fr(1001, 1000))


@mut("rho_prime = rho (no room for the snap)")
def m_rho_prime(case, cert, tree):
    cert["rho_prime"] = cert["rho"]


@mut("r_t above r* (p - q S2_LO)")
def m_r_t(case, cert, tree):
    p, q = Fr(cert["r_star"]["p"]), Fr(cert["r_star"]["q"])
    cert["r_t"] = str(p - q * Fr(cert["S2_LO"])); tree["r_t"] = cert["r_t"]


@mut("Minv entry + 1/10^70")
def m_minv(case, cert, tree):
    a, b = cert["Minv"][0][0]; cert["Minv"][0][0] = [str(Fr(a) + Fr(1, 10 ** 70)), b]


@mut("normMinv_hi too small")
def m_norm(case, cert, tree):
    cert["normMinv_hi"] = str(Fr(cert["normMinv_hi"]) * Fr(99, 100))


@mut("basis B with a repeated row")
def m_basis(case, cert, tree):
    cert["B"][1] = cert["B"][0]


@mut("c* centre moved by 1e-60")
def m_cstar(case, cert, tree):
    a, b = cert["c_star"][1][0]; cert["c_star"][1][0] = [str(Fr(a) + Fr(1, 10 ** 60)), b]


@mut("tight set T drops a contact")
def m_T(case, cert, tree):
    cert["T"] = cert["T"][:-1]


@mut("slack_lo above the exact slack")
def m_slack(case, cert, tree):
    k = next(iter(cert["slack_lo"])); cert["slack_lo"][k] = str(Fr(cert["slack_lo"][k]) + Fr(1, 10 ** 20))


@mut("snapped image moved by 1e-45")
def m_image(case, cert, tree):
    X, Y = cert["images_point_frame"][0]["points"][0]
    cert["images_point_frame"][0]["points"][0] = [str(Fr(X) + Fr(1, 10 ** 45)), Y]


@mut("S2_HI below sqrt2")
def m_s2(case, cert, tree):
    cert["S2_HI"] = cert["S2_LO"]


@mut("certificate edited without re-binding the tree (sha256)")
def m_bind(case, cert, tree):
    cert["authors"] = "x"
    return "nohash"


# ---- tree (section 5)
@mut("grid G + 1 (odd)")
def m_G(case, cert, tree):
    tree["G"] += 1; cert["grid_G"] += 1


@mut("root: one cell removed from the cover")
def m_cell(case, cert, tree):
    rf = tree["nodes"][0]["fate"]; rf["cells"] = rf["cells"][:-1]


@mut("root: gap in the cover (x-intervals start at 1)")
def m_gap(case, cert, tree):
    for cc in tree["nodes"][0]["fate"]["cells"]:
        if cc[0] == 0:
            cc[0] = 1


@mut("root: a cell stretched to diameter >= d")
def m_bigcell(case, cert, tree):
    rf = tree["nodes"][0]["fate"]; top = max(cc[3] for cc in rf["cells"])
    for cc in rf["cells"]:
        if cc[3] == top:
            cc[2] = 0


@mut("root: a kept start assignment deleted")
def m_startdel(case, cert, tree):
    rf = tree["nodes"][0]["fate"]; rf["children"] = rf["children"][:-1]; rf["assign"] = rf["assign"][:-1]


@mut("root: a kept start relabelled as rejected with a false witness")
def m_startwit(case, cert, tree):
    rf = tree["nodes"][0]["fate"]
    a = rf["assign"].pop(0); rf["children"].pop(0)
    rf["discards"].insert(0, [a, ["xorder", 0]])


@mut("pair witness names a pair that is not too close")
def m_pairw(case, cert, tree):
    n = tree["N"]
    rec = first(tree, lambda r: r["fate"]["t"] == "discard" and r["fate"]["w"][0] == "pair")
    w = rec["fate"]["w"]
    for a in range(n):
        for b in range(a + 1, n):
            if (a, b) != (w[1], w[2]):
                rec["fate"]["w"] = ["pair", a, b]
                return


@mut("empty witness on a non-empty interval")
def m_emptyw(case, cert, tree):
    rec = first(tree, lambda r: r["fate"]["t"] == "discard" and r["fate"]["w"][0] == "empty")
    w = rec["fate"]["w"]; rec["fate"]["w"] = ["empty", w[1], 1 - w[2]]


@mut("strip removal pushed one grid step too far")
def m_strip(case, cert, tree):
    rec = first(tree, lambda r: any(op[0] == "sh" for op in r["ops"]))
    for op in rec["ops"]:
        if op[0] == "sh":
            op[4] += 1 if op[3] in ("x0", "y0") else -1
            return


@mut("non-existent cut operation (e.g. a strict y-cut)")
def m_fakeop(case, cert, tree):
    rec = first(tree, lambda r: len(r["ops"]) > 0)
    rec["ops"].append(["y1_strict"])


@mut("order clip with an out-of-range index")
def m_ohi(case, cert, tree):
    rec = first(tree, lambda r: len(r["ops"]) > 0)
    rec["ops"].insert(0, ["o_hi", tree["N"] - 1])


@mut("acceptance with a wrong permutation")
def m_perm(case, cert, tree):
    rec = first(tree, lambda r: r["fate"]["t"] == "accept")
    p = rec["fate"]["perm"]; p[0], p[1] = p[1], p[0]


@mut("acceptance with a wrong mirror")
def m_mirror(case, cert, tree):
    rec = first(tree, lambda r: r["fate"]["t"] == "accept")
    rec["fate"]["mirror"] = (rec["fate"]["mirror"] + 1) % 4


@mut("split drops a child")
def m_child(case, cert, tree):
    rec = first(tree, lambda r: r["fate"]["t"] == "split")
    rec["fate"]["children"] = rec["fate"]["children"][:1]


@mut("bisection point shifted")
def m_mid(case, cert, tree):
    rec = first(tree, lambda r: r["fate"]["t"] == "split")
    rec["fate"]["how"][3] += 1


@mut("leaf marked undecided")
def m_undecided(case, cert, tree):
    rec = first(tree, lambda r: r["fate"]["t"] == "discard")
    rec["fate"] = {"t": "undecided"}


@mut("node box edited")
def m_box(case, cert, tree):
    tree["nodes"][5]["box"][0][1] -= 1


@mut("tightened box edited")
def m_tbox(case, cert, tree):
    rec = first(tree, lambda r: "tbox" in r)
    rec["tbox"][0][0] += 1


@mut("orphan node appended")
def m_orphan(case, cert, tree):
    rec = copy.deepcopy(tree["nodes"][-1]); rec["id"] = len(tree["nodes"]); tree["nodes"].append(rec)


def main(only):
    tmp = tempfile.mkdtemp(prefix="mutate_check_")
    R.DIR = tmp
    rejected = total = 0; identity_ok = True; bad = []
    try:
        for case in CASES:
            for name, f in MUTS:
                if only and name != IDENTITY and not any(o in name for o in only):
                    continue
                cert, tree = load(case)
                res = f(case, cert, tree)
                write(tmp, case, cert, tree, rehash=(res != "nohash"))
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    passed = R.run(case)
                last = buf.getvalue().strip().splitlines()[-1]
                if name == IDENTITY:
                    good = passed; identity_ok &= passed
                else:
                    good = not passed; total += 1; rejected += (not passed)
                if not good:
                    bad.append(f"{case}: {name}")
                print(f"{'ok ' if good else 'BAD'} {case}  {name:58s} -> {last[:140]}", flush=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"SUMMARY: {rejected}/{total} corruptions rejected ({total // len(CASES)} per case); unmodified files "
          f"{'PASS' if identity_ok else 'FAIL'} on both cases" + ("" if not bad else f"; NOT OK: {bad}"))
    return 0 if (not bad and identity_ok) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
