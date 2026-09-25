"""Grok's independent csq constructions (grok/csq_constructions.py: 'stretched' and 'mixed' row lattices, written by Grok 09-24) on
EVERY listed csq N > 1000: keep a construction only if it beats our current best for N (cand_big_hp) AND passes the csq gate
(Packomania's monotone envelope + Amore). Certificate as in hats_csq_polish (radius = clearance * (1 - 1e-12), rounded down) and
BOTH O(N) exact checkers. -> cand_big_hp/csq/csq_<N>.txt + .json (tag csq_grok)   usage: py -3.11 hats_csq_grok.py
"""
import os, sys, json, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "grok"))
import numpy as np, mpmath as mp
from decimal import Decimal, ROUND_FLOOR, getcontext
import csq_constructions as G, lai_compare as L, certify_big, hats_flagship as HF
getcontext().prec = 60
T = HF.table("csq"); T2, T7 = L.load(); out = open(os.path.join(HERE, "out", "hats_csq_grok.jsonl"), "a"); kept = 0
for n in sorted((k for k in T if k > 1000), reverse=True):          # descending: hats_csq_polish owns N <= 2000 until it finishes
    best = (0.0, None, None)
    try:
        for name, c, r in G.candidates(n):
            if r > best[0]: best = (r, name, c)
    except Exception as e:
        out.write(json.dumps({"N": n, "error": repr(e)}) + "\n"); continue
    if best[2] is None: continue
    r_claim = (Decimal(repr(float(best[0]))) * (1 - Decimal(10) ** -12)).quantize(Decimal(10) ** -25, rounding=ROUND_FLOOR)
    f = os.path.join(HERE, "cand_big_hp", "csq", f"csq_{n}.txt")
    ours = Decimal(open(f).readline().split()[1]) if os.path.exists(f) else Decimal(0)
    ok, why = L.gate("csq", n, mp.mpf(str(r_claim)), T2, T7)
    row = {"N": n, "name": best[1], "r_claim": str(r_claim), "vs_table": float(r_claim / Decimal(T[n]) - 1), "gate": why or "claimable",
           "vs_ours": float(r_claim / ours - 1) if ours else None}
    if ok and r_claim > ours:
        tmp = f + ".grok.tmp"; c = best[2]
        open(tmp, "w", newline="\n").write(f"r {r_claim}\n" + "".join(f"{format(Decimal(repr(float(x))), 'f')} {format(Decimal(repr(float(y))), 'f')}\n" for x, y in c))
        a = certify_big.check("square", tmp, n, T[n])
        o = subprocess.run([sys.executable, os.path.join(HERE, "grok", "verify_exact_big.py"), "square", tmp, str(n), T[n]], capture_output=True, text=True).stdout
        b = [l for l in o.splitlines() if l.startswith("VERDICT")]; b = b[0].split(":")[1].strip() if b else "ERROR"
        row.update(claude=a, grok=b)
        if a == "IMPROVES" and b == "IMPROVES":
            os.replace(tmp, f); kept += 1; row["kept"] = True
            json.dump({"shelf": "csq", "N": n, "tag": "csq_grok", "how": f"Grok's construction {best[1]} (grok/csq_constructions.py)",
                       "r_new": str(r_claim), "claude": a, "grok": b, "lopt": "NOT_ATTEMPTED"}, open(f[:-4] + ".json", "w"), indent=1)
        else: os.remove(tmp)
    out.write(json.dumps(row) + "\n"); out.flush()
print("DONE kept", kept, flush=True)
