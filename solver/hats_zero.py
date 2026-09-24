"""Four hats, zero-cost checks (vision/CLAUDE-SEALED-HATS.md): pure arithmetic on best-known radii, no search.
best(s, N) = max(Packomania snapshot, our certified records, cand_hp). Float screen; anything within 1e-12 is a tie, and any
real violation is re-checked exactly before it counts.
 TAO   T0  monotone r(N) >= r(N+k); height-monotone r_crc_b(N) <= r_crc_k(N) for b < k; stacking crc_a + crc_b -> crc_(a+b);
           semicircle = 2 quarter discs; triangle = 2 half triangles (altitude, scale 1/sqrt2) and corner square + 2 triangles (1/2)
 EINSTEIN E0 outbound mirrors: ccq(N) -> cci(4N); csc(N) -> cci(2N); two semicircles -> cci; crt(N) -> csq(2N);
           crc_k(N) -> crc_mk(mN); crc_a + crc_b -> csq (a + b = 1000)
Writes out/hats_zero.json; prints violations and the closest misses per check.
"""
import os, re, json, glob, math, collections
H = os.path.dirname(os.path.abspath(__file__)); P = lambda *a: os.path.join(H, *a)
TABLES = ["crt", "ccq", "csc"] + [f"crc_{k}" for k in range(100, 900, 100)]
def snap(s):
    p = P("data", "radius.txt") if s == "crt" else P("shelves", s, "radius.txt")
    return {int(l.split()[0]): float(l.split()[1]) for l in open(p) if l.strip()}
def live(s):
    t = open(P("data", "refs", f"packomania_{s}_2026-09-24.html"), encoding="utf-8", errors="ignore").read()
    return {int(n): float(v) for n, v in re.findall(r'name="[a-z]+(\d+)">\s*<strong>\s*\d+</strong></a></td>\s*<td>(?:<strong>)?\s*([0-9.]+)', t)}
pub = {s: snap(s) for s in TABLES}; pub["cci"] = live("cci"); pub["csq"] = live("csq")
ours = collections.defaultdict(dict)
for r in json.load(open(P("out", "final_records.json"))):
    ours[r["shelf"]][r["n"]] = max(ours[r["shelf"]].get(r["n"], 0.0), float(r["r_new"]))
for f in glob.glob(P("cand_hp", "*", "*.txt")):
    s = os.path.basename(os.path.dirname(f)); n = int(os.path.basename(f)[:-4].split("_")[-1])
    try: v = float(open(f).readline().split()[1])
    except Exception: continue
    ours[s][n] = max(ours[s].get(n, 0.0), v)
best = {s: {n: max(pub[s].get(n, 0.0), ours[s].get(n, 0.0)) for n in set(pub[s]) | set(ours.get(s, {}))} for s in pub}
INF = float("inf")
def b(s, n): return INF if n == 0 else best[s].get(n, 0.0)
TOL = 1e-12
out = collections.defaultdict(list)          # check -> [(ratio constructed/best, table, N, how)]
def test(check, s, n, val, how):
    ref = best[s].get(n)
    if ref and val < INF: out[check].append((val / ref - 1, s, n, how))

# TAO: monotone in N (delete k circles from N+k)
for s in TABLES + ["cci", "csq"]:
    ns = sorted(best[s]); run = 0.0
    for n in reversed(ns):                       # run = max best(m) over m > n
        if run > 0: test("T0 monotone N", s, n, run, "delete circles from a larger N")
        run = max(run, best[s][n])
# TAO: height-monotone + stacking of rectangles (N1 = 0 is height monotonicity)
K = list(range(100, 900, 100))
for k in K:
    tgt = f"crc_{k}"
    for n in sorted(best[tgt]):
        cand = []
        for a in K:
            if a >= k: continue
            bb = k - a
            if bb in K:
                for n1 in range(0, n + 1):
                    v = min(b(f"crc_{a}", n1), b(f"crc_{bb}", n - n1))
                    cand.append((v, f"crc_{a}({n1}) + crc_{bb}({n - n1})"))
            cand.append((b(f"crc_{a}", n), f"crc_{a}({n}) inside"))
        if cand:
            v, how = max(cand); test("T0 rect stack/height", tgt, n, v, how)
# TAO: semicircle = two quarter discs
Hs = {}
for n in sorted(best["csc"]):
    if n < 2: continue
    v, how = max((min(b("ccq", n1), b("ccq", n - n1)), f"ccq({n1}) + ccq({n - n1})") for n1 in range(1, n // 2 + 1))
    test("T0 csc = 2 ccq", "csc", n, v, how)
# TAO: triangle = 2 half triangles (altitude) ; corner square + 2 triangles at scale 1/2
for n in sorted(best["crt"]):
    if n < 2: continue
    v, how = max((min(b("crt", n1), b("crt", n - n1)) / math.sqrt(2), f"crt({n1}) + crt({n - n1}) at 1/sqrt2") for n1 in range(1, n // 2 + 1))
    test("T0 crt = 2 crt", "crt", n, v, how)
    cq = []
    for a in range(0, n + 1):
        sq = b("csq", a) / 2
        for n1 in range(0, (n - a) // 2 + 1):
            cq.append((min(sq, b("crt", n1) / 2, b("crt", n - a - n1) / 2), f"csq({a}) + crt({n1}) + crt({n - a - n1}) at 1/2"))
    v, how = max(cq); test("T0 crt = sq + 2 crt", "crt", n, v, how)
# EINSTEIN: outbound mirrors into other tables
for n, v in best["ccq"].items(): test("E0 ccq -> cci(4N)", "cci", 4 * n, v, f"ccq({n}) x4 mirror")
semi = {n: max(best["csc"].get(n, 0.0), max((min(b("ccq", n1), b("ccq", n - n1)) for n1 in range(1, n // 2 + 1)), default=0.0))
        for n in set(best["csc"]) | {n1 + n2 for n1 in best["ccq"] for n2 in best["ccq"] if n1 + n2 <= 900}}
for n, v in best["csc"].items(): test("E0 csc -> cci(2N)", "cci", 2 * n, v, f"csc({n}) x2 mirror")
for N in sorted(best["cci"]):
    if N > 1800: break
    c = [(min(semi.get(n1, 0.0), semi.get(N - n1, 0.0)), f"semi({n1}) + semi({N - n1})") for n1 in range(1, N // 2 + 1)]
    if c: v, how = max(c); test("E0 2 semicircles -> cci", "cci", N, v, how)
for n, v in best["crt"].items(): test("E0 crt -> csq(2N)", "csq", 2 * n, v, f"crt({n}) mirrored in the hypotenuse")
for k in K:
    for m in range(2, 9):
        if k * m in K:
            for n, v in best[f"crc_{k}"].items(): test("E0 crc_k -> crc_mk", f"crc_{k * m}", m * n, v, f"crc_{k}({n}) x{m}")
for a in K:
    bb = 1000 - a
    if bb in K and a <= bb:
        for N in sorted(best["csq"]):
            if N > 1200: break
            c = [(min(b(f"crc_{a}", n1), b(f"crc_{bb}", N - n1)), f"crc_{a}({n1}) + crc_{bb}({N - n1})") for n1 in range(1, N)]
            if c: v, how = max(c); test("E0 crc_a + crc_b -> csq", "csq", N, v, how)
res = {}
for check, L in out.items():
    L.sort(reverse=True); viol = [x for x in L if x[0] > TOL]
    res[check] = {"tested": len(L), "violations": viol[:50], "n_violations": len(viol), "closest": L[:8]}
    print(f"\n{check}: tested {len(L)}, violations {len(viol)}")
    for x in (viol[:10] if viol else L[:5]): print(f"   {x[0]:+.3e}  {x[1]} N={x[2]}  <- {x[3]}")
json.dump(res, open(P("out", "hats_zero.json"), "w"), indent=1)
