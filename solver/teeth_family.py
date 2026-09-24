"""recon4 (Moki&Julio, 2026-09-24): Molnar 'teeth' = Furedi 1991 (DCG 6:95-106) Example 1.1.
(1) C2 as a universal statement is false: explicit teeth packings beat two straight rows for large N.
(2) Teeth as a finite family for 1 x h rectangles: ties and wins vs best known (table + our RECORDS).
Units in (1): D = 1 (min distance), strip height tau.  Float checks only; run the exact checkers before claiming.
Standard library only."""
import math, re, os, itertools
S3 = math.sqrt(3); W = 8
ROOT = '.'   # run from the folder holding shelves/

# ---------- (1) crossover: two straight rows vs teeth (n = 1 blocks), even N ----------
def place(ys):                      # exact min x-extent for a fixed height sequence in x-order
    xs = [0.0]
    for j in range(1, len(ys)):
        xj = xs[-1]
        for i in range(max(0, j - W), j):
            dy = ys[j] - ys[i]
            if abs(dy) < 1: xj = max(xj, xs[i] + math.sqrt(1 - dy * dy))
        xs.append(xj)
    return xs
def mind(xs, ys):
    return min(math.hypot(xs[i] - xs[j], ys[i] - ys[j])
               for i in range(len(xs)) for j in range(i + 1, min(len(xs), i + 2 * W)))
def word(N, tau, st):
    out, b = [], st
    while len(out) < N:
        out += [0, S3 / 2, 0] if b == 0 else [tau, tau - S3 / 2, tau]; b ^= 1
    return out[:N]
def alt(n, tau, p): return [(0 if (i + p) % 2 == 0 else tau) for i in range(n)]
def best_teeth_extent(N, tau):
    best = 1e9
    for pre in range(5):
        for suf in range(5):
            if N - pre - suf < 0: continue
            for st in (0, 1):
                for p0 in (0, 1):
                    for s0 in (0, 1):
                        ys = alt(pre, tau, p0) + word(N - pre - suf, tau, st) + alt(suf, tau, s0)
                        xs = place(ys)
                        if xs[-1] < best and mind(xs, ys) > 1 - 1e-12: best = xs[-1]
    return best
def crossover(tau):
    s = math.sqrt(1 - tau * tau)
    for m in range(2, 400):
        e2, et = (m - 1) + s, best_teeth_extent(2 * m, tau)
        if et < e2 - 1e-12: return 2 * m, e2, et
    return None

# ---------- (2) teeth (general n) as a family in 1 x h ----------
def teeth_sel(N, tau):
    n = int(math.floor(2 * tau / S3 + 1e-12))
    if n < 1: return None
    a = (n + math.sqrt(max(0.0, 4 - (2 * tau - n * S3) ** 2))) / 2
    blk = [(j / 2 + i, j * S3 / 2) for j in range(n + 1) for i in range(n + 1 - j)]
    pts, k = [], 0
    while len(pts) < N + len(blk):
        for (x, y) in blk: pts.append((k * a + x, y if k % 2 == 0 else tau - y))
        k += 1
    pts.sort(); best = None
    for st in range(len(blk)):
        sel = pts[st:st + N]
        if len(sel) == N and (best is None or sel[-1][0] - sel[0][0] < best[0]):
            best = (sel[-1][0] - sel[0][0], sel)
    e, sel = best
    md = min(math.hypot(p[0] - q[0], p[1] - q[1]) for p, q in itertools.combinations(sel, 2))
    return e, md, sel
def teeth_r(N, h):
    lo, hi, ok = 1e-6, h / 2, None
    for _ in range(80):
        r = (lo + hi) / 2; D = 2 * r; t = teeth_sel(N, (h - D) / D)
        if t is None or t[1] < 1 - 1e-9: hi = r; continue
        if (t[0] + 1) * D <= 1: lo = ok = r
        else: hi = r
    return ok
def best_known(t):
    d = {}
    for l in open(f'{ROOT}/shelves/{t}/radius.txt'):
        p = l.split()
        if len(p) >= 2: d[int(p[0])] = float(p[1])
    fn = f'{ROOT}/RECORDS_{t}.md'
    if os.path.exists(fn):
        for l in open(fn, encoding='utf8'):
            m = re.match(r'\|\s*(\d+)\s*\|\s*([0-9.]+)\s*\|', l)
            if m: d[int(m.group(1))] = max(d.get(int(m.group(1)), 0), float(m.group(2)))
    return d

if __name__ == '__main__':
    for name, tau in [('crc_200 N18', 0.891), ('crc_800 N4', 0.918), ('crc_300 N12', 0.919),
                      ('crc_100 N38', 0.935), ('crc_600 N6', 0.963), ('tau=1 edge cases', 1.0)]:
        print(f'{name:18s} tau={tau}: first even N where teeth < two rows ->', crossover(tau))
    ties, wins = {}, []
    for h in (100, 200, 300, 400, 500, 600, 700, 800):
        t = f'crc_{h}'; H = h / 1000; bk = best_known(t)
        for N in sorted(bk):
            r0 = bk[N]; tau0 = (H - 2 * r0) / (2 * r0)
            if not (0.87 < tau0 < 3 * S3 / 2 - 0.01): continue
            r = teeth_r(N, H)
            if r is None: continue
            if r > r0 * (1 + 1e-9): wins.append((t, N, round(tau0, 3), r, r0, (r - r0) / r0))
            elif abs(r - r0) / r0 < 1e-9: ties.setdefault(t, []).append(N)
    print('entries equal to teeth (rel 1e-9):', ties)
    print('teeth beats best known:', wins)
