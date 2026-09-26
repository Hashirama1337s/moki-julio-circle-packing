# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Recover the missing (x5, x6) of Packomania's TRUNCATED hsp6 files by RIGID GROWTH on contact coincidences (2026-09-26;
Moki&Julio). Float only; the deliverable is the completed packing polished by slpd.polish.

Setting (recon56.py): the file gives a_i = (x1..x4) of every centre; the missing y_i = (x5, x6) must satisfy
  |y_i| <= m_i        (m_i^2 = (1 - r)^2 - |a_i|^2; equality = the ball touches the wall: all but the 0-2 'core' balls)
  |y_i - y_j| >= g_ij (g_ij^2 = 4 r^2 - |a_i - a_j|^2, pairs with g_ij^2 > 0).
A CONTACT of the 6-D packing is an equality in the second family, a wall contact an equality in the first. The packings are
jammed (~10 contacts per ball), so in the (x5, x6) plane every non-rattler ball sits on >= 3 of its 'circles' (the wall
circle |y| = m_i and the circles |y - y_j| = g_ij around its placed neighbours) AT ONCE. Two circles through a point is
what any intersection gives; a THIRD circle through the same point (to ~1e-9) is a coincidence that essentially never happens by
chance (~1e-8 per candidate). So:
  seed   two balls a, b: y_a = (m_a, 0) (fixes the rotation of the plane), y_b on its wall circle at the angle theta_ab from a
         (fixes the reflection) - seeds ranked by how many common neighbours close a triangle a-b-c in angle;
  grow   repeatedly: for every unplaced ball u with placed neighbours, candidates = all intersections of its circles (wall x
         neighbour and neighbour x neighbour, the latter for an interior ball); support = number of circles through the point
         (tol 1e-8 in distance, circles at the EXACT printed radius); feasible = every constraint at the loosened radius
         (recon56.needs: r - 1e-9, 12-decimal allowance) holds to 1e-8. Place the ball whose best feasible point has the
         highest support >= 3 and is unambiguous (no second distinct feasible point with the same support); the point is
         refined by Gauss-Newton on its supporting circles (so errors do not accumulate along the growth);
  rest   balls never reached with support >= 3 are the rattlers (and a central core ball): each goes, most-constrained first,
         to the candidate point (intersections + 2,000 random disc points) that maximises its worst slack;
  check  L-BFGS on the Cartesian penalty (recon56.penalty) with the four known columns fixed; worst violation, geomd.rmin.
If the growth from the best seed leaves a violation, the next of 8 ranked seeds is tried; if none succeeds, the LADDER:
GrowerBT (the same growth with backtracking on ambiguous balls and a dead-end test) over up to 300 / 300 / 1000 seeds at
coincidence tolerance 1e-8 / 1e-7 / 1e-6. The completion is then released to slpd.polish (all 6 coordinates free): a coarse
polish (60 s cpu) and a fine one (delta0 = 1e-6, 30 s cpu; an exact completion at printed - 4e-12 reads as stationary to the
coarse LP), the better kept.

usage: python recover6.py --cells 30,66,100 | --all [--nmin 2] [--polish 60] [--workers 3]
Writes out/recon56/hsp6_<N>.npy (the POLISHED packing) and out/recon56/recover6.json (per N: printed, recovered radius,
polished radius, method, seconds; log rows in recover6.jsonl).
RESULT (2026-09-26, all of N = 2..250; N = 1 is trivial): 249 / 249 recovered and validated - completed radius >= printed - 1e-8
(violation 0; worst completed gap -1.01e-9 = the loosening) AND polished >= printed - 1e-12 (worst -5.3e-13), re-checked from
the saved .npy files. 246 cells by the plain growth from the top 8 seeds (median 1.9 s cpu, max 24 s), N = 68, 70 by GrowerBT
at ts = 1e-8 (20-28 seeds), N = 67 at ts = 1e-6 (near-contacts with ~1e-9 gaps, 40 rattlers). All 95 bar_source = 'page'
cells (N >= 2) recovered. 112 polished packings exceed printed + 2e-12 (39 of them page cells; largest N = 139 +3.3e-6,
N = 158 +9.9e-7, N = 250 +1.3e-7) - candidates for the sweep / certification, NOT claims.
"""
import os, sys, time, json, argparse
import numpy as np
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prio                                                      # noqa: E402
import geomd, slpd                                               # noqa: E402
from recon56 import trunc, needs, penalty, viol_of               # noqa: E402
OUT = os.path.join(geomd.OUT, 'recon56')
SUMMARY = os.path.join(OUT, 'recover6.json')
LOG = os.path.join(OUT, 'recover6.jsonl')


class Prob:
    """The circles of every ball: exact values (printed r, no allowance) for coincidences, loosened ones for feasibility."""
    def __init__(self, A, r):
        n = len(A); self.n = n; self.r = r
        self.m_ex = np.sqrt(np.maximum((1 - r) ** 2 - (A * A).sum(1), 0.0))
        self.m_lo = needs(A, r)[0]
        iu = np.triu_indices(n, 1); d2 = ((A[iu[0]] - A[iu[1]]) ** 2).sum(1)
        g2 = 4 * r * r - d2; k = g2 > 0
        I, J = iu[0][k], iu[1][k]; gex = np.sqrt(g2[k])
        re_ = r - 1e-9; glo = np.sqrt(np.maximum(4 * re_ * re_ - d2[k] - 1e-11, 0.0))
        self.I, self.J, self.gex, self.glo = I, J, gex, glo
        nb = [[] for _ in range(n)]
        for p, (i, j) in enumerate(zip(I, J)): nb[i].append((j, p)); nb[j].append((i, p))
        self.nbj = [np.array([j for j, _ in L], dtype=int) for L in nb]
        self.nbp = [np.array([p for _, p in L], dtype=int) for L in nb]
        # angular data for seeding (wall-wall pairs)
        mi, mj = self.m_ex[I], self.m_ex[J]
        kap = (mi ** 2 + mj ** 2 - gex ** 2) / np.maximum(2 * mi * mj, 1e-300)
        self.kap = kap; self.th = np.arccos(np.clip(kap, -1, 1))
        # loosened constraint set for the final penalty / violation
        self.L = needs(A, r)


def circ_int(P1, r1, P2, r2):
    """Intersections of circles (P1, r1) and (P2, r2), vectorised over rows; near-tangent pairs give the tangent point.
    Returns points (2k, 2) and a validity mask (2k,)."""
    P1 = np.broadcast_to(P1, P2.shape) if np.ndim(P1) == 1 else P1
    D = P2 - P1; d0 = np.sqrt((D * D).sum(1)); d = np.maximum(d0, 1e-12)
    a = (r1 * r1 - r2 * r2 + d * d) / (2 * d); h2 = r1 * r1 - a * a
    ok = (h2 > -1e-14 * np.maximum(r1 * r1, 1e-30) - 1e-18) & (d0 > 1e-12)
    h = np.sqrt(np.maximum(h2, 0.0))
    ex = D / d[:, None]; B = P1 + a[:, None] * ex; perp = np.stack([-ex[:, 1], ex[:, 0]], 1)
    pts = np.concatenate([B + h[:, None] * perp, B - h[:, None] * perp])
    return pts, np.r_[ok, ok]


class Grower:
    def __init__(self, P, ts=1e-8, tv=1e-8, distinct=1e-6):
        self.P = P; self.ts, self.tv, self.distinct = ts, tv, distinct
        n = P.n; self.Y = np.full((n, 2), np.nan); self.placed = np.zeros(n, bool); self.cache = {}
        self.order = []

    def place(self, u, y):
        self.Y[u] = y; self.placed[u] = True; self.order.append(u)
        self.cache.pop(u, None)
        for j in self.P.nbj[u]:
            if not self.placed[j]: self.cache.pop(int(j), None)

    def eval(self, u, interior=True):
        """(best support, margin, point, n_supporting) of ball u's feasible candidate points; None if < 1 placed neighbour."""
        if u in self.cache: return self.cache[u]
        P = self.P; nj = P.nbj[u]; pp = P.nbp[u]; mk = self.placed[nj]
        if mk.sum() < 1: self.cache[u] = None; return None
        pj = nj[mk]; Q = self.Y[pj]; gx = P.gex[pp[mk]]; gl = P.glo[pp[mk]]; k = len(pj)
        pts, ok = circ_int(np.zeros(2), np.full(k, P.m_ex[u]), Q, gx)
        pts = pts[ok]
        if interior and k >= 2:
            a, b = np.triu_indices(k, 1)
            p2, ok2 = circ_int(Q[a], gx[a], Q[b], gx[b]); pts = np.concatenate([pts, p2[ok2]])
        if not len(pts): self.cache[u] = (0, 0, None, 0); return self.cache[u]
        nr = np.sqrt((pts * pts).sum(1))
        dist = np.sqrt(((pts[:, None, :] - Q[None]) ** 2).sum(-1))
        res_n = np.abs(dist - gx[None]) < self.ts
        sup = res_n.sum(1) + (np.abs(nr - P.m_ex[u]) < self.ts)
        slack = np.minimum(P.m_lo[u] - nr, (dist - gl[None]).min(1))
        feas = slack >= -self.tv
        if not feas.any(): self.cache[u] = (0, 0, None, 0); return self.cache[u]
        pts, sup, res_n = pts[feas], sup[feas], res_n[feas]
        o = np.argsort(-sup); b0 = o[0]; best = int(sup[b0])
        far = np.sqrt(((pts - pts[b0]) ** 2).sum(1)) > self.distinct
        second = int(sup[far].max()) if far.any() else 2
        y = self.refine(u, pts[b0], Q[res_n[b0]], gx[res_n[b0]], abs(np.sqrt((pts[b0] ** 2).sum()) - P.m_ex[u]) < self.ts)
        self.cache[u] = (best, best - second, y, best)
        return self.cache[u]

    def refine(self, u, y, Q, g, wall, it=3):
        """Gauss-Newton on the supporting circles (least squares of |y - q| - g and |y| - m)."""
        m = self.P.m_ex[u]
        for _ in range(it):
            D = y - Q; d = np.sqrt((D * D).sum(1)); Jm = D / d[:, None]; res = d - g
            if wall:
                ny = np.sqrt((y * y).sum()); Jm = np.vstack([Jm, y / ny]); res = np.r_[res, ny - m]
            if len(res) < 2: break
            step, *_ = np.linalg.lstsq(Jm, -res, rcond=None); y = y + step
            if np.abs(step).max() < 1e-15: break
        return y

    def grow(self, t_cap=600.0, t0=None):
        t0 = time.process_time() if t0 is None else t0
        n = self.P.n
        while True:
            best = None
            for u in range(n):
                if self.placed[u]: continue
                e = self.eval(u)
                if e is None or e[2] is None or e[0] < 3: continue
                key = (e[1] > 0, e[0], e[1])
                if best is None or key > best[0]: best = (key, u, e)
            if best is None or not best[0][0]: return best
            self.place(best[1], best[2][2])
            if time.process_time() - t0 > t_cap: return 'time'


class GrowerBT(Grower):
    """Rigid growth WITH BACKTRACKING (fallback for the symmetric cells near the E6 kissing configuration, N = 67-72, where
    the plain growth either stalls on mirror-ambiguous balls or commits to a structural (non-chance) coincidence). Each step:
    every unplaced ball's DISTINCT feasible points with support >= 3 (clusters at 1e-6); a ball with >= 2 placed neighbours
    and NO feasible intersection point at all is a contradiction -> backtrack. Otherwise branch on the ball with the fewest
    alternatives (then the highest support), alternatives in support order."""

    def alts(self, u):
        if u in self.cache: return self.cache[u]
        P = self.P; nj = P.nbj[u]; pp = P.nbp[u]; mk = self.placed[nj]; k = int(mk.sum())
        if k < 1: self.cache[u] = ('free', []); return self.cache[u]
        pj = nj[mk]; Q = self.Y[pj]; gx = P.gex[pp[mk]]; gl = P.glo[pp[mk]]
        pts, ok = circ_int(np.zeros(2), np.full(k, P.m_ex[u]), Q, gx); pts = pts[ok]
        if k >= 2:
            a, b = np.triu_indices(k, 1)
            p2, ok2 = circ_int(Q[a], gx[a], Q[b], gx[b]); pts = np.concatenate([pts, p2[ok2]])
        th = np.linspace(0, 2 * np.pi, 64, endpoint=False)
        pts = np.concatenate([pts, P.m_lo[u] * np.stack([np.cos(th), np.sin(th)], 1), np.zeros((1, 2))])
        nr = np.sqrt((pts * pts).sum(1)); dist = np.sqrt(((pts[:, None, :] - Q[None]) ** 2).sum(-1))
        res_n = np.abs(dist - gx[None]) < self.ts
        wall = np.abs(nr - P.m_ex[u]) < self.ts
        sup = res_n.sum(1) + wall
        slack = np.minimum(P.m_lo[u] - nr, (dist - gl[None]).min(1))
        feas = slack >= -self.tv
        if not feas.any():
            self.cache[u] = ('dead' if k >= 2 else 'free', []); return self.cache[u]
        f = np.flatnonzero(feas & (sup >= 3))
        f = f[np.argsort(-sup[f], kind='stable')]
        out = []
        for i in f:
            if any(np.sqrt(((pts[i] - y) ** 2).sum()) < self.distinct for _, y in out): continue
            out.append((int(sup[i]), self.refine(u, pts[i], Q[res_n[i]], gx[res_n[i]], bool(wall[i]))))
            if len(out) >= 6: break
        self.cache[u] = ('ok', out); return self.cache[u]

    def unplace_to(self, L):
        for u in self.order[L:]:
            self.placed[u] = False; self.Y[u] = np.nan
        del self.order[L:]; self.cache = {}

    def grow_bt(self, t_cap=300.0, node_cap=20000):
        t0 = time.process_time(); n = self.P.n; stack = []; nodes = 0; best = (-1, None); dead_ends = 0
        while nodes < node_cap and time.process_time() - t0 < t_cap:
            nodes += 1; choice = None; dead = False
            for u in range(n):
                if self.placed[u]: continue
                st, out = self.alts(u)
                if st == 'dead': dead = True; break
                if not out: continue
                key = (-len(out), out[0][0])
                if choice is None or key > choice[0]: choice = (key, u, out)
            if not dead and choice is None:
                return dict(done=True, nodes=nodes, dead_ends=dead_ends)
            if not dead:
                L = len(self.order); stack.append((choice[1], choice[2], 0, L)); self.place(choice[1], choice[2][0][1])
                if len(self.order) > best[0]: best = (len(self.order), list(self.order), self.Y.copy())
                continue
            dead_ends += 1                                        # backtrack to the last decision with an untried alternative
            while stack:
                u, out, i, L = stack.pop(); self.unplace_to(L)
                if i + 1 < len(out): stack.append((u, out, i + 1, L)); self.place(u, out[i + 1][1]); break
            if not stack: break
        if best[1] is not None:                                   # budget out / exhausted: keep the deepest state
            self.unplace_to(0)
            for u in best[1]: self.place(u, best[2][u])
        return dict(done=False, nodes=nodes, dead_ends=dead_ends)


def seeds(P, top=12, tol=1e-7):
    """Seed pairs (a, b), ranked by the number of common neighbours c closing the triangle in angle (some signs s, t with
    theta_ab + t theta_bc - s theta_ac = 0 mod 2 pi to tol)."""
    n = P.n; T = {}
    for p, (i, j) in enumerate(zip(P.I, P.J)): T[(i, j)] = p; T[(j, i)] = p
    act = (P.kap < 1) & (P.kap > -1)
    nbs = [set() for _ in range(n)]
    for p, (i, j) in enumerate(zip(P.I, P.J)):
        if act[p]: nbs[i].add(j); nbs[j].add(i)
    sc = []
    for p, (a, b) in enumerate(zip(P.I, P.J)):
        if not act[p]: continue
        tab = P.th[p]; cnt = 0
        for c in nbs[a] & nbs[b]:
            tac = P.th[T[(a, c)]]; tbc = P.th[T[(b, c)]]
            for s in (1, -1):
                for t in (1, -1):
                    v = (tab + t * tbc - s * tac) % (2 * np.pi)
                    if min(v, 2 * np.pi - v) < tol: cnt += 1
        sc.append((cnt, p))
    sc.sort(reverse=True)
    return [(int(P.I[p]), int(P.J[p]), c) for c, p in sc[:top]]


def place_rest(G, rng, nrand=2000):
    """Unplaced balls (rattlers / a free core ball), most placed neighbours first: the candidate point (all circle
    intersections + random disc points) maximising the worst slack at the loosened radius."""
    P = G.P; n = P.n
    while not G.placed.all():
        U = np.flatnonzero(~G.placed)
        u = int(U[np.argmax([G.placed[P.nbj[x]].sum() + 1e-3 * P.m_ex[x] for x in U])])
        nj = P.nbj[u]; mk = G.placed[nj]; pj = nj[mk]; Q = G.Y[pj]; gl = P.glo[P.nbp[u][mk]]; gx = P.gex[P.nbp[u][mk]]
        k = len(pj); m = P.m_lo[u]
        c = rng.uniform(-1, 1, (nrand, 2)); c = c[(c * c).sum(1) <= 1] * m
        cands = [c, np.zeros((1, 2))]
        if k:
            pts, ok = circ_int(np.zeros(2), np.full(k, m), Q, gl); cands.append(pts[ok])
            if k >= 2:
                a, b = np.triu_indices(k, 1); p2, ok2 = circ_int(Q[a], gl[a], Q[b], gl[b]); cands.append(p2[ok2])
        C = np.concatenate(cands)
        sl = m - np.sqrt((C * C).sum(1))
        if k: sl = np.minimum(sl, (np.sqrt(((C[:, None] - Q[None]) ** 2).sum(-1)) - gl[None]).min(1))
        G.place(u, C[int(np.argmax(sl))])
    return G


def finish(P, Y):
    """L-BFGS on the Cartesian penalty (loosened constraints), first 4 columns fixed; returns (Y, violation)."""
    m, I, J, G = P.L
    res = minimize(penalty, Y.ravel(), args=(m, I, J, G, 2), jac=True, method='L-BFGS-B',
                   options=dict(maxiter=20000, ftol=1e-30, gtol=1e-18))
    Y2 = res.x.reshape(-1, 2)
    v0, v2 = viol_of(Y, m, I, J, G), viol_of(Y2, m, I, J, G)
    return (Y2, v2) if v2 <= v0 else (Y, v0)


# fallback ladder after the plain growth (stage 0: 8 ranked seeds, ts = 1e-8): (method, coincidence tol ts, seeds, per-seed
# budget). Measured 2026-09-26 on the three cells stage 0 missed (all at r ~ 1/3, subsets of / near the E6 kissing
# configuration): N = 70 needs a seed outside the top 8 (the 9th, 10-45, grows all 70; the top ones close triangles
# STRUCTURALLY without being contacts), N = 68 needs ts = 1e-7, N = 67 (r = 1/3 + 1.4e-9, 40 rattlers, 53 contacts: its
# near-contacts have ~1e-9 gaps) needs ts = 1e-6.
LADDER = (('grow6bt', 1e-8, 300, 5.0), ('grow6bt', 1e-7, 300, 5.0), ('grow6bt', 1e-6, 1000, 5.0))


def recover6(A, r, rng, n_seeds=8, t_cap=900.0, ts=1e-8, tv=1e-8, ladder=LADDER):
    """Rigid growth from the ranked seeds (stage 0 = Grower, ts = tv = 1e-8, n_seeds); if no seed gives a completion with
    violation < 1e-12, the LADDER stages (GrowerBT with backtracking, looser coincidence tolerance, many more seeds) until
    one does or t_cap. Returns (C (n, 6), violation, info) of the best completion."""
    t0 = time.process_time(); P = Prob(A, r); best = None; tried = []
    stages = [('grow6', ts, n_seeds, t_cap)] + list(ladder)
    for method, ts_, ns, tseed in stages:
        for (a, b, cnt) in (seeds(P, top=ns, tol=max(1e-7, 10 * ts_)) or [(None, None, 0)]):
            G = Grower(P, ts=ts_, tv=ts_) if method == 'grow6' else GrowerBT(P, ts=ts_, tv=ts_)
            if a is not None:                                     # (no wall-wall pair constrains anything: all 'rest')
                G.place(a, np.array([P.m_ex[a], 0.0]))
                t = th_of(P, a, b); G.place(b, P.m_ex[b] * np.array([np.cos(t), np.sin(t)]))
                if method == 'grow6': G.grow(t_cap=t_cap, t0=t0)
                else: G.grow_bt(t_cap=tseed, node_cap=2000)
            grown = int(G.placed.sum())
            if method != 'grow6' and grown < 0.5 * P.n and len(tried) > 0: # a failed seed: skip the (slow) rest + finish
                tried.append(dict(method=method, ts=ts_, seed=(a, b), grown=grown)); continue
            place_rest(G, rng)
            Y, v = finish(P, G.Y.copy())
            tried.append(dict(method=method, ts=ts_, seed=(a, b), tri=cnt, grown=grown, viol=v))
            if best is None or v < best[1]: best = (Y, v, grown, (a, b), method, ts_)
            if v < 1e-12 or time.process_time() - t0 > t_cap: break
        if best[1] < 1e-12 or time.process_time() - t0 > t_cap: break
    Y, v, grown, sd, method, ts_ = best
    return np.hstack([A, Y]), v, dict(method=method if ts_ == 1e-8 else f'{method}@ts={ts_:.0e}', grown=grown, seed=sd,
                                      tried=tried, secs=round(time.process_time() - t0, 1))


def th_of(P, a, b):
    i, j = min(a, b), max(a, b)
    p = np.flatnonzero((P.I == i) & (P.J == j))[0]
    return P.th[p]


def run(args):
    try:
        return run1(args)
    except Exception as e:                                        # a failed cell is a row, not a dead pool
        import traceback
        return dict(N=args[0], ok=False, error=repr(e), tb=traceback.format_exc()[-800:])


def run1(args):
    n, polish_cap, seed = args
    prio.lower()
    rng = np.random.default_rng(seed); t0 = time.process_time()
    r_pr, A = trunc(6, n)
    C, viol, info = recover6(A, r_pr, rng)
    r_c = geomd.rmin(C); t1 = time.process_time()
    out = dict(N=n, printed=r_pr, viol=viol, r_completed=r_c, gap_completed=r_c - r_pr, method=info['method'],
               grown=info['grown'], seeds_tried=len(info['tried']), t_complete=round(t1 - t0, 1))
    if polish_cap > 0:
        # coarse polish (default trust region 0.1 r), then a FINE one (delta0 = 1e-6): a completion that is already exact
        # (printed - ~4e-12, the files' 12-decimal rounding) reads as 'stationary' to the coarse LP (its row perturbation
        # 1e-9 x delta is larger than the gap); the fine trust region resolves it (measured N = 20: -4.4e-12 -> -1.7e-13).
        s1, s2 = {}, {}
        q, r_p = slpd.polish(C, t_cap=polish_cap, seed=seed, stats=s1)
        q2, r_p2 = slpd.polish(q, t_cap=min(30.0, polish_cap), seed=seed + 1, delta0=1e-6, stats=s2)
        if r_p2 > r_p: q, r_p = q2, r_p2
        os.makedirs(OUT, exist_ok=True); np.save(os.path.join(OUT, f'hsp6_{n}.npy'), q)
        out.update(r_polished=r_p, gap_polished=r_p - r_pr, t_polish=round(time.process_time() - t1, 1),
                   polish_stop=[s1.get('stop'), s2.get('stop')])
        out['ok'] = bool(r_c >= r_pr - 1e-8 and r_p >= r_pr - 1e-12)
    out['seconds'] = round(time.process_time() - t0, 1)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--cells', default='30,66,100,160,200,250')
    ap.add_argument('--polish', type=float, default=60.0); ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--all', action='store_true'); ap.add_argument('--nmin', type=int, default=20)
    a = ap.parse_args(); prio.lower()
    cells = list(range(250, a.nmin - 1, -1)) if a.all else [int(x) for x in a.cells.split(',')]
    os.makedirs(OUT, exist_ok=True)
    S = json.load(open(SUMMARY)) if os.path.exists(SUMMARY) else {}
    from multiprocessing import Pool
    with Pool(a.workers) as pool:
        for x in pool.imap_unordered(run, [(n, a.polish, 20260926 + n) for n in cells]):
            print(json.dumps(x), flush=True)
            with open(LOG, 'a') as f: f.write(json.dumps(x) + '\n')
            S[str(x['N'])] = x
            S = {k: S[k] for k in sorted(S, key=int)}
            with open(SUMMARY + '.tmp', 'w') as f: json.dump(S, f, indent=1)
            os.replace(SUMMARY + '.tmp', SUMMARY)


if __name__ == '__main__':
    main()
