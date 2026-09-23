"""LOCAL-OPTIMALITY certificate for an equal-circle packing (circle form; containers tri / rect:h / quad as in certify_circ.py).

Given a certificate z* = (c_1..c_n, r*) that is exactly feasible, this proves, rigorously:

  THEOREM. Let B be the "backbone" circles (those carrying contact stress; the others are rattlers). Every feasible packing whose
  backbone centres AND radius are within rho (max-norm) of z* -- rattlers anywhere -- has radius <= r* + Delta; and any feasible
  packing in that box with radius >= r* is within t0 of z*. Hence the box contains a genuine LOCAL MAXIMUM of the packing problem,
  within t0 of z*, whose radius lies in [r*, r* + Delta].   (Typically Delta ~ 1e-44 and t0 ~ 1e-42: z* IS the local optimum.)

PROOF (all constants computed below; exact rationals except K, which is a rigorous floating-point upper bound).
Write every kept constraint as g_a(z) >= 0 with g_a a polynomial of degree <= 2 (pairs |ci-cj|^2 - 4r^2, straight walls linear,
quadrant arc (1-r)^2 - |ci|^2). Let S = kept constraints, eps_a = g_a(z*) (>= 0), G_a = grad g_a(z*), q_a(d) = g_a(z*+d) -
eps_a - G_a.d (exact, quadratic), |q_a(d)| <= kappa_a t^2 with t = |d|_inf over the backbone coordinates and r
(kappa: pair 8, arc 2, straight wall 0).  Multipliers lambda_a > 0 (from an LP, then taken as EXACT rationals) define the exact
residual  res = sum_a lambda_a G_a + e_r.  Put E = sum lambda eps, Q = sum lambda kappa, lmin = min lambda.
 (1) For feasible z = z* + d:  0 <= sum lambda_a g_a(z) = E - dr + res.d + sum lambda_a q_a  =>  dr <= E + |res|_1 t + Q t^2.
 (2) If dr >= 0: each u_a = g_a(z) >= 0 has lambda_a u_a <= sum lambda u <= E + |res|_1 t + Q t^2, and G_a.d = u_a - eps_a - q_a,
     so |G d|_inf <= A0 + A1 t + A2 t^2 with A0 = E/lmin + max eps, A1 = |res|_1/lmin, A2 = Q/lmin + max kappa.
 (3) Rigidity: for ANY matrix Y, d = Y G d + (I - Y G) d, so t <= K |G d|_inf with K = |Y|_inf / (1 - eta), eta >= |I - Y G|_inf.
     K is bounded rigorously in IEEE double: Higham's bound |fl(A B) - A B| <= gamma_k |A||B| (any summation order, FMA or not),
     the entrywise rounding of G, and 1e-9 relative slack on the final float sums; requires eta < 1/2.
 (4) So t <= K A0 + K A1 t + K A2 t^2. With rho = (1/2 - K A1) / (K A2): t <= rho  =>  t <= t0 = 2 K A0, and then by (1)
     dr <= Delta = E + |res|_1 t0 + Q t0^2.  The box {t <= rho} meets the (compact) feasible set, which contains z*, so r attains a
     maximum there, >= r*; by the above it sits at t <= t0 < rho, inside the box: a local maximum of the full problem.  QED
S and lambda are chosen by a max-support LP (every stressable near-active constraint), then lambda maximises min lambda.  The theorem
holds for ANY choice; a poor choice only weakens or fails the certificate, never makes it wrong.  The only irrational number,
sqrt2 (triangle hypotenuse), is enclosed in [S2_LO, S2_HI] (60 digits, checked exactly).
Verdicts: CERTIFIED | NOT_RIGID (backbone has a first-order flex: slides or buckling, see 'flex') | NOT_STATIONARY (no KKT
multipliers: an improving direction exists) | INFEASIBLE (certificate not feasible) | K_FAIL (eta >= 1/2 or K A1 >= 1/2).
usage: py -3.11 lopt.py <container> <certificate.txt>
"""
import sys, os, math, json, numpy as np
from fractions import Fraction as F
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, eye as speye, hstack, vstack, csr_matrix

U = 2.0 ** -53                                   # IEEE double unit roundoff (round to nearest)
S2_LO = F(math.isqrt(2 * 10 ** 120), 10 ** 60); S2_HI = S2_LO + F(1, 10 ** 60)
assert S2_LO ** 2 < 2 < S2_HI ** 2
NEAR = F(1, 10 ** 30)                            # kept if the exact normalised slack is below this
PRE = 1e-9                                       # float prefilter for candidate contacts (selection only, never a proof step)

# ---- numbers a + b*sqrt2 (b != 0 only on the triangle hypotenuse) --------------------------------------------------------------
def hi(v): a, b = v; return a + b * (S2_HI if b > 0 else S2_LO)
def lo(v): a, b = v; return a + b * (S2_LO if b > 0 else S2_HI)
def absup(v): return max(abs(hi(v)), abs(lo(v)))
def tofloat(v):
    a, b = v
    if b == 0: return float(a)                   # Fraction -> float is correctly rounded
    assert a == 0 and b == -1; return -math.sqrt(2.0)   # IEEE sqrt is correctly rounded

def load(path):
    L = [l.split() for l in open(path) if l.strip()]
    assert L[0][0] == "r"
    return F(L[0][1]), [(F(a), F(b)) for a, b in L[1:]]

def feasible(cont, C, r):
    """Exact, same rules as certify_circ.py; pairs in integers at a common power-of-ten scale (all inputs are finite decimals)."""
    import certify_circ
    if not (r > 0 and all(certify_circ.fits(cont, x, y, r) for x, y in C)): return False
    den = 1
    for v in [r] + [w for xy in C for w in xy]:
        den = math.lcm(den, v.denominator)          # 09-23: no 10^80 cap (a centre on an axis converges to ~1e-41 -> 85+ decimals)
    X = [(int(x * den), int(y * den)) for x, y in C]; R4 = 4 * int(r * den) ** 2
    assert all(F(a, den) == x and F(b, den) == y for (a, b), (x, y) in zip(X, C)) and F(int(r * den), den) == r
    for i, (a0, b0) in enumerate(X):
        for a1, b1 in X[i + 1:]:
            if (a0 - a1) ** 2 + (b0 - b1) ** 2 < R4: return False
    return True

def constraints(cont, C, r):
    """Candidate constraints near contact. Each: (name, slack (a,b), normalised slack (Fraction upper bound),
    gradient {col: (a,b)}, kappa).  Column 2i, 2i+1 = centre i; column R = 2n = radius."""
    n = len(C); R = 2 * n; Z = F(0); out = []
    cf = np.array([[float(x), float(y)] for x, y in C]); rf = float(r)
    D = np.sqrt(((cf[:, None] - cf[None]) ** 2).sum(-1)); iu = np.triu_indices(n, 1)
    for i, j in zip(*iu):
        if D[i, j] <= 2 * rf * (1 + PRE):
            dx, dy = C[i][0] - C[j][0], C[i][1] - C[j][1]; g = dx * dx + dy * dy - 4 * r * r
            out.append((("p", int(i), int(j)), (g, Z), g / (4 * r * r),
                        {2 * i: (2 * dx, Z), 2 * i + 1: (2 * dy, Z), 2 * j: (-2 * dx, Z), 2 * j + 1: (-2 * dy, Z), R: (-8 * r, Z)}, 8))
    for i, (x, y) in enumerate(C):
        xf, yf = cf[i]
        if cont == "tri":
            W = [(xf - rf, (x - r, Z), {2 * i: (F(1), Z), R: (F(-1), Z)}),
                 (yf - rf, (y - r, Z), {2 * i + 1: (F(1), Z), R: (F(-1), Z)}),
                 ((1 - xf - yf) / math.sqrt(2) - rf, (1 - x - y, -r), {2 * i: (F(-1), Z), 2 * i + 1: (F(-1), Z), R: (Z, F(-1))})]
            kap = [0, 0, 0]
        elif cont.startswith("rect:"):
            h = F(cont.split(":")[1]); hf = float(h)
            W = [(xf + 0.5 - rf, (x + F(1, 2) - r, Z), {2 * i: (F(1), Z), R: (F(-1), Z)}),
                 (0.5 - xf - rf, (F(1, 2) - x - r, Z), {2 * i: (F(-1), Z), R: (F(-1), Z)}),
                 (yf + hf / 2 - rf, (y + h / 2 - r, Z), {2 * i + 1: (F(1), Z), R: (F(-1), Z)}),
                 (hf / 2 - yf - rf, (h / 2 - y - r, Z), {2 * i + 1: (F(-1), Z), R: (F(-1), Z)})]
            kap = [0, 0, 0, 0]
        elif cont == "quad":
            W = [(xf - rf, (x - r, Z), {2 * i: (F(1), Z), R: (F(-1), Z)}),
                 (yf - rf, (y - r, Z), {2 * i + 1: (F(1), Z), R: (F(-1), Z)}),
                 (1 - math.hypot(xf, yf) - rf, ((1 - r) ** 2 - x * x - y * y, Z),
                  {2 * i: (-2 * x, Z), 2 * i + 1: (-2 * y, Z), R: (-2 * (1 - r), Z)})]
            kap = [0, 0, 2]
        else: raise ValueError(cont)
        for k, ((sf, g, grad), kp) in enumerate(zip(W, kap)):
            if sf <= PRE * rf:
                norm = hi(g) / (r if not (cont == "quad" and k == 2) else 2 * r)   # arc slack ~ 2 r * distance slack
                out.append((("w", i, k), g, norm, grad, kp))
    return out

def lp(c, A_eq, b_eq, A_ub, b_ub, bounds):
    """HiGHS; on an unrecognised status (09-23: 'model_status Unknown' on large ccq) retry with dual simplex, then IPM."""
    for method in ("highs", "highs-ds", "highs-ipm"):
        res = linprog(c, A_eq=A_eq, b_eq=b_eq, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method=method,
                      options={"primal_feasibility_tolerance": 1e-10, "dual_feasibility_tolerance": 1e-10, "time_limit": 120.0})
        if res.status in (0, 2): return res
    return res

def certify(cont, path, want_flex=True):
    """Try every plausible contact cut (widest multiplicative gaps in the sorted normalised slacks, <= 1e-12, gap >= 1e3),
    tightest first; return the first CERTIFIED result, else the result at the widest gap.  Any cut is valid for the theorem
    (it holds for any kept set S); a looser cut only means a larger Delta, which is reported exactly."""
    r, C = load(path); n = len(C)
    base = {"container": cont, "file": os.path.basename(path), "n": n, "r": str(r)}
    if not feasible(cont, C, r):
        base["verdict"] = "INFEASIBLE"; return base
    cand = constraints(cont, C, r)
    sl = sorted(max(float(c[2]), 1e-46) for c in cand) + [PRE]          # exact zeros join the 45-digit cluster
    gaps = [(sl[k + 1] / sl[k], F(math.sqrt(sl[k] * sl[k + 1]))) for k in range(len(sl) - 1) if sl[k] <= 1e-12]
    widest = max(gaps, key=lambda g: g[0]) if gaps else (0.0, NEAR)
    cuts = sorted({g[1] for g in gaps if g[0] >= 1e3} | {widest[1]})
    first = None
    for near in cuts:
        out = dict(base); out["contact_cut"] = float(near)
        out["cut_gap_ratio"] = next(g for g in gaps if g[1] == near)[0] if gaps else 0.0
        out["cuts_tried"] = cuts.index(near) + 1
        A = [c for c in cand if c[2] < near]
        rest = [float(c[2]) for c in cand if c[2] >= near]
        out["near_active"] = len(A); out["gap_min_excluded_slack"] = min(rest) if rest else PRE
        out["max_kept_slack"] = float(max(c[2] for c in A)) if A else None
        if not A: out["verdict"] = "NOT_STATIONARY"
        else: out = core(cont, C, r, A, out, want_flex)
        if out["verdict"] == "CERTIFIED": return out
        if near == widest[1]: first = out
    return first

def core(cont, C, r, A, out, want_flex):
    n = len(C); R = 2 * n
    m = len(A)
    # float gradient matrix over ALL columns
    rows, cols, vals = [], [], []
    for a, c in enumerate(A):
        for k, v in c[3].items(): rows.append(a); cols.append(k); vals.append(tofloat(v))
    Gall = coo_matrix((vals, (rows, cols)), shape=(m, 2 * n + 1)).tocsr()
    def circles(a): nm = A[a][0]; return nm[1:3] if nm[0] == "p" else nm[1:2]

    def support_l1(allowed):
        """Positive multipliers with the smallest force-balance residual (the theorem needs lam > 0 and a SMALL residual,
        not an exact zero: it enters only through K |res|_1 / lmin < 1/2 and Delta).  ell_1 residual e = G^T lam + e_R.
          phase 1: min |e|_1 over lam >= 0                                  -> res1
          phase 2: max sum y, y <= lam, y <= cap, |e|_1 <= max(res1 (1+1e-3), 1e-9)   (max support)
          phase 3: on the support, max tau, lam >= tau, same residual bound; drop numerically-zero stresses and repeat.
        Returns (S, lam) or (None, msg)."""
        Ga = Gall[allowed]; cu = np.unique(Ga.nonzero()[1].tolist() + [R]); pu = len(cu); ma = len(allowed)
        iR = int(np.nonzero(cu == R)[0][0]); eR = np.zeros(pu); eR[iR] = 1.0; Gt = Ga[:, cu].T.tocsr()
        I_p = speye(pu, format="csr")
        # vars [lam (ma), ep (pu), em (pu)]:  Gt lam - ep + em = -eR
        Aeq = hstack([Gt, -I_p, I_p]).tocsr()
        r1 = lp(np.r_[np.zeros(ma), np.ones(2 * pu)], Aeq, -eR, None, None, [(0, None)] * (ma + 2 * pu))
        if r1.status != 0: return None, "LP_FAIL: phase1 " + r1.message
        res1 = float(r1.fun); lam1 = r1.x[:ma]; out.setdefault("kkt_residual_l1", res1)
        if res1 > 1e-6: return None, f"NOT_STATIONARY (min force-balance residual {res1:.2e})"
        bound = max(res1 * (1 + 1e-3), 1e-9); cap = 1e-4 * max(lam1.max(), 1e-300)   # 1e-9 ~ LP resolution; exact res later
        # vars [lam (ma), ep, em, y (ma)]
        Aeq2 = hstack([Gt, -I_p, I_p, csr_matrix((pu, ma))]).tocsr()
        Aub2 = vstack([hstack([-speye(ma), csr_matrix((ma, 2 * pu)), speye(ma)]),
                       hstack([csr_matrix((1, ma)), csr_matrix(np.ones((1, 2 * pu))), csr_matrix((1, ma))])]).tocsr()
        r2 = lp(np.r_[np.zeros(ma + 2 * pu), -np.ones(ma) / cap], Aeq2, -eR, Aub2, np.r_[np.zeros(ma), bound],
                [(0, None)] * (ma + 2 * pu) + [(0, cap)] * ma)
        if r2.status != 0: return None, "LP_FAIL: phase2 " + r2.message
        lam2 = r2.x[:ma]; y = r2.x[ma + 2 * pu:]
        # support = clearly stressable (y at its cap) + small stresses only as far as balance needs them: escalate the
        # threshold 1e-7 -> 1e-10 -> any, first that balances wins (09-23: 1e-12 alone let LP noise in; y-only cut out
        # required small stresses)
        for thr in (1e-7, 1e-10, 0.0):
            S = [allowed[q] for q in range(ma) if y[q] > 0.5 * cap or lam2[q] > thr * lam2.max()]
            for _ in range(20):
                cs_ = np.unique(Gall[S].nonzero()[1].tolist() + [R]); pS = len(cs_); mS = len(S)
                e = np.zeros(pS); e[int(np.nonzero(cs_ == R)[0][0])] = 1.0; GtS = Gall[S][:, cs_].T.tocsr(); I_s = speye(pS, format="csr")
                # vars [lam (mS), ep, em, tau]: GtS lam - ep + em = -e ; -lam + tau <= 0 ; sum(ep+em) <= bound ; max tau
                Aeq3 = hstack([GtS, -I_s, I_s, csr_matrix((pS, 1))]).tocsr()
                Aub3 = vstack([hstack([-speye(mS), csr_matrix((mS, 2 * pS)), csr_matrix(np.ones((mS, 1)))]),
                               hstack([csr_matrix((1, mS)), csr_matrix(np.ones((1, 2 * pS))), csr_matrix((1, 1))])]).tocsr()
                r3 = lp(np.r_[np.zeros(mS + 2 * pS), -1.0], Aeq3, -e, Aub3, np.r_[np.zeros(mS), bound],
                        [(0, None)] * (mS + 2 * pS) + [(None, None)])
                if r3.status != 0: break
                lam3 = r3.x[:mS]; tau = r3.x[-1]; big = lam3.max()
                if tau > 1e-9 * big: return S, lam3
                S = [a for a, l in zip(S, lam3) if l > 1e-9 * big]
        return None, "NOT_STATIONARY (no balanced positive support)"

    def support_cone(allowed):
        """Original method (certified the first 740): max-support KKT cone LP, vars [lam, s, y, ys], G^T lam + s e_R = 0,
        y <= lam, ys <= s, y, ys in [0,1]; then normalised max-min LP on the support, dropping stresses <= 1e-9 of the max."""
        Ga = Gall[allowed]; cu = np.unique(Ga.nonzero()[1].tolist() + [R]); Gu = Ga[:, cu]; pu = len(cu); ma = len(allowed)
        eR = np.zeros(pu); eR[int(np.nonzero(cu == R)[0][0])] = 1.0
        Aeq = hstack([Gu.T, csr_matrix(eR[:, None]), csr_matrix((pu, ma + 1))]).tocsr()
        Aub = vstack([hstack([-speye(ma), csr_matrix((ma, 1)), speye(ma), csr_matrix((ma, 1))]),
                      hstack([csr_matrix((1, ma)), csr_matrix([[-1.0]]), csr_matrix((1, ma)), csr_matrix([[1.0]])])]).tocsr()
        res = lp(np.r_[np.zeros(ma + 1), -np.ones(ma + 1)], Aeq, np.zeros(pu), Aub, np.zeros(ma + 1),
                 [(0, None)] * (ma + 1) + [(0, 1)] * (ma + 1))
        if res.status != 0: return None, "LP_FAIL: cone " + res.message
        y = res.x[ma + 1:]
        if y[-1] < 0.5: return None, "NOT_STATIONARY (cone)"
        S = [allowed[q] for q in range(ma) if y[q] > 0.5]
        for _ in range(20):
            lam, tau = exact_maxmin(S)
            if lam is None: return None, "LP_FAIL: cone maxmin"
            big = lam.max()
            if tau > 1e-9 * big: return S, lam
            S = [a for a, l in zip(S, lam) if l > 1e-9 * big]
        return None, "LP_FAIL: cone support did not settle"

    def exact_maxmin(S):
        """max tau s.t. G_S^T lam = -e_R (equality), lam >= tau.  Returns (lam, tau) or (None, None)."""
        cs_ = np.unique(Gall[S].nonzero()[1].tolist() + [R]); pS = len(cs_); mS = len(S)
        e = np.zeros(pS); e[int(np.nonzero(cs_ == R)[0][0])] = 1.0; GS = Gall[S][:, cs_]
        r2 = lp(np.r_[np.zeros(mS), -1.0], hstack([GS.T, csr_matrix((pS, 1))]).tocsr(), -e,
                hstack([-speye(mS), csr_matrix(np.ones((mS, 1)))]).tocsr(), np.zeros(mS), [(0, None)] * mS + [(None, None)])
        if r2.status != 0: return None, None
        return r2.x[:mS], r2.x[-1]

    def support(allowed):
        S, lam = support_cone(allowed)
        if S is not None: return S, lam
        msg1 = lam; S, lam = support_l1(allowed)
        if S is None: return None, lam if lam.startswith("NOT_STATIONARY") and msg1.startswith("NOT_STATIONARY") else (lam + " | " + msg1)
        lam_e, tau = exact_maxmin(S)                                            # re-balance exactly on the l1 support
        if lam_e is not None and tau > 1e-9 * lam_e.max(): lam = lam_e
        return S, lam

    # Choose the stressed set S.  With redundant contacts the multipliers are not unique; a circle that can move at first order
    # (a straight chain, a slide) may carry OPTIONAL stress.  Loop: find first-order flex circles of the current backbone, drop
    # their constraints, and ask the LP again.  If multipliers still exist, those circles are rattlers for this certificate;
    # if not, the stress through them is REQUIRED and the packing has a genuine first-order flex -> NOT_RIGID (buckling).
    allowed = list(range(m)); dropped = set(); flex = None
    for loop in range(50):
        S, lam0 = support(allowed)
        if S is None:
            if not dropped:
                out["verdict"] = "NOT_STATIONARY" if lam0.startswith("NOT_STATIONARY") else "LP_FAIL"; out["lp_msg"] = lam0; return out
            out["verdict"] = "NOT_RIGID"; out["flex"] = flex; out["stress_required_at"] = sorted(flex_c)
            if want_flex == "second": out["_so"] = (S_prev, lam_prev, bcols_prev, G_prev, ker_prev)
            return out
        # backbone = every circle in a stressed constraint; its columns = BOTH coordinates of each such circle + r, whether or
        # not a gradient entry is zero (a zero first-order entry still moves the constraint at second order: 09-23 control)
        backbone = sorted(set(i for a in S for i in circles(a)))
        bcols = np.array(sorted([2 * i for i in backbone] + [2 * i + 1 for i in backbone] + [R])); p = len(bcols)
        jR = int(np.nonzero(bcols == R)[0][0])
        G = Gall[S][:, bcols].toarray()                                       # |S| x p, float (entrywise correctly rounded)
        M = G.T @ G
        # rank test on G itself (SVD error ~ u*sigma_max); eigenvalues of G^T G square the error (09-23: a zero column read as
        # 1.8e-16).  This test only SELECTS; a wrong "rigid" here cannot certify anything: eta/K below are rigorous.
        _, sv, Vt = np.linalg.svd(G, full_matrices=True); sv = np.r_[sv, np.zeros(p - len(sv))]
        ev = sv[::-1] ** 2
        if len(S) >= p and sv[-1] > 1e-10 * sv[0]: break
        ker = Vt[sv <= 1e-10 * sv[0]].T; mv = {}
        S_prev, lam_prev, bcols_prev, G_prev, ker_prev = S, lam0, bcols, G, ker
        for q, k in enumerate(bcols):
            if k != R: mv[int(k) // 2] = max(mv.get(int(k) // 2, 0.0), float(np.abs(ker[q]).max()))
        flex_c = {i for i, w in mv.items() if w > 1e-8}
        flex = {"dr": float(np.abs(ker[jR]).max()), "moving": [[i, round(w, 4)] for i, w in sorted(mv.items(), key=lambda kv: -kv[1])[:6]]}
        k0 = ker[:, 0] / np.abs(ker[:, 0]).max()
        flex["vec"] = [[int(bcols[q]) // 2, int(bcols[q]) % 2, float(k0[q])] for q in range(p) if bcols[q] != R and abs(k0[q]) > 1e-9]
        if flex["dr"] > 1e-8 or not flex_c:
            out["verdict"] = "NOT_RIGID"; out["flex"] = flex; return out
        dropped |= flex_c; allowed = [a for a in allowed if not (set(circles(a)) & flex_c)]
    else:
        out["verdict"] = "NOT_RIGID"; out["flex"] = flex; return out
    out.update({"stressed": len(S), "p": p, "redundancy": len(S) - p, "backbone": len(backbone), "rattlers": n - len(backbone),
                "optional_stress_dropped": sorted(dropped)})
    out["sigma_min"] = float(math.sqrt(max(ev[0], 0.0))); out["sigma_max"] = float(math.sqrt(ev[-1]))
    ms = len(S); eRp = np.zeros(p); eRp[jR] = 1.0
    lam_lp = lam0                                                             # from support(): G_S^T lam = -e_R, min lam > 0
    try: lam = lam_lp - G @ np.linalg.solve(M, G.T @ lam_lp + eRp)          # least-squares polish of G^T lam = -e_R
    except np.linalg.LinAlgError: out["verdict"] = "K_FAIL"; out["why"] = "singular G^T G"; return out
    out["lambda_polished"] = bool(np.all(lam > 0))
    if not np.all(lam > 0): lam = lam_lp                                      # polish flipped a tiny stress: use the LP's
    if not np.all(lam > 0): out["verdict"] = "LAMBDA_FAIL"; return out        # (the residual is computed EXACTLY either way)
    # ---- rigorous K: Y = (G^T G)^-1 G^T; eta >= |I - Y G|_inf; K = |Y|_inf / (1 - eta) ----------------------------------
    Y = np.linalg.solve(M, G.T); P = Y @ G; Wm = np.abs(Y) @ np.abs(G)
    if not (np.all(np.isfinite(Y)) and np.all(np.isfinite(P))): out["verdict"] = "K_FAIL"; out["why"] = "non-finite Y"; return out
    if not np.all((np.diag(P) >= 0.5) & (np.diag(P) <= 2.0)): out["verdict"] = "K_FAIL"; out["why"] = "diag(YG) off"; return out
    k_in, g_k = ms, ms * U / (1 - ms * U); g_p = p * U / (1 - p * U)
    dg = np.diag(P); assert np.all((dg >= 0.5) & (dg <= 2.0)), "diag(YG) outside [1/2,2]: Sterbenz step not exact"
    IP = np.abs(P); np.fill_diagonal(IP, 1.0 - dg); IP = np.abs(IP)
    c1 = (g_k + 2 * U) / (1 - g_k)
    bracket = IP.sum(1) + c1 * Wm.sum(1)
    eta = F(float(bracket.max())) * (1 + F(1, 10 ** 9)) / (1 - F(g_p))
    Yinf = F(float(np.abs(Y).sum(1).max())) * (1 + F(1, 10 ** 9)) / (1 - F(g_k))
    if eta >= F(1, 2): out["verdict"] = "K_FAIL"; out["eta"] = float(eta); return out
    K = Yinf / (1 - eta)
    # ---- exact constants ------------------------------------------------------------------------------------------------------
    lamF = [F(float(v)) for v in lam]
    E = F(0); emax = F(0); Q = F(0); kmax = 0; resid = {}
    for lf, a in zip(lamF, S):
        name, g, _, grad, kp = A[a]
        E += lf * hi(g); emax = max(emax, absup(g)); Q += lf * kp; kmax = max(kmax, kp)
        for k, (va, vb) in grad.items():
            ra, rb = resid.get(k, (F(0), F(0))); resid[k] = (ra + lf * va, rb + lf * vb)
    ra, rb = resid[R]; resid[R] = (ra + 1, rb)
    res1 = sum(absup(v) for v in resid.values())
    lmin = min(lamF)
    A0 = max(E, F(0)) / lmin + emax; A1 = res1 / lmin; A2 = Q / lmin + kmax
    if K * A1 >= F(1, 2): out["verdict"] = "K_FAIL"; out["KA1"] = float(K * A1); return out
    rho = (F(1, 2) - K * A1) / (K * A2) if A2 > 0 else F(1)
    t0 = 2 * K * A0
    Delta = max(E, F(0)) + res1 * t0 + Q * t0 * t0
    if not t0 < rho:                                    # 09-23: was an assert (crashed crc_600 211/245); it is simply no certificate
        out["verdict"] = "K_FAIL"; out["why"] = "t0 >= rho"; out["t0"] = float(t0); out["rho"] = float(rho); return out
    out.update({"verdict": "CERTIFIED", "rho": float(rho), "rho_rel": float(rho / r), "t0": float(t0), "Delta": float(Delta),
                "Delta_rel": float(Delta / r), "K": float(K), "eta": float(eta), "lambda_min": float(lmin),
                "lambda_ratio": float(sum(lamF) / lmin), "residual_l1": float(res1),
                "S": [list(A[a][0]) for a in S], "lambda": [v.hex() for v in lam]})
    return out

if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    o = certify(sys.argv[1], sys.argv[2])
    print(json.dumps({k: v for k, v in o.items() if k not in ("S", "lambda")}, indent=1))
