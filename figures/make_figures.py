"""Figures for the circle-packing repository, drawn only from files the repository ships.

usage:  python figures/make_figures.py [repo_root]      (repo_root defaults to the folder above this script)
writes: figures/packing_gallery.svg, figures/scoreboard.svg, figures/triangle16_two_optima.svg  (next to this script)

Sources: MANIFEST.csv, SUPERSEDED.csv, RESULTS_TABLE.md, certificates/<table>.zip, proofs/triangle16/out/local_cert16.json.
Every packing drawn is first re-checked in floating point (inside the container, no overlaps); the script stops if one fails.
"""
import csv, io, json, math, os, re, sys, zipfile
from fractions import Fraction

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Rectangle, Wedge, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(HERE)
OUT = HERE if os.path.basename(HERE) == "figures" else os.path.join(HERE, "figures")   # in the repo: figures/ itself
os.makedirs(OUT, exist_ok=True)

# ---- lab style (STYLE.md) -------------------------------------------------------------------------------------------
BG, BORDER, GRID, AX, TXT, ACC = "#1c2026", "#2f353d", "#22272e", "#9aa4ae", "#e8edf2", "#39ff14"
CAT = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"]
NEUTRAL = "#5b636c"
plt.rcParams.update({
    "svg.fonttype": "none", "svg.hashsalt": "moki-julio-circle-packing",
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 14, "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "axes.edgecolor": AX, "axes.labelcolor": TXT, "text.color": TXT, "xtick.color": AX, "ytick.color": AX,
    "axes.titlecolor": TXT, "grid.color": GRID, "grid.linewidth": 0.8, "lines.linewidth": 2,
})
W_IN = 12.5  # 900 px viewBox at 72 px per inch


def frame(fig):
    fig.patches.append(FancyBboxPatch((0.004, 0.004), 0.992, 0.992, boxstyle="round,pad=0,rounding_size=0.012",
                                      transform=fig.transFigure, fill=False, ec=BORDER, lw=1.5, zorder=-1))


def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, format="svg", metadata={"Date": None, "Creator": None})
    plt.close(fig)
    print("wrote figures/" + name)


def fmt_int(n):
    return f"{n:,}"


# ---- shipped data -----------------------------------------------------------------------------------------------------
released = list(csv.DictReader(open(os.path.join(ROOT, "MANIFEST.csv"), newline="")))
superseded = list(csv.DictReader(open(os.path.join(ROOT, "SUPERSEDED.csv"), newline="")))
gone = {(r["table"], r["N"]): r["status"] for r in superseded}

# container names: the summary table at the top of RESULTS_TABLE.md ("| shelf | container | records | ...")
NAMES = {}
for line in open(os.path.join(ROOT, "RESULTS_TABLE.md"), encoding="utf-8"):
    m = re.match(r"\|\s*([a-z0-9_]+)\s*\|\s*([^|]+?)\s*\|\s*(\d+)\s*\|", line)
    if m and m.group(1) != "shelf":
        NAMES.setdefault(m.group(1), m.group(2))


def load_cert(table, n):
    """certificate text from certificates/<table>.zip (or its split parts)."""
    d = os.path.join(ROOT, "certificates")
    zips = sorted(z for z in os.listdir(d) if z == table + ".zip" or z.startswith(table + "_part"))
    for z in zips:
        with zipfile.ZipFile(os.path.join(d, z)) as zf:
            nm = f"{table}_{n}.txt"
            if nm in zf.namelist():
                lines = [l.split() for l in io.TextIOWrapper(zf.open(nm), encoding="ascii") if l.strip()]
                assert lines[0][0] == "r"
                return float(lines[0][1]), [(float(a), float(b)) for a, b in lines[1:]]
    raise FileNotFoundError(f"{table}_{n}.txt")


# ---- containers in Packomania's frames (TECHNICAL.md, "same coordinate frame as Packomania's own files") -------------
def regular(k):
    """circumradius 1, centred, one vertex at (0, 1) (pentagon, 15-gon)."""
    return [(math.cos(math.radians(90 + 360 * j / k)), math.sin(math.radians(90 + 360 * j / k))) for j in range(k)]


def slack_polygon(P, x, y, r):
    s = math.inf
    for (x1, y1), (x2, y2) in zip(P, P[1:] + P[:1]):
        nx, ny = y2 - y1, -(x2 - x1)                     # outward normal for counter-clockwise vertices
        L = math.hypot(nx, ny)
        s = min(s, -((x - x1) * nx + (y - y1) * ny) / L - r)
    return s


def container(table):
    """returns (draw(ax), slack(x, y, r)) for one 2-D table."""
    if table.startswith("crc_"):
        h = int(table[4:]) / 1000
        return (lambda ax: ax.add_patch(Rectangle((-0.5, -h / 2), 1, h, fill=False, ec=AX, lw=1.5)),
                lambda x, y, r: min(0.5 - abs(x), h / 2 - abs(y)) - r)
    if table == "crt":
        P = [(0, 0), (1, 0), (0, 1)]
        return (lambda ax: ax.add_patch(Polygon(P, closed=True, fill=False, ec=AX, lw=1.5)),
                lambda x, y, r: slack_polygon(P, x, y, r))
    if table == "ccq":
        return (lambda ax: ax.add_patch(Wedge((0, 0), 1, 0, 90, fill=False, ec=AX, lw=1.5)),
                lambda x, y, r: min(x - r, y - r, 1 - math.hypot(x, y) - r))
    if table == "csc":
        return (lambda ax: ax.add_patch(Wedge((0, 0), 1, 0, 180, fill=False, ec=AX, lw=1.5)),
                lambda x, y, r: min(y - r, 1 - math.hypot(x, y) - r))
    if table in ("cpt", "cpd"):
        P = regular(5 if table == "cpt" else 15)
        return (lambda ax: ax.add_patch(Polygon(P, closed=True, fill=False, ec=AX, lw=1.5)),
                lambda x, y, r: slack_polygon(P, x, y, r))
    raise ValueError(table)


def check_packing(table, r, C):
    _, slack = container(table)
    worst_wall = min(slack(x, y, r) for x, y in C)
    worst_pair = min(math.hypot(a[0] - b[0], a[1] - b[1]) - 2 * r for i, a in enumerate(C) for b in C[i + 1:])
    if worst_wall < -1e-9 or worst_pair < -1e-9:
        raise SystemExit(f"{table} N={len(C)}: packing does not fit in the assumed frame (wall {worst_wall:.2e}, pair {worst_pair:.2e})")


# ---- figure 1: gallery of still-standing records ------------------------------------------------------------------------
GALLERY = ["crc_800", "crt", "ccq", "csc", "cpt", "cpd"]
GAIN_MIN = 1e-4          # selection rule: per container, the smallest N whose record still stands and gained >= 0.01 %


def pick(table):
    rows = [r for r in released if r["shelf"] == table and (table, r["N"]) not in gone and float(r["relative_gain"]) >= GAIN_MIN]
    return min(rows, key=lambda r: int(r["N"]))


def short_name(table):
    n = NAMES[table]
    return {"crc_800": "rectangle 1 x 0.8", "crt": "right triangle", "ccq": "quarter circle", "csc": "half circle",
            "cpt": "regular pentagon", "cpd": "regular 15-gon"}.get(table, n)


def fig_gallery():
    fig, axes = plt.subplots(2, 3, figsize=(W_IN, 9.2))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.80, bottom=0.03, wspace=0.06, hspace=0.30)
    frame(fig)
    fig.text(0.03, 0.945, "Six record packings that were still the best known on 2026-10-18", fontsize=21, weight="bold")
    fig.text(0.03, 0.895, "Equal circles, as large as possible, without overlapping each other or the edge "
                          "(drawn from the shipped certificates).", fontsize=14, color=AX)
    picks = []
    for ax, table in zip(axes.flat, GALLERY):
        row = pick(table)
        n = int(row["N"])
        r, C = load_cert(table, n)
        assert len(C) == n and abs(r - float(row["new_radius"])) <= 1e-12 * r
        check_packing(table, r, C)
        draw, _ = container(table)
        for x, y in C:
            ax.add_patch(Circle((x, y), r, fc=ACC, ec=BG, lw=0.6, alpha=0.9))
        draw(ax)
        ax.set_aspect("equal")
        ax.autoscale_view()
        xs = [x for x, _ in C]; ys = [y for _, y in C]
        if table.startswith("crc_"):
            ax.set_xlim(-0.53, 0.53); ax.set_ylim(-0.43, 0.43)
        elif table == "crt" or table == "ccq":
            ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
        elif table == "csc":
            ax.set_xlim(-1.03, 1.03); ax.set_ylim(-0.03, 1.03)
        else:
            ax.set_xlim(-1.03, 1.03); ax.set_ylim(-1.03, 1.03)
        ax.axis("off")
        gain = float(row["relative_gain"]) * 100
        ax.set_title(f"{short_name(table)}, {n} circles\nradius +{gain:.3g} % over the record it beat", fontsize=14.5,
                     color=TXT, pad=6)
        picks.append((table, n, row["new_radius"], row["published_radius"], row["relative_gain"], row["kind"]))
    save(fig, "packing_gallery.svg")
    return picks


# ---- figure 2: the scoreboard ---------------------------------------------------------------------------------------------
def fig_scoreboard():
    tables = sorted({r["shelf"] for r in released})
    cnt = {t: dict(rel=0, beaten=0, tied=0) for t in tables}
    for r in released:
        cnt[r["shelf"]]["rel"] += 1
    for s in superseded:
        cnt[s["table"]]["tied" if s["status"] == "TIED" else "beaten"] += 1
    order = sorted(tables, key=lambda t: (cnt[t]["rel"], t))
    fig, ax = plt.subplots(figsize=(W_IN, 9.6))
    fig.subplots_adjust(left=0.30, right=0.97, top=0.86, bottom=0.08)
    frame(fig)
    tot = {k: sum(c[k] for c in cnt.values()) for k in ("rel", "beaten", "tied")}
    still = tot["rel"] - tot["beaten"] - tot["tied"]
    fig.text(0.03, 0.945, f"{fmt_int(tot['rel'])} records released; {fmt_int(still)} still the best known on 2026-10-18",
             fontsize=21, weight="bold")
    fig.text(0.03, 0.905, "Per container. Records others have since beaten or tied are no longer claimed. "
                                                  "Number at each bar: still best known (of released).",
             fontsize=14, color=AX)
    y = range(len(order))
    s_ = [cnt[t]["rel"] - cnt[t]["beaten"] - cnt[t]["tied"] for t in order]
    b_ = [cnt[t]["beaten"] for t in order]
    t_ = [cnt[t]["tied"] for t in order]
    ax.barh(y, s_, height=0.72, color=CAT[0], label="still best known", ec=BG, lw=1)
    ax.barh(y, b_, left=s_, height=0.72, color=CAT[1], label="since beaten by others", ec=BG, lw=1)
    ax.barh(y, t_, left=[a + b for a, b in zip(s_, b_)], height=0.72, color=CAT[2], label="since tied", ec=BG, lw=1)
    ax.set_yticks(list(y))
    ax.set_yticklabels([NAMES[t].replace(" (equal 4-D balls)", "").replace(" (equal 5-D balls)", "")
                        .replace(" (equal 6-D balls)", "").replace(" (spheres)", ", spheres")
                        .replace("isosceles right triangle (legs 1)", "right triangle (legs 1)") for t in order], fontsize=12.5)
    for yi, t in zip(y, order):
        c = cnt[t]; st = c["rel"] - c["beaten"] - c["tied"]
        lab = f"{fmt_int(st)} of {fmt_int(c['rel'])}" if st != c["rel"] else fmt_int(c["rel"])
        ax.text(c["rel"] + 18, yi, lab, va="center", fontsize=12, color=TXT)
    ax.set_xlim(0, max(c["rel"] for c in cnt.values()) * 1.16)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("number of record packings released")
    ax.grid(axis="x"); ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower right", frameon=False, fontsize=13)
    save(fig, "scoreboard.svg")
    return cnt, tot


# ---- figure 3: the two optimal arrangements of 16 circles in a triangle ------------------------------------------------
SQ3 = math.sqrt(3)


def bary_xy(a, b):
    """barycentric (a, b, c = 1 - a - b) in the unit-side triangle with vertices V_a = (0, 0), V_b = (1, 0), V_c = (1/2, sqrt3/2)."""
    c = 1 - a - b
    return (b * 1 + c * 0.5, c * SQ3 / 2)


def fig_triangle16():
    cert = json.load(open(os.path.join(ROOT, "proofs", "triangle16", "out", "local_cert16.json")))
    m = math.sqrt(float(Fraction(cert["data"][0]["m_lo_sq"])))           # point form: min distance of the 16 points
    r16 = m / (2 * (1 + SQ3 * m))                                       # circle form (proofs/triangle16/README.md, Section 1)

    def frame_pts(q, swap=False):
        D = cert["data"][q]
        P = [(float(Fraction(u)), float(Fraction(v))) for u, v in D["frame_uv"]]
        R = [(float(Fraction(u)), float(Fraction(v))) for u, v in D["rattlers_uv"].values()]
        if swap:
            P = [(v, u) for u, v in P]; R = [(v, u) for u, v in R]
        return P, R, D["pairs"]

    PA, RA, pairsA = frame_pts(0, swap=True)     # frame A, reflected (a <-> b) so its shared part lines up with frame B
    PB, RB, pairsB = frame_pts(1)
    keyB = [min(range(15), key=lambda j: math.hypot(p[0] - PB[j][0], p[1] - PB[j][1])) for p in PA]
    sameA = [math.hypot(p[0] - PB[k][0], p[1] - PB[k][1]) < 1e-9 for p, k in zip(PA, keyB)]
    sameB = [any(math.hypot(p[0] - q[0], p[1] - q[1]) < 1e-9 for q in PA) for p in PB]
    n_same = sum(sameA)
    assert n_same == 11 and sum(sameB) == 11, n_same      # README: 11 of 15 frame points coincide, 4 do not

    def check(P, R, pairs):
        X = [bary_xy(*p) for p in P + R]
        dmin = min(math.dist(a, b) for i, a in enumerate(X) for b in X[i + 1:])
        assert dmin >= m * (1 - 1e-12), dmin
        for i, j in pairs:
            assert abs(math.dist(X[i], X[j]) - m) < 1e-12
        for a, b in P + R:
            assert min(a, b, 1 - a - b) >= -1e-15
        return X

    XA, XB = check(PA, RA, pairsA), check(PB, RB, pairsB)
    assert len(pairsA) == 20 and len(pairsB) == 20

    fig, axes = plt.subplots(1, 2, figsize=(W_IN, 7.6))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.76, bottom=0.10, wspace=0.04)
    frame(fig)
    fig.text(0.03, 0.935, "16 circles in a triangle: proved optimal, two arrangements", fontsize=21, weight="bold")
    fig.text(0.03, 0.885, f"Largest radius for a triangle of side 1: r = {r16:.12f}...  "
                          "Two arrangements are known;", fontsize=13.5, color=AX)
    fig.text(0.03, 0.848, "the proof does not say which of them is optimal, or whether both are.", fontsize=13.5, color=AX)
    g = 1 / (2 * SQ3)                                       # inradius of the unit triangle
    s = (g + m / 2) / g                                     # container = centre triangle pushed out by one circle radius
    cx, cy = 0.5, SQ3 / 6
    T = [(cx + s * (x - cx), cy + s * (y - cy)) for x, y in [(0, 0), (1, 0), (0.5, SQ3 / 2)]]
    for ax, X, same, pairs, title in ((axes[0], XA, sameA, pairsA, "arrangement A (shown mirrored)"),
                                      (axes[1], XB, sameB, pairsB, "arrangement B")):
        ax.add_patch(Polygon(T, closed=True, fill=False, ec=AX, lw=1.5))
        for i, j in pairs:
            ax.plot([X[i][0], X[j][0]], [X[i][1], X[j][1]], color=TXT, lw=1.0, alpha=0.85, zorder=3)
        for k, (x, y) in enumerate(X[:15]):
            ax.add_patch(Circle((x, y), m / 2, fc=CAT[0] if same[k] else CAT[1], ec=BG, lw=1.0, zorder=2))
        x, y = X[15]
        ax.add_patch(Circle((x, y), m / 2, fc="none", ec=TXT, lw=1.6, ls=(0, (4, 3)), zorder=2))
        ax.set_aspect("equal"); ax.set_xlim(T[0][0] - 0.03, T[1][0] + 0.03); ax.set_ylim(T[0][1] - 0.03, T[2][1] + 0.03)
        ax.axis("off")
        ax.set_title(title, fontsize=15, pad=4)
    hs = [Circle((0, 0), 1, fc=CAT[0]), Circle((0, 0), 1, fc=CAT[1]), Circle((0, 0), 1, fc="none", ec=TXT, ls="--"),
          matplotlib.lines.Line2D([], [], color=TXT, lw=1.2)]
    fig.legend(hs, ["same place in both (11)", "placed differently (4)", "loose circle (one valid spot shown)",
                    "contacts the proof uses (20)"],
               loc="lower center", ncol=4, frameon=False, fontsize=13, handlelength=1.2,
               handler_map={Circle: matplotlib.legend_handler.HandlerPatch()})
    save(fig, "triangle16_two_optima.svg")
    return r16, m, n_same


if __name__ == "__main__":
    import matplotlib.legend_handler  # noqa: F401
    picks = fig_gallery()
    for p in picks:
        print("  gallery:", p)
    cnt, tot = fig_scoreboard()
    print("  scoreboard totals:", tot)
    r16, m, n_same = fig_triangle16()
    print(f"  triangle16: m = {m:.15f}, r16 = {r16:.15f}, shared frame points {n_same}/15")
