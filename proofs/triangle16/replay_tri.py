#!/usr/bin/env python3
"""Independent replay of a triangle branch-and-bound tree (Moki&Julio).

Standard library only.  Does NOT import tri_engine.py and shares none of its geometry code: regions are kept as raw
integer bounds (never "tightened"), the maximum squared distance between two regions is computed by enumerating the
integer vertices of each region and taking the maximum over vertex pairs (a convex function over a product of convex
polygons attains its maximum at a pair of vertices), and every threshold is recomputed from the header's rational m_t.

Checks
  header   k^2 m_t^2 > 1 (tile diameter 1/k < m_t), G % k == 0, Tc == ceil(2 m_t^2 G^2), the tile list
  cover    pigeonhole (k^2 < n), or: every n-subset of the k^2 tiles has an image under the 6 permutations of the tile
           triples (the triangle's symmetry group) that is a root of this log and that root's tree is fully discarded
  tree     every node: each reduction strip lies within distance < m_t of every point of the other region
           (max vertex-pair distance^2 < Tc), each discard pair has max vertex-pair distance^2 < Tc (or a region is
           empty), each split's two children are the two closed halves of the parent, every child appears exactly once
  witness  FEASIBLE: all points in the triangle (A, B, C >= 0, A + B + C = G) and every pair with
           dA^2 + dB^2 + dC^2 >= Tc, i.e. distance >= m_t exactly

Usage: python replay_tri.py <tree.jsonl.gz> [more trees ...]      exit code 0 iff every file PASSES
Several files may form one run (parts); pass them together with --join.
"""
import gzip
import itertools
import json
import sys
from fractions import Fraction

PERM6 = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]


def vertices(b, G):
    """Integer vertices of {A in [b0,b1], B in [b2,b3], C in [b4,b5], A+B+C = G}; empty list if the set is empty."""
    a0, a1, b0, b1, c0, c1 = b
    out = set()
    for A in (a0, a1):
        for B in (b0, b1):
            C = G - A - B
            if c0 <= C <= c1:
                out.add((A, B, C))
        for C in (c0, c1):
            B = G - A - C
            if b0 <= B <= b1:
                out.add((A, B, C))
    for B in (b0, b1):
        for C in (c0, c1):
            A = G - B - C
            if a0 <= A <= a1:
                out.add((A, B, C))
    return sorted(out)


def max_d2(P, Q):
    best = -1
    for p in P:
        for q in Q:
            v = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2
            if v > best:
                best = v
    return best


class Fail(Exception):
    pass


class Stream:
    """One log file, read lazily (constant memory): R[0] is the header record, R.rest() streams the others."""

    def __init__(self, path):
        self.path = path
        with gzip.open(path, 'rt') as fh:
            self.head = json.loads(fh.readline())

    def __getitem__(self, i):
        if i != 0:
            raise IndexError(i)
        return self.head

    def rest(self):
        with gzip.open(self.path, 'rt') as fh:
            fh.readline()
            for line in fh:
                yield json.loads(line)


def check(files):
    recs = [Stream(f) for f in files]
    hdr0 = recs[0][0]['header']
    for R in recs:
        h = R[0]['header']
        for key in ('n', 'm_t', 'k', 'G', 'Tc'):
            if h[key] != hdr0[key]:
                raise Fail(f'parts disagree on {key}')
    n = hdr0['n']; mt = Fraction(hdr0['m_t']); k = hdr0['k']; G = hdr0['G']; Tc = hdr0['Tc']
    if not (k * k * mt * mt > 1):
        raise Fail('tile diameter 1/k is not < m_t')
    if G % k:
        raise Fail('G not divisible by k')
    t = 2 * mt * mt * G * G
    if Tc != -((-t.numerator) // t.denominator):
        raise Fail('Tc != ceil(2 m_t^2 G^2)')
    g = G // k
    tiles = sorted((i, j, l) for i in range(k) for j in range(k - i) for l in (k - 1 - i - j, k - 2 - i - j) if l >= 0)
    if len(tiles) != k * k:
        raise Fail('tile count')
    for R in recs:
        if [tuple(x) for x in R[0]['tiles']] != tiles:
            raise Fail('tile list differs')
    stats = dict(n=n, m_t=str(mt), k=k, nodes=0, ops=0, discards=0, splits=0, roots=0)
    if k * k < n:
        for R in recs:
            if not any('pigeonhole' in r for r in R.rest()):
                raise Fail('pigeonhole case without pigeonhole record')
        stats['result'] = 'PROVED (pigeonhole: %d tiles < %d points)' % (k * k, n)
        return stats
    tile_box = [(i * g, (i + 1) * g, j * g, (j + 1) * g, l * g, (l + 1) * g) for (i, j, l) in tiles]
    proved_roots = set()
    witness_ok = False
    owners = []
    clusters = {}
    for R in recs:
        for cl in R[0].get('accept', []):
            if len(cl['points']) != n or len(cl['frame']) != n or cl['rho'] < 0:
                raise Fail('bad cluster record')
            if cl['id'] in clusters and clusters[cl['id']] != cl:
                raise Fail('parts disagree on a cluster')
            clusters[cl['id']] = cl
    for R in recs:
        pending = {}
        root_of = {}
        for r in R.rest():
            if 'root' in r:
                combo = tuple(r['root'])
                if list(combo) != sorted(set(combo)) or len(combo) != n or not all(0 <= x < k * k for x in combo):
                    raise Fail('bad root combo')
                owner = {'combo': combo, 'ok': True, 'open': 1, 'accepted': 0}
                owners.append(owner)
                if r['id'] in pending:
                    raise Fail('duplicate root id')
                pending[r['id']] = [list(tile_box[x]) for x in combo]
                root_of[r['id']] = owner
                stats['roots'] += 1
                continue
            if 'fate' not in r:
                continue
            nid = r['id']
            if nid not in pending:
                raise Fail(f'node {nid} has no pending parent region')
            reg = pending.pop(nid)
            owner = root_of.pop(nid)
            owner['open'] -= 1
            stats['nodes'] += 1
            V = [vertices(b, G) for b in reg]
            empty = [i for i in range(n) if not V[i]]
            fate = r['fate']
            if empty and fate[0] not in ('E', 'D'):
                raise Fail(f'node {nid}: empty region but fate {fate[0]}')
            for op in r['ops']:
                if empty:
                    break
                i, side, tv, j = op
                if not (0 <= i < n and 0 <= j < n and i != j and side in (0, 1, 2, 3, 4, 5) and isinstance(tv, int)):
                    raise Fail(f'node {nid}: malformed op {op}')
                q = side // 2
                low = side % 2 == 0
                b = reg[i]
                if not (b[2 * q] <= tv <= b[2 * q + 1]):
                    raise Fail(f'node {nid}: op {op} outside the current bound range')
                strip = list(b)
                if low:
                    strip[2 * q + 1] = tv
                else:
                    strip[2 * q] = tv
                Vs = vertices(strip, G)
                if Vs and not (max_d2(Vs, V[j]) < Tc):
                    raise Fail(f'node {nid}: op {op}: strip not within m_t of every point of region {j}')
                b[side] = tv
                V[i] = vertices(b, G)
                stats['ops'] += 1
                if not V[i]:
                    empty = [i]
            if fate[0] == 'D':
                i, j = fate[1], fate[2]
                if not (0 <= i < n and 0 <= j < n and i != j):
                    raise Fail(f'node {nid}: malformed discard pair ({i},{j})')
                if not empty and not (max_d2(V[i], V[j]) < Tc):
                    raise Fail(f'node {nid}: discard pair ({i},{j}) not too close everywhere')
                stats['discards'] += 1
            elif fate[0] == 'E':
                if not empty:
                    raise Fail(f'node {nid}: E fate but no empty region')
                stats['discards'] += 1
            elif fate[0] == 'S':
                a, qq, c, lo_id, hi_id = fate[1:]
                if not (0 <= a < n and qq in (0, 1, 2) and isinstance(c, int)):
                    raise Fail(f'node {nid}: malformed split')
                b = reg[a]
                if not (b[2 * qq] <= c <= b[2 * qq + 1]):
                    raise Fail(f'node {nid}: split value outside range')
                lo = [list(x) for x in reg]; lo[a][2 * qq + 1] = c
                hi = [list(x) for x in reg]; hi[a][2 * qq] = c
                if lo_id in pending or hi_id in pending:
                    raise Fail('duplicate child id')
                pending[lo_id] = lo; pending[hi_id] = hi
                root_of[lo_id] = owner; root_of[hi_id] = owner
                owner['open'] += 2
                stats['splits'] += 1
            elif fate[0] == 'A':
                cl = clusters.get(fate[1])
                if cl is None or tuple(cl['combo']) != owner['combo']:
                    raise Fail(f'node {nid}: accept cites an unknown cluster or a cluster of another combination')
                if empty:
                    raise Fail(f'node {nid}: accept on an empty region')
                rho = cl['rho']
                for i in range(n):
                    if not cl['frame'][i]:
                        continue
                    c = cl['points'][i]
                    for v in V[i]:
                        if any(abs(v[u] - c[u]) > rho for u in range(3)):
                            raise Fail(f'node {nid}: accept, region {i} leaves the cluster box')
                owner['accepted'] += 1
                stats['accepted'] = stats.get('accepted', 0) + 1
            elif fate[0] == 'F':
                pts = [tuple(p) for p in fate[1]]
                if len(pts) != n:
                    raise Fail('witness size')
                for p in pts:
                    if min(p) < 0 or sum(p) != G:
                        raise Fail('witness point outside the triangle')
                for p, s2 in itertools.combinations(pts, 2):
                    if (p[0] - s2[0]) ** 2 + (p[1] - s2[1]) ** 2 + (p[2] - s2[2]) ** 2 < Tc:
                        raise Fail('witness pair closer than m_t')
                witness_ok = True
                owner['ok'] = False
            else:  # 'U' or unknown: this root is not proved
                owner['ok'] = False
    for o in owners:
        if o['ok'] and o['open'] == 0:
            proved_roots.add(o['combo'])
    stats['witness'] = witness_ok
    if witness_ok:
        stats['result'] = 'FEASIBLE (exact witness at m_t)'
        return stats
    # cover: every n-subset maps to a proved root
    T = k * k
    tindex = {tt: q for q, tt in enumerate(tiles)}
    tperm = [[tindex[(tt[p[0]], tt[p[1]], tt[p[2]])] for tt in tiles] for p in PERM6]
    missing = 0
    total = 0
    for c in itertools.combinations(range(T), n):
        total += 1
        if not any(tuple(sorted(p[x] for x in c)) in proved_roots for p in tperm):
            missing += 1
    stats['subsets'] = total
    stats['uncovered_subsets'] = missing
    acc = sum(o['accepted'] for o in owners)
    if missing == 0 and acc:
        stats['result'] = ('PROVED MODULO CLUSTERS: every configuration with all distances >= m_t has a D3 image whose '
                           'frame points lie within rho (barycentric sup-norm, grid units) of a listed cluster '
                           f'configuration ({acc} accepted leaves in {len(clusters)} clusters)')
    else:
        stats['result'] = 'PROVED' if missing == 0 else f'NOT PROVED ({missing} of {total} tile subsets uncovered)'
    return stats


def main():
    files = [a for a in sys.argv[1:] if not a.startswith('--')]
    join = '--join' in sys.argv
    groups = [files] if join else [[f] for f in files]
    ok = True
    for grp in groups:
        try:
            st = check(grp)
            good = st['result'].startswith('PROVED') or st['result'].startswith('FEASIBLE')
            print(('PASS ' if good else 'FAIL ') + ' + '.join(grp) + ' ' + json.dumps(st))
            ok &= good
        except Fail as e:
            print('FAIL ' + ' + '.join(grp) + ': ' + str(e))
            ok = False
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
