"""Exhaustive verification of the closed forms in
"Reachability, exact diameters, and label-erasure heuristics in two-row
sliding-domino puzzles".

The program enumerates every legal configuration for small boards, builds the
unit, straight and finger move graphs, and compares breadth-first-search
results with the formulas of the paper.  It uses no randomness.

Usage:  python verify_two_row.py [--nmax-q1 N] [--nmax-gen N] [--kmax-lab K]
        python verify_two_row.py --heuristic-only --out heuristic_log.txt
"""
import argparse
import itertools
import math
import sys
import time
from collections import deque

# ---------------------------------------------------------------------------
# Generic engine.  A configuration is a sorted tuple of pieces (sid, r, c),
# where sid indexes SHAPES (height, width) for unlabeled pieces, or a label
# whose shape is looked up in a dict for labeled pieces.
# ---------------------------------------------------------------------------


class Puzzle:
    def __init__(self, n, shape_of, labeled):
        self.n = n
        self.shape_of = shape_of      # sid -> (height, width)
        self.labeled = labeled

    def cells(self, sid, r, c):
        hgt, wid = self.shape_of[sid]
        return [(r + i, c + j) for i in range(hgt) for j in range(wid)]

    def legal(self, sid, r, c, occ):
        hgt, wid = self.shape_of[sid]
        if r < 0 or c < 0 or r + hgt > 2 or c + wid > self.n:
            return False
        return all(cell not in occ for cell in self.cells(sid, r, c))

    def occupied_except(self, state, idx):
        occ = set()
        for j, (sid, r, c) in enumerate(state):
            if j != idx:
                occ.update(self.cells(sid, r, c))
        return occ

    def canon(self, pieces):
        return tuple(sorted(pieces))

    def neighbors(self, state, mode):
        out = set()
        dirs = ((1, 0), (-1, 0), (0, 1), (0, -1))
        for idx, (sid, r, c) in enumerate(state):
            occ = self.occupied_except(state, idx)
            rest = state[:idx] + state[idx + 1:]
            if mode == 'U':
                targets = [(r + dr, c + dc) for dr, dc in dirs
                           if self.legal(sid, r + dr, c + dc, occ)]
            elif mode == 'S':
                targets = []
                for dr, dc in dirs:
                    rr, cc = r + dr, c + dc
                    while self.legal(sid, rr, cc, occ):
                        targets.append((rr, cc))
                        rr, cc = rr + dr, cc + dc
            else:  # finger
                seen = {(r, c)}
                dq = deque([(r, c)])
                while dq:
                    a, b = dq.popleft()
                    for dr, dc in dirs:
                        t = (a + dr, b + dc)
                        if t not in seen and self.legal(sid, t[0], t[1], occ):
                            seen.add(t)
                            dq.append(t)
                seen.discard((r, c))
                targets = list(seen)
            for rr, cc in targets:
                out.add(self.canon(rest + ((sid, rr, cc),)))
        return out


def enumerate_states(n, multiset, shape_of, zero_vac=None):
    """All placements of the pieces in `multiset` (list of sids, repeated for
    indistinguishable pieces) on the 2 x n board."""
    cells_order = [(r, c) for c in range(n) for r in range(2)]
    total_area = sum(shape_of[s][0] * shape_of[s][1] for s in multiset)
    vac = 2 * n - total_area
    remaining = {}
    for s in multiset:
        remaining[s] = remaining.get(s, 0) + 1
    results = []
    occ = set()
    placed = []

    def rec(pos, vac_left):
        while pos < len(cells_order) and cells_order[pos] in occ:
            pos += 1
        if pos == len(cells_order):
            if vac_left == 0 and all(v == 0 for v in remaining.values()):
                results.append(tuple(sorted(placed)))
            return
        r, c = cells_order[pos]
        if vac_left > 0:
            occ.add((r, c))
            rec(pos + 1, vac_left - 1)
            occ.discard((r, c))
        for s in list(remaining):
            if remaining[s] == 0:
                continue
            hgt, wid = shape_of[s]
            cl = [(r + i, c + j) for i in range(hgt) for j in range(wid)]
            if r + hgt > 2 or c + wid > n or any(x in occ for x in cl):
                continue
            remaining[s] -= 1
            occ.update(cl)
            placed.append((s, r, c))
            rec(pos + 1, vac_left)
            placed.pop()
            for x in cl:
                occ.discard(x)
            remaining[s] += 1

    if vac >= 0:
        rec(0, vac)
    return results


def build_graph(pz, states, mode):
    idx = {s: i for i, s in enumerate(states)}
    adj = [[] for _ in states]
    for i, s in enumerate(states):
        for t in pz.neighbors(s, mode):
            adj[i].append(idx[t])
    return adj


def components(adj):
    comp = [-1] * len(adj)
    comps = []
    for s in range(len(adj)):
        if comp[s] >= 0:
            continue
        cid = len(comps)
        comp[s] = cid
        lst = [s]
        dq = deque([s])
        while dq:
            u = dq.popleft()
            for w in adj[u]:
                if comp[w] < 0:
                    comp[w] = cid
                    lst.append(w)
                    dq.append(w)
        comps.append(lst)
    return comp, comps


def bfs(adj, s):
    dist = {s: 0}
    dq = deque([s])
    while dq:
        u = dq.popleft()
        for w in adj[u]:
            if w not in dist:
                dist[w] = dist[u] + 1
                dq.append(w)
    return dist


H, V = 0, 1
DOMINO = {H: (1, 2), V: (2, 1)}


def domino_states(n, h, v):
    return enumerate_states(n, [H] * h + [V] * v, DOMINO)


def signature(state):
    """Horizontal counts between consecutive vertical dominoes."""
    vcols = sorted(c for s, r, c in state if s == V)
    sig = [0] * (len(vcols) + 1)
    for s, r, c in state:
        if s == H:
            sig[sum(1 for x in vcols if x < c)] += 1
    return tuple(sig)


def binom(a, b):
    if b < 0 or a < 0 or a < b:
        return 0
    return math.comb(a, b)


# ---------------------------------------------------------------------------
# Formulas from the paper
# ---------------------------------------------------------------------------

def c_general(h, v, Q):
    tot = 0
    for j in range(0, min(v + 1, h, Q) + 1):
        if (h - j) % 2 == 0:
            tot += binom(v + 1, j) * binom((h - j) // 2 + v, v)
    return tot


def counts_q1(h, v):
    if h % 2 == 0:
        K = h // 2
        return (binom(K + v, v),
                (v + 1) * (binom(K + v + 2, v + 2) + binom(K + v + 1, v + 2)))
    K = (h + 1) // 2
    return (v + 1) * binom(K + v - 1, v), 2 * (v + 1) * binom(K + v + 1, v + 2)


def q1_coordinates(state, n):
    """Return the coordinates used in Sections 4-5 for a two-vacancy state."""
    occ = {}
    for s, r, c in state:
        for cell in [(r, c), (r, c + 1)] if s == H else [(0, c), (1, c)]:
            occ[cell] = (s, r, c)
    vac = [(r, c) for r in range(2) for c in range(n) if (r, c) not in occ]
    vcols = sorted(c for s, r, c in state if s == V)
    bounds = [-1] + vcols + [n]
    h = sum(1 for s, _, _ in state if s == H)
    if h % 2 == 0:
        # chamber index of the vacancies and gap indices in each row
        vc = vac[0][1]
        i = sum(1 for x in vcols if x < vc)
        lo, hi = bounds[i], bounds[i + 1]
        ks = []
        for j in range(len(bounds) - 1):
            a, b = bounds[j], bounds[j + 1]
            ks.append(sum(1 for s, r, c in state if s == H and r == 0 and a < c < b))
        gap = []
        for row in (0, 1):
            vcol = [c for r, c in vac if r == row][0]
            gap.append(sum(1 for s, r, c in state
                           if s == H and r == row and lo < c < vcol))
        return ('chain', tuple(ks), i, gap[0], gap[1])
    r = vac[0][0]
    cols = sorted(c for _, c in vac)
    i = sum(1 for x in vcols if x < cols[0])
    start = bounds[i] + 1          # first column (0-based) of active chamber
    x, y = cols[0] - start + 1, cols[1] - start + 1
    assert x % 2 == 1 and y % 2 == 0, (state, x, y)
    return ('tri', tuple(vcols), i, r, (x - 1) // 2, (y - 2) // 2)


def q1_formula(X, Y, mode):
    if X[0] == 'chain':
        _, ks, s, x, y = X
        _, ks2, t, u, w = Y
        if s == t:
            return abs(x - u) + abs(y - w)
        if s > t:
            s, t, x, y, u, w = t, s, u, w, x, y
        return (2 * ks[s] - x - y + (t - s) + 2 * sum(ks[s + 1:t]) + u + w)
    _, _, _, r, p, q = X
    _, _, _, s2, p2, q2 = Y
    a, b = q - p, q2 - p2
    g = max(0, p2 - q, p - q2)
    M = abs(p - p2) + abs(q - q2)
    if mode == 'U':
        return M if r == s2 else a + b + 2 * g + 1
    return min(M, a + b + g) if r == s2 else a + b + g + 1


def is_partial_cube(adj, nodes):
    """Bipartite + Djokovic-Winkler relation transitive.  Returns
    (is_partial_cube, number_of_Theta_classes or None)."""
    sub = set(nodes)
    dist = {u: bfs(adj, u) for u in nodes}
    colour = {}
    for u in nodes:
        colour[u] = dist[nodes[0]][u] % 2
    for u in nodes:
        for w in adj[u]:
            if colour[u] == colour[w]:
                return False, None
    edges = sorted({(min(u, w), max(u, w)) for u in nodes for w in adj[u] if w in sub})
    m = len(edges)

    def theta(e, f):
        a, b = e
        x, y = f
        return dist[a][x] + dist[b][y] != dist[a][y] + dist[b][x]

    cls = [-1] * m
    ncls = 0
    for i in range(m):
        if cls[i] >= 0:
            continue
        members = [j for j in range(m) if theta(edges[i], edges[j])]
        for j in members:
            if cls[j] >= 0 and cls[j] != ncls:
                return False, None
            cls[j] = ncls
        # transitivity: every member must relate to the same set
        mset = set(members)
        for j in members:
            if {t for t in range(m) if theta(edges[j], edges[t])} != mset:
                return False, None
        ncls += 1
    return True, ncls


def diameters(adj, comps):
    best = 0
    for lst in comps:
        for u in lst:
            d = bfs(adj, u)
            best = max(best, max(d.values()))
    return best


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_two_vacancies(nmax, log):
    log("== Two vacancies (Sections 4-7) ==")
    log(" n  C_n  S_n  max|C|  diamU  diamS  diamF  pairs-checked")
    rows = []
    for n in range(1, nmax + 1):
        Cn = Sn = maxc = dU = dS = dF = 0
        pairs = 0
        for h in range(0, n):
            v = n - 1 - h
            states = domino_states(n, h, v)
            pz = Puzzle(n, DOMINO, False)
            gU = build_graph(pz, states, 'U')
            gS = build_graph(pz, states, 'S')
            gF = build_graph(pz, states, 'F')
            assert gS == gF or all(set(a) == set(b) for a, b in zip(gS, gF)), \
                "straight and finger graphs differ"
            _, comps = components(gU)
            c_f, s_f = counts_q1(h, v)
            assert len(comps) == c_f, (n, h, v, len(comps), c_f)
            assert len(states) == s_f, (n, h, v, len(states), s_f)
            Cn += len(comps)
            Sn += len(states)
            coords = [q1_coordinates(s, n) for s in states]
            for lst in comps:
                maxc = max(maxc, len(lst))
                ok, dim = is_partial_cube(gU, lst)
                assert ok, ("not a partial cube", n, h, v)
                ecc = 0
                for u in lst:
                    du = bfs(gU, u)
                    ds = bfs(gS, u)
                    df = bfs(gF, u)
                    for w in lst:
                        assert du[w] == q1_formula(coords[u], coords[w], 'U')
                        assert ds[w] == q1_formula(coords[u], coords[w], 'S')
                        assert df[w] == ds[w]
                        pairs += 1
                    ecc = max(ecc, max(du.values()))
                    dU = max(dU, max(du.values()))
                    dS = max(dS, max(ds.values()))
                    dF = max(dF, max(df.values()))
                assert dim == ecc, ("isometric dimension != diameter", n, h, v)
        assert maxc == (n + 1) ** 2 // 4
        assert dU == dS == dF == n - 1
        rows.append((n, Cn, Sn, maxc, dU, dS, dF, pairs))
        log(f"{n:2d} {Cn:4d} {Sn:5d} {maxc:6d} {dU:6d} {dS:6d} {dF:6d} {pairs:9d}")
    # generating functions
    c = [0] * (nmax + 5)
    f = [0] * (nmax + 5)
    f[0], f[1] = 1, 1
    for i in range(2, nmax + 5):
        f[i] = f[i - 1] + f[i - 2]
    sq = [sum(f[i] * f[m - i] for i in range(m + 1)) for m in range(nmax + 5)]
    for n, Cn, Sn, *_ in rows:
        gfC = sq[n - 1] - (sq[n - 3] if n >= 3 else 0) - (sq[n - 4] if n >= 4 else 0)
        gfS = sum(sq[m] + (sq[m - 1] if m >= 1 else 0) for m in range(n))
        assert Cn == gfC and Sn == gfS, (n, Cn, gfC, Sn, gfS)
    log("generating functions for C_n and S_n: OK")
    return rows


def check_general(nmax, log):
    log("== Arbitrary vacancy counts (Section 3, Theorems 3.4, 3.7, 3.8) ==")
    log(" n  Q  N  D_U   QN  D_S  min(n-1,2N)  D_F   N")
    rows = []
    for n in range(1, nmax + 1):
        for Q in range(0, n + 1):
            N = n - Q
            DU = DS = DF = 0
            for h in range(0, N + 1):
                v = N - h
                states = domino_states(n, h, v)
                if not states:
                    continue
                pz = Puzzle(n, DOMINO, False)
                gU = build_graph(pz, states, 'U')
                comp, comps = components(gU)
                # signature is a complete invariant, count formula
                sig_of_comp = {}
                for lst in comps:
                    sigs = {signature(states[u]) for u in lst}
                    assert len(sigs) == 1, (n, Q, h, v)
                    sg = sigs.pop()
                    assert sg not in sig_of_comp, (n, Q, h, v, sg)
                    sig_of_comp[sg] = True
                    assert sum(x % 2 for x in sg) <= Q
                assert len(comps) == c_general(h, v, Q), (n, Q, h, v)
                if Q == 0:
                    assert all(len(lst) == 1 for lst in comps)
                    continue
                gS = build_graph(pz, states, 'S')
                gF = build_graph(pz, states, 'F')
                DU = max(DU, diameters(gU, comps))
                DS = max(DS, diameters(gS, comps))
                DF = max(DF, diameters(gF, comps))
            if 1 <= Q < n:
                assert DU == Q * N, (n, Q, DU)
                assert DS == min(n - 1, 2 * N), (n, Q, DS)
                assert DF == N, (n, Q, DF)
                rows.append((n, Q, N, DU, DS, DF))
                log(f"{n:2d} {Q:2d} {N:2d} {DU:4d} {Q*N:4d} {DS:4d} {min(n-1, 2*N):12d} {DF:4d} {N:3d}")
    log("component counts, chamber-signature invariant and exact diameters: OK")
    return rows


def check_four_vacancy_example(log):
    n = 4
    states = domino_states(n, 2, 0)
    pz = Puzzle(n, DOMINO, False)
    gU = build_graph(pz, states, 'U')
    _, comps = components(gU)
    assert len(states) == 11 and len(comps) == 1
    ok, _ = is_partial_cube(gU, comps[0])
    assert not ok
    log("four-vacancy example (2 horizontal dominoes on 2x4): 11 states, not a partial cube: OK")


def labeled_chamber(k):
    n = 2 * k
    shape_of = {i: (1, 2) for i in range(2 * k - 1)}
    states = enumerate_states(n, list(range(2 * k - 1)), shape_of)
    return Puzzle(n, shape_of, True), states


def vacancy_info(state, n):
    occ = set()
    for _, r, c in state:
        occ.update([(r, c), (r, c + 1)])
    vac = sorted((r, c) for r in range(2) for c in range(n) if (r, c) not in occ)
    return vac


def check_labeled(kmax, log):
    log("== Labeled chamber (Section 8) ==")
    log(" k  states  comps  comp-size  Delta_U  Delta_S  Delta_F")
    rows = []
    for k in range(2, kmax + 1):
        n = 2 * k
        pz, states = labeled_chamber(k)
        assert len(states) == k * (k + 1) * math.factorial(2 * k - 1)
        idx = {s: i for i, s in enumerate(states)}
        gU = build_graph(pz, states, 'U')
        gS = build_graph(pz, states, 'S')
        gF = build_graph(pz, states, 'F')
        assert all(set(a) == set(b) for a, b in zip(gS, gF))
        _, comps = components(gU)
        expect = k * (k + 1) * math.factorial(2 * k - 1) // 2
        assert len(comps) == 2 and all(len(c) == expect for c in comps)
        if k <= 3:
            DU = diameters(gU, comps)
            DS = diameters(gS, comps)
            DF = DS
            assert 2 * k * (k - 1) <= DU and k * (k - 1) <= DS
        else:
            DU = DS = DF = None
        # reflection pair of Section 10 (unlabeled endpoints adjacent)
        dX, dY, dXu, dYu = distortion_pair(k, pz, states, idx, gU, gS)
        assert dX >= 2 * (k - 1) ** 2 and dY >= (k - 1) ** 2
        rows.append((k, len(states), len(comps), expect, DU, DS, DF, dX, dY))
        log(f"{k:2d} {len(states):7d} {len(comps):5d} {expect:9d} "
            f"{'-' if DU is None else DU:>8} {'-' if DS is None else DS:>8} {'-' if DF is None else DF:>8}"
            f"   pair: d_U={dX} d_S={dY} (unlabeled distance 1)")
    return rows


def tile_from_coarse(k, coarse):
    """coarse: dict (row, pos) -> label or None (blank).  Return the aligned
    domino state."""
    pieces = []
    for (row, pos), lab in coarse.items():
        if lab is not None:
            pieces.append((lab, row, 2 * pos))
    return tuple(sorted(pieces))


def distortion_pair(k, pz, states, idx, gU, gS):
    labels = list(range(2 * k - 1))
    # X: blank at upper left, labels in reading order
    cells = [(0, j) for j in range(k)] + [(1, j) for j in range(k)]
    X = {}
    it = iter(labels)
    for cell in cells:
        X[cell] = None if cell == (0, 0) else next(it)
    # reflected target
    T = {(r, k - 1 - j): X[(r, j)] for (r, j) in X}
    if k % 2 == 0:
        T[(0, 0)], T[(1, 0)] = T[(1, 0)], T[(0, 0)]
    # move blank down at right boundary, then left along lower row
    br, bc = (0, k - 1)
    assert T[(br, bc)] is None
    T[(0, k - 1)], T[(1, k - 1)] = T[(1, k - 1)], None
    for j in range(k - 1, 0, -1):
        T[(1, j)], T[(1, j - 1)] = T[(1, j - 1)], None
    sx = tile_from_coarse(k, X)
    sy = tile_from_coarse(k, T)
    iX, iY = idx[sx], idx[sy]
    dU = bfs(gU, iX)
    dS = bfs(gS, iX)
    assert iY in dU, "pair not reachable"
    return dU[iY], dS[iY], 1, 1


def check_tile_reduction(kmax, log):
    """d_U on aligned states equals weighted (2,1) tile distance and
    d_S equals unweighted tile distance (Theorem 6.1)."""
    import heapq
    for k in range(2, kmax + 1):
        pz, states = labeled_chamber(k)
        n = 2 * k
        idx = {s: i for i, s in enumerate(states)}
        gU = build_graph(pz, states, 'U')
        gS = build_graph(pz, states, 'S')
        aligned = []
        for i, s in enumerate(states):
            vac = vacancy_info(s, n)
            if vac[0][0] == vac[1][0] and vac[1][1] == vac[0][1] + 1 and vac[0][1] % 2 == 0:
                aligned.append(i)

        def coarse_of(s):
            g = {}
            for lab, r, c in s:
                g[(r, c // 2)] = lab
            return tuple(g.get((r, j)) for r in range(2) for j in range(k))

        cid = {coarse_of(states[i]): i for i in aligned}
        assert len(cid) == len(aligned)

        def tile_nb(t):
            b = t.index(None)
            r, j = divmod(b, k)
            for rr, jj, w in ((1 - r, j, 1), (r, j - 1, 2), (r, j + 1, 2)):
                if 0 <= jj < k:
                    o = rr * k + jj
                    lst = list(t)
                    lst[b], lst[o] = lst[o], None
                    yield tuple(lst), w

        sources = aligned if k <= 3 else aligned[:3]
        for src in sources:
            t0 = coarse_of(states[src])
            for weighted, g in ((True, gU), (False, gS)):
                dist = {t0: 0}
                pq = [(0, 0, t0)]
                cnt = 0
                while pq:
                    d, _, t = heapq.heappop(pq)
                    if d > dist[t]:
                        continue
                    for t2, w in tile_nb(t):
                        nd = d + (w if weighted else 1)
                        if nd < dist.get(t2, 10 ** 9):
                            dist[t2] = nd
                            cnt += 1
                            heapq.heappush(pq, (nd, cnt, t2))
                fine = bfs(g, src)
                for t, d in dist.items():
                    assert fine[cid[t]] == d, (k, weighted)
                assert len(dist) == sum(1 for i in aligned if i in fine)
        log(f"tile-puzzle reduction, k={k}: d_U = d_(2,1) and d_S = d_tile on aligned states "
            f"({len(sources)} sources x {len(aligned)} targets): OK")


def check_mixed_example(log):
    # shapes: id -> (height, width); labeled pieces in Example 9.3
    shape_of = {0: (1, 2), 1: (1, 4), 2: (1, 3), 3: (1, 3), 4: (2, 2),
                5: (1, 6), 6: (1, 2), 7: (1, 2), 8: (1, 2)}
    start = ((0, 0, 1), (1, 0, 3), (2, 1, 1), (3, 1, 4), (4, 0, 7),
             (5, 0, 9), (6, 1, 9), (7, 1, 11), (8, 1, 13))
    pz = Puzzle(15, shape_of, True)
    for mode in 'USF':
        seen = {start}
        dq = deque([start])
        adj = {}
        while dq:
            s = dq.popleft()
            adj[s] = list(pz.neighbors(s, mode))
            for t in adj[s]:
                if t not in seen:
                    seen.add(t)
                    dq.append(t)
        nodes = list(seen)
        ix = {s: i for i, s in enumerate(nodes)}
        g = [[ix[t] for t in adj[s]] for s in nodes]
        diam = diameters(g, [list(range(len(nodes)))])
        assert len(nodes) == 17 and diam == 9, (mode, len(nodes), diam)
    log("mixed-rectangle example (2x15 board): 17 states, diameter 9 under U, S, F: OK")


def unlabel(state):
    return tuple(sorted((H, r, c) for _, r, c in state))


def heuristic_quality(kmax, log, nsrc_large=24):
    """Quality of the label-erasure heuristic h(X,Y) = d(rho X, rho Y) in the
    labeled horizontal chamber.  For k <= 3 all ordered pairs inside one
    component are used; for larger k the first `nsrc_large` states of the
    sorted state list serve as sources (deterministic, no random sampling)."""
    log("== Label-erasure heuristic in the labeled chamber ==")
    log(" k  metric  pairs      mean d   mean h   mean h/d  max d/h   h=d share")
    rows = []
    for k in range(2, kmax + 1):
        n = 2 * k
        pz, states = labeled_chamber(k)
        uz = Puzzle(n, DOMINO, False)
        ustates = domino_states(n, 2 * k - 1, 0)
        uidx = {s: i for i, s in enumerate(ustates)}
        rho = [uidx[unlabel(s)] for s in states]
        for mode in 'US':
            g = build_graph(pz, states, mode)
            ug = build_graph(uz, ustates, mode)
            udist = [bfs(ug, i) for i in range(len(ustates))]
            _, comps = components(g)
            comp0 = sorted(comps[0])
            sources = comp0 if k <= 3 else comp0[:nsrc_large]
            tot_d = tot_h = tot_r = 0.0
            cnt = eq = 0
            worst = 0.0
            for src in sources:
                d = bfs(g, src)
                for t, dt in d.items():
                    if t == src:
                        continue
                    hv = udist[rho[src]][rho[t]]
                    tot_d += dt
                    tot_h += hv
                    tot_r += hv / dt
                    cnt += 1
                    eq += (hv == dt)
                    if hv > 0:
                        worst = max(worst, dt / hv)
            row = (k, mode, cnt, tot_d / cnt, tot_h / cnt, tot_r / cnt, worst, eq / cnt)
            rows.append(row)
            log(f"{k:2d}  {mode:>6}  {cnt:8d}  {row[3]:8.2f} {row[4]:8.2f} {row[5]:9.3f} {row[6]:8.1f} {row[7]:10.3f}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nmax-q1', type=int, default=10)
    ap.add_argument('--nmax-gen', type=int, default=8)
    ap.add_argument('--kmax-lab', type=int, default=4)
    ap.add_argument('--out', default='verification_log.txt')
    ap.add_argument('--heuristic-only', action='store_true')
    args = ap.parse_args()
    lines = []

    def log(s):
        print(s)
        sys.stdout.flush()
        lines.append(s)

    t0 = time.time()
    if args.heuristic_only:
        heuristic_quality(args.kmax_lab, log)
        log(f"done in {time.time() - t0:.1f} s")
        with open(args.out, 'w', encoding='utf-8') as fh:
            fh.write('\n'.join(lines) + '\n')
        return
    check_two_vacancies(args.nmax_q1, log)
    check_general(args.nmax_gen, log)
    check_four_vacancy_example(log)
    check_tile_reduction(min(args.kmax_lab, 4), log)
    check_labeled(args.kmax_lab, log)
    check_mixed_example(log)
    log(f"all checks passed in {time.time() - t0:.1f} s")
    with open(args.out, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
