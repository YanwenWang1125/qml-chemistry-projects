"""Composite-reduction search (Jeon et al. Sec. IV) with a checker on every child.

A composite reduction = zero or one enabling rewrite (restricted commutation or target swap)
that exposes a new Bell pair, an optional basis swap on that pair, then one Bell reduction.
Enabling rewrites are used only when they immediately expose a Bell pair (Jeon's setting, "L=1").
Every child is checked against the root (verify.check); failing children are counted and dropped.
Unlike Jeon, we also deduplicate with an approximate canonical form.
"""
from collections import Counter, deque

from .ftlib import resources
from .rules import basis_swap, bell_pairs, bell_reduction, restricted_commutation, target_swap
from .verify import check


def canonical(g):
    """Approximate canonical form: ASAP layers, sort within layers, relabel ancillas by first use."""
    def layered(ops):
        last, keyed = {}, []
        for op in ops:
            qs = op[1:]
            layer = max((last.get(q, -1) for q in qs), default=-1) + 1
            for q in qs:
                last[q] = layer
            keyed.append((layer, op))
        return [op for _, op in sorted(keyed, key=lambda t: (t[0], t[1]))]

    ops = layered(g["ops"])
    mapping, nxt = {}, 7
    for op in ops:
        for q in op[1:]:
            if q >= 7 and q not in mapping:
                mapping[q] = nxt
                nxt += 1
    rl = lambda op: tuple([op[0]] + [mapping.get(x, x) for x in op[1:]])
    ops = layered([rl(op) for op in ops])
    syn = tuple(sorted(mapping[q] for q in g["syndrome"]))
    flags = tuple(sorted(mapping[q] for q in g["flags"]))
    return tuple(ops), syn, flags


def _wkey(g):
    w = {}
    for op in g["ops"]:
        for q in op[1:]:
            w.setdefault(q, []).append(op)
    return tuple(sorted((q, tuple(s)) for q, s in w.items()))


def composite_reductions(g, L=1):
    """Yield (child, steps) for every composite reduction of g.

    L=1 is Jeon's setting: at most one enabling rewrite, and only if it exposes a new Bell pair.
    L>1 allows up to L enabling rewrites of any kind before the reduction (App. C needs 2: one of
    them reorders a flag and exposes nothing).
    """
    seen_pairs = {(a, b) for a, b, _ in bell_pairs(g)}
    starts, layer, seen = [(g, [])], [(g, [])], {_wkey(g)}
    for _ in range(L):
        nxt = []
        for h0, s0 in layer:
            for r in (restricted_commutation, target_swap):
                for h, desc in r(h0):
                    k = _wkey(h)
                    if k in seen:
                        continue
                    if L == 1 and not ({(a, b) for a, b, _ in bell_pairs(h)} - seen_pairs):
                        continue
                    seen.add(k)
                    nxt.append((h, s0 + [desc]))
        starts += nxt
        layer = nxt
    for h, steps in starts:
        for a, b, i in bell_pairs(h):
            if steps == [] or L > 1 or (a, b) not in seen_pairs:
                for child, d in bell_reduction(h, (a, b, i)):
                    yield child, steps + [d]
                for sw, ds in basis_swap(h, (a, b, i)):
                    pb = [p for p in bell_pairs(sw) if {p[0], p[1]} == {a, b}]
                    for child, d in bell_reduction(sw, pb[0]) if pb else []:
                        yield child, steps + [ds, d]


def bfs(root, dedup=True, max_nodes=100000, check=check, L=1):
    """Breadth-first search over composite reductions.

    Returns a dict with all accepted circuits (by canonical key), BFS node counts, rejected-rewrite
    rewrites and their reasons, and the terminal circuits (no accepted child).
    """
    q = deque([(root, [])])
    accepted = {canonical(root): (root, [])}
    nodes, rejected, terminals, generated, verdict = 0, Counter(), [], 0, {}
    while q and nodes < max_nodes:
        g, path = q.popleft()
        nodes += 1
        has_child = False
        for child, steps in composite_reductions(g, L):
            generated += 1
            key = canonical(child)
            if key not in verdict:  # the checker only depends on the circuit
                verdict[key] = check(root, child)
            ok, why = verdict[key]
            if not ok:
                rejected[why] += 1  # counts only: storing paths ran out of memory (E1.3-fix)
                continue
            has_child = True
            if dedup and key in accepted:
                continue
            accepted.setdefault(key, (child, path + steps))
            q.append((child, path + steps))
        if not has_child:
            terminals.append((g, path))
    return {"accepted": accepted, "nodes": nodes, "rejected": rejected, "terminals": terminals,
            "generated": generated, "checked": len(verdict)}


def summary_rows(result):
    rows = []
    for key, (g, path) in result["accepted"].items():
        r = resources(g["ops"])
        rows.append({**r, "steps": len([s for s in path if s.startswith("Bell")]), "path": path, "gadget": g})
    rows.sort(key=lambda r: (r["ancillas"], r["cnot_depth"], r["cnots"]))
    return rows
