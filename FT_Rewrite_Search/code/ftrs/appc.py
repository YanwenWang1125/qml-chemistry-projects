"""Jeon et al. Appendix C (PDF page 13): the five circuits of the Steane-gadget derivation.

Read gate by gate from a 400 dpi render. Data q0..q6 are the top seven wires; ancillas start at q7.
Columns where one vertical line carries two CNOTs are split by the pattern of the neighbouring
figures; the one ambiguous column (F3, x~323) is given in both readings, F3a and F3b.
"""


def _gadget(name, anc, plus, flag, cx):
    ops = [("R", q) for q in anc] + [("H", q) for q in plus] + [("CX", c, t) for c, t in cx]
    ops += [("H", flag)] + [("M", q) for q in anc]
    return {"name": name, "kind": "Z", "targets": [(0, 1, 2, 3), (1, 2, 4, 5), (2, 3, 5, 6)],
            "syndrome": [q for q in anc if q != flag], "flags": [flag], "ops": ops}


F0 = _gadget("AppC F0", range(7, 15), (8, 9, 10, 12, 14), 14, [
    (8, 7), (12, 11), (9, 11), (8, 11), (4, 11), (10, 13), (10, 7), (12, 13), (0, 7), (9, 10),
    (2, 9), (6, 13), (14, 8), (1, 8), (14, 10), (3, 10), (14, 12), (5, 12)])

F1 = _gadget("AppC F1", range(7, 14), (8, 9, 11, 13), 13, [
    (0, 7), (8, 10), (11, 10), (9, 12), (9, 7), (11, 12), (7, 10), (4, 10), (8, 9), (2, 8),
    (13, 9), (6, 12), (3, 9), (13, 7), (1, 7), (13, 11), (5, 11)])

F2 = _gadget("AppC F2", range(7, 13), (8, 10, 12), 12, [
    (0, 7), (4, 9), (10, 9), (8, 11), (8, 7), (10, 11), (7, 9), (9, 8), (2, 9), (6, 11),
    (12, 7), (1, 7), (12, 10), (5, 10), (12, 8), (3, 8)])

_F3_head = [(0, 7), (5, 9), (8, 10)]
_F3_tail = [(9, 10), (4, 9), (6, 10), (7, 9), (9, 8), (2, 9), (11, 7), (1, 7), (11, 8), (3, 8)]
F3a = _gadget("AppC F3a", range(7, 12), (8, 11), 11, _F3_head + [(8, 7), (11, 9)] + _F3_tail)
F3b = _gadget("AppC F3b", range(7, 12), (8, 11), 11, _F3_head + [(8, 9), (11, 7)] + _F3_tail)

F4 = _gadget("AppC F4", range(7, 11), (10,), 10, [
    (0, 7), (5, 8), (6, 9), (10, 8), (8, 9), (4, 8), (9, 7), (3, 9), (7, 8), (8, 9), (2, 8),
    (10, 7), (1, 7), (10, 9)])


def wires(g):
    """Per-qubit op sequences: two op lists with equal wires are the same circuit (same DAG)."""
    w = {}
    for op in g["ops"]:
        for q in op[1:]:
            w.setdefault(q, []).append(op)
    return w


def isomorphic(g, h, n_data=7):
    """Same circuit up to commuting gates and relabelling ancillas (data qubits fixed)."""
    wg, wh = wires(g), wires(h)
    ag = sorted(q for q in wg if q >= n_data)
    ah = sorted(q for q in wh if q >= n_data)
    if len(ag) != len(ah) or len(g["ops"]) != len(h["ops"]):
        return False
    sig = lambda w, q: (len(w[q]), tuple(op[0] for op in w[q]))

    def rl(op, m):
        return tuple([op[0]] + [m.get(x, x) for x in op[1:]])

    def go(m, rest):
        if not rest:
            return all([rl(op, m) for op in wg[q]] == wh[m.get(q, q)] for q in wg) \
                and {m[q] for q in g["flags"]} == set(h["flags"])
        q = rest[0]
        for r in ah:
            if r in m.values() or sig(wg, q) != sig(wh, r):
                continue
            m[q] = r
            # prune: ops on q whose partners are already mapped must match
            if all(rl(op, m) == wh[r][k] for k, op in enumerate(wg[q])
                   if all(x < n_data or x in m for x in op[1:])) and go(m, rest[1:]):
                return True
            del m[q]
        return False

    return go({}, ag)
