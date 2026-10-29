"""Jeon et al.'s four fault-equivalent rewrite rules, written on the op-list circuits of ftlib.

Source: arXiv:2610.09749, Sec. III, Eqs. (5)-(9) and the figure on PDF page 4.
Each rule function returns a list of (new_gadget, description). A rule only fires when its
preconditions hold; whether the result is really fault-equivalent is checked separately
(verify.py), as the experiment plan requires for every rewrite.

Conventions (as in ftlib): ("R", q) prepares |0>, ("H", q), ("CX", c, t), ("M", q) measures Z.
A qubit prepared in |+> is R then H; an X measurement is H then M. Flag qubits never take part
in a Bell reduction or target swap (Jeon: "a verification flag is never merged into the
preparation circuit"); restricted commutation may act on a flag's first two CNOTs (App. C, F0->F1).
"""


def timeline(ops, q):
    return [i for i, op in enumerate(ops) if q in op[1:]]


def is_ancilla(g, q):
    return q in g["syndrome"]  # rules act only on syndrome (cat) ancillas, never on flags


def fresh_plus_first_gate(ops, q):
    """Index of q's first gate after an R, H preparation, or None."""
    tl = timeline(ops, q)
    if len(tl) >= 3 and ops[tl[0]] == ("R", q) and ops[tl[1]] == ("H", q):
        return tl[2]
    return None


def fresh_zero_first_gate(ops, q):
    """Index of q's first gate after an R preparation not followed by H, or None."""
    tl = timeline(ops, q)
    if len(tl) >= 2 and ops[tl[0]] == ("R", q) and ops[tl[1]][0] != "H":
        return tl[1]
    return None


def bell_pairs(g):
    """All (a, b, i): op i is CX a->b, the first gate of a fresh |+> qubit a and of a fresh |0> qubit b."""
    ops = g["ops"]
    out = []
    for i, op in enumerate(ops):
        if op[0] != "CX":
            continue
        a, b = op[1], op[2]
        if not (is_ancilla(g, a) and is_ancilla(g, b)):
            continue
        if fresh_plus_first_gate(ops, a) == i and fresh_zero_first_gate(ops, b) == i:
            out.append((a, b, i))
    return out


def _touches_between(ops, qubits, i, j):
    return any(set(op[1:]) & set(qubits) for op in ops[i + 1:j])


def _with(g, ops, syndrome=None, note=""):
    return {**g, "ops": ops, "syndrome": list(syndrome if syndrome is not None else g["syndrome"]),
            "name": note}


# --- Rule 1: restricted commutation, Eqs. (5) and (6) -------------------------------------

def restricted_commutation(g):
    """Eqs. (5)/(6), read as spider fusion: a fresh |+> qubit and the run of CNOTs it controls
    right after preparation are one Z spider (a fresh |0> and the CNOTs targeting it, one X
    spider), so the run may be put in any order. Jeon's two-gate swap is the case of a run of two;
    App. C (F1 -> F2) reorders a flag's run of three."""
    ops = g["ops"]
    out = []
    for q in g["syndrome"] + g["flags"]:  # flags may commute; they are only kept out of reductions
        tl = timeline(ops, q)
        for eq, i, role in (("5", fresh_plus_first_gate(ops, q), "c"),
                            ("6", fresh_zero_first_gate(ops, q), "t")):
            if i is None:
                continue
            k0 = tl.index(i)
            run = []
            for k in tl[k0:]:
                if _role(ops[k], q) != role:
                    break
                run.append(k)
            if len(run) < 2:
                continue
            from itertools import islice, permutations
            for perm in islice(permutations(run), 1, MAX_RUN_ORDERS):  # skip the identity
                new = _reorder_on_wire(ops, q, run, list(perm))
                if new is not None:
                    desc = " ".join(f"{ops[k][1]}>{ops[k][2]}" for k in perm)
                    out.append((_with(g, new), f"commute ({eq}) on q{q}: {desc}"))
    return out


def _reorder_on_wire(ops, q, run, perm):
    """Put q's run of ops in the order perm; other wires keep theirs. None if that makes a cycle."""
    order_on = {}
    for k, op in enumerate(ops):
        for x in op[1:]:
            order_on.setdefault(x, []).append(k)
    seq = order_on[q]
    s = seq.index(run[0])
    seq[s:s + len(run)] = perm
    return _topo(dict(enumerate(ops)), order_on)


# --- Rule 2: basis swap, Eq. (7) ------------------------------------------------------------

def basis_swap(g, pair=None):
    """CX a->b from |+>(a) to |0>(b) becomes CX b->a from |+>(b) to |0>(a)."""
    ops = g["ops"]
    out = []
    for a, b, i in ([pair] if pair else bell_pairs(g)):
        h_a = timeline(ops, a)[1]
        new = []
        for k, op in enumerate(ops):
            if k == h_a:
                continue
            if k == i:
                new.append(("H", b))
                new.append(("CX", b, a))
                continue
            new.append(op)
        out.append((_with(g, new), f"basis swap on Bell pair (q{a}, q{b})"))
    return out


# --- Rule 3: target swap, Eq. (8) -----------------------------------------------------------

def target_swap(g):
    """q2 |+>, q3 |0>, CX q2->q3; then CX q1->q2 from a fresh |+> q1 may be redirected to q3."""
    ops = g["ops"]
    out = []
    for q2, q3, i in bell_pairs(g):
        tl2 = timeline(ops, q2)
        k = tl2.index(i)
        if k + 1 >= len(tl2):
            continue
        j = tl2[k + 1]
        op = ops[j]
        if op[0] != "CX" or op[2] != q2:
            continue
        q1 = op[1]
        if not is_ancilla(g, q1) or fresh_plus_first_gate(ops, q1) != j:
            continue
        if _touches_between(ops, [q3], i, j):
            continue
        new = list(ops)
        new[j] = ("CX", q1, q3)
        out.append((_with(g, new), f"target swap: CX{q1}->{q2} becomes CX{q1}->{q3}"))
    return out


# --- Rule 4: Bell reduction, Eq. (9) --------------------------------------------------------

def bell_reduction(g, pair=None):
    """Remove b of a Bell pair (a |+>, b |0>): one wire from |0>, then V^T (b's ops reversed), then U.

    Preconditions (Jeon Sec. III D): after the Bell CNOT no gate acts on both a and b; b's
    measurement only enters the syndrome parity (b is a syndrome qubit, not a flag).
    For CNOT and H, the transpose on the eliminated wire leaves each gate unchanged, so V^T is
    b's gate list in reverse order, relabelled to a.
    """
    ops = g["ops"]
    out = []
    for a, b, i in ([pair] if pair else bell_pairs(g)):
        tl_a, tl_b = timeline(ops, a), timeline(ops, b)
        if any(k > i and a in ops[k][1:] and b in ops[k][1:] for k in tl_a):
            continue
        if ops[tl_b[-1]] != ("M", b):
            continue
        v = [k for k in tl_b if k > i and k != tl_b[-1]]  # b's gates after the Bell CNOT
        relabel = lambda op: tuple([op[0]] + [a if x == b else x for x in op[1:]])
        drop = {tl_b[0], tl_a[1], tl_b[-1], *v}  # R(b), H(a), M(b), V
        syn = [q for q in g["syndrome"] if q != b]
        # Circuit extraction from the diagram: every other wire keeps its gate order; wire a
        # carries R, V^T, then U. Topologically sort; a cycle means the diagram is not extractable.
        keep = [k for k in range(len(ops)) if k not in drop - set(v) and k != i]
        node = {k: relabel(ops[k]) for k in keep}
        order_on = {}
        for k in keep:
            for q in node[k][1:]:
                if q != a:
                    order_on.setdefault(q, []).append(k)
        wire_a = list(reversed(v)) + [k for k in tl_a if k > i]
        seen = set()
        # On wire a, a run of consecutive CNOTs that all target a (or all control from a) is one
        # fused spider in the diagram, so any order of the run is a valid extraction (App. C uses this).
        for perm in _run_orders(wire_a, lambda k: _role(node[k], a)):
            order_on[a] = [tl_a[0]] + perm
            new = _topo(node, order_on)
            if new is None:
                continue
            key = tuple(sorted(_wire_key(new).items()))
            if key in seen:
                continue
            seen.add(key)
            out.append((_with(g, new, syn), f"Bell reduction ({a},{b})"))
    return out


MAX_RUN_ORDERS = 720


def _role(op, q):
    if op[0] != "CX":
        return None
    return "t" if op[2] == q else "c"


def _run_orders(seq, role):
    """All orders of seq that permute only within maximal runs of equal (non-None) role."""
    from itertools import permutations, product
    runs, cur = [], []
    for k in seq:
        if cur and role(k) is not None and role(k) == role(cur[-1]):
            cur.append(k)
        else:
            if cur:
                runs.append(cur)
            cur = [k]
    if cur:
        runs.append(cur)
    choices = [list(permutations(r)) if role(r[0]) is not None else [tuple(r)] for r in runs]
    for n, combo in enumerate(product(*choices)):
        if n >= MAX_RUN_ORDERS:
            return
        yield [k for part in combo for k in part]


def _wire_key(ops):
    w = {}
    for op in ops:
        for q in op[1:]:
            w.setdefault(q, []).append(op)
    return {q: tuple(s) for q, s in w.items()}


def _topo(node, order_on):
    """Kahn sort of the gate DAG given per-wire orders; ties broken by original position."""
    import heapq
    succ, indeg = {n: set() for n in node}, {n: 0 for n in node}
    for seq in order_on.values():
        for x, y in zip(seq, seq[1:]):
            if y not in succ[x]:
                succ[x].add(y)
                indeg[y] += 1
    heap = [n for n in node if indeg[n] == 0]
    heapq.heapify(heap)
    out = []
    while heap:
        n = heapq.heappop(heap)
        out.append(node[n])
        for m in succ[n]:
            indeg[m] -= 1
            if indeg[m] == 0:
                heapq.heappush(heap, m)
    return out if len(out) == len(node) else None
