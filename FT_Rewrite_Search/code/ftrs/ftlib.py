"""Small fault-propagation toolkit shared by the FT_Rewrite_Search lessons.

Data qubits are 0..n_data-1 of a CSS code (a `Code` object; default STEANE, the [[7,1,3]] code,
so qubits 0-6). Ancillas start at n_data, matching the q7, q8, ... labels in Jeon et al.
(arXiv:2610.09749, Fig. 1 and App. B). A gadget may carry its code under the key "code".

A circuit is a list of ops:
    ("R", q)      prepare |0>
    ("H", q)      Hadamard
    ("CX", c, t)  CNOT, control c, target t
    ("M", q)      measure Z (all measurements are at the end)
A |+> preparation is written R then H, and an X measurement is H then M, as in the paper.

A Pauli on n qubits is a pair of int bitmasks (x, z): bit q of x is set when the Pauli has
an X or Y on qubit q, bit q of z when it has a Z or Y.
"""
import itertools
from collections import Counter

import numpy as np
import stim

STEANE_SUPPORTS = [(0, 1, 2, 3), (1, 2, 4, 5), (2, 3, 5, 6)]
N_DATA = 7


def mask(qubits):
    m = 0
    for q in qubits:
        m |= 1 << q
    return m


def popcount(m):
    return bin(m).count("1")


def stabilizer_group(x_supports, z_supports):
    """All products of the given X-type and Z-type generators, as (x, z) masks."""
    gens = [(mask(s), 0) for s in x_supports] + [(0, mask(s)) for s in z_supports]
    group = []
    for bits in itertools.product([0, 1], repeat=len(gens)):
        x = z = 0
        for b, (gx, gz) in zip(bits, gens):
            if b:
                x ^= gx
                z ^= gz
        group.append((x, z))
    return group


STEANE_GROUP = stabilizer_group(STEANE_SUPPORTS, STEANE_SUPPORTS)


def _span(masks):
    out = {0}
    for m in masks:
        out |= {v ^ m for v in out}
    return sorted(out)


class Code:
    """A CSS code on data qubits 0..n_data-1: X- and Z-type stabilizer supports, logical supports.

    weight_mod and canonical_mod give the same values as the module functions with the full
    group list, but run over the X and Z halves of the group with numpy (2^16 elements for
    [[17,1,5]] would be too slow as a Python loop per fault) and cache their results.
    """

    def __init__(self, name, n_data, x_supports, z_supports, logical_x=(), logical_z=()):
        self.name, self.n_data = name, n_data
        self.x_supports, self.z_supports = [tuple(s) for s in x_supports], [tuple(s) for s in z_supports]
        self.logical_x, self.logical_z = [tuple(s) for s in logical_x], [tuple(s) for s in logical_z]
        self.data_mask = (1 << n_data) - 1
        self._sx = np.array(_span([mask(s) for s in self.x_supports]), dtype=np.int64)
        self._sz = np.array(_span([mask(s) for s in self.z_supports]), dtype=np.int64)
        v = np.arange(1 << n_data, dtype=np.int64)
        self._pc = sum((v >> b) & 1 for b in range(n_data)).astype(np.uint8)
        self._cache = {}

    @property
    def group(self):
        return [(int(a), int(b)) for a in self._sx for b in self._sz]

    def _weights(self, x, z):
        return self._pc[(x ^ self._sx)[:, None] | (z ^ self._sz)[None, :]]

    def _rep(self, x, z):
        """(weight, x', z') of the coset representative with the smallest (weight, x', z')."""
        key = (x, z)
        if key not in self._cache:
            w = self._weights(x, z)
            i, j = np.nonzero(w == w.min())
            rx, rz = min(zip((x ^ self._sx[i]).tolist(), (z ^ self._sz[j]).tolist()))
            self._cache[key] = (int(w.min()), rx, rz)
        return self._cache[key]

    def weight_mod(self, x, z):
        return self._rep(x, z)[0]

    def canonical_mod(self, x, z):
        return self._rep(x, z)[1:]


STEANE = Code("[[7,1,3]] Steane", N_DATA, STEANE_SUPPORTS, STEANE_SUPPORTS,
              logical_x=[range(7)], logical_z=[range(7)])


def code_of(gadget, code=None):
    return code or gadget.get("code") or STEANE


def weight_mod(x, z, group):
    """Minimum weight of the Pauli (x, z) times any element of the stabilizer group."""
    return min(popcount((x ^ sx) | (z ^ sz)) for sx, sz in group)


def canonical_mod(x, z, group):
    """A fixed representative of the coset (x, z)*group: the smallest (weight, x, z)."""
    return min((popcount((x ^ sx) | (z ^ sz)), x ^ sx, z ^ sz) for sx, sz in group)[1:]


def pauli_str(x, z, qubits):
    out = []
    for q in qubits:
        bx, bz = (x >> q) & 1, (z >> q) & 1
        if bx or bz:
            out.append({(1, 0): "X", (0, 1): "Z", (1, 1): "Y"}[(bx, bz)] + str(q))
    return "".join(out) if out else "I"


def propagate(x, z, ops):
    """Push a Pauli forward through the Clifford ops (R and M are skipped)."""
    for op in ops:
        if op[0] == "H":
            q = op[1]
            bx, bz = (x >> q) & 1, (z >> q) & 1
            if bx != bz:  # H swaps X and Z on q
                x ^= 1 << q
                z ^= 1 << q
        elif op[0] == "CX":
            c, t = op[1], op[2]
            if (x >> c) & 1:  # X on control spreads to target
                x ^= 1 << t
            if (z >> t) & 1:  # Z on target spreads to control
                z ^= 1 << c
    return x, z


TWO_QUBIT_PAULIS = [p for p in itertools.product("IXYZ", repeat=2) if p != ("I", "I")]
ONE_QUBIT_BITS = {"I": (0, 0), "X": (1, 0), "Y": (1, 1), "Z": (0, 1)}


def located_faults(ops):
    """Every single fault of the circuit-level model used in these lessons.

    after R q: X on q (preparation error)
    after H q: X, Y or Z on q
    after CX c t: any of the 15 non-identity two-qubit Paulis on (c, t)
    before M q: X on q (measurement error)
    Yields (location = index of the faulty op, label, index of the first op the fault passes
    through, x, z). Faults with the same location are alternatives, not a multi-fault event.
    """
    for i, op in enumerate(ops):
        if op[0] == "R":
            yield i, f"prep X on q{op[1]}", i + 1, 1 << op[1], 0
        elif op[0] == "H":
            for p in "XYZ":
                bx, bz = ONE_QUBIT_BITS[p]
                yield i, f"{p} on q{op[1]} after H", i + 1, bx << op[1], bz << op[1]
        elif op[0] == "CX":
            c, t = op[1], op[2]
            for pc, pt in TWO_QUBIT_PAULIS:
                xc, zc = ONE_QUBIT_BITS[pc]
                xt, zt = ONE_QUBIT_BITS[pt]
                yield (i, f"{pc}{pt} on (q{c},q{t}) after CX{c}->{t}", i + 1,
                       (xc << c) | (xt << t), (zc << c) | (zt << t))
        elif op[0] == "M":
            yield i, f"measurement flip on q{op[1]}", i, 1 << op[1], 0


def single_faults(ops):
    """located_faults without the location: (label, first op index, x, z)."""
    for _, label, start, x, z in located_faults(ops):
        yield label, start, x, z


def fault_columns(ops):
    """Propagation matrix of plan Sec. 2.4: each elementary fault pushed once to the end.

    Returns {location: [(label, x, z)]} with x, z the end-of-circuit Pauli on all qubits.
    """
    cols = {}
    for loc, label, start, x, z in located_faults(ops):
        cols.setdefault(loc, []).append((label, *propagate(x, z, ops[start:])))
    return cols


def faults_upto(ops, t):
    """Every set of 1..t faults at distinct locations, in order of size.

    Propagation is linear, so a set's end-of-circuit Pauli is the XOR of its columns.
    Yields (s, labels, x, z).
    """
    cols = fault_columns(ops)
    locs = sorted(cols)
    for s in range(1, t + 1):
        for combo in itertools.combinations(locs, s):
            for choice in itertools.product(*(cols[l] for l in combo)):
                x = z = 0
                for _, cx, cz in choice:
                    x ^= cx
                    z ^= cz
                yield s, tuple(c[0] for c in choice), x, z


def ft_report(gadget, t=1, code=None):
    """FT_t (plan Sec. 0.3): every unflagged set of s <= t faults leaves weight <= s mod S."""
    code = code_of(gadget, code)
    flags = mask(gadget["flags"])
    sets, undetected, hist, bad = Counter(), Counter(), Counter(), []
    for s, labels, x, z in faults_upto(gadget["ops"], t):
        sets[s] += 1
        if x & flags:  # a Z measurement of a flag flips on X or Y
            continue
        undetected[s] += 1
        dx, dz = x & code.data_mask, z & code.data_mask
        w = code.weight_mod(dx, dz)
        hist[(s, w)] += 1
        if w > s:
            bad.append({"faults": labels, "data": pauli_str(dx, dz, range(code.n_data)),
                        "weight_mod_S": w})
    return {"t": t, "fault_sets": dict(sets), "undetected": dict(undetected),
            "undetected_weight_hist": dict(sorted(hist.items())), "violations": bad, "FT": not bad}


def fe_signature_t(gadget, t=1, code=None):
    """{(data coset representative, syndrome flip): smallest number of faults producing it}
    over the unflagged sets of at most t faults (the C2 effect set of plan Sec. 0.3)."""
    code = code_of(gadget, code)
    flags, syn = mask(gadget["flags"]), mask(gadget["syndrome"])
    sig = {}
    for s, labels, x, z in faults_upto(gadget["ops"], t):
        if x & flags:
            continue
        key = (code.canonical_mod(x & code.data_mask, z & code.data_mask), popcount(x & syn) % 2)
        sig.setdefault(key, s)
    return sig


def fe_compare_t(a, b, t=1, code=None):
    """Order-t necessary condition for fault equivalence: both gadgets must give the same
    unflagged effects with the same minimal fault count (the no-fault effect is ignored).
    Returns ({effect: count} only or differently in a, same for b)."""
    code = code_of(a, code)
    ident = (code.canonical_mod(0, 0), 0)
    sa, sb = fe_signature_t(a, t, code), fe_signature_t(b, t, code)
    only_a = {k: v for k, v in sa.items() if k != ident and sb.get(k) != v}
    only_b = {k: v for k, v in sb.items() if k != ident and sa.get(k) != v}
    return only_a, only_b


# --- vectorised versions for the search checkers (t <= 2) ----------------------------------

def _parity(v):
    for sh in (32, 16, 8, 4, 2, 1):
        v = v ^ (v >> sh)
    return v & 1


def effect_arrays(ops, t, flags):
    """Unflagged sets of 1..t faults (t <= 2) as numpy arrays (s, x, z) of end-of-circuit Paulis;
    the same sets as faults_upto, without labels."""
    if t > 2:
        raise NotImplementedError("the vectorised path covers t <= 2; use faults_upto")
    loc, xs, zs = [], [], []
    for l, _, start, x, z in located_faults(ops):
        x, z = propagate(x, z, ops[start:])
        loc.append(l)
        xs.append(x)
        zs.append(z)
    # int64 masks hold qubits 0-62; larger circuits (e.g. a 47-ancilla [[17,1,5]] round) use Python ints
    wide = max((q for op in ops for q in op[1:]), default=0) >= 62
    dt = object if wide else np.int64
    loc, X, Z = np.array(loc), np.array(xs, dtype=dt), np.array(zs, dtype=dt)
    S, XS, ZS = [np.ones(len(X), np.int64)], [X], [Z]
    if t == 2:
        i, j = np.triu_indices(len(X), 1)
        keep = loc[i] != loc[j]
        i, j = i[keep], j[keep]
        S.append(np.full(len(i), 2, np.int64))
        XS.append(X[i] ^ X[j])
        ZS.append(Z[i] ^ Z[j])
    s, x, z = (np.concatenate(a) for a in (S, XS, ZS))
    keep = (x & mask(flags)) == 0  # a Z measurement of a flag flips on X or Y
    return s[keep], x[keep], z[keep]


def effect_signature(gadget, t, readouts, code=None):
    """Vectorised C2 effect set and FT_t test.

    readouts: lists of measured qubits; readout r flips when the fault set leaves an odd number
    of X/Y on its qubits (a one-bit syndrome gadget has readouts = [gadget["syndrome"]]).
    Returns ({((x', z') data coset representative, flips): smallest fault count}, number of
    unflagged sets whose residual weight mod S exceeds their fault count).
    """
    code = code_of(gadget, code)
    s, x, z = effect_arrays(gadget["ops"], t, gadget["flags"])
    n, dm = code.n_data, code.data_mask
    pair = ((x & dm) << n) | (z & dm)
    uniq, inv = np.unique(pair, return_inverse=True)
    reps = np.array([code._rep(int(p) >> n, int(p) & dm) for p in uniq], dtype=np.int64).reshape(-1, 3)
    w, rep = reps[inv, 0], (reps[inv, 1] << n) | reps[inv, 2]
    key = rep
    for r in readouts:
        key = (key << 1) | _parity(x & mask(r))
    order = np.argsort(s, kind="stable")
    keys, first = np.unique(key[order], return_index=True)
    mins = s[order][first]
    R = len(readouts)
    sig = {}
    for k, m in zip(keys.tolist(), mins.tolist()):
        flips = tuple((k >> (R - 1 - b)) & 1 for b in range(R))
        r = k >> R
        sig[((r >> n, r & dm), flips)] = m
    return sig, int(np.count_nonzero(w > s))


def effect(ops, start, x, z, syndrome_qubits, flag_qubits, n_data=N_DATA):
    """Propagate a fault to the end; return (data x, data z, syndrome flip, flag flips)."""
    x, z = propagate(x, z, ops[start:])
    data = (1 << n_data) - 1
    flips = lambda qs: [(x >> q) & 1 for q in qs]  # a Z measurement flips on X or Y
    syndrome_flip = sum(flips(syndrome_qubits)) % 2
    return x & data, z & data, syndrome_flip, tuple(flips(flag_qubits))


def analyse(gadget, code=None):
    """Run every single fault through the gadget and classify it."""
    code = code_of(gadget, code)
    rows = []
    for label, start, x, z in single_faults(gadget["ops"]):
        dx, dz, sflip, fflips = effect(gadget["ops"], start, x, z,
                                       gadget["syndrome"], gadget["flags"], code.n_data)
        rows.append({
            "fault": label,
            "data": pauli_str(dx, dz, range(code.n_data)),
            "dx": dx, "dz": dz,
            "weight_mod_S": code.weight_mod(dx, dz),
            "syndrome_flip": sflip,
            "detected": any(fflips),
        })
    return rows


def ft1_report(gadget, code=None):
    """ft_report at t = 1, in the field names the t = 1 scripts use."""
    r = ft_report(gadget, 1, code)
    return {
        "faults": r["fault_sets"].get(1, 0),
        "detected": r["fault_sets"].get(1, 0) - r["undetected"].get(1, 0),
        "undetected": r["undetected"].get(1, 0),
        "undetected_weight_hist": {w: n for (_, w), n in r["undetected_weight_hist"].items()},
        "violations": r["violations"],
        "FT1": r["FT"],
    }


def fe_signature(gadget, code=None):
    """fe_signature_t at t = 1: {(data coset representative, syndrome flip): 1}."""
    return fe_signature_t(gadget, 1, code)


def fe_compare(a, b, code=None):
    """Single-fault check of 2-fault-equivalence between gadgets a and b (fe_compare_t, t = 1).

    Every undetected single-fault effect of one gadget must be produced by a fault of weight
    at most 1 in the other (weight 0 = no fault = identity effect with no syndrome flip).
    """
    return fe_compare_t(a, b, 1, code)


def resources(ops):
    """Ancilla count (qubits prepared by R), CNOT count and CNOT depth (greedy layering)."""
    ancillas = sorted({op[1] for op in ops if op[0] == "R"})
    layer = {}
    depth = 0
    for op in ops:
        if op[0] == "CX":
            d = max(layer.get(op[1], 0), layer.get(op[2], 0)) + 1
            layer[op[1]] = layer[op[2]] = d
            depth = max(depth, d)
    return {"ancillas": len(ancillas), "cnots": sum(op[0] == "CX" for op in ops),
            "cnot_depth": depth}


def to_stim(ops, n_data=N_DATA, data_init="plus", z_on=(), p=None, x_on=()):
    """Build a stim.Circuit: prepare data, run the gadget, then measure the data.

    data_init: "plus" (|+>^n, measure data in X at the end) or "zero" (|0>^n, Z at the end).
    z_on: data qubits that get a Z before the gadget (flips X-stabilizer eigenvalues).
    x_on: data qubits that get an X before the gadget (flips Z-stabilizer eigenvalues).
    p: if given, add noise with the ratios of Jeon et al. Sec. V A: two-qubit gates p,
       single-qubit gates 3p/100, preparation and measurement p. Idle noise is left out.
       The data preparation and the final data measurement stay noiseless.
    Measurement record order: gadget measurements in op order, then data 0..n-1.
    """
    c = stim.Circuit()
    data = list(range(n_data))
    c.append("RX" if data_init == "plus" else "R", data)
    for q in z_on:
        c.append("Z", [q])
    for q in x_on:
        c.append("X", [q])
    for op in ops:
        if op[0] == "R":
            c.append("R", [op[1]])
            if p:
                c.append("X_ERROR", [op[1]], p)
        elif op[0] == "H":
            c.append("H", [op[1]])
            if p:
                c.append("DEPOLARIZE1", [op[1]], 3 * p / 100)
        elif op[0] == "CX":
            c.append("CX", [op[1], op[2]])
            if p:
                c.append("DEPOLARIZE2", [op[1], op[2]], p)
        elif op[0] == "M":
            if p:
                c.append("X_ERROR", [op[1]], p)
            c.append("M", [op[1]])
    c.append("MX" if data_init == "plus" else "M", data)
    return c


def measured_qubits(ops):
    return [op[1] for op in ops if op[0] == "M"]


# ---------------------------------------------------------------------------
# Gadgets measuring the X-type generator X0 X1 X2 X3 of the Steane code.
# Read off Jeon et al., Fig. 1 (page 6). Check against the PDF figure.
# ---------------------------------------------------------------------------

def _x_measure_tail(cat):
    return [("H", q) for q in cat]


BASELINE = {
    "name": "baseline (Fig. 1a)",
    "ops": [("R", q) for q in (7, 8, 9, 10, 11)] + [
        ("H", 7),
        ("CX", 7, 8), ("CX", 7, 9), ("CX", 8, 10),   # cat state on q7..q10
        ("CX", 7, 11),                                 # verification: parity Z7 Z8 -> flag q11
        ("CX", 7, 0), ("CX", 8, 11),
        ("CX", 8, 1), ("CX", 9, 2), ("CX", 10, 3),    # each cat qubit touches one data qubit
    ] + _x_measure_tail((7, 8, 9, 10)) + [("M", q) for q in (7, 8, 9, 10, 11)],
    "syndrome": [7, 8, 9, 10],
    "flags": [11],
}

NO_FLAG = {
    "name": "baseline without verification",
    "ops": [op for op in BASELINE["ops"] if 11 not in op[1:]],
    "syndrome": [7, 8, 9, 10],
    "flags": [],
}

REDUCED = {
    "name": "reduced (Fig. 1b)",
    "ops": [("R", q) for q in (7, 8, 9)] + [
        ("H", 7),
        ("CX", 7, 2),
        ("CX", 7, 9),          # first flag CNOT
        ("CX", 7, 0),
        ("CX", 7, 8),
        ("CX", 7, 9),          # second flag CNOT
        ("CX", 8, 3),
        ("H", 8),
        ("CX", 7, 1),
        ("H", 7),
    ] + [("M", q) for q in (7, 8, 9)],
    "syndrome": [7, 8],
    "flags": [9],
}


def dual(gadget):
    """Z-type twin of an X-type gadget: conjugate the whole circuit by H (CX reverses,
    |0> <-> |+>, Z <-> X measurement), then cancel the H H pairs this leaves on a wire."""
    ops = []
    for op in gadget["ops"]:
        if op[0] == "R":
            ops += [op, ("H", op[1])]
        elif op[0] == "M":
            ops += [("H", op[1]), op]
        elif op[0] == "CX":
            ops.append(("CX", op[2], op[1]))
        else:
            ops.append(op)
    out = []
    for op in ops:
        if op[0] == "H":
            k = max((i for i, o in enumerate(out) if op[1] in o[1:]), default=None)
            if k is not None and out[k] == op:
                del out[k]
                continue
        out.append(op)
    kind = "Z" if gadget.get("kind", "X") == "X" else "X"
    return {**gadget, "ops": out, "kind": kind, "name": f"dual of {gadget.get('name', '')}"}


def ideal_syndrome_check(gadget, shots=200, code=None):
    """Noiseless: does the syndrome parity equal the eigenvalue of the measured stabilizer?

    gadget["kind"] is "X" (default: data |+>, stabilizer X on gadget["target"], default the
    Steane generator X0X1X2X3) or "Z" (data |0>, Z stabilizer). The -1 case flips the first
    target qubit with Z (or X).
    Returns {case: (set of syndrome parities seen, set of flag values seen)}.
    """
    code = code_of(gadget, code)
    kind, target = gadget.get("kind", "X"), gadget.get("target", STEANE_SUPPORTS[0])
    meas = measured_qubits(gadget["ops"])
    out = {}
    flip = "z_on" if kind == "X" else "x_on"
    for label, kw in ((f"{kind}{target} = +1", {}), (f"{kind}{target} = -1", {flip: (target[0],)})):
        c = to_stim(gadget["ops"], code.n_data, "plus" if kind == "X" else "zero", **kw)
        s = c.compile_sampler().sample(shots).astype(int)
        syn = s[:, [meas.index(q) for q in gadget["syndrome"]]].sum(axis=1) % 2
        flg = (s[:, [meas.index(q) for q in gadget["flags"]]].sum(axis=1) % 2
               if gadget["flags"] else np.zeros(shots, int))
        out[label] = (sorted(set(syn.tolist())), sorted(set(flg.tolist())))
    return out
