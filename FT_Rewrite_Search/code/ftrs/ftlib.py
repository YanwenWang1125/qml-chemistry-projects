"""Small fault-propagation toolkit shared by the FT_Rewrite_Search lessons.

Qubits 0-6 are the data qubits of the [[7,1,3]] Steane code. Ancillas start at 7,
matching the q7, q8, ... labels in Jeon et al. (arXiv:2610.09749, Fig. 1 and App. B).

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


def single_faults(ops):
    """Every single fault of the circuit-level model used in these lessons.

    after R q: X on q (preparation error)
    after H q: X, Y or Z on q
    after CX c t: any of the 15 non-identity two-qubit Paulis on (c, t)
    before M q: X on q (measurement error)
    Yields (label, index of the first op the fault passes through, x, z).
    """
    for i, op in enumerate(ops):
        if op[0] == "R":
            yield f"prep X on q{op[1]}", i + 1, 1 << op[1], 0
        elif op[0] == "H":
            for p in "XYZ":
                bx, bz = ONE_QUBIT_BITS[p]
                yield f"{p} on q{op[1]} after H", i + 1, bx << op[1], bz << op[1]
        elif op[0] == "CX":
            c, t = op[1], op[2]
            for pc, pt in TWO_QUBIT_PAULIS:
                xc, zc = ONE_QUBIT_BITS[pc]
                xt, zt = ONE_QUBIT_BITS[pt]
                yield (f"{pc}{pt} on (q{c},q{t}) after CX{c}->{t}", i + 1,
                       (xc << c) | (xt << t), (zc << c) | (zt << t))
        elif op[0] == "M":
            yield f"measurement flip on q{op[1]}", i, 1 << op[1], 0


def effect(ops, start, x, z, syndrome_qubits, flag_qubits, n_data=N_DATA):
    """Propagate a fault to the end; return (data x, data z, syndrome flip, flag flips)."""
    x, z = propagate(x, z, ops[start:])
    data = (1 << n_data) - 1
    flips = lambda qs: [(x >> q) & 1 for q in qs]  # a Z measurement flips on X or Y
    syndrome_flip = sum(flips(syndrome_qubits)) % 2
    return x & data, z & data, syndrome_flip, tuple(flips(flag_qubits))


def analyse(gadget, group=STEANE_GROUP, n_data=N_DATA):
    """Run every single fault through the gadget and classify it."""
    rows = []
    for label, start, x, z in single_faults(gadget["ops"]):
        dx, dz, sflip, fflips = effect(gadget["ops"], start, x, z,
                                       gadget["syndrome"], gadget["flags"], n_data)
        rows.append({
            "fault": label,
            "data": pauli_str(dx, dz, range(n_data)),
            "dx": dx, "dz": dz,
            "weight_mod_S": weight_mod(dx, dz, group),
            "syndrome_flip": sflip,
            "detected": any(fflips),
        })
    return rows


def ft1_report(gadget, group=STEANE_GROUP, n_data=N_DATA):
    rows = analyse(gadget, group, n_data)
    undetected = [r for r in rows if not r["detected"]]
    bad = [r for r in undetected if r["weight_mod_S"] > 1]
    hist = {}
    for r in undetected:
        hist[r["weight_mod_S"]] = hist.get(r["weight_mod_S"], 0) + 1
    return {
        "faults": len(rows),
        "detected": len(rows) - len(undetected),
        "undetected": len(undetected),
        "undetected_weight_hist": dict(sorted(hist.items())),
        "violations": bad,
        "FT1": not bad,
    }


def fe_signature(gadget, group=STEANE_GROUP, n_data=N_DATA):
    """Effects of the undetected single faults: {(data coset representative, syndrome flip)}.

    This is a simplified, single-fault version of fault equivalence: it compares only what
    an undetected fault does to the data (modulo the stabilizer group) and to the syndrome bit.
    """
    sig = {}
    for label, start, x, z in single_faults(gadget["ops"]):
        dx, dz, sflip, fflips = effect(gadget["ops"], start, x, z,
                                       gadget["syndrome"], gadget["flags"], n_data)
        if any(fflips):
            continue
        key = (canonical_mod(dx, dz, group), sflip)
        sig.setdefault(key, label)
    return sig


def fe_compare(a, b, group=STEANE_GROUP, n_data=N_DATA):
    """Single-fault check of 2-fault-equivalence between gadgets a and b.

    Every undetected single-fault effect of one gadget must be produced by a fault of weight
    at most 1 in the other (weight 0 = no fault = identity effect with no syndrome flip).
    Returns the effects found only in a and only in b, with one example fault each.
    """
    ident = (canonical_mod(0, 0, group), 0)
    sa, sb = fe_signature(a, group, n_data), fe_signature(b, group, n_data)
    only_a = {k: v for k, v in sa.items() if k not in sb and k != ident}
    only_b = {k: v for k, v in sb.items() if k not in sa and k != ident}
    return only_a, only_b


def resources(ops):
    """Ancilla count (qubits >= 7), CNOT count and CNOT depth (greedy layering)."""
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


def to_stim(ops, n_data=N_DATA, data_init="plus", z_on=(), p=None):
    """Build a stim.Circuit: prepare data, run the gadget, then measure the data.

    data_init: "plus" (|+>^n, measure data in X at the end) or "zero" (|0>^n, Z at the end).
    z_on: data qubits that get a Z before the gadget (flips X-stabilizer eigenvalues).
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


def ideal_syndrome_check(gadget, shots=200):
    """Noiseless: does the syndrome parity equal the eigenvalue of X0X1X2X3?

    Returns {case: (set of syndrome parities seen, set of flag values seen)}.
    """
    meas = measured_qubits(gadget["ops"])
    out = {}
    for label, z_on in (("data |+>^7 (X0X1X2X3 = +1)", ()),
                        ("data Z0|+>^7 (X0X1X2X3 = -1)", (0,))):
        c = to_stim(gadget["ops"], data_init="plus", z_on=z_on)
        s = c.compile_sampler().sample(shots).astype(int)
        syn = s[:, [meas.index(q) for q in gadget["syndrome"]]].sum(axis=1) % 2
        flg = (s[:, [meas.index(q) for q in gadget["flags"]]].sum(axis=1) % 2
               if gadget["flags"] else np.zeros(shots, int))
        out[label] = (sorted(set(syn.tolist())), sorted(set(flg.tolist())))
    return out
