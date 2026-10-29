"""[[7,1,3]] memory experiment with Shor-style gadgets (Jeon et al., Sec. V B).

A cycle is three rounds of syndrome measurement followed by one decoding step (Shor's rule is a
property of the decoder, not of the circuit).  Each round measures the three X-type and three
Z-type generators with a cat-state gadget (ftlib.BASELINE or ftlib.REDUCED); every gadget carries
its own verification flag.

Noise (plan Sec. 0.4 / Jeon Sec. V A): two-qubit depolarizing p after every CX, single-qubit
depolarizing 3p/100 after every H, X flip with probability p after every preparation and before
every measurement, single-qubit depolarizing p/10 on every live qubit idle in a layer.  Layers
come from a greedy ASAP schedule inside each round; rounds run one after another.

The circuit is emitted with DETECTOR annotations (one per gadget parity, one per flag) and one
logical observable, so that any stim-compatible decoder can be attached.  Two decoders are
provided in ftrs/decoders.py: a flag-conditioned lookup table (Jeon Sec. V B 3) and BP+OSD on
the detector error model.
"""
import numpy as np
import stim

from .ftlib import BASELINE, REDUCED, STEANE_SUPPORTS

N_DATA = 7
NGEN = 6


def template_spec(template_ops):
    """Ancillas of a template gadget: the cat qubits and the verification flag (last ancilla)."""
    anc = sorted({op[1] for op in template_ops if op[0] == "R"})
    return {"anc": anc, "cat": anc[:-1], "flag": anc[-1]}


def relabel(template_ops, support, anc0, spec):
    m = {0: support[0], 1: support[1], 2: support[2], 3: support[3]}
    for i, a in enumerate(spec["anc"]):
        m[a] = anc0 + i
    return [tuple([op[0]] + [m[q] for q in op[1:]]) for op in template_ops]


def round_gadgets(anc_start, spec):
    """The six gadgets of one round: three Z-type (conjugated by H on data), three X-type."""
    gs, anc, n, ncat = [], anc_start, len(spec["anc"]), len(spec["cat"])
    for kind in ("Z", "X"):
        for sup in STEANE_SUPPORTS:
            gs.append({"type": kind, "support": sup, "anc": list(range(anc, anc + ncat)),
                       "flag": anc + n - 1})
            anc += n
    return gs


def round_ops(template_ops, anc_start):
    """Op list of one round (no noise): H-wrap the Z-type sub-round, then the X-type one."""
    spec = template_spec(template_ops)
    gs = round_gadgets(anc_start, spec)
    ops = [("H", q) for q in range(N_DATA)]
    for g in gs[:3]:
        ops += relabel(template_ops, g["support"], g["anc"][0], spec)
    ops += [("H", q) for q in range(N_DATA)]
    for g in gs[3:]:
        ops += relabel(template_ops, g["support"], g["anc"][0], spec)
    return ops, gs


def schedule(ops):
    """Greedy ASAP layering; returns (n_layers, layer_of_op)."""
    free, layers = {}, []
    for op in ops:
        qs = op[1:]
        m = max((free.get(q, 0) for q in qs), default=0)
        layers.append(m)
        for q in qs:
            free[q] = m + 1
    return (max(free.values()) if free else 0), layers


def encode_circuit(memory="Z"):
    """Ideal preparation of |0>_L (or |+>_L) as a Clifford circuit on qubits 0..6."""
    def P(sup, k):
        return stim.PauliString({q: k for q in sup})
    stabs = [P(s, "X") for s in STEANE_SUPPORTS] + [P(s, "Z") for s in STEANE_SUPPORTS]
    stabs.append(P(range(N_DATA), "Z"))          # Z_L: |0>_L is a +1 eigenstate
    c = stim.Tableau.from_stabilizers(stabs).to_circuit("elimination")
    if memory == "X":
        c.append("H", list(range(N_DATA)))        # logical H
    return c


def schedule_serial(ops, block_len):
    """Layers when each block (H wrap or one gadget) starts after the previous block ends."""
    lay, base, i = [], 0, 0
    for n in block_len:
        nb, lb = schedule(ops[i:i + n])
        lay += [base + m for m in lb]
        base += nb
        i += n
    return base, lay


def build_circuit(template_ops, n_cycles=10, rounds_per_cycle=3, memory="Z", p=1e-3,
                  noisy=True, detectors=True, perfect_final_round=True,
                  schedule_mode="asap", idle_factor=0.1):
    """Build the memory circuit.  Returns (stim.Circuit, info dict).

    schedule_mode "asap" layers a whole round greedily (every cat state is prepared in the first
    layers and idles until its gadget reaches the data); "serial" runs the H wraps and the six
    gadgets one after another.  Idle qubits get DEPOLARIZE1(idle_factor * p) per layer.
    """
    spec = template_spec(template_ops)
    nanc = len(spec["anc"])
    n_rounds = n_cycles * rounds_per_cycle
    c = stim.Circuit()
    c += encode_circuit(memory)

    live = set(range(N_DATA))
    nrec = 0
    rec_anc, rec_flag = {}, {}
    n_det = 0

    for r in range(n_rounds):
        ops, gs = round_ops(template_ops, N_DATA + NGEN * nanc * r)
        if schedule_mode == "serial":
            L = len(template_ops)
            n_layers, lay = schedule_serial(ops, [N_DATA, L, L, L, N_DATA, L, L, L])
        else:
            n_layers, lay = schedule(ops)
        by_layer = [[] for _ in range(n_layers)]
        for op, m in zip(ops, lay):
            by_layer[m].append(op)
        round_rec = {}
        noisy_here = noisy and not (perfect_final_round and r == n_rounds - 1)
        for m in range(n_layers):
            touched = set()
            for op in by_layer[m]:
                name, qs = op[0], op[1:]
                touched |= set(qs)
                if name == "R":
                    c.append("R", [qs[0]])
                    live.add(qs[0])
                    if noisy_here:
                        c.append("X_ERROR", [qs[0]], p)
                elif name == "H":
                    c.append("H", [qs[0]])
                    if noisy_here:
                        c.append("DEPOLARIZE1", [qs[0]], 3 * p / 100)
                elif name == "CX":
                    c.append("CX", [qs[0], qs[1]])
                    if noisy_here:
                        c.append("DEPOLARIZE2", [qs[0], qs[1]], p)
                elif name == "M":
                    if noisy_here:
                        c.append("X_ERROR", [qs[0]], p)
                    c.append("M", [qs[0]])
                    live.discard(qs[0])
                    nrec += 1
                    round_rec[qs[0]] = nrec - 1
            idle = sorted(live - touched)
            if noisy_here and idle and idle_factor:
                c.append("DEPOLARIZE1", idle, idle_factor * p)
        # detectors: differences between consecutive rounds (the standard memory-experiment
        # detector), so every error mechanism touches only a local window of detectors.
        for gi, g in enumerate(gs):
            for k, q in enumerate(g["anc"]):
                rec_anc[(r, gi, k)] = round_rec[q]
            rec_flag[(r, gi)] = round_rec[g["flag"]]
            if detectors:
                offs = [stim.target_rec(round_rec[q] - nrec) for q in g["anc"]]
                if r > 0:
                    offs += [stim.target_rec(prev_rec[(gi, k)] - nrec)
                             for k, q in enumerate(g["anc"])]
                c.append("DETECTOR", offs)
                n_det += 1
                offs = [stim.target_rec(round_rec[g["flag"]] - nrec)]
                if r > 0:
                    offs.append(stim.target_rec(prev_rec[(gi, len(g["anc"]))] - nrec))
                c.append("DETECTOR", offs)
                n_det += 1
        prev_rec = {}
        for gi, g in enumerate(gs):
            for k, q in enumerate(g["anc"]):
                prev_rec[(gi, k)] = round_rec[q]
            prev_rec[(gi, len(g["anc"]))] = round_rec[g["flag"]]

    rec_data = nrec
    data = list(range(N_DATA))
    if memory == "Z":
        c.append("M", data)
    else:
        c.append("MX", data)
    if detectors:
        c.append("OBSERVABLE_INCLUDE", [stim.target_rec(i - N_DATA)
                                        for i in range(N_DATA)], 0)
    info = {"n_rounds": n_rounds, "rec_anc": rec_anc, "rec_flag": rec_flag,
            "rec_data": rec_data, "n_detectors": n_det, "spec": spec}
    return c, info


# ---------------------------------------------------------------- lookup-table decoder

def single_fault_table(template_ops):
    """Lookup table from one round's single faults: 6-bit syndrome -> (dx, dz).

    Faults that raise a flag are left out (the protocol discards a flagged round).  Ambiguous
    syndromes keep the correction of minimal data weight.
    """
    from .ftlib import STEANE_GROUP, canonical_mod, effect, single_faults, weight_mod
    from .ftlib import propagate
    ops, gs = round_ops(template_ops, N_DATA)
    syndrome_q = [q for g in gs for q in g["anc"]]
    flag_q = [g["flag"] for g in gs]
    table, conflicts = {}, 0
    for label, start, x, z in single_faults(ops):
        dx, dz, _, fflips = effect(ops, start, x, z, syndrome_q, flag_q, N_DATA)
        if any(fflips):
            continue
        xx, zz = propagate(x, z, ops[start:])
        key = tuple(int(sum((xx >> q) & 1 for q in g["anc"]) % 2) for g in gs)
        if key in table:
            if canonical_mod(*table[key], STEANE_GROUP) != canonical_mod(dx, dz, STEANE_GROUP):
                conflicts += 1
            if weight_mod(dx, dz, STEANE_GROUP) < weight_mod(*table[key], STEANE_GROUP):
                table[key] = (dx, dz)
        else:
            table[key] = (dx, dz)
    return table, {"entries": len(table), "conflicts": conflicts}


def syndromes_and_flags(samples, info):
    """(syn[n, n_rounds, 6], flag[n, n_rounds, 6]) from a raw measurement sample."""
    rec, n_rounds = info["rec_anc"], info["n_rounds"]
    n = samples.shape[0]
    n_cat = len(info["spec"]["cat"])
    syn = np.zeros((n, n_rounds, NGEN), dtype=np.uint8)
    flg = np.zeros((n, n_rounds, NGEN), dtype=np.uint8)
    for r in range(n_rounds):
        for gi in range(NGEN):
            s = np.zeros(n, dtype=np.uint8)
            for k in range(n_cat):
                s ^= samples[:, rec[(r, gi, k)]].astype(np.uint8)
            syn[:, r, gi] = s
            flg[:, r, gi] = samples[:, info["rec_flag"][(r, gi)]].astype(np.uint8)
    return syn, flg


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (0.0, 0.0)
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))
