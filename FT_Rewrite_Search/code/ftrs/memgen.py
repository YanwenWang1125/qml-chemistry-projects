"""Memory experiment for any CSS code and any multi-bit syndrome gadgets (E2.0 step 7, setting B).

A round runs each half-round gadget in turn (e.g. the Z-stabilizer gadget, then its Hadamard
dual for the X stabilizers), each ASAP-scheduled on its own. Ancillas are reset and reused
every round. Readouts are the measurement parities found by multi.derive. Detectors: each
readout XOR the same readout one round earlier (in round 0 only readouts of the memory basis,
which are deterministic), and each flag measurement on its own (flags are deterministic 0).
After the rounds the data are measured in the memory basis; final detectors compare the data
parities of the memory-basis stabilizers with the last round, and the observable is the
logical operator of that basis.

Noise (plan Sec. 0.4, Jeon's setting): DEPOLARIZE2(p) after CX, DEPOLARIZE1(3p/100) after H,
X_ERROR(p) after R and before M, DEPOLARIZE1(idle_factor * p) on every live idle qubit per
ASAP layer. The encoding and the final data measurement are noiseless; with
perfect_final_round the last round is noiseless too.
"""
import stim

from .multi import derive


def encode(code, memory="Z"):
    n = code.n_data
    ps = lambda sup, k: stim.PauliString("".join(k if q in sup else "_" for q in range(n)))
    stabs = [ps(s, "X") for s in code.x_supports] + [ps(s, "Z") for s in code.z_supports]
    stabs.append(ps((code.logical_z if memory == "Z" else code.logical_x)[0], memory))
    return stim.Tableau.from_stabilizers(stabs).to_circuit("elimination")


def _layers(ops):
    free, lay = {}, []
    for op in ops:
        m = max((free.get(q, 0) for q in op[1:]), default=0)
        lay.append(m)
        for q in op[1:]:
            free[q] = m + 1
    n = max(free.values()) if free else 0
    out = [[] for _ in range(n)]
    for op, m in zip(ops, lay):
        out[m].append(op)
    return out


def build(code, halves, rounds, memory="Z", p=1e-3, idle_factor=0.1, perfect_final_round=True,
          h_factor=0.03, prep_factor=1.0, meas_factor=1.0):
    """halves: list of multi-bit gadgets (kind, targets, syndrome, flags, ops) run each round.
    Error rates relative to p: H h_factor (DEPOLARIZE1), preparation prep_factor and measurement
    meas_factor (X_ERROR), idle idle_factor; the defaults are plan Sec. 0.4. Liou 2024 Sec. II B
    (gamma = beta = 1) is h_factor 1, prep_factor 2/3, meas_factor 2/3, idle_factor 1.
    Returns (stim.Circuit, info)."""
    n = code.n_data
    readouts = []
    for g in halves:
        par, ok = derive(g, shots=2000)
        assert ok, f"{g.get('name')} fails its noiseless check"
        readouts.append(par)
    c = stim.Circuit()
    c += encode(code, memory)
    live, nrec, prev = set(range(n)), 0, None
    n_det = 0
    for r in range(rounds):
        noisy = p > 0 and not (perfect_final_round and r == rounds - 1)
        recs = []  # one map per half: the halves reuse ancilla labels
        for g in halves:
            rec = {}
            recs.append(rec)
            for layer in _layers(g["ops"]):
                touched = set()
                for op in layer:
                    name, qs = op[0], list(op[1:])
                    touched |= set(qs)
                    if name == "R":
                        c.append("R", qs)
                        live.add(qs[0])
                        if noisy:
                            c.append("X_ERROR", qs, prep_factor * p)
                    elif name == "H":
                        c.append("H", qs)
                        if noisy:
                            c.append("DEPOLARIZE1", qs, h_factor * p)
                    elif name == "CX":
                        c.append("CX", qs)
                        if noisy:
                            c.append("DEPOLARIZE2", qs, p)
                    elif name == "M":
                        if noisy:
                            c.append("X_ERROR", qs, meas_factor * p)
                        c.append("M", qs)
                        live.discard(qs[0])
                        rec[qs[0]] = nrec
                        nrec += 1
                idle = sorted(live - touched)
                if noisy and idle and idle_factor:
                    c.append("DEPOLARIZE1", idle, idle_factor * p)
        cur = []
        for g, par, rec in zip(halves, readouts, recs):
            for k, qs in enumerate(par):
                cur.append((g["kind"], [rec[q] for q in qs]))
            for f in g["flags"]:
                c.append("DETECTOR", [stim.target_rec(rec[f] - nrec)])
                n_det += 1
        for i, (kind, recs) in enumerate(cur):
            targets = [stim.target_rec(x - nrec) for x in recs]
            if prev is not None:
                targets += [stim.target_rec(x - nrec) for x in prev[i][1]]
            elif kind != memory:
                continue  # first round: the other basis is random
            c.append("DETECTOR", targets)
            n_det += 1
        prev = cur
    data = list(range(n))
    c.append("M" if memory == "Z" else "MX", data)
    nrec_final = nrec + n
    sups = code.z_supports if memory == "Z" else code.x_supports
    i0 = 0
    for g, par in zip(halves, readouts):
        if g["kind"] == memory:
            for k, sup in enumerate(g["targets"]):
                last = prev[i0 + k][1]
                c.append("DETECTOR", [stim.target_rec(nrec + q - nrec_final) for q in sup] +
                         [stim.target_rec(x - nrec_final) for x in last])
                n_det += 1
        i0 += len(par)
    logical = (code.logical_z if memory == "Z" else code.logical_x)[0]
    c.append("OBSERVABLE_INCLUDE", [stim.target_rec(nrec + q - nrec_final) for q in logical], 0)
    return c, {"detectors": n_det, "readouts": readouts}
