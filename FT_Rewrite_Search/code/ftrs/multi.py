"""Checker for gadgets whose syndrome has several bits (e.g. Steane-style extraction).

After a rewrite the meaning of each measurement can change, so the syndrome bits are re-derived
for every circuit: sample the noiseless circuit with random data errors, then solve over GF(2)
for the combinations of measurement bits that reproduce each target stabilizer bit.
Gadget dict keys: ops, syndrome (ancillas allowed in rules), flags, and
    kind: "Z" (data in |0>, random X errors, targets = Z stabilizers) or "X" (the reverse),
    targets: list of data supports.
"""
import numpy as np
import stim

from .ftlib import code_of, effect, effect_signature, propagate, single_faults


def _circuit(g, shots_err=True):
    c = stim.Circuit()
    d = list(range(code_of(g).n_data))
    if g["kind"] == "Z":
        c.append("R", d)
        if shots_err:
            c.append("X_ERROR", d, 0.5)
    else:
        c.append("RX", d)
        if shots_err:
            c.append("Z_ERROR", d, 0.5)
    for op in g["ops"]:
        c.append({"R": "R", "H": "H", "CX": "CX", "M": "M"}[op[0]], list(op[1:]))
    c.append("M" if g["kind"] == "Z" else "MX", d)
    return c


def _solve(A, b):
    """Solve A x = b over GF(2); return x or None."""
    A = A.copy() % 2
    b = b.copy() % 2
    n_rows, n_cols = A.shape
    piv_cols, r = [], 0
    for c in range(n_cols):
        rows = np.nonzero(A[r:, c])[0]
        if len(rows) == 0:
            continue
        p = r + rows[0]
        A[[r, p]] = A[[p, r]]
        b[[r, p]] = b[[p, r]]
        for k in range(n_rows):
            if k != r and A[k, c]:
                A[k] ^= A[r]
                b[k] ^= b[r]
        piv_cols.append(c)
        r += 1
        if r == n_rows:
            break
    if b[r:].any():
        return None
    x = np.zeros(n_cols, dtype=np.uint8)
    for i, c in enumerate(piv_cols):
        x[c] = b[i]
    return x


def derive(g, shots=600, seed=5):
    """Return (parities, ok): parities[i] = list of measured qubits whose XOR gives target i."""
    meas = [op[1] for op in g["ops"] if op[0] == "M"]
    s = _circuit(g).compile_sampler(seed=seed).sample(shots).astype(np.uint8)
    m, data = s[:, :len(meas)], s[:, len(meas):]
    A = np.hstack([m, np.ones((shots, 1), np.uint8)])
    parities = []
    for sup in g["targets"]:
        x = _solve(A, data[:, list(sup)].sum(axis=1) % 2)
        if x is None or x[-1]:  # unsolvable, or needs a constant offset
            return None, False
        parities.append([meas[k] for k in np.nonzero(x[:-1])[0]])
    for f in g["flags"]:  # flags must be deterministically 0 without noise
        if m[:, meas.index(f)].any():
            return None, False
    return parities, True


def effects(g, parities, code=None):
    """{(data coset, syndrome flip vector)} for undetected single faults, plus FT_1 violations."""
    code = code_of(g, code)
    sig, bad = {}, []
    for label, start, x, z in single_faults(g["ops"]):
        dx, dz, _, fflips = effect(g["ops"], start, x, z, [], g["flags"], code.n_data)
        if any(fflips):
            continue
        xx, _ = propagate(x, z, g["ops"][start:])
        flips = tuple(sum((xx >> q) & 1 for q in par) % 2 for par in parities)
        sig.setdefault((code.canonical_mod(dx, dz), flips), label)
        if code.weight_mod(dx, dz) > 1:
            bad.append(label)
    return sig, bad


def check_factory(root, t=1):
    rp, ok = derive(root)
    assert ok, "root gadget fails its own ideal check"
    rsig, rbad = effect_signature(root, t, rp)
    assert not rbad, f"root gadget is not FT_{t}"
    # Compare only the error type this gadget handles (X data errors for kind "Z", Z for "X")
    # plus syndrome flips; the other type only has to keep weight <= s (the FT_t test below).
    # Jeon's App. C steps differ from the baseline only in residual weight-1 Z errors (E1.3-fix).
    part = 0 if root["kind"] == "Z" else 1

    def proj(sig):  # keep the smallest fault count per projected effect
        out = {}
        for k, s in sig.items():
            pk = (k[0][part], k[1])
            out[pk] = min(out.get(pk, s), s)
        return out

    rsig = proj(rsig)
    ident = (code_of(root).canonical_mod(0, 0)[part], tuple(0 for _ in rp))

    def check(_root, g):
        p, ok = derive(g)
        if not ok:
            return False, "ideal syndrome or flag wrong"
        sig, bad = effect_signature(g, t, p)
        sig = proj(sig)
        a = [k for k, v in rsig.items() if k != ident and sig.get(k) != v]
        b = [k for k, v in sig.items() if k != ident and rsig.get(k) != v]
        if a or b:
            order = "single-fault" if t == 1 else f"{t}-fault"
            return False, f"not {order} equivalent ({len(a)} only in root, {len(b)} only in child)"
        if bad:
            return False, f"not FT_{t}"
        return True, "ok"
    return check


# Goto's |+>_L ancilla block for Steane-style Z-syndrome extraction (Jeon Fig. 3a, PDF page 9).
# Block qubits q7..q13 pair with data q0..q6; q14 is the verification qubit (flag).
GOTO = {
    "name": "Goto baseline (Fig. 3a)",
    "kind": "Z",
    "targets": [(0, 1, 2, 3), (1, 2, 4, 5), (2, 3, 5, 6)],
    "syndrome": list(range(7, 14)),
    "flags": [14],
    "ops": [("R", q) for q in range(7, 15)]
    + [("H", q) for q in (8, 9, 10, 12, 14)]
    + [("CX", 8, 7), ("CX", 12, 11),
       ("CX", 10, 13),
       ("CX", 9, 11), ("CX", 12, 13),
       ("CX", 10, 7),
       ("CX", 8, 11),
       ("CX", 9, 10),
       ("CX", 0, 7), ("CX", 14, 8),
       ("CX", 1, 8), ("CX", 14, 10),
       ("CX", 2, 9), ("CX", 14, 12),
       ("CX", 3, 10), ("H", 14),
       ("CX", 4, 11),
       ("CX", 5, 12),
       ("CX", 6, 13)]
    + [("M", q) for q in range(7, 15)],
}
