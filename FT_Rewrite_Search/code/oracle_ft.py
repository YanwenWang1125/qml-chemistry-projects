"""Independent FT_t oracle for E1.5b (supplementary plan): does not import ftlib.propagate.

Each elementary fault (same circuit-level model as ftlib: X after R; X/Y/Z after H; the 15
two-qubit Paulis after CX; X before M) is pushed to the end with stim: the H/CX ops after the
fault location are compiled to a stim.Tableau, which conjugates the fault's PauliString.
A set of faults at distinct locations has the XOR of their end Paulis (linearity).
Residual weight modulo the stabilizer group is a brute-force minimum over the full group,
with its own popcount table (not ftlib.Code).

FT_t (plan Sec. 0.3): every unflagged set of s <= t faults leaves weight <= s on the data
modulo S. Flags are measured in Z at the end (X-basis flags carry their own H in the op list),
so a set is flagged when its end Pauli has X or Y on a flag qubit.
"""
import itertools

import numpy as np
import stim

_PAULI1 = {"X": (1, 0), "Y": (1, 1), "Z": (0, 1)}


def _elementary(ops):
    """[(location, label, start index, {qubit: 'X'|'Y'|'Z'})]."""
    out = []
    for i, op in enumerate(ops):
        if op[0] == "R":
            out.append((i, f"prep X q{op[1]}", i + 1, {op[1]: "X"}))
        elif op[0] == "H":
            for p in "XYZ":
                out.append((i, f"{p} q{op[1]} after H", i + 1, {op[1]: p}))
        elif op[0] == "CX":
            for a, b in itertools.product("IXYZ", repeat=2):
                if a == b == "I":
                    continue
                d = {}
                if a != "I":
                    d[op[1]] = a
                if b != "I":
                    d[op[2]] = b
                out.append((i, f"{a}{b} on CX{op[1]}->{op[2]}", i + 1, d))
        elif op[0] == "M":
            out.append((i, f"meas flip q{op[1]}", i, {op[1]: "X"}))
    return out


def end_paulis(ops):
    """[(location, label, x mask, z mask)] of every elementary fault at the end of the circuit."""
    nq = 1 + max(q for op in ops for q in op[1:])
    tabs = {}
    out = []
    for loc, label, start, paulis in _elementary(ops):
        if start not in tabs:
            c = stim.Circuit()
            c.append("I", list(range(nq)))  # fixes the tableau size to all qubits
            for op in ops[start:]:
                if op[0] in ("H", "CX"):
                    c.append(op[0], list(op[1:]))
            tabs[start] = c.to_tableau()
        ps = stim.PauliString(nq)
        for q, p in paulis.items():
            ps[q] = p
        xs, zs = tabs[start](ps).to_numpy()
        x = sum(1 << int(q) for q in np.nonzero(xs)[0])
        z = sum(1 << int(q) for q in np.nonzero(zs)[0])
        out.append((loc, label, x, z))
    return out


class Weight:
    """Minimum weight of a data Pauli times any stabilizer, by brute force over the group."""

    def __init__(self, n_data, x_supports, z_supports):
        def span(sups):
            els = [0]
            for s in sups:
                m = sum(1 << q for q in s)
                els = els + [e ^ m for e in els]
            return np.array(sorted(set(els)), dtype=np.int64)
        self.sx, self.sz = span(x_supports), span(z_supports)
        self.pc = np.array([bin(v).count("1") for v in range(1 << n_data)], dtype=np.int64)
        self.mask = (1 << n_data) - 1
        self.cache = {}

    def __call__(self, x, z):
        key = (x & self.mask, z & self.mask)
        if key not in self.cache:
            self.cache[key] = int(self.pc[(key[0] ^ self.sx)[:, None] | (key[1] ^ self.sz)[None, :]].min())
        return self.cache[key]


def ft_oracle(ops, flags, weight, t):
    """(FT_t holds, first violating set or None, number of sets checked)."""
    cols = end_paulis(ops)
    fmask = sum(1 << q for q in flags)
    by_loc = {}
    for loc, label, x, z in cols:
        by_loc.setdefault(loc, []).append((label, x, z))
    locs = sorted(by_loc)
    checked = 0
    for s in range(1, t + 1):
        for combo in itertools.combinations(locs, s):
            for choice in itertools.product(*(by_loc[l] for l in combo)):
                x = z = 0
                for _, cx, cz in choice:
                    x ^= cx
                    z ^= cz
                checked += 1
                if x & fmask:
                    continue
                if weight(x, z) > s:
                    return False, [c[0] for c in choice], checked
    return True, None, checked
