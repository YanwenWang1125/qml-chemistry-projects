"""Jeon et al. Sec. V B 3: repeated-syndrome protocol with flag-conditioned lookup tables.

Protocol (paper text, Sec. V B 3):
  * A cycle is three rounds followed by one decoding step, using Shor's rule: take the first
    syndrome that repeats in two consecutive rounds, otherwise the last one measured.
  * Corrections accumulate in a Pauli frame; the decoder acts on the difference between the
    reported syndrome and what the frame already accounts for.
  * When a flag is raised in a round, that round's syndrome is discarded entirely.  The next
    round's syndrome is looked up in the table of the gadget that raised the flag, and the
    correction goes into the frame.  (The last round of ftrs.mem is noiseless, so a flag in the
    last noisy round is resolved by it; the perfect data readout is decoded at the end.)
  * Each gadget's table: enumerate all single faults of a round, keep those that raise that
    gadget's flag, map the syndrome of the residual (projected on the memory's sensitive basis,
    modulo stabilizers) to the residual itself.

Choices the paper does not fix (ours):
  * several flags in one round: the first flagged gadget in time order selects the table;
  * "consecutive rounds" in Shor's rule means consecutive among the cycle's unflagged rounds;
  * each round's syndrome is compared with the frame as it stood at that round.

Only the error type the memory is sensitive to is tracked: X errors (Z-type syndromes) for
|0>_L, Z errors (X-type syndromes) for |+>_L.  Everything is vectorised over shots, with the
7-qubit error as a bitmask.
"""
import numpy as np

from .ftlib import STEANE_SUPPORTS, effect, single_faults
from .mem import N_DATA, NGEN, round_ops

SUP_MASKS = [sum(1 << q for q in s) for s in STEANE_SUPPORTS]


def _parity(x):
    return bin(x).count("1") & 1


# syndrome (3 bits, bit k = generator k) of every 7-bit error of one type
SYN = np.array([sum(_parity(e & m) << k for k, m in enumerate(SUP_MASKS))
                for e in range(1 << N_DATA)], dtype=np.int64)
# stabilizers of one type (the 8 products of the three generators)
STABS = sorted({a ^ b ^ c for a in (0, SUP_MASKS[0]) for b in (0, SUP_MASKS[1])
                for c in (0, SUP_MASKS[2])})


def min_rep(e):
    """Minimum-weight element of e + stabilizers."""
    return min((e ^ s for s in STABS), key=lambda v: (bin(v).count("1"), v))


# Hamming decoder: syndrome -> weight-<=1 correction
HAM = np.zeros(8, dtype=np.int64)
for _q in range(N_DATA):
    HAM[SYN[1 << _q]] = 1 << _q


def flag_tables(template_ops, memory="Z"):
    """Per-gadget tables, shape (6, 8): syndrome -> correction mask; plus a summary."""
    ops, gs = round_ops(template_ops, N_DATA)
    flags = [g["flag"] for g in gs]
    entries = [dict() for _ in gs]
    n_faults = [0] * len(gs)
    conflicts = 0
    for _label, start, x, z in single_faults(ops):
        dx, dz, _s, fflips = effect(ops, start, x, z, [], flags)
        for gi, f in enumerate(fflips):
            if not f:
                continue
            n_faults[gi] += 1
            e = min_rep(dx if memory == "Z" else dz)
            s = int(SYN[e])
            if s in entries[gi] and min_rep(entries[gi][s] ^ e) != 0:
                conflicts += 1
                if bin(e).count("1") < bin(entries[gi][s]).count("1"):
                    entries[gi][s] = e
            else:
                entries[gi].setdefault(s, e)
    T = np.zeros((len(gs), 8), dtype=np.int64)
    for gi, ent in enumerate(entries):
        for s, e in ent.items():
            if s:
                T[gi, s] = e
    summary = [{"gadget": f"{g['type']}{list(g['support'])}", "faults_raising_flag": n,
                "table": {format(s, "03b"): format(e, "07b")[::-1] for s, e in ent.items()
                          if s and e}}
               for g, n, ent in zip(gs, n_faults, entries)]
    return T, {"gadgets": summary, "conflicts": conflicts}


class FlagDecoder:
    def __init__(self, template_ops, info, memory="Z", rounds_per_cycle=3):
        self.T, self.summary = flag_tables(template_ops, memory)
        self.info, self.memory, self.rpc = info, memory, rounds_per_cycle
        self.n_cat = len(info["spec"]["cat"])
        # Z-type gadgets come first in a round (ftrs.mem.round_gadgets)
        self.gsel = slice(0, 3) if memory == "Z" else slice(3, 6)

    def _syndromes(self, m):
        info, n = self.info, m.shape[0]
        R = info["n_rounds"]
        syn = np.zeros((n, R), dtype=np.int64)
        flg = np.zeros((n, R, NGEN), dtype=bool)
        gids = list(range(NGEN))[self.gsel]
        for r in range(R):
            for k, gi in enumerate(gids):
                b = np.zeros(n, dtype=np.int64)
                for j in range(self.n_cat):
                    b ^= m[:, info["rec_anc"][(r, gi, j)]].astype(np.int64)
                syn[:, r] |= b << k
            for gi in range(NGEN):
                flg[:, r, gi] = m[:, info["rec_flag"][(r, gi)]]
        return syn, flg

    def decode(self, m):
        """Logical outcome (0 = no logical error) and flag statistics from raw measurements."""
        n = m.shape[0]
        syn, flg = self._syndromes(m)
        R = syn.shape[1]
        frame = np.zeros(n, dtype=np.int64)
        pend = np.full(n, -1, dtype=np.int64)
        rows = np.arange(n)
        for c0 in range(0, R, self.rpc):
            found = np.zeros(n, bool)
            have_prev = np.zeros(n, bool)
            any_u = np.zeros(n, bool)
            prev = np.zeros(n, np.int64)
            last = np.zeros(n, np.int64)
            chosen = np.zeros(n, np.int64)
            for r in range(c0, min(c0 + self.rpc, R)):
                s = syn[:, r]
                hit = pend >= 0
                if hit.any():
                    d = s ^ SYN[frame]
                    corr = self.T[np.where(hit, pend, 0), d]
                    frame ^= np.where(hit, corr, 0)
                    pend[:] = -1
                det = s ^ SYN[frame]
                fl = flg[:, r, :]
                flagged = fl.any(axis=1)
                pend = np.where(flagged, fl.argmax(axis=1), -1)
                u = ~flagged
                match = u & have_prev & (det == prev) & ~found
                chosen = np.where(match, prev, chosen)
                found |= match
                prev = np.where(u, det, prev)
                last = np.where(u, det, last)
                have_prev |= u
                any_u |= u
            pick = np.where(found, chosen, last)
            frame ^= np.where(any_u, HAM[pick], 0)
        # perfect readout of the data, decoded with the Hamming table
        rd = self.info["rec_data"]
        data = np.zeros(n, dtype=np.int64)
        for q in range(N_DATA):
            data |= m[:, rd + q].astype(np.int64) << q
        data ^= frame
        data ^= HAM[SYN[data]]
        logical = np.array([_parity(int(v)) for v in data], dtype=np.uint8)
        flags_per_cycle = flg.any(axis=2).sum(axis=1) / (R / self.rpc)
        return logical, {"flag_rounds_per_cycle": float(flags_per_cycle.mean()),
                         "shots_with_flag": float(flg.any(axis=(1, 2)).mean())}
