"""Decoders for the flag/cat-type memory circuits of ftrs.mem.

1. `bposd_decoder(circuit)` -- BP+OSD on the detector error model (plan Sec. 0.7: graphlike
   MWPM is not allowed on flag/cat circuits, BP+OSD is the recommended replacement).
2. `LookupDecoder` -- the flag-conditioned lookup table of Jeon et al. Sec. V B 3, decoding one
   cycle at a time on the raw syndrome history.

Both decoders are scored identically: a shot is a logical error when the predicted logical
observable differs from the one actually measured.
"""
import numpy as np
import stim


# ------------------------------------------------------------------ BP+OSD on the DEM

def dem_matrices(dem):
    """(H, L, priors) of a detector error model: detectors x errors, observables x errors."""
    rows, obs, probs = [], [], []
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        p = inst.args_copy()[0]
        dets, ob = [], []
        for t in inst.targets_copy():
            if t.is_relative_detector_id():
                dets.append(t.val)
            elif t.is_logical_observable_id():
                ob.append(t.val)
        rows.append(sorted(set(dets)))
        obs.append(sorted(set(ob)))
        probs.append(p)
    n_det = dem.num_detectors
    H = np.zeros((n_det, len(rows)), dtype=np.uint8)
    L = np.zeros((dem.num_observables, len(rows)), dtype=np.uint8)
    for j, (dets, ob) in enumerate(zip(rows, obs)):
        H[dets, j] = 1
        L[ob, j] = 1
    return H, L, np.array(probs, dtype=np.float64)


class BpOsd:
    def __init__(self, circuit, max_iter=30, osd_order=4):
        import ldpc
        dem = circuit.detector_error_model(decompose_errors=False, allow_gauge_detectors=True)
        self.H, self.L, self.priors = dem_matrices(dem)
        self.dec = ldpc.BpOsdDecoder(
            self.H,
            error_channel=list(self.priors),
            bp_method="product_sum", max_iter=max_iter,
            osd_method="osd_cs", osd_order=osd_order)

    def decode_batch(self, dets):
        """Predicted observable flips, one row per shot."""
        out = np.zeros((dets.shape[0], self.L.shape[0]), dtype=np.uint8)
        for i in range(dets.shape[0]):
            e = np.asarray(self.dec.decode(dets[i].astype(np.uint8)), dtype=np.uint8)
            out[i] = (self.L @ e) % 2
        return out.astype(bool)


class DemLookup:
    """Exact most-likely decoding for syndromes explained by <= 2 DEM errors, else BP+OSD.

    BP+OSD on these DEMs misdecodes some single faults (returns 4-9 errors when one explains the
    syndrome), which leaves an O(p) LER term.  A pair (a, b) explaining syndrome s must contain
    an error touching the lowest fired detector, so pairs are found by scanning those columns
    and looking up s ^ col in the single-error table.
    """

    def __init__(self, circuit, **bposd_kw):
        self.fallback = BpOsd(circuit, **bposd_kw)
        H, L, pr = self.fallback.H, self.fallback.L, self.fallback.priors
        self.lp = np.log(pr / (1 - pr))
        self.Hp = np.packbits(H.T, axis=1)          # errors x packed detectors
        self.L = L
        self.single = {}                             # syndrome -> best single error
        for j in range(H.shape[1]):
            k = self.Hp[j].tobytes()
            if k not in self.single or self.lp[j] > self.lp[self.single[k]]:
                self.single[k] = j
        self.by_det = [np.flatnonzero(H[d]) for d in range(H.shape[0])]
        self.stats = {"empty": 0, "w1": 0, "w2": 0, "bposd": 0}

    def decode_one(self, s):
        if not s.any():
            self.stats["empty"] += 1
            return np.zeros(self.L.shape[0], dtype=np.uint8)
        sp = np.packbits(s)
        best, cols = -np.inf, None
        j = self.single.get(sp.tobytes())
        if j is not None:
            best, cols = self.lp[j], (j,)
        for a in self.by_det[int(np.flatnonzero(s)[0])]:
            b = self.single.get((sp ^ self.Hp[a]).tobytes())
            if b is not None and b != a and self.lp[a] + self.lp[b] > best:
                best, cols = self.lp[a] + self.lp[b], (a, b)
        if cols is None:
            self.stats["bposd"] += 1
            e = np.asarray(self.fallback.dec.decode(s.astype(np.uint8)), dtype=np.uint8)
            return (self.L @ e) % 2
        self.stats["w%d" % len(cols)] += 1
        return self.L[:, list(cols)].sum(axis=1) % 2

    def decode_batch(self, dets):
        """Predicted observable flips, one row per shot."""
        return np.array([self.decode_one(d) for d in dets], dtype=bool).reshape(
            dets.shape[0], self.L.shape[0])


def run_bposd(circuit, shots, seed=1, **kw):
    sampler = circuit.compile_detector_sampler(seed=seed)
    dets, obs = sampler.sample(shots, separate_observables=True)
    dec = BpOsd(circuit, **kw)
    pred = dec.decode_batch(dets)
    err = (pred != obs).any(axis=1)
    return int(err.sum()), int(shots), dets, obs, pred


# ------------------------------------------------------------------ flag lookup table

class LookupDecoder:
    """Jeon Sec. V B 3: one correction per cycle from a flag-conditioned single-fault table.

    The table maps a 6-bit syndrome to the data residual of the single fault that produces it.
    Rounds whose flag fired are dropped; Shor's rule (first syndrome repeated in two consecutive
    rounds, else the last round) chooses the round a cycle is decoded from.
    """

    def __init__(self, table, frame_syndrome, n_cycles, rounds_per_cycle=3, ngen=6):
        self.keys = np.array(list(table.keys()), dtype=np.uint8).reshape(-1, ngen)
        self.cx = np.array([table[k][0] for k in table], dtype=np.int64)
        self.cz = np.array([table[k][1] for k in table], dtype=np.int64)
        enc = np.zeros(len(self.keys), dtype=np.int64)
        for s in range(ngen):
            enc |= self.keys[:, s].astype(np.int64) << s
        self.order = np.argsort(enc)
        self.sorted_keys = enc[self.order]
        self.frame_syndrome = frame_syndrome
        self.n_cycles = n_cycles
        self.rpc = rounds_per_cycle
        self.ngen = ngen

    def decode(self, syn, flg):
        n = syn.shape[0]
        fx = np.zeros(n, dtype=np.int64)
        fz = np.zeros(n, dtype=np.int64)
        rows = np.arange(n)
        for c in range(self.n_cycles):
            r0 = c * self.rpc
            rounds = syn[:, r0:r0 + self.rpc]
            usable = flg[:, r0:r0 + self.rpc].sum(axis=2) == 0
            idx = np.full(n, self.rpc - 1, dtype=np.int64)
            for r in range(self.rpc - 2, -1, -1):
                same = (rounds[:, r] == rounds[:, r + 1]).all(axis=1)
                idx[usable[:, r] & same] = r
            chosen = rounds[rows, idx]
            det = chosen ^ self.frame_syndrome(fx, fz)
            enc = np.zeros(n, dtype=np.int64)
            for s in range(self.ngen):
                enc |= det[:, s].astype(np.int64) << s
            pos = np.clip(np.searchsorted(self.sorted_keys, enc), 0, len(self.order) - 1)
            ok = self.sorted_keys[pos] == enc
            src = self.order[pos]
            fx ^= np.where(ok, self.cx[src], 0)
            fz ^= np.where(ok, self.cz[src], 0)
        return fx, fz
