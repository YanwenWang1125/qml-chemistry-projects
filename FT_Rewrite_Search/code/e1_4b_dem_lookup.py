"""E1.4b: rerun E1.4 with an exact <=2-error DEM decoder (ftrs.decoders.DemLookup).

ClawsGO's E1.4 (BP+OSD only) failed the slope check; local diagnosis showed BP+OSD misdecodes
single faults although every single-fault syndrome is unambiguous.  This script
  0. checks DemLookup on every injected single DEM error and on random injected pairs (N=1),
  1. measures LER vs p on the baseline gadget (Z memory) at N=1 and N=10,
  2. measures Jeon's ratios (reduced / baseline) at N=10, p=1e-3, both memories.

Outputs results/e1_4b_dem_lookup.json and _output.txt.
"""
import json
import os
import sys
import time

import numpy as np

from ftrs.decoders import DemLookup
from ftrs.ftlib import BASELINE, REDUCED
from ftrs.mem import build_circuit, wilson

T0 = time.time()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
SHOTS_SCALE = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
log = open(os.path.join(OUT, "e1_4b_dem_lookup_output.txt"), "w", encoding="utf-8")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.write(s + "\n")


def measure(ops, n_cycles, memory, p, shots, seed):
    circ, _ = build_circuit(ops, n_cycles=n_cycles, memory=memory, p=p)
    dets, obs = circ.compile_detector_sampler(seed=seed).sample(shots, separate_observables=True)
    dec = DemLookup(circ, max_iter=50, osd_order=6)
    t = time.time()
    k = int((dec.decode_batch(dets) != obs).any(axis=1).sum())
    lo, hi = wilson(k, shots)
    return {"p": p, "n_cycles": n_cycles, "memory": memory, "shots": shots, "failures": k,
            "ler": k / shots, "ci95": [lo, hi], "path_counts": dec.stats,
            "decode_seconds": round(time.time() - t, 1)}


rows = {"inject": {}, "slope": [], "ratio": []}

# ---- 0. injected faults, N = 1
circ, _ = build_circuit(BASELINE["ops"], n_cycles=1, memory="Z", p=1e-3)
dec = DemLookup(circ, max_iter=50, osd_order=6)
H, L = dec.fallback.H, dec.L
n = H.shape[1]
single_fail = sum(bool(((dec.decode_one(H[:, j]) ^ L[:, j]) % 2).any()) for j in range(n))
rng = np.random.default_rng(5)
pairs = rng.choice(n, size=(4000, 2))
pairs = pairs[pairs[:, 0] != pairs[:, 1]]
pair_fail = sum(bool(((dec.decode_one((H[:, a] ^ H[:, b]).astype(np.uint8))
                       ^ (L[:, a] ^ L[:, b])) % 2).any()) for a, b in pairs)
rows["inject"] = {"columns": n, "single_fail": single_fail, "pairs": len(pairs),
                  "pair_fail": pair_fail, "paths": dict(dec.stats)}
say("inject", rows["inject"])

# ---- 1. LER vs p, baseline, Z memory
for N, plist in ((1, (1e-4, 3e-4, 1e-3)), (10, (3e-5, 1e-4, 3e-4, 1e-3))):
    for p in plist:
        shots = int(SHOTS_SCALE * {3e-5: 40000, 1e-4: 40000, 3e-4: 20000, 1e-3: 10000}[p]
                    * (4 if N == 1 else 1))
        r = measure(BASELINE["ops"], N, "Z", p, shots, seed=7)
        rows["slope"].append(r)
        say("slope", r)


def slopes(N):
    pts = [(r["p"], r["ler"]) for r in rows["slope"] if r["n_cycles"] == N and r["ler"] > 0]
    return {f"{a[0]:g}->{b[0]:g}": round(float(np.log(b[1] / a[1]) / np.log(b[0] / a[0])), 2)
            for a, b in zip(pts, pts[1:])}


rows["slopes"] = {"N=1": slopes(1), "N=10": slopes(10)}
say("slopes", rows["slopes"])

# ---- 2. Jeon's ratios, N = 10, p = 1e-3
for memory in ("Z", "X"):
    pair = {}
    for name, ops in (("baseline", BASELINE["ops"]), ("reduced", REDUCED["ops"])):
        r = measure(ops, 10, memory, 1e-3, int(SHOTS_SCALE * 10000), seed=11)
        r["gadget"] = name
        rows["ratio"].append(r)
        pair[name] = r["ler"]
        say("ratio", r)
    say(f"{memory} ratio reduced/baseline = {pair['reduced'] / pair['baseline']:.3f}")

rows["seconds"] = round(time.time() - T0, 1)
json.dump(rows, open(os.path.join(OUT, "e1_4b_dem_lookup.json"), "w"), indent=1)
say("done in", rows["seconds"], "s")
