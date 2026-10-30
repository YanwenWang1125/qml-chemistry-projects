"""E2.0 step 7, setting B (AUTOLOOP B4): LER of the Liou [[17,1,5]] round under our default noise.

One round = Liou Fig. 20 on the Z stabilizers, then its Hadamard dual (ftlib.dual) on the X
stabilizers. d = 5 rounds, the last one noiseless, then a Z (or X) data measurement.
Decoder: decoders.DemLookup (exact for <= 2 DEM errors, BP+OSD beyond; qualified 2026-10-29:
0 wrong on all single faults and 20000 random pairs, both bases). Samples until max_errors logical errors or
max_shots. Writes results/e2_0_ler_b_p<p>_<memory>.json.
Run:  python e2_0_ler_b.py <p> [Z|X] [max_errors] [max_shots]
"""
import json
import os
import sys
import time

import numpy as np

from ftrs.codes import COLOR_17
from ftrs.decoders import DemLookup
from ftrs.ftlib import dual
from ftrs.mem import wilson
from ftrs.memgen import build

p = float(sys.argv[1])
memory = sys.argv[2] if len(sys.argv) > 2 else "Z"
max_errors = int(sys.argv[3]) if len(sys.argv) > 3 else 400
max_shots = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 7
here = os.path.dirname(os.path.abspath(__file__))
L = json.load(open(os.path.join(here, "..", "data", "circuits", "liou_17_1_5.json")))
Z = {**L, "ops": [tuple(o) for o in L["ops"]], "targets": [tuple(t) for t in L["targets"]], "code": COLOR_17}
X = {**dual(Z), "targets": [tuple(s) for s in COLOR_17.x_supports]}

t0 = time.time()
circ, info = build(COLOR_17, [Z, X], rounds=5, memory=memory, p=p)
noiseless, _ = build(COLOR_17, [Z, X], rounds=5, memory=memory, p=0)
d, o = noiseless.compile_detector_sampler(seed=3).sample(200, separate_observables=True)
assert not d.any() and not o.any(), "noiseless memory circuit fires detectors or flips the observable"
dec = DemLookup(circ)
sampler = circ.compile_detector_sampler(seed=11)
errors = shots = 0
while errors < max_errors and shots < max_shots:
    dets, obs = sampler.sample(20000, separate_observables=True)
    pred = dec.decode_batch(dets)
    errors += int((pred != obs).any(axis=1).sum())
    shots += len(dets)
    print(f"p={p} {memory}: {shots} shots, {errors} logical errors, {time.time() - t0:.0f} s", flush=True)
lo, hi = wilson(errors, shots)
out = {"p": p, "memory": memory, "rounds": 5, "shots": shots, "errors": errors, "ler": errors / shots,
       "ci95": [lo, hi], "detectors": info["detectors"], "dem_errors": int(dec.fallback.H.shape[1]), "decoder_stats": dec.stats,
       "seconds": round(time.time() - t0, 1)}
print(json.dumps(out))
json.dump(out, open(os.path.join(here, "..", "results", f"e2_0_ler_b_p{p:g}_{memory}.json"), "w"), indent=1)
