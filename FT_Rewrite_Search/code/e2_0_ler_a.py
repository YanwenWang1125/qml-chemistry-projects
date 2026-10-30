"""E2.0 step 7, setting A (simplified, AUTOLOOP B4): the Liou round under Liou's noise model.

Noise as in Liou & Lai 2024 Sec. II B with gamma = beta = 1 (their "all components have the same
error rate" case): two-qubit gates DEPOLARIZE2(p), single-qubit gates DEPOLARIZE1(p), preparation
and measurement flips 2p/3, idle DEPOLARIZE1(p) per layer. Not their protocol: we run 5 rounds
(4 noisy + 1 perfect) of Z half + X half and decode with DemLookup, then convert to a per-round
logical error rate 1 - (1 - LER)^(1/4). Their reported pseudothreshold for [[17,1,5]] is about
6e-5 (p.2, p.13); agreement within the same order of magnitude is the check.
Writes results/e2_0_ler_a_p<p>_<memory>.json.
Run:  python e2_0_ler_a.py <p> [Z|X] [max_errors]
"""
import json
import os
import sys
import time

from ftrs.codes import COLOR_17
from ftrs.decoders import DemLookup
from ftrs.ftlib import dual
from ftrs.mem import wilson
from ftrs.memgen import build

p = float(sys.argv[1])
memory = sys.argv[2] if len(sys.argv) > 2 else "Z"
max_errors = int(sys.argv[3]) if len(sys.argv) > 3 else 400
here = os.path.dirname(os.path.abspath(__file__))
L = json.load(open(os.path.join(here, "..", "data", "circuits", "liou_17_1_5.json")))
Z = {**L, "ops": [tuple(o) for o in L["ops"]], "targets": [tuple(t) for t in L["targets"]], "code": COLOR_17}
X = {**dual(Z), "targets": [tuple(s) for s in COLOR_17.x_supports]}
noise = dict(h_factor=1.0, prep_factor=2 / 3, meas_factor=2 / 3, idle_factor=1.0)

t0 = time.time()
circ, _ = build(COLOR_17, [Z, X], rounds=5, memory=memory, p=p, **noise)
dec = DemLookup(circ)
sampler = circ.compile_detector_sampler(seed=13)
errors = shots = 0
while errors < max_errors and shots < 5 * 10 ** 8:
    d, o = sampler.sample(50000, separate_observables=True)
    errors += int((dec.decode_batch(d) != o).any(axis=1).sum())
    shots += len(d)
lo, hi = wilson(errors, shots)
per_round = lambda x: 1 - (1 - x) ** 0.25
out = {"p": p, "memory": memory, "noise": noise, "rounds": 5, "noisy_rounds": 4, "shots": shots, "errors": errors,
       "ler": errors / shots, "ci95": [lo, hi], "ler_per_round": per_round(errors / shots),
       "ci95_per_round": [per_round(lo), per_round(hi)], "decoder_stats": dec.stats, "seconds": round(time.time() - t0, 1)}
print(json.dumps(out))
json.dump(out, open(os.path.join(here, "..", "results", f"e2_0_ler_a_p{p:g}_{memory}.json"), "w"), indent=1)
