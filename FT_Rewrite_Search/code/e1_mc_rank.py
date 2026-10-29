"""Monte Carlo of the circuits reachable from the baseline (E1.2) plus Jeon Fig. 1b.

Metric (as in L02): single gadget, data in |0>^7, flag silent and X error of weight >= 2 modulo
the X stabilizers. Usage: python e1_mc_rank.py <p> <shots> <out.json>
"""
import json
import sys

import numpy as np

from ftrs.ftlib import REDUCED, STEANE_SUPPORTS, measured_qubits, to_stim

p, shots, out = float(sys.argv[1]), int(float(sys.argv[2])), sys.argv[3]
data = json.load(open("e1_2_shor_search.json"))
circuits = [{"name": f"reachable #{k} ({c['ancillas']}/{c['cnots']}/{c['cnot_depth']})",
             "ops": [tuple(o) for o in c["ops"]], "syndrome": c["syndrome"], "flags": c["flags"]}
            for k, c in enumerate(data["circuits"])] + [{**REDUCED, "name": "Jeon Fig. 1b (3/7/6)"}]

xg = np.zeros((3, 7), bool)
for i, s in enumerate(STEANE_SUPPORTS):
    xg[i, list(s)] = True
group = np.array([(a * xg[0]) ^ (b * xg[1]) ^ (c * xg[2]) for a in (0, 1) for b in (0, 1) for c in (0, 1)])

rows = []
for g in circuits:
    sampler = to_stim(g["ops"], data_init="zero", p=p).compile_sampler(seed=11)
    meas = measured_qubits(g["ops"])
    fi = [meas.index(q) for q in g["flags"]]
    bad = flagged = 0
    left = shots
    while left:
        n = min(left, 2_000_000)
        s = sampler.sample(n)
        flag = s[:, fi].any(axis=1)
        d = s[:, len(meas):len(meas) + 7]
        w = np.min([(d ^ x).sum(axis=1) for x in group], axis=0)
        bad += int(((~flag) & (w >= 2)).sum())
        flagged += int(flag.sum())
        left -= n
    rows.append({"name": g["name"], "p": p, "shots": shots, "malignant": bad,
                 "rate": bad / shots, "flag_rate": flagged / shots})
    print(rows[-1], flush=True)
json.dump(rows, open(out, "w"), indent=1)
