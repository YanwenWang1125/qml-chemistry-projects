"""E1.4c: why is the reduced gadget (Fig. 1b) worse than the baseline under DemLookup?

E1.4b reproduced the baseline LER but gave ratio reduced/baseline = 1.33 (Z) / 1.18 (X) where
Jeon reports 0.82 / 0.84.  Two hypotheses:
  H1 decoder: at N=10, p=1e-3 most shots exceed 2 errors and fall to BP+OSD, which may be worse
     on the reduced gadget's DEM.  Test: per-path failure counts, and ratios at N=1 (where the
     exact <=2-error lookup covers almost every shot).
  H2 circuit: the reduced gadget has more weight-2 fault pairs that flip the logical (Jeon notes
     4/33 single faults of the reduced Z gadget leave a weight-1 residual).  Test: structure task.

Usage (one Slurm array task each):
  python e1_4c_reduced_diag.py list               -> number of tasks
  python e1_4c_reduced_diag.py <task_id> <outdir>
Task ids 0..len(JOBS)-1 sample one chunk; the last id runs the structure analysis.
"""
import itertools
import json
import os
import sys
import time

import numpy as np

from ftrs.decoders import DemLookup
from ftrs.ftlib import BASELINE, REDUCED
from ftrs.mem import build_circuit

GADGETS = {"baseline": BASELINE["ops"], "reduced": REDUCED["ops"]}
CHUNKS = 10
# (n_cycles, p, shots per chunk)
SETTINGS = [(1, 1e-3, 100000), (1, 3e-4, 100000), (10, 1e-3, 10000)]
JOBS = [(g, m, n, p, s, c) for (n, p, s) in SETTINGS for g in GADGETS for m in ("Z", "X")
        for c in range(CHUNKS)]


def sample_task(gadget, memory, n_cycles, p, shots, chunk):
    circ, _ = build_circuit(GADGETS[gadget], n_cycles=n_cycles, memory=memory, p=p)
    seed = 1000 * chunk + 17
    dets, obs = circ.compile_detector_sampler(seed=seed).sample(shots, separate_observables=True)
    dec = DemLookup(circ, max_iter=50, osd_order=6)
    paths = {"empty": [0, 0], "w1": [0, 0], "w2": [0, 0], "bposd": [0, 0]}  # [shots, failures]
    for d, o in zip(dets, obs):
        before = dict(dec.stats)
        pred = dec.decode_one(d)
        path = next(k for k in dec.stats if dec.stats[k] != before[k])
        paths[path][0] += 1
        paths[path][1] += int(bool((pred != o).any()))
    return {"gadget": gadget, "memory": memory, "n_cycles": n_cycles, "p": p, "shots": shots,
            "chunk": chunk, "seed": seed, "paths": paths,
            "failures": sum(v[1] for v in paths.values())}


def structure():
    """Single- and two-error structure of each gadget's DEM at N=1, p=1e-3."""
    out = {}
    for (g, ops), m in itertools.product(GADGETS.items(), ("Z", "X")):
        circ, _ = build_circuit(ops, n_cycles=1, memory=m, p=1e-3)
        dec = DemLookup(circ, max_iter=50, osd_order=6)
        H, L, pr = dec.fallback.H, dec.L, dec.fallback.priors
        n = H.shape[1]
        single_fail = sum(bool(((dec.decode_one(H[:, j]) ^ L[:, j]) % 2).any()) for j in range(n))
        # every pair: probability mass of pairs the decoder gets wrong (leading p^2 LER term)
        bad_mass, bad_pairs, n_pairs = 0.0, 0, 0
        for a in range(n):
            for b in range(a + 1, n):
                s = (H[:, a] ^ H[:, b]).astype(np.uint8)
                n_pairs += 1
                if (((dec.decode_one(s) ^ L[:, a] ^ L[:, b]) % 2).any()):
                    bad_pairs += 1
                    bad_mass += pr[a] * pr[b]
        undetectable = len(circ.search_for_undetectable_logical_errors(
            dont_explore_detection_event_sets_with_size_above=4,
            dont_explore_edges_with_degree_above=99,
            dont_explore_edges_increasing_symptom_degree=False))
        out[f"{g}/{m}"] = {"columns": n, "total_prior": float(pr.sum()),
                           "single_fail": single_fail, "pairs": n_pairs,
                           "bad_pairs": bad_pairs, "bad_pair_mass": bad_mass,
                           "min_undetectable_logical": undetectable,
                           "paths": dict(dec.stats)}
        print(g, m, out[f"{g}/{m}"], flush=True)
    return out


if __name__ == "__main__":
    if sys.argv[1] == "list":
        print(len(JOBS) + 1)
        sys.exit()
    tid, outdir = int(sys.argv[1]), sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    t = time.time()
    res = structure() if tid == len(JOBS) else sample_task(*JOBS[tid])
    res["task"], res["seconds"] = tid, round(time.time() - t, 1)
    json.dump(res, open(os.path.join(outdir, f"task_{tid:03d}.json"), "w"), indent=1)
    print("done", res["task"], res["seconds"], flush=True)
