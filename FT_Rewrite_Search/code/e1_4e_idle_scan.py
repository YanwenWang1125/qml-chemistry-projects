"""E1.4e: idle-noise strength x schedule scan under Jeon's flag protocol (ftrs.flagdec).

E1.4d left one question: with no idle noise the baseline matches Jeon's absolute LER (4.2e-2)
but reduced/baseline = 0.56; serial scheduling with idle p/10 matches the ratio (0.84) but the
baseline is 2.5x too high.  Scan idle_factor x schedule ("asap", "jit", "serial") for a point
that matches both (baseline ~4.23e-2 Z / 4.03e-2 X, ratio ~0.82 / 0.84).

N = 10, p = 1e-3, both memories, 2e5 shots per circuit.

Usage:  python e1_4e_idle_scan.py list | <task_id> <outdir>
"""
import json
import os
import sys
import time

import numpy as np

from ftrs.flagdec import FlagDecoder
from ftrs.ftlib import BASELINE, REDUCED
from ftrs.mem import build_circuit, wilson

GADGETS = {"baseline": BASELINE["ops"], "reduced": REDUCED["ops"]}
IDLE = (0.0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1)
JOBS = [(s, f) for s in ("asap", "jit", "serial") for f in IDLE]
SHOTS = 200000


def task(mode, idle):
    out = []
    for m in "ZX":
        for g, ops in GADGETS.items():
            c, info = build_circuit(ops, n_cycles=10, memory=m, p=1e-3,
                                    schedule_mode=mode, idle_factor=idle)
            lg, st = FlagDecoder(ops, info, m).decode(
                c.compile_sampler(seed=31).sample(SHOTS))
            k = int(lg.sum())
            out.append({"schedule": mode, "idle_factor": idle, "memory": m, "gadget": g,
                        "shots": SHOTS, "failures": k, "ler": k / SHOTS,
                        "ci95": [float(x) for x in wilson(k, SHOTS)], **st})
    return out


if __name__ == "__main__":
    if sys.argv[1] == "list":
        print(len(JOBS))
        sys.exit()
    tid, outdir = int(sys.argv[1]), sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    t = time.time()
    rows = task(*JOBS[tid])
    json.dump({"task": tid, "rows": rows, "seconds": round(time.time() - t, 1)},
              open(os.path.join(outdir, f"task_{tid:02d}.json"), "w"), indent=1)
    print("done", tid, round(time.time() - t, 1), flush=True)
