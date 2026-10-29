"""E1.4d: Jeon's flag-conditioned protocol (ftrs.flagdec, Sec. V B 3) and the effect of scheduling.

E1.4c showed that under near-optimal DEM decoding the reduced gadget is 13-17% worse than the
baseline, while Jeon reports ratios 0.82 / 0.84.  This script decodes with the paper's own
protocol and varies what the paper does not pin down for the Shor case: how a round is scheduled
("asap": all cat states prepared up front and idling; "serial": one gadget after another) and
whether idle noise is applied.

  1. tables vs Table II (printed), p = 0 sanity check;
  2. N = 10, p = 1e-3: LER, flag rate and ratio for every (schedule, idle) setting;
  3. N = 1, p sweep: log-log slopes (paper: 1.94-1.97 for p <= 2e-3).

Outputs results/e1_4d_flag_protocol.json and _output.txt.
"""
import json
import os
import time

import numpy as np

from ftrs.flagdec import FlagDecoder, flag_tables
from ftrs.ftlib import BASELINE, REDUCED
from ftrs.mem import build_circuit, wilson

T0 = time.time()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
log = open(os.path.join(OUT, "e1_4d_flag_protocol_output.txt"), "w", encoding="utf-8")
GADGETS = {"baseline": BASELINE["ops"], "reduced": REDUCED["ops"]}
SETTINGS = [("asap", 0.1), ("serial", 0.1), ("asap", 0.0), ("serial", 0.0)]
PAPER = {"Z": (4.23e-2, 3.45e-2), "X": (4.03e-2, 3.37e-2)}


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.write(s + "\n")


def run(ops, n_cycles, memory, p, shots, mode, idle, seed):
    c, info = build_circuit(ops, n_cycles=n_cycles, memory=memory, p=p,
                            schedule_mode=mode, idle_factor=idle)
    lg, st = FlagDecoder(ops, info, memory).decode(c.compile_sampler(seed=seed).sample(shots))
    k = int(lg.sum())
    return {"n_cycles": n_cycles, "memory": memory, "p": p, "shots": shots, "schedule": mode,
            "idle_factor": idle, "failures": k, "ler": k / shots,
            "ci95": [float(x) for x in wilson(k, shots)], **st}


rows = {"tables": {}, "p0": [], "n10": [], "sweep": []}

# ---- 1. tables and p = 0
for g, ops in GADGETS.items():
    for m in "ZX":
        rows["tables"][f"{g}/{m}"] = flag_tables(ops, m)[1]
        for mode, idle in SETTINGS[:2]:
            r = run(ops, 10, m, 0.0, 2000, mode, idle, 1)
            rows["p0"].append({"gadget": g, **r})
            assert r["failures"] == 0, r
say("tables (Z memory):", json.dumps({k: v for k, v in rows["tables"].items() if k.endswith("Z")},
                                     indent=1))
say("p=0: all settings 0 failures")

# ---- 2. N = 10, p = 1e-3
for mode, idle in SETTINGS:
    for m in "ZX":
        res = {g: run(ops, 10, m, 1e-3, 200000, mode, idle, 11) for g, ops in GADGETS.items()}
        for g, r in res.items():
            rows["n10"].append({"gadget": g, **r})
        b, rd = res["baseline"], res["reduced"]
        ratio = rd["ler"] / b["ler"]
        se = np.sqrt(1 / rd["failures"] + 1 / b["failures"])
        say(f"{mode:6s} idle={idle:<4} {m}: baseline {b['ler']:.3e} "
            f"[{b['ci95'][0]:.2e},{b['ci95'][1]:.2e}] flags/cycle {b['flag_rounds_per_cycle']:.3f} | "
            f"reduced {rd['ler']:.3e} [{rd['ci95'][0]:.2e},{rd['ci95'][1]:.2e}] "
            f"flags/cycle {rd['flag_rounds_per_cycle']:.3f} | ratio {ratio:.3f} "
            f"[{ratio * np.exp(-1.96 * se):.3f},{ratio * np.exp(1.96 * se):.3f}]  "
            f"(paper {PAPER[m][0]:.2e} -> {PAPER[m][1]:.2e}, {PAPER[m][1] / PAPER[m][0]:.2f})")

# ---- 3. N = 1 slopes
PS = (2.5e-4, 5e-4, 1e-3, 2e-3)
for mode, idle in SETTINGS[:2]:
    for g, ops in GADGETS.items():
        pts = []
        for p in PS:
            r = run(ops, 1, "Z", p, 2000000, mode, idle, 21)
            rows["sweep"].append({"gadget": g, **r})
            pts.append((p, r["ler"]))
        x, y = np.log([q for q, _ in pts]), np.log([v for _, v in pts])
        slope = float(np.polyfit(x, y, 1)[0])
        say(f"slope N=1 {mode} idle={idle} {g}: {slope:.2f}  LER " +
            ", ".join(f"{q:g}:{v:.2e}" for q, v in pts))

rows["seconds"] = round(time.time() - T0, 1)
json.dump(rows, open(os.path.join(OUT, "e1_4d_flag_protocol.json"), "w"), indent=1)
say("done in", rows["seconds"], "s")
