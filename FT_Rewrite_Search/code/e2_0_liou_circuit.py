"""E2.0 step 3-6: transcribe Liou Fig. 20 ([2;2;2;1;1] flag sharing, Z stabilizers of [[17,1,5]])
from results/liou_fig20_cnots.txt into an op list, then check it.

Qubits: data 0-16 (paper 1-17), measurement qubits m1..m8 -> 17..24, flags f1..f7 -> 25..31.
Measurement qubits start in |0> and are measured in Z; flags start in |+> (R, H) and are
measured in X (H, M), as in the paper's Fig. 3 [inferred: Fig. 20 omits them]; the noiseless
check below confirms the choice.
Writes data/circuits/liou_17_1_5.json and results/e2_0_liou_circuit.json.
Run:  python e2_0_liou_circuit.py
"""
import json
import os
import re
import time
from collections import Counter

from ftrs.codes import COLOR_17, COLOR_17_SUPPORTS
from ftrs.ftlib import effect_signature, resources
from ftrs.multi import derive

here = os.path.dirname(os.path.abspath(__file__))
lines = open(os.path.join(here, "..", "results", "liou_fig20_cnots.txt")).read().splitlines()
cnots = [(int(s), c, t) for s, c, t in
         (re.match(r"step\s+(\d+)\s+x=\s*[\d.]+\s+CX\s+(\S+) -> (\S+)", l).groups() for l in lines
          if l.startswith("step"))]


def q(label):
    if label.startswith("m"):
        return 16 + int(label[1:])
    if label.startswith("f"):
        return 24 + int(label[1:])
    return int(label) - 1


meas, flags = list(range(17, 25)), list(range(25, 32))
ops = [("R", m) for m in meas] + [("R", f) for f in flags] + [("H", f) for f in flags]
# Gates are kept in the figure's left-to-right order. The dashed separators group gates but are
# not strict time steps (group 6 has two CNOTs onto m8), so the depth comes from ASAP layering.
steps = {}
for s, c, t in cnots:
    steps.setdefault(s, []).append(("CX", q(c), q(t)))
clashes = {s: g for s, g in steps.items()
           if len([x for op in g for x in op[1:]]) != len({x for op in g for x in op[1:]})}
print("dashed groups that reuse a qubit:", {s: [op for op in g] for s, g in clashes.items()})
for s in sorted(steps):
    ops += steps[s]
ops += [("H", f) for f in flags] + [("M", m) for m in meas + flags]
gadget = {"name": "Liou 2024 Fig. 20, [2;2;2;1;1] flag sharing (Z stabilizers)", "kind": "Z",
          "targets": COLOR_17_SUPPORTS, "syndrome": meas, "flags": flags, "ops": ops, "code": COLOR_17}

out = {"cnots": len(cnots), "steps": len(steps), "ancillas": len(meas) + len(flags)}
res = resources(ops)
out.update({"resources": res})
print(f"CNOTs {len(cnots)}, time steps {len(steps)}, ancillas {len(meas)} + {len(flags)}; resources {res}")
print("paper Table I [2;2;2;1;1]: 56 CNOTs, depth 14, 17+15 qubits")

parities, ok = derive(gadget, shots=2000)
out["ideal_ok"] = ok
out["parities"] = parities
want = [[17 + i] for i in range(8)]
print("noiseless check (syndrome bit i = measurement qubit m_i alone, flags always 0):", ok,
      "| each stabilizer read by its own m:", parities == want)
out["each_stabilizer_own_m"] = parities == want

t0 = time.time()
for t in (1, 2):
    sig, bad = effect_signature(gadget, t, parities)
    out[f"single_round_FT{t}_violations"] = bad
    print(f"single-round FT_{t} (plan Sec. 0.3; only reported, not the pass criterion): "
          f"{bad} violating fault sets, {len(sig)} distinct unflagged effects, {time.time() - t0:.1f} s")

os.makedirs(os.path.join(here, "..", "data", "circuits"), exist_ok=True)
json.dump({k: v for k, v in gadget.items() if k != "code"} | {"code": COLOR_17.name},
          open(os.path.join(here, "..", "data", "circuits", "liou_17_1_5.json"), "w"), indent=1)
json.dump(out, open(os.path.join(here, "..", "results", "e2_0_liou_circuit.json"), "w"), indent=1)
print("saved")
