"""Item 2, next step: weight-4 gadgets of A-2 (one measurement qubit, one flag bracketing data
CNOTs 1-3) plus Liou's own g8 gadget (m8 with flags f4-f6, CNOT sequence of Fig. 20), then the
corrected greedy flag merge and interleaving, accepted on the memory circuit distance.

Why: Liou's g8 alone violates the per-gadget condition K of d1 (8774 fault sets) although Liou's
round has circuit distance 5, so K is sufficient but not necessary and made A-2 put 5 flags on
weight 8. Here acceptance is only noiseless + single-round FT_2 + memory distance 5.
Writes results/a5_liou_g8.json, data/circuits/a5_*.json.
Run:  python a5_liou_g8.py [max_errors]
"""
import itertools
import json
import os
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
max_errors = int(sys.argv[1]) if len(sys.argv) > 1 else 400
src = open(os.path.join(here, "a2_flag_start.py"), encoding="utf-8").read()
ns = {"__file__": os.path.join(here, "a2_flag_start.py")}
argv, sys.argv = sys.argv, ["x", str(max_errors)]
exec(src[:src.index("LIOU_G8 = (13, 9, 14, 5, 2, 6, 3, 10)\nper = {}")], ns)
exec(src[src.index("def halves("):src.index("ok, why = accepted(Zh)")], ns)
isrc = open(os.path.join(here, "a3_interleave.py"), encoding="utf-8").read()
exec(isrc[isrc.index("def interleave("):isrc.index("d = json.load(open(path))")], ns)
sys.argv = argv
COLOR_17, gadget, accepted, merge, ler, resources, interleave = (ns[k] for k in (
    "COLOR_17", "gadget", "accepted", "merge", "ler", "resources", "interleave"))
t0 = time.time()

L = json.load(open(os.path.join(here, "..", "data", "circuits", "liou_17_1_5.json")))
lops = [tuple(o) for o in L["ops"]]
ops, syn, flags, base = [], [], [], 17
for sup in COLOR_17.z_supports[:7]:
    g = gadget(sup, [(0, 2)], base)
    ops += g["ops"]
    syn += g["syndrome"]
    flags += g["flags"]
    base += 2
m8 = {24: base, 28: base + 1, 29: base + 2, 30: base + 3}
g8 = [tuple([o[0]] + [m8.get(q, q) for q in o[1:]]) for o in lops
      if (o[0] == "CX" and o[2] == 24) or (o[0] != "CX" and o[1] in m8)]
ops += g8
syn.append(base)
flags += [base + 1, base + 2, base + 3]
Z0 = {"name": "w4 CR gadgets + Liou g8", "kind": "Z", "targets": [tuple(s) for s in COLOR_17.z_supports],
      "syndrome": syn, "flags": flags, "ops": ops, "code": COLOR_17}
ok, why = accepted(Z0)
out = {"start": resources(ops), "start_flags": len(flags), "start_accepted": ok, "start_reason": why}
print(json.dumps(out), flush=True)
cur, steps = Z0, []
while ok:
    best, verdicts = None, {}
    for fa, fb in itertools.combinations(cur["flags"], 2):
        cand = merge(cur, fa, fb)
        good, w = accepted(cand)
        verdicts[w] = verdicts.get(w, 0) + 1
        if good and (best is None or len(cand["ops"]) <= len(best[0]["ops"])):
            best = (cand, (fa, fb))
    print(f"merge step {len(steps) + 1}: {verdicts} ({time.time() - t0:.0f} s)", flush=True)
    if best is None:
        break
    cur = best[0]
    steps.append({"merge": best[1], "resources": resources(cur["ops"]), "flags": len(cur["flags"])})
for name, g in (("start", Z0), ("final", cur)):
    if not ok:
        break
    iops, depth = interleave(g["ops"])
    gi = {**g, "ops": iops}
    good, w = accepted(gi)
    out[name + "_interleaved"] = {"resources": resources(iops), "flags": len(g["flags"]), "accepted": good, "reason": w}
    if good:
        out[name + "_interleaved"]["ler"] = [ler(gi, p, m) for p in (1e-3, 2e-3) for m in ("Z", "X")]
    json.dump({k: v for k, v in gi.items() if k != "code"} | {"code": COLOR_17.name},
              open(os.path.join(here, "..", "data", "circuits", f"a5_{name}.json"), "w"), indent=1)
out["merge_steps"] = steps
out["seconds"] = round(time.time() - t0, 1)
print(json.dumps(out))
json.dump(out, open(os.path.join(here, "..", "results", "a5_liou_g8.json"), "w"), indent=1)
