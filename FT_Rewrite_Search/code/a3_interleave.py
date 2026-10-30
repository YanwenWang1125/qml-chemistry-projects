"""Follow-up item 2a: interleave a Z half (Liou Lemma 6) and re-check it.

In a Z half every CNOT targets a measurement qubit m: data -> m or flag -> m. CNOTs that share
a control commute, so only the order of the CNOTs on each m matters for the gadget's behaviour.
interleave() keeps each m's CNOT sequence and list-schedules all CNOTs into layers (each layer
uses a qubit once; the m with the most CNOTs left goes first). The result is checked like A-2
(noiseless, single-round FT_2, 5-round memory distance 5) and its setting-B LER is measured.
Writes results/a3_interleave_<name>.json and data/circuits/a3_<name>.json.
Run:  python a3_interleave.py <data/circuits/... .json> <name> [max_errors]
"""
import json
import os
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
path, name = sys.argv[1], sys.argv[2]
max_errors = int(sys.argv[3]) if len(sys.argv) > 3 else 400
src = open(os.path.join(here, "a2_flag_start.py"), encoding="utf-8").read()
ns = {"__file__": os.path.join(here, "a2_flag_start.py")}
argv, sys.argv = sys.argv, ["x", str(max_errors)]  # a2_flag_start reads max_errors from argv
exec(src[:src.index("per = {}")], ns)
sys.argv = argv
exec(src[src.index("def halves("):src.index("ok, why = accepted(Zh)")], ns)
COLOR_17, resources, accepted, ler = ns["COLOR_17"], ns["resources"], ns["accepted"], ns["ler"]


def interleave(ops):
    cx = [op for op in ops if op[0] == "CX"]
    first_cx, last_cx = {}, {}
    for i, op in enumerate(ops):
        if op[0] == "CX":
            for q in op[1:]:
                first_cx.setdefault(q, i)
                last_cx[q] = i
    # each ancilla's R/H before its first CNOT go in front, its H/M after its last CNOT at the end
    pre = [op for i, op in enumerate(ops) if op[0] != "CX" and i < first_cx.get(op[1], -1)]
    post = [op for i, op in enumerate(ops) if op[0] != "CX" and i > last_cx.get(op[1], 10 ** 9)]
    assert len(pre) + len(post) + len(cx) == len(ops), "an ancilla op sits between its CNOTs"
    seq = {}
    for op in cx:
        seq.setdefault(op[2], []).append(op)
    layers = []
    while any(seq.values()):
        used, layer = set(), []
        for m in sorted(seq, key=lambda q: -len(seq[q])):
            if seq[m] and m not in used and seq[m][0][1] not in used:
                op = seq[m].pop(0)
                layer.append(op)
                used |= {op[1], op[2]}
        layers.append(layer)
    return pre + [op for l in layers for op in l] + post, len(layers)


d = json.load(open(path))
Zh = {**d, "ops": [tuple(o) for o in d["ops"]], "targets": [tuple(t) for t in d["targets"]], "code": COLOR_17}
t0 = time.time()
new_ops, depth = interleave(Zh["ops"])
Zi = {**Zh, "ops": new_ops}
ok, why = accepted(Zi)
out = {"name": name, "before": resources(Zh["ops"]), "after": resources(new_ops), "layers": depth,
       "flags": len(Zh["flags"]), "accepted": ok, "reason": why}
print(json.dumps(out), flush=True)
if ok:
    out["ler"] = [ler(Zi, p, m) for p in (1e-3, 2e-3) for m in ("Z", "X")]
out["seconds"] = round(time.time() - t0, 1)
print(json.dumps(out))
json.dump({k: v for k, v in Zi.items() if k != "code"} | {"code": COLOR_17.name},
          open(os.path.join(here, "..", "data", "circuits", f"a3_{name}.json"), "w"), indent=1)
json.dump(out, open(os.path.join(here, "..", "results", f"a3_interleave_{name}.json"), "w"), indent=1)
