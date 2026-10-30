"""E2.7 (early, AUTOLOOP): LER of full [[17,1,5]] rounds built from rewrite-reduced Shor-type
gadgets, against the Liou baseline under the same setting B.

A variant names, per generator weight, which per-gadget circuit to use (the root or the best
circuit of an e2_1_per_gadget run). Each X generator gets that circuit with its data qubits
mapped position by position onto the generator's support and fresh ancillas; the eight gadgets
form the X half, its Hadamard dual the Z half. Checks before sampling: noiseless correctness,
single-round FT_2 of both halves, circuit distance 5 of the memory circuit, and DemLookup
correct on every single DEM error. Then the same memory experiment as e2_0_ler_b.py.
Writes results/e2_7_round_ler_<variant>_p<p>_<memory>.json.
Run:  python e2_7_round_ler.py <variant> <p> [Z|X] [max_errors]
  variant: <w4 source>-<w8 source>, each "root", "jeon" or "fusion" (best of that BFS);
  with suffix "_d5" the distance-5 start of D1 ("root_d5", "jeon_d5", "fusion_d5")
"""
import json
import os
import sys
import time

from ftrs.codes import COLOR_17, COLOR_17_SUPPORTS
from ftrs.decoders import DemLookup
from ftrs.ftlib import dual, effect_signature, resources
from ftrs.mem import wilson
from ftrs.memgen import build
from ftrs.multi import derive

variant, p = sys.argv[1], float(sys.argv[2])
memory = sys.argv[3] if len(sys.argv) > 3 else "Z"
max_errors = int(sys.argv[4]) if len(sys.argv) > 4 else 400
here = os.path.dirname(os.path.abspath(__file__))
res = os.path.join(here, "..", "results")
src = dict(zip((4, 8), variant.split("-")))

base = {}
for w in (4, 8):
    d = json.load(open(os.path.join(res, f"e2_1_per_gadget_w{w}_{'jeon' if src[w] in ('root', 'root_d5') else src[w]}.json")))
    if src[w] == "root_d5":
        sh = json.load(open(os.path.join(here, "..", "data", "circuits", "shor_17_1_5_d5.json")))
        ns = {"__file__": os.path.join(here, "e2_1_shor_start.py")}
        code_src = open(os.path.join(here, "e2_1_shor_start.py"), encoding="utf-8").read()
        sys_argv, sys.argv = sys.argv, ["x", "jeon"]
        exec(code_src[:code_src.index("def ft2(")], ns)
        sys.argv = sys_argv
        g = ns["gadget"](tuple(d["support"]), [tuple(x) for x in sh["verification_pairs"][str(w)]], 17)
        base[w] = (tuple(d["support"]), g["ops"])
        continue
    if src[w] == "root":
        sh = json.load(open(os.path.join(here, "..", "data", "circuits", "shor_17_1_5.json")))
        k = next(i for i, s in enumerate(sh["targets"]) if len(s) == w)
        # the root gadget for this weight is the k-th block of the Shor round: rebuild it instead
        ns = {"__file__": os.path.join(here, "e2_1_shor_start.py")}
        code_src = open(os.path.join(here, "e2_1_shor_start.py"), encoding="utf-8").read()
        sys_argv, sys.argv = sys.argv, ["x", "jeon"]
        exec(code_src[:code_src.index("def ft2(")], ns)
        sys.argv = sys_argv
        g = ns["gadget"](tuple(d["support"]), [tuple(x) for x in sh["verification_pairs"][str(w)]], 17)
        base[w] = (tuple(d["support"]), g["ops"])
    elif src[w].endswith("_d5"):  # distance-5 start (D1): cheapest circuit passing the d1 condition
        base[w] = (tuple(d["support"]), [tuple(o) for o in d["best_d5_ops"]])
    else:
        base[w] = (tuple(d["support"]), [tuple(o) for o in d["best_ops"]])


def place(w, support, nxt):
    sup0, ops = base[w]
    m = dict(zip(sup0, support))
    anc = sorted({q for op in ops for q in op[1:] if q >= 17})
    m.update({a: nxt + i for i, a in enumerate(anc)})
    new = [tuple([op[0]] + [m[q] for q in op[1:]]) for op in ops]
    meas = [op[1] for op in new if op[0] == "M"]
    return new, meas, nxt + len(anc)


ops, meas, nxt = [], [], 17
for sup in COLOR_17_SUPPORTS:
    o, m, nxt = place(len(sup), sup, nxt)
    ops += o
    meas += m
Xh = {"name": f"X half {variant}", "kind": "X", "targets": [tuple(s) for s in COLOR_17_SUPPORTS],
      "syndrome": meas, "flags": [], "ops": ops, "code": COLOR_17}
par, ok = derive(Xh, shots=2000)
used = {q for pr in par for q in pr}
Xh["flags"] = [q for q in meas if q not in used]  # measurements outside every readout are flags
par, ok = derive(Xh, shots=2000)
assert ok, "X half fails its noiseless check"
Zh = {**dual(Xh), "targets": [tuple(s) for s in COLOR_17.z_supports]}
pz, okz = derive(Zh, shots=2000)
ft = (effect_signature(Xh, 2, par)[1], effect_signature(Zh, 2, pz)[1])
print(f"{variant}: X half {resources(ops)}, flags {len(Xh['flags'])}, ideal {ok}/{okz}, single-round FT_2 violations {ft}", flush=True)
assert okz and ft == (0, 0)

t0 = time.time()
circ, info = build(COLOR_17, [Zh, Xh], rounds=5, memory=memory, p=p)
err = circ.search_for_undetectable_logical_errors(dont_explore_detection_event_sets_with_size_above=4,
      dont_explore_edges_with_degree_above=4, dont_explore_edges_increasing_symptom_degree=False,
      canonicalize_circuit_errors=True)
dec = DemLookup(circ)
H, Lm = dec.fallback.H, dec.L
bad1 = sum((dec.decode_one(H[:, j]) != Lm[:, j]).any() for j in range(H.shape[1]))
print(f"circuit distance {len(err)}, DemLookup wrong on {bad1}/{H.shape[1]} single errors", flush=True)
sampler = circ.compile_detector_sampler(seed=11)
errors = shots = 0
while errors < max_errors and shots < 2 * 10 ** 8:
    dets, obs = sampler.sample(20000, separate_observables=True)
    errors += int((dec.decode_batch(dets) != obs).any(axis=1).sum())
    shots += len(dets)
lo, hi = wilson(errors, shots)
out = {"variant": variant, "p": p, "memory": memory, "x_half": resources(ops), "flags": len(Xh["flags"]),
       "circuit_distance": len(err), "dem_single_wrong": int(bad1), "shots": shots, "errors": errors,
       "ler": errors / shots, "ci95": [lo, hi], "seconds": round(time.time() - t0, 1)}
print(json.dumps(out))
json.dump(out, open(os.path.join(res, f"e2_7_round_ler_{variant}_p{p:g}_{memory}.json"), "w"), indent=1)
