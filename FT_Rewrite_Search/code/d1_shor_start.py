"""D1: rebuild the Shor-type [[17,1,5]] start so that the 5-round memory circuit has distance 5.

Why: the first start (e2_1_shor_start.py, single-round FT_2 per gadget) has circuit distance 4.
Its minimum logical error: one X fault early in the weight-8 cat tree leaves a weight-4 data
error that fires two flags; two flag measurement errors hide it, and one more data error makes
a logical (4 faults). Distance 5 needs, for every set of s faults that fires f flags and leaves
residual weight w (mod S) on the data, s + f >= w: then at least w - s - f more faults are
needed to hide it, and 5 - w more to complete a logical. Checked here for s <= 3 per gadget
(ftlib.faults_upto), flags chosen greedily (most violations removed first, then pruning), and
then confirmed on the whole round by stim's undetectable-logical-error search (Z and X memory).
Writes data/circuits/shor_17_1_5_d5.json and results/d1_shor_start.json.
Run:  python d1_shor_start.py
"""
import itertools
import json
import os
import sys
import time

from ftrs.codes import COLOR_17, COLOR_17_SUPPORTS
from ftrs.ftlib import dual, faults_upto, mask, resources
from ftrs.memgen import build
from ftrs.multi import derive

here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "e2_1_shor_start.py"), encoding="utf-8").read()
ns = {"__file__": os.path.join(here, "e2_1_shor_start.py")}
argv, sys.argv = sys.argv, ["x", "jeon"]
exec(src[:src.index("def ft2(")], ns)
sys.argv = argv
gadget = ns["gadget"]
t_start = time.time()


def violations(g, smax=3):
    fl, dm = mask(g["flags"]), COLOR_17.data_mask
    bad = 0
    for s, _, x, z in faults_upto(g["ops"], smax):
        f = bin(x & fl).count("1")
        if s + f < 4 and COLOR_17.weight_mod(x & dm, z & dm) > s + f:
            bad += 1
    return bad


def choose(support):
    allp = list(itertools.combinations(range(len(support)), 2))
    v = lambda ps: violations(gadget(support, ps, 17))
    pairs, cur = [], v([])
    while cur:
        scores = {p: v(pairs + [p]) for p in allp if p not in pairs}
        best = min(scores, key=scores.get)
        pairs.append(best)
        cur = scores[best]
        print(f"  weight {len(support)}: +{best} -> {cur} violations ({time.time() - t_start:.0f} s)", flush=True)
        if len(pairs) > 12:
            raise RuntimeError("no flag set found")
    for p in pairs[:]:
        if not v([q for q in pairs if q != p]):
            pairs.remove(p)
    return pairs


per = {}
for sup in sorted({tuple(s) for s in COLOR_17_SUPPORTS}, key=len):
    if len(sup) not in per:
        per[len(sup)] = choose(sup)
        print(f"weight {len(sup)}: verification pairs {per[len(sup)]}", flush=True)

ops, syn, flags, base = [], [], [], 17
for sup in COLOR_17_SUPPORTS:
    g = gadget(sup, per[len(sup)], base)
    ops += g["ops"]
    syn += g["syndrome"]
    flags += g["flags"]
    base += len(g["syndrome"]) + len(g["flags"])
Xh = {"name": "Shor-type X half, distance-5 flags", "kind": "X", "targets": [tuple(s) for s in COLOR_17_SUPPORTS],
      "syndrome": syn, "flags": flags, "ops": ops, "code": COLOR_17}
Zh = {**dual(Xh), "targets": [tuple(s) for s in COLOR_17.z_supports]}
assert derive(Xh)[1] and derive(Zh)[1]
dist = {}
for mem in ("Z", "X"):
    c, _ = build(COLOR_17, [Zh, Xh], rounds=5, memory=mem, p=1e-3)
    err = c.search_for_undetectable_logical_errors(dont_explore_detection_event_sets_with_size_above=4,
          dont_explore_edges_with_degree_above=4, dont_explore_edges_increasing_symptom_degree=False,
          canonicalize_circuit_errors=True)
    dist[mem] = len(err)
out = {"verification_pairs": {str(k): v for k, v in per.items()}, "x_half": resources(ops),
       "circuit_distance": dist, "seconds": round(time.time() - t_start, 1)}
print(json.dumps(out))
json.dump({k: v for k, v in Xh.items() if k != "code"} | {"code": COLOR_17.name, "verification_pairs": per},
          open(os.path.join(here, "..", "data", "circuits", "shor_17_1_5_d5.json"), "w"), indent=1)
json.dump(out, open(os.path.join(here, "..", "results", "d1_shor_start.json"), "w"), indent=1)
