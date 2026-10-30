"""E2.1 part 3: full reachable sets of the per-generator Shor-type gadgets of [[17,1,5]].

The round built by e2_1_shor_start.py is eight gadgets on disjoint ancillas; FT_2 and the
unflagged effect sets compose across them (a fault touches one gadget), so the best round is
the best gadget per generator. Here: BFS (t = 2 checker, verify.check) from the weight-4 gadget
on g9 = X1X2X3X4 and the weight-8 gadget on g16, in the given commutation mode.
Writes results/e2_1_per_gadget_w<4|8>_<mode>.json.
Run:  python e2_1_per_gadget.py <4|8> [jeon|fusion] [max BFS nodes]
"""
import json
import os
import sys
import time
from collections import Counter

from ftrs import rules
from ftrs.codes import COLOR_17_SUPPORTS
from ftrs.ftlib import resources
from ftrs.search import bfs
from ftrs.verify import check

W = int(sys.argv[1])
MODE = sys.argv[2] if len(sys.argv) > 2 else "fusion"
rules.COMMUTATION_MODE = MODE
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "e2_1_shor_start.py"), encoding="utf-8").read()
ns = {"__file__": os.path.join(here, "e2_1_shor_start.py")}
exec(src[:src.index("def ft2(")], ns)
rules.COMMUTATION_MODE = MODE  # the exec above resets it from its own sys.argv
START = sys.argv[4] if len(sys.argv) > 4 else "old"
pairs = ({4: [(0, 1)], 8: [(2, 3), (0, 1), (0, 3), (1, 2)]} if START == "d5"  # d1_shor_start (job 1239453)
         else {4: [(0, 1)], 8: [(2, 3), (0, 1), (0, 2)]})[W]  # e2_1_shor_start (job 1230344): distance 4
support = next(s for s in COLOR_17_SUPPORTS if len(s) == W)
root = ns["gadget"](support, pairs, 17)

t0 = time.time()
r = bfs(root, dedup=True, max_nodes=int(sys.argv[3]) if len(sys.argv) > 3 else 100000, check=lambda a, g: check(a, g, 2), L=1)
rows = [(g, path) for g, path in r["accepted"].values()]
cost = lambda g: tuple(resources(g["ops"]).values())
costs = Counter(cost(g) for g, _ in rows)
nd = [c for c in costs if not any(d != c and all(x <= y for x, y in zip(d, c)) for d in costs)]
best = min(rows, key=lambda gp: cost(gp[0]))
out = {"weight": W, "mode": MODE, "support": support, "root": cost(root), "bfs_nodes": r["nodes"],
       "accepted_asap": len(rows), "generated": r["generated"], "rejected": dict(r["rejected"]),
       "cost_classes": {str(k): v for k, v in sorted(costs.items())}, "non_dominated": sorted(nd),
       "best": cost(best[0]), "best_ops": best[0]["ops"], "best_path": best[1],
       "seconds": round(time.time() - t0, 1)}
if START == "d5":  # post-filter with the distance-5 condition of d1_shor_start (s + flags >= weight, s <= 3)
    dsrc = open(os.path.join(here, "d1_shor_start.py"), encoding="utf-8").read()
    dns = {"__file__": os.path.join(here, "d1_shor_start.py")}
    argv0, sys.argv = sys.argv, ["x"]
    exec(dsrc[:dsrc.index("def choose(")], dns)
    sys.argv = argv0
    passing, tried = [], 0
    for g, path in sorted(rows, key=lambda gp: cost(gp[0])):
        tried += 1
        if dns["violations"](g) == 0:
            passing.append((cost(g), g["ops"], path))
            if len(passing) >= 5:
                break
        if tried >= 400:
            break
    out.update({"d5_tried": tried, "d5_passing": [c for c, _, _ in passing],
                "best_d5": passing[0][0] if passing else None, "best_d5_ops": passing[0][1] if passing else None,
                "best_d5_path": passing[0][2] if passing else None})
print(json.dumps({k: v for k, v in out.items() if not k.endswith("_ops")}, indent=1))
json.dump(out, open(os.path.join(here, "..", "results", f"e2_1_per_gadget_w{W}_{MODE}" + ("_d5" if START == "d5" else "") + ".json"), "w"), indent=1)
