"""E2.1 (AUTOLOOP C2): single-step reachable set of the [[17,1,5]] Liou Fig. 20 circuit.

On this circuit no Bell pair and no target swap exist, so composite reductions (which remove
an ancilla) do not apply; the only rewrites are restricted commutations. This script checks
every one-rewrite child with the t = 2 checker and counts them in the four ways of A3:
ASAP canonical form, structure + cost classes, cost classes, non-dominated cost vectors
(cost = ancillas, CNOTs, CNOT depth, total depth).
Writes results/e2_1_single_step_<mode>.json.
Run:  python e2_1_single_step.py [jeon|fusion]
"""
import json
import os
import sys
import time
from collections import Counter

from ftrs import rules
from ftrs.codes import COLOR_17
from ftrs.ftlib import resources
from ftrs.multi import check_factory
from ftrs.rules import bell_pairs, restricted_commutation, target_swap
from ftrs.search import canonical

MODE = sys.argv[1] if len(sys.argv) > 1 else "fusion"
rules.COMMUTATION_MODE = MODE
here = os.path.dirname(os.path.abspath(__file__))
L = json.load(open(os.path.join(here, "..", "data", "circuits", "liou_17_1_5.json")))
L = {**L, "ops": [tuple(o) for o in L["ops"]], "targets": [tuple(t) for t in L["targets"]], "code": COLOR_17}

src = open(os.path.join(here, "e1_3_recount.py"), encoding="utf-8").read()
ns = {}
exec(src[src.index("def ancilla_labels"):src.index("keys = {")].replace("DATA", "17"), ns)


def total_depth(ops):
    last, d = {}, 0
    for op in ops:
        if op[0] in ("H", "CX"):
            k = max(last.get(q, 0) for q in op[1:]) + 1
            for q in op[1:]:
                last[q] = k
            d = max(d, k)
    return d


def cost(g):
    r = resources(g["ops"])
    return (r["ancillas"], r["cnots"], r["cnot_depth"], total_depth(g["ops"]))


t0 = time.time()
check = check_factory(L, 2)
children = restricted_commutation(L)
print(f"mode {MODE}: Bell pairs {len(bell_pairs(L))}, target swaps {len(target_swap(L))}, "
      f"restricted commutations {len(children)}; root cost {cost(L)}", flush=True)
accepted, reasons = {}, Counter()
for child, desc in children:
    key = canonical(child)
    if key in accepted:
        continue
    ok, why = check(L, child)
    if ok:
        accepted[key] = child
    else:
        reasons[why] += 1
costs = Counter(cost(g) for g in accepted.values())
struct = {(ns["cx_commute_key"](g["ops"]), cost(g)) for g in accepted.values()}
nd = [c for c in costs if not any(d != c and all(x <= y for x, y in zip(d, c)) for d in costs)]
out = {"mode": MODE, "children": len(children), "accepted_asap": len(accepted), "rejected": dict(reasons),
       "structure_cost_classes": len(struct), "cost_classes": {str(k): v for k, v in costs.items()},
       "non_dominated": sorted(nd), "root_cost": cost(L), "seconds": round(time.time() - t0, 1)}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(here, "..", "results", f"e2_1_single_step_{MODE}.json"), "w"), indent=1)
