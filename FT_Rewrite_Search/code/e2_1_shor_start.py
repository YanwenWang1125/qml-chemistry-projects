"""E2.1 part 2: a Shor-type starting circuit for [[17,1,5]] (plan Sec. 5 table 10) and its
single-step reachable set.

Built on the X-type generators (self-dual code; Liou measure Z and note X is the Hadamard
image), so the rewrite rules run in Jeon's setting: cat qubits in |+> control the data.
Per generator of weight w: cat qubits c_0..c_{w-1} (R, H on c_0, a CNOT tree from c_0),
verification flags (each R f, CX c_a->f, CX c_b->f, measured in Z: checks Z_a Z_b), then
CX c_i -> data_i, H and M on every cat qubit. Flags are added greedily (fewest pairs) until
the single-generator gadget is FT_2 on COLOR_17 (ftlib.ft_report, t = 2). The eight gadgets
are then run one after another as one multi-bit gadget (kind "X") and checked again with the
t = 2 multi checker; then every composite reduction one step from it is checked.
Writes data/circuits/shor_17_1_5.json and results/e2_1_shor_start_<mode>.json.
Run:  python e2_1_shor_start.py [jeon|fusion]
"""
import itertools
import json
import os
import sys
import time
from collections import Counter

from ftrs import rules
from ftrs.codes import COLOR_17, COLOR_17_SUPPORTS
from ftrs.ftlib import effect_signature, ft_report, resources
from ftrs.multi import check_factory, derive
from ftrs.search import canonical, composite_reductions
from ftrs.verify import ideal_ok

MODE = sys.argv[1] if len(sys.argv) > 1 else "fusion"
rules.COMMUTATION_MODE = MODE
here = os.path.dirname(os.path.abspath(__file__))
t_start = time.time()


def tree(w):
    """CNOT tree on cat positions: each new qubit copies from one already in the cat."""
    edges, have = [], [0]
    while len(have) < w:
        new = []
        for a in have:
            if len(have) + len(new) < w:
                b = len(have) + len(new)
                edges.append((a, b))
                new.append(b)
        have += new
    return edges


def gadget(support, pairs, base):
    w = len(support)
    cat = list(range(base, base + w))
    fl = list(range(base + w, base + w + len(pairs)))
    ops = [("R", q) for q in cat + fl] + [("H", cat[0])]
    ops += [("CX", cat[a], cat[b]) for a, b in tree(w)]
    for f, (a, b) in zip(fl, pairs):
        ops += [("CX", cat[a], f), ("CX", cat[b], f)]
    ops += [("CX", c, d) for c, d in zip(cat, support)]
    ops += [("H", c) for c in cat] + [("M", q) for q in cat + fl]
    return {"ops": ops, "syndrome": cat, "flags": fl, "kind": "X", "target": tuple(support),
            "code": COLOR_17}


def ft2(g):
    return ideal_ok(g) and ft_report(g, 2)["FT"]


def fewest_pairs(support):
    """Exhaustive for up to 2 flags, then greedy: add the pair that leaves the fewest FT_2
    violations, then drop any pair that is not needed. Not guaranteed minimal for weight 8."""
    w = len(support)
    allp = list(itertools.combinations(range(w), 2))
    for k in range(0, 3):
        for pairs in itertools.combinations(allp, k):
            if ft2(gadget(support, list(pairs), 17)):
                return list(pairs)
    viol = lambda ps: len(ft_report(gadget(support, ps, 17), 2)["violations"])
    pairs = []
    while viol(pairs):
        pairs.append(min((p for p in allp if p not in pairs), key=lambda p: viol(pairs + [p])))
        if len(pairs) > 10:
            raise RuntimeError(f"no FT_2 verification with <= 10 flags for {support}")
    for p in pairs[:]:
        if not viol([q for q in pairs if q != p]):
            pairs.remove(p)
    assert ft2(gadget(support, pairs, 17))
    return pairs


per = {}
for sup in sorted({tuple(s) for s in COLOR_17_SUPPORTS}, key=len):
    if len(sup) not in per:
        per[len(sup)] = fewest_pairs(sup)
        print(f"weight {len(sup)}: FT_2 with verification pairs {per[len(sup)]} "
              f"({time.time() - t_start:.0f} s)", flush=True)

ops, syndrome, flags, base = [], [], [], 17
for sup in COLOR_17_SUPPORTS:
    g = gadget(sup, per[len(sup)], base)
    ops += g["ops"]
    syndrome += g["syndrome"]
    flags += g["flags"]
    base += len(g["syndrome"]) + len(g["flags"])
ROOT = {"name": "Shor-type X-stabilizer round for [[17,1,5]]", "kind": "X", "targets": COLOR_17_SUPPORTS,
        "syndrome": syndrome, "flags": flags, "ops": ops, "code": COLOR_17}
parities, ok = derive(ROOT, shots=2000)
sig, bad = effect_signature(ROOT, 2, parities)
print(f"root: {resources(ops)}, ideal {ok}, single-round FT_2 violations {bad}", flush=True)
assert ok and not bad
os.makedirs(os.path.join(here, "..", "data", "circuits"), exist_ok=True)
json.dump({k: v for k, v in ROOT.items() if k != "code"} | {"code": COLOR_17.name, "verification_pairs": per},
          open(os.path.join(here, "..", "data", "circuits", "shor_17_1_5.json"), "w"), indent=1)

check = check_factory(ROOT, 2)
accepted, reasons, n = {}, Counter(), 0
for child, steps in composite_reductions(ROOT, 1):
    n += 1
    key = canonical(child)
    if key in accepted:
        continue
    good, why = check(ROOT, child)
    if good:
        accepted[key] = (child, steps)
    else:
        reasons[why] += 1
costs = Counter(tuple(resources(g["ops"]).values()) for g, _ in accepted.values())
nd = [c for c in costs if not any(d != c and all(x <= y for x, y in zip(d, c)) for d in costs)]
out = {"mode": MODE, "root_resources": resources(ops), "verification_pairs": {str(k): v for k, v in per.items()},
       "composite_reductions": n, "accepted_asap": len(accepted), "rejected": dict(reasons),
       "cost_classes": {str(k): v for k, v in costs.items()}, "non_dominated": sorted(nd),
       "seconds": round(time.time() - t_start, 1)}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(here, "..", "results", f"e2_1_shor_start_{MODE}.json"), "w"), indent=1)
