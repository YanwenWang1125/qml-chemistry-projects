"""Follow-up item 3: does a local condition predict which flag merges keep fault tolerance?

Condition (notes/flag_merge_theory.md): after merging, every unflagged-or-flagged set of s <= 2
faults must satisfy s + (number of flag qubits it flips) >= residual data weight mod S. A fault
set that flipped both merged flags flips the shared one twice, i.e. loses 2 from that count.
Vectorised with ftlib.effect_arrays-style enumeration over the Z half and its dual X half.
Compared with the full check used in A-2 (noiseless, single-round FT_2, memory distance 5) on
every first-step merge of the given start. Writes results/a4_merge_predict_<tag>.json.
Run:  python a4_merge_predict.py <start json> <tag>
"""
import itertools
import json
import os
import sys
import time

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
path, tag = sys.argv[1], sys.argv[2]
src = open(os.path.join(here, "a2_flag_start.py"), encoding="utf-8").read()
ns = {"__file__": os.path.join(here, "a2_flag_start.py")}
argv, sys.argv = sys.argv, ["x", "1"]
exec(src[:src.index("LIOU_G8 = (13, 9, 14, 5, 2, 6, 3, 10)\nper = {}")], ns)
exec(src[src.index("def halves("):src.index("ok, why = accepted(Zh)")], ns)
sys.argv = argv
from ftrs.ftlib import located_faults, propagate  # noqa: E402

COLOR_17, accepted, merge, halves = ns["COLOR_17"], ns["accepted"], ns["merge"], ns["halves"]


def local_ok(g):
    for h in halves(g):
        ops, fl = h["ops"], h["flags"]
        loc, X, Z = [], [], []
        for l, _, start, x, z in located_faults(ops):
            x, z = propagate(x, z, ops[start:])
            loc.append(l), X.append(x), Z.append(z)
        loc, X, Z = np.array(loc), np.array(X, dtype=object), np.array(Z, dtype=object)
        i, j = np.triu_indices(len(X), 1)
        keep = loc[i] != loc[j]
        s = np.concatenate([np.ones(len(X), int), np.full(keep.sum(), 2)])
        xs = np.concatenate([X, X[i[keep]] ^ X[j[keep]]])
        zs = np.concatenate([Z, Z[i[keep]] ^ Z[j[keep]]])
        fired = sum(((xs >> f) & 1).astype(int) for f in fl)
        dm = COLOR_17.data_mask
        uniq = {}
        for k in range(len(xs)):
            if s[k] + fired[k] >= 4:
                continue
            key = (int(xs[k]) & dm, int(zs[k]) & dm)
            if key not in uniq:
                uniq[key] = COLOR_17.weight_mod(*key)
            if uniq[key] > s[k] + fired[k]:
                return False
    return True


d = json.load(open(path))
Z0 = {**d, "ops": [tuple(o) for o in d["ops"]], "targets": [tuple(t) for t in d["targets"]], "code": COLOR_17}
t0 = time.time()
rows = []
for fa, fb in itertools.combinations(Z0["flags"], 2):
    g = merge(Z0, fa, fb)
    pred = local_ok(g)
    full, why = accepted(g)
    rows.append({"merge": [fa, fb], "local": pred, "full": full, "reason": why})
from collections import Counter  # noqa: E402
conf = Counter((r["local"], r["full"]) for r in rows)
out = {"tag": tag, "pairs": len(rows), "confusion_local_vs_full": {f"local {a} / full {b}": n for (a, b), n in conf.items()},
       "seconds": round(time.time() - t0, 1), "rows": rows}
print(json.dumps({k: v for k, v in out.items() if k != "rows"}))
json.dump(out, open(os.path.join(here, "..", "results", f"a4_merge_predict_{tag}.json"), "w"), indent=1)
