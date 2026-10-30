"""E1.3: reproduce Jeon et al.'s search on the Steane-style gadget (Goto baseline, Fig. 3a).

Targets (Sec. V C, Table III): baseline 8 ancillas / 18 CNOT / 6 H / depth 8 / CNOT depth 7;
the search generates 391 candidates; best has 4 ancillas / 14 CNOT / 2 H / depth 8 / CNOT depth 7.
Run:  python e1_3_steane_search.py [t] [fusion|jeon]   (checker order, default 1; commutation mode, default fusion)
"""
import json
import os
import sys
import time

from ftrs.ftlib import resources
from ftrs.multi import GOTO, check_factory, derive, effects
from ftrs.search import bfs, summary_rows
from ftrs import rules


def total_depth(ops):
    last, d = {}, 0
    for op in ops:
        if op[0] in ("H", "CX"):
            k = max(last.get(q, 0) for q in op[1:]) + 1
            for q in op[1:]:
                last[q] = k
            d = max(d, k)
    return d


def describe(g):
    r = resources(g["ops"])
    return {**r, "H": sum(op[0] == "H" for op in g["ops"]), "depth": total_depth(g["ops"])}


t0 = time.time()
p, ok = derive(GOTO)
sig, bad = effects(GOTO, p)
print("baseline:", describe(GOTO), "| ideal ok:", ok, "| parities:", p, "| FT_1:", not bad)

L = 2  # E1.3-fix: App. C needs two enabling rewrites per step
T = int(sys.argv[1]) if len(sys.argv) > 1 else 1  # checker order; t = 2 since the E1.2 regression
MODE = sys.argv[2] if len(sys.argv) > 2 else "fusion"  # restricted commutation: "fusion", "jeon" or "jeon_noflag"
rules.COMMUTATION_MODE = MODE
print("checker order t =", T, "| commutation mode", MODE)
check = check_factory(GOTO, T)
res = {}
for dedup in (True,):  # without dedup the L=2 tree is exponential; report generated children instead
    r = bfs(GOTO, dedup=dedup, check=check, L=L)
    res[dedup] = r
    reasons = dict(r["rejected"])
    print(f"dedup={dedup}: BFS nodes {r['nodes']}, distinct accepted {len(r['accepted'])}, "
          f"children generated {r['generated']} (distinct checked {r['checked']}), "
          f"rejected rewrites {sum(r['rejected'].values())} {reasons}, terminals {len(r['terminals'])}")

rows = summary_rows(res[True])
for row in rows:
    row.update(describe(row["gadget"]))
rows.sort(key=lambda r: (r["ancillas"], r["cnots"], r["cnot_depth"], r["H"]))
print("\nbest 8 distinct circuits:")
for row in rows[:8]:
    print(f"  ancillas {row['ancillas']}, CNOT {row['cnots']}, H {row['H']}, depth {row['depth']}, "
          f"CNOT depth {row['cnot_depth']}, Bell reductions {row['steps']}")
hist = {}
for row in rows:
    hist[row["ancillas"]] = hist.get(row["ancillas"], 0) + 1
print("distinct circuits by ancilla count:", dict(sorted(hist.items())))
print("best ops:", rows[0]["gadget"]["ops"])
print(f"wall time {time.time() - t0:.1f} s")

out = os.path.join(os.path.dirname(__file__), "..", "results", "e1_3_steane_search" + ("_t2" if T == 2 else "") + ("_" + MODE if MODE != "fusion" else "") + ".json")
json.dump({"L": L, "children_generated": res[True]["generated"], "distinct_checked": res[True]["checked"],
           "rejected": dict(res[True]["rejected"]), "bfs_nodes_dedup": res[True]["nodes"],
           "distinct": len(res[True]["accepted"]),
           "circuits": [{k: v for k, v in row.items() if k != "gadget"} | {"ops": row["gadget"]["ops"],
                        "syndrome": row["gadget"]["syndrome"], "flags": row["gadget"]["flags"]} for row in rows]},
          open(out, "w"), indent=1)
print("saved", os.path.normpath(out))
