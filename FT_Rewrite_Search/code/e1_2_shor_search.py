"""E1.2 (with E1.1 audit): reproduce Jeon et al.'s search on the Shor-style X0X1X2X3 gadget.

Target (Jeon Sec. V B 2, Table I): 7 reachable circuits per gadget; best 3 ancillas, 7 CNOTs,
CNOT depth 6. Every rewrite is checked (ideal, single-fault equivalence to the baseline, FT_1).
Writes results/e1_2_shor_search.json.
Run:  C:/Users/Administrator/venvs/qml/Scripts/python.exe e1_2_shor_search.py
"""
import json
import os
import time

from ftrs.ftlib import BASELINE, REDUCED, resources
from ftrs.search import bfs, canonical, summary_rows

t0 = time.time()
res = {}
for dedup in (False, True):
    r = bfs(BASELINE, dedup=dedup)
    res[dedup] = r
    print(f"dedup={dedup}: BFS nodes {r['nodes']}, distinct accepted circuits {len(r['accepted'])}, "
          f"rejected rewrites {sum(r['rejected'].values())}, terminal nodes {len(r['terminals'])}")

r = res[True]
reasons = dict(r["rejected"])
print("rejected rewrites by reason:", reasons or "none")

rows = summary_rows(r)
print("\ndistinct circuits reachable from the baseline (sorted by ancillas, CNOT depth, CNOTs):")
for k, row in enumerate(rows):
    print(f"  #{k}: ancillas {row['ancillas']}, CNOTs {row['cnots']}, CNOT depth {row['cnot_depth']}, "
          f"Bell reductions {row['steps']}; path: {' | '.join(row['path']) or '(baseline)'}")

best = rows[0]
same_as_paper = canonical(best["gadget"]) == canonical(REDUCED)
print("\nbest circuit:", resources(best["gadget"]["ops"]))
print("identical (canonical form) to Jeon Fig. 1b:", same_as_paper)
print("Jeon Fig. 1b among the reachable circuits:", canonical(REDUCED) in r["accepted"])
print("ops of the best circuit:", best["gadget"]["ops"])
print(f"\nwall time {time.time() - t0:.1f} s")

out = os.path.join(os.path.dirname(__file__), "..", "results", "e1_2_shor_search.json")
json.dump({
    "bfs_nodes_no_dedup": res[False]["nodes"], "bfs_nodes_dedup": r["nodes"],
    "distinct_circuits": len(r["accepted"]), "rejected": reasons,
    "circuits": [{k: v for k, v in row.items() if k != "gadget"} | {"ops": row["gadget"]["ops"],
                 "syndrome": row["gadget"]["syndrome"], "flags": row["gadget"]["flags"]} for row in rows],
    "best_equals_fig1b": same_as_paper,
}, open(out, "w"), indent=1)
print("saved", os.path.normpath(out))
