"""E1.2 (with E1.1 audit): reproduce Jeon et al.'s search on the Shor-style X0X1X2X3 gadget.

Target (Jeon Sec. V B 2, Table I): 7 reachable circuits per gadget; best 3 ancillas, 7 CNOTs,
CNOT depth 6. Every rewrite is checked (ideal, single-fault equivalence to the baseline, FT_1).
Writes results/e1_2_shor_search.json (or e1_2_shor_search_t2.json with checker order t = 2).
Run:  python e1_2_shor_search.py [t] [fusion|jeon]   (checker order, default 1; commutation mode, default fusion)
"""
import json
import os
import sys
import time

from ftrs.ftlib import BASELINE, REDUCED, resources
from ftrs.search import bfs, canonical, summary_rows
from ftrs.verify import check
from ftrs import rules

T = int(sys.argv[1]) if len(sys.argv) > 1 else 1
MODE = sys.argv[2] if len(sys.argv) > 2 else "fusion"  # restricted commutation: "fusion", "jeon" or "jeon_noflag"
rules.COMMUTATION_MODE = MODE
print("checker order t =", T, "| commutation mode", MODE)

t0 = time.time()
res = {}
for dedup in (True,):  # the no-dedup tree OOMs at 32G under the new rules (Fir job 1158940_1)
    r = bfs(BASELINE, dedup=dedup, check=lambda root, g: check(root, g, T))
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

out = os.path.join(os.path.dirname(__file__), "..", "results", "e1_2_shor_search" + ("_t2" if T == 2 else "") + ("_" + MODE if MODE != "fusion" else "") + ".json")
json.dump({
    "bfs_nodes_dedup": r["nodes"],
    "distinct_circuits": len(r["accepted"]), "rejected": reasons,
    "circuits": [{k: v for k, v in row.items() if k != "gadget"} | {"ops": row["gadget"]["ops"],
                 "syndrome": row["gadget"]["syndrome"], "flags": row["gadget"]["flags"]} for row in rows],
    "best_equals_fig1b": same_as_paper,
}, open(out, "w"), indent=1)
print("saved", os.path.normpath(out))
