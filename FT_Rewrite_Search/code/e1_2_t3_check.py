"""A2: order-3 effect-set check of the E1.2 (t=2 checker) circuits against the Shor baseline.

The t=2 search still reaches 2 ancillas / 6 CNOTs, below Jeon's best 3/7/6. If those circuits
also agree with the baseline on all unflagged sets of <= 3 faults, they stay "candidates, FE
undecided"; if not, the widened restricted-commutation rule is the cause and must be tightened.
Checks every distinct circuit of results/e1_2_shor_search_t2.json plus Jeon Fig. 1b (control).
Writes results/e1_2_t3_check.json.
Run on Fir:  sbatch e1_2_t3_check.sh
"""
import json
import os
import time
from collections import Counter

from ftrs.ftlib import BASELINE, REDUCED, fe_compare_t, ft_report, resources

here = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(here, "..", "results", "e1_2_shor_search_t2.json")))["circuits"]
gads = [("Jeon Fig. 1b", REDUCED)] + [
    (f"#{k}", {"ops": [tuple(o) for o in r["ops"]], "syndrome": r["syndrome"], "flags": r["flags"]})
    for k, r in enumerate(rows)]

out, summary = [], Counter()
for name, g in gads:
    t0 = time.time()
    a, b = fe_compare_t(BASELINE, g, 3)
    ft3 = ft_report(g, 3)["FT"]
    res = resources(g["ops"])
    key = (res["ancillas"], res["cnots"], res["cnot_depth"])
    ok = not (a or b)
    summary[(key, ok)] += 1
    out.append({"name": name, "resources": key, "t3_consistent": ok, "only_root": len(a),
                "only_child": len(b), "FT3": ft3, "seconds": round(time.time() - t0, 1)})
    print(f"{name}: {key} t=3 consistent {ok} ({len(a)}/{len(b)}), FT_3 {ft3}, {out[-1]['seconds']} s", flush=True)

print("\nby resources (anc, CNOT, CNOT depth) and t=3 consistency:")
for (key, ok), n in sorted(summary.items()):
    print(f"  {key} consistent={ok}: {n}")
json.dump(out, open(os.path.join(here, "..", "results", "e1_2_t3_check.json"), "w"), indent=1)
print("saved")
