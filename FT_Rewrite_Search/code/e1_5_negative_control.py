"""E1.5: negative control. Does the checker catch changes that are NOT fault-equivalent?

Mutants of the baseline Shor gadget, in two groups:
  (a) "same ideal map, broken fault behaviour" (the dangerous kind): the checker must reject
      them through the fault checks, not the ideal check.
      - remove the verification flag (both flag CNOTs and the flag qubit)
      - move a data CNOT to another cat qubit (cat qubits share their Z value, so the ideal map
        is unchanged, but one cat qubit now touches two data qubits)
  (b) random edits (delete a CNOT, delete an H, retarget a CNOT): most change the ideal map.
Every mutant is classified by which check rejects it; mutants that pass are listed for review
(a pass is not automatically a miss: some edits are genuinely harmless).
Run:  C:/Users/Administrator/venvs/qml/Scripts/python.exe e1_5_negative_control.py
"""
import json
import os
import random

from ftrs.ftlib import BASELINE, fe_compare, ft1_report
from ftrs.verify import ideal_ok

random.seed(2026)
root = BASELINE


def classify(g):
    if not ideal_ok(g):
        return "rejected: ideal"
    a, b = fe_compare(root, g)
    if a or b:
        return "rejected: fault equivalence"
    if not ft1_report(g)["FT1"]:
        return "rejected: FT_1"
    return "PASSED"


mutants = []
# (a1) remove the flag
mutants.append(("a", "remove flag q11", {**root, "ops": [op for op in root["ops"] if 11 not in op[1:]], "flags": []}))
# (a2) move each data CNOT to each other cat qubit
cat = [7, 8, 9, 10]
for k, op in enumerate(root["ops"]):
    if op[0] == "CX" and op[2] < 7:
        for c in cat:
            if c != op[1]:
                ops = list(root["ops"])
                ops[k] = ("CX", c, op[2])
                mutants.append(("a", f"data CNOT {op} -> control q{c}", {**root, "ops": ops}))
# (b) random edits
cx_idx = [k for k, op in enumerate(root["ops"]) if op[0] == "CX"]
h_idx = [k for k, op in enumerate(root["ops"]) if op[0] == "H"]
for k in cx_idx:
    mutants.append(("b", f"delete {root['ops'][k]}", {**root, "ops": [o for i, o in enumerate(root["ops"]) if i != k]}))
for k in h_idx:
    mutants.append(("b", f"delete {root['ops'][k]}", {**root, "ops": [o for i, o in enumerate(root["ops"]) if i != k]}))
for _ in range(20):
    k = random.choice(cx_idx)
    c, t = root["ops"][k][1:]
    new_t = random.choice([q for q in range(12) if q not in (c, t)])
    ops = list(root["ops"])
    ops[k] = ("CX", c, new_t)
    mutants.append(("b", f"retarget {root['ops'][k]} -> CX{c}->{new_t}", {**root, "ops": ops}))

counts, passed = {}, []
for grp, desc, g in mutants:
    verdict = classify(g)
    counts.setdefault(grp, {}).setdefault(verdict, 0)
    counts[grp][verdict] += 1
    if verdict == "PASSED":
        passed.append((grp, desc))

for grp, name in (("a", "same ideal map, broken fault behaviour"), ("b", "random edits")):
    n = sum(counts.get(grp, {}).values())
    print(f"group {grp} ({name}): {n} mutants -> {counts.get(grp, {})}")
print("mutants that passed every check:")
for grp, desc in passed:
    print(f"  [{grp}] {desc}")
if not passed:
    print("  none")

out = os.path.join(os.path.dirname(__file__), "..", "results", "e1_5_negative_control.json")
json.dump({"counts": counts, "passed": passed}, open(out, "w"), indent=1)
print("saved", os.path.normpath(out))
