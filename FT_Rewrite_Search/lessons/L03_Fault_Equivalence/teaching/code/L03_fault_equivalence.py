"""L03 verification: a simplified, single-fault check of fault equivalence.

For each pair of gadgets with the same ideal behaviour, list the effects (data error modulo
the stabilizers, syndrome flip) that an undetected single fault produces in one gadget but
that no fault of weight <= 1 produces in the other. Both lists empty = the pair passes this
single-fault check. This is NOT the full definition (see the lesson): it uses one noise model,
compares only weight-1 faults, and summarises the faulty map by the data coset and syndrome bit.
Run from this folder:  C:/Users/Administrator/venvs/qml/Scripts/python.exe L03_fault_equivalence.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                "L01_FT_and_Shor_Gadget", "teaching", "code"))
from ftlib import (BASELINE, NO_FLAG, REDUCED, fe_compare, ft1_report,  # noqa: E402
                   ideal_syndrome_check, pauli_str, resources)


def swapped(gadget, op_a, op_b, name):
    ops = list(gadget["ops"])
    i, j = ops.index(op_a), ops.index(op_b)
    ops[i], ops[j] = ops[j], ops[i]
    return {**gadget, "ops": ops, "name": name}


# Two reorderings of CNOTs that share the control q7. Ideally both commute.
REORDER_FRESH = swapped(BASELINE, ("CX", 7, 8), ("CX", 7, 9),
                        "baseline, CX7->8 and CX7->9 swapped (q7 freshly prepared in |+>)")
REORDER_LATE = swapped(BASELINE, ("CX", 7, 11), ("CX", 7, 0),
                       "baseline, CX7->11 (flag) and CX7->0 (data) swapped")

pairs = [
    (BASELINE, NO_FLAG),
    (BASELINE, REORDER_FRESH),
    (BASELINE, REORDER_LATE),
    (BASELINE, REDUCED),
]

for a, b in pairs:
    print("=" * 72)
    print(f"A = {a['name']}\nB = {b['name']}")
    for g in (a, b):
        same = all(v[0] == ([0] if "+1" in k else [1]) for k, v in ideal_syndrome_check(g).items())
        print(f"  {g['name'][:45]:<45} ideal syndrome correct: {same}; FT_1: {ft1_report(g)['FT1']}; "
              f"{resources(g['ops'])}")
    only_a, only_b = fe_compare(a, b)
    for side, only in (("A", only_a), ("B", only_b)):
        if not only:
            print(f"  effects only in {side}: none")
        for ((dx, dz), sflip), fault in sorted(only.items(), key=lambda kv: str(kv[0])):
            print(f"  effect only in {side}: data {pauli_str(dx, dz, range(7)):<10} syndrome flip {sflip}"
                  f"   e.g. {fault}")
    print("  passes the single-fault check:", not only_a and not only_b)
