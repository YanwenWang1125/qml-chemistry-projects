"""L04 verification: Bell reduction on the smallest example, then the reduced Shor gadget.

Part 1. Measure X0 X1 on two data qubits in two ways:
  A: a Bell pair of ancillas (q2 in |+>, q3 in |0>, CX2->3), q2 controls d0, q3 controls d1,
     both ancillas measured in X, only the parity m2 xor m3 is used.
  B: one ancilla q2 in |+> controls d0 and d1, measured in X.
  Bell reduction says B is fault-equivalent to A with one ancilla and one CNOT fewer.
  We check: same ideal syndrome, the single-fault effect check of L03, and the fact that the
  two circuits nevertheless have different numbers of fault locations and error rates.
Part 2. Resource counts and FT_1 of the gadgets before and after Jeon et al.'s reduction.
Run from this folder:  C:/Users/Administrator/venvs/qml/Scripts/python.exe L04_bell_reduction.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                "L01_FT_and_Shor_Gadget", "teaching", "code"))
from ftlib import (BASELINE, REDUCED, fe_compare, ft1_report, measured_qubits,  # noqa: E402
                   pauli_str, resources, single_faults, stabilizer_group, to_stim)

N2 = 2
GROUP2 = stabilizer_group([(0, 1)], [])  # {I, X0X1}: the operator being measured

A = {"name": "A: Bell pair of ancillas",
     "ops": [("R", 2), ("R", 3), ("H", 2), ("CX", 2, 3), ("CX", 2, 0), ("CX", 3, 1),
             ("H", 2), ("H", 3), ("M", 2), ("M", 3)],
     "syndrome": [2, 3], "flags": []}
B = {"name": "B: one ancilla",
     "ops": [("R", 2), ("H", 2), ("CX", 2, 0), ("CX", 2, 1), ("H", 2), ("M", 2)],
     "syndrome": [2], "flags": []}

print("Part 1: Bell reduction on the smallest example")
for g in (A, B):
    meas = measured_qubits(g["ops"])
    seen = {}
    for label, z_on in (("X0X1=+1", ()), ("X0X1=-1", (0,))):
        s = to_stim(g["ops"], n_data=N2, data_init="plus", z_on=z_on).compile_sampler().sample(200)
        syn = s[:, [meas.index(q) for q in g["syndrome"]]].sum(axis=1) % 2
        seen[label] = sorted(set(syn.astype(int).tolist()))
    print(f"  {g['name']:<26} {resources(g['ops'])}, single-fault locations: "
          f"{sum(1 for _ in single_faults(g['ops']))}, noiseless syndrome {seen}")

only_a, only_b = fe_compare(A, B, group=GROUP2, n_data=N2)
print("  effects only in A:", {pauli_str(*k[0], range(N2)) + f" flip {k[1]}": v for k, v in only_a.items()} or "none")
print("  effects only in B:", {pauli_str(*k[0], range(N2)) + f" flip {k[1]}": v for k, v in only_b.items()} or "none")
print("  passes the single-fault check:", not only_a and not only_b)

# Same fault-equivalence class, different error rates: Monte Carlo of a wrong syndrome bit
# or a data X error (data in |0>^2, measured in Z, so only X errors on data are seen).
p, shots = 1e-2, 2_000_000
print(f"\n  Monte Carlo at p = {p}, {shots} shots, data in |0>^2:")
for g in (A, B):
    s = to_stim(g["ops"], n_data=N2, data_init="zero", p=p).compile_sampler(seed=7).sample(shots)
    n_g = len(measured_qubits(g["ops"]))
    data = s[:, n_g:n_g + N2]
    # modulo X0X1: an X error counts only if exactly one data qubit flipped
    x_err = data.sum(axis=1) == 1
    print(f"    {g['name']:<26} P(data X error mod X0X1) = {x_err.mean():.5f}")

print("\nPart 2: Jeon et al.'s Shor gadget before and after two Bell reductions")
for g in (BASELINE, REDUCED):
    rep = ft1_report(g)
    print(f"  {g['name']:<20} {resources(g['ops'])}, FT_1: {rep['FT1']}, "
          f"undetected single faults: {rep['undetected']}")
