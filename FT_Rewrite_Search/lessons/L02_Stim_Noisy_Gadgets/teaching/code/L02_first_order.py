"""L02 hand-check: predict the first-order malignant rate of the gadget without verification.

Each single fault that leaves an undetected X error of weight >= 2 (mod X stabilizers) adds its
probability to the coefficient c1 in P(malignant) ~ c1 * p. Probabilities in the noise model of
ftlib.to_stim: preparation and measurement flips p, each two-qubit Pauli after a CNOT p/15,
each single-qubit Pauli after H (3p/100)/3 = p/100.
Run from this folder:  C:/Users/Administrator/venvs/qml/Scripts/python.exe L02_first_order.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                "L01_FT_and_Shor_Gadget", "teaching", "code"))
from ftlib import NO_FLAG, STEANE_SUPPORTS, analyse, stabilizer_group, weight_mod  # noqa: E402

x_group = stabilizer_group(STEANE_SUPPORTS, [])
c1 = 0.0
for r in analyse(NO_FLAG):
    if r["detected"] or weight_mod(r["dx"], 0, x_group) < 2:
        continue
    f = r["fault"]
    prob = 1.0 if f.startswith(("prep", "measurement")) else (1 / 15 if "after CX" in f else 0.01)
    c1 += prob
    print(f"{f:<42} data {r['data']:<8} probability / p = {prob:.4f}")
print(f"first-order coefficient c1 = {c1:.4f}")
print(f"predicted P(malignant) at p = 2e-3: {c1 * 2e-3:.2e};  at p = 5e-3: {c1 * 5e-3:.2e}")
