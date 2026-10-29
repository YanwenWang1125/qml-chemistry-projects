"""L01 verification: the three Shor-style gadgets for X0X1X2X3 on the Steane code.

1. Noiseless: the syndrome parity equals the stabilizer eigenvalue, and the flag stays 0.
2. Resource counts (ancillas, CNOTs, CNOT depth) against Jeon et al. Table I.
3. Every single fault: is the gadget FT_1 (no undetected single fault leaves a data error of
   weight > 1 modulo the stabilizer group)?
Run from this folder:  C:/Users/Administrator/venvs/qml/Scripts/python.exe L01_check_gadgets.py
"""
from ftlib import BASELINE, NO_FLAG, REDUCED, ft1_report, ideal_syndrome_check, resources

for g in (NO_FLAG, BASELINE, REDUCED):
    print("=" * 72)
    print(g["name"])
    print("  resources:", resources(g["ops"]))
    for case, (syn, flg) in ideal_syndrome_check(g).items():
        print(f"  noiseless, {case}: syndrome parity seen {syn}, flag seen {flg}")
    rep = ft1_report(g)
    print(f"  single faults: {rep['faults']}, raise the flag: {rep['detected']}, "
          f"undetected: {rep['undetected']}")
    print("  undetected faults by data weight mod S:", rep["undetected_weight_hist"])
    print("  FT_1:", rep["FT1"])
    for r in rep["violations"][:6]:
        print(f"    violation: {r['fault']:<40} data error {r['data']:<12} "
              f"weight mod S = {r['weight_mod_S']}")
    if len(rep["violations"]) > 6:
        print(f"    ... {len(rep['violations'])} violations in total")
