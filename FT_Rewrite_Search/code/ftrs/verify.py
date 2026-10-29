"""Checker run on every rewrite: ideal behaviour, single-fault equivalence to the root, FT_1.

This is check chain C1/C2 of the experiment plan at t = 1, restricted to X-type gadgets on the
Steane code (the only ones built so far). A pass is evidence, not a proof (see L03, Sec. 4).
"""
from .ftlib import fe_compare, ft1_report, ideal_syndrome_check


def ideal_ok(g):
    res = ideal_syndrome_check(g, shots=64)
    for case, (syn, flg) in res.items():
        want = [0] if "= +1" in case else [1]
        if syn != want or (g["flags"] and flg != [0]):
            return False
    return True


def check(root, g):
    if not ideal_ok(g):
        return False, "ideal syndrome or flag wrong"
    only_root, only_g = fe_compare(root, g)
    if only_root or only_g:
        return False, f"not single-fault equivalent ({len(only_root)} effects only in root, {len(only_g)} only in child)"
    if not ft1_report(g)["FT1"]:
        return False, "not FT_1"
    return True, "ok"
