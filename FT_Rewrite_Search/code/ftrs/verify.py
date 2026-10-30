"""Checker run on every rewrite: ideal behaviour, fault equivalence to the root, FT_t.

This is check chain C1/C2 of the experiment plan for gadgets with a one-bit syndrome: X-type
(default) or Z-type (gadget["kind"] = "Z", see ftlib.dual), on any CSS code carried in
gadget["code"] (default Steane). t = 1 by default; t = 2 catches rewrites that keep the
single-fault effects but change the two-fault ones (E1.2 regression, 2026-10-29).
A pass is evidence, not a proof (see L03, Sec. 4).
"""
from .ftlib import code_of, effect_signature, ideal_syndrome_check

_root_sigs = {}


def ideal_ok(g):
    res = ideal_syndrome_check(g, shots=64)
    for case, (syn, flg) in res.items():
        want = [0] if "= +1" in case else [1]
        if syn != want or (g["flags"] and flg != [0]):
            return False
    return True


def check(root, g, t=1):
    if not ideal_ok(g):
        return False, "ideal syndrome or flag wrong"
    if (id(root), t) not in _root_sigs:
        _root_sigs[(id(root), t)] = effect_signature(root, t, [root["syndrome"]])[0]
    rsig = _root_sigs[(id(root), t)]
    sig, bad = effect_signature(g, t, [g["syndrome"]])
    ident = (code_of(root).canonical_mod(0, 0), (0,))
    only_root = [k for k, v in rsig.items() if k != ident and sig.get(k) != v]
    only_g = [k for k, v in sig.items() if k != ident and rsig.get(k) != v]
    if only_root or only_g:
        order = "single-fault" if t == 1 else f"{t}-fault"
        return False, f"not {order} equivalent ({len(only_root)} effects only in root, {len(only_g)} only in child)"
    if bad:
        return False, f"not FT_{t}"
    return True, "ok"
