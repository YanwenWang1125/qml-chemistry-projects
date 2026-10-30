"""Checks for the ftrs generalisation (plan Sec. 0.2): Code object, faults_upto(t), Z-type gadgets.

1. Code.weight_mod / canonical_mod equal the group-list functions on every 7-qubit Pauli.
2. fault_columns equal effect() of every single fault.
3. Linearity: a 2-fault set injected at its two real times equals the XOR of the columns.
4. ft_report(t=1) equals the old single-fault analysis on the Shor gadgets.
5. Z-type: dual gadgets pass/fail exactly as their X-type originals.
6. t = 2 report on the Steane baseline, and weight_mod speed on a 17-qubit code.
Run:  C:/Users/Administrator/venvs/qml/Scripts/python.exe test_generalize.py
"""
import random
import time

import numpy as np

from ftrs.ftlib import (BASELINE, NO_FLAG, REDUCED, STEANE, STEANE_GROUP, Code, canonical_mod, dual,
                        effect, fault_columns, fe_compare, fe_compare_t, ft1_report, ft_report,
                        located_faults, propagate, single_faults, weight_mod)
from ftrs.verify import check, ideal_ok

random.seed(7)
ok = True


def report(name, cond):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")


# 1
bad = sum(STEANE.weight_mod(x, z) != weight_mod(x, z, STEANE_GROUP)
          or STEANE.canonical_mod(x, z) != canonical_mod(x, z, STEANE_GROUP)
          for x in range(128) for z in range(128))
report(f"Code coset functions = group-list functions on 16384 Paulis ({bad} differ)", bad == 0)

gadgets = [BASELINE, REDUCED, NO_FLAG]
# 2
for g in gadgets:
    cols = fault_columns(g["ops"])
    singles = list(single_faults(g["ops"]))
    flat = [c for loc in sorted(cols) for c in cols[loc]]
    same = len(flat) == len(singles) and all(
        (lab == l2) and propagate(x, z, g["ops"][st:]) == (cx, cz)
        for (l2, cx, cz), (lab, st, x, z) in zip(flat, singles))
    report(f"fault columns = single faults on {g['name']} ({len(flat)} faults)", same)

# 3
for g in gadgets:
    faults = list(located_faults(g["ops"]))
    cols = {(f[0], f[1]): propagate(f[3], f[4], g["ops"][f[2]:]) for f in faults}
    worst = 0
    for _ in range(2000):
        a, b = random.sample(faults, 2)
        if a[0] == b[0]:
            continue
        a, b = sorted((a, b), key=lambda f: f[2])
        x, z = propagate(a[3], a[4], g["ops"][a[2]:b[2]])
        x, z = propagate(x ^ b[3], z ^ b[4], g["ops"][b[2]:])
        ca, cb = cols[(a[0], a[1])], cols[(b[0], b[1])]
        worst += (x, z) != (ca[0] ^ cb[0], ca[1] ^ cb[1])
    report(f"2-fault linearity on {g['name']} ({worst} of 2000 differ)", worst == 0)

# 4
for g in gadgets:
    old_bad, old_und = 0, 0
    for label, start, x, z in single_faults(g["ops"]):
        dx, dz, _, ff = effect(g["ops"], start, x, z, g["syndrome"], g["flags"])
        if not any(ff):
            old_und += 1
            old_bad += weight_mod(dx, dz, STEANE_GROUP) > 1
    r = ft1_report(g)
    report(f"ft_report(t=1) = old analysis on {g['name']} (FT1 {r['FT1']}, undetected {r['undetected']}, "
           f"violations {len(r['violations'])})",
           r["undetected"] == old_und and len(r["violations"]) == old_bad)

# 5
dB, dR, dN = dual(BASELINE), dual(REDUCED), dual(NO_FLAG)
report("dual gadgets pass the ideal check (Z-type)", ideal_ok(dB) and ideal_ok(dR) and ideal_ok(dN))
report("an X-type gadget fails the ideal check when relabelled Z-type", not ideal_ok({**BASELINE, "kind": "Z"}))
report("FT_1: dual(baseline) and dual(reduced) pass, dual(no flag) fails",
       ft1_report(dB)["FT1"] and ft1_report(dR)["FT1"] and not ft1_report(dN)["FT1"])
xa, xb = fe_compare(BASELINE, REDUCED)
za, zb = fe_compare(dB, dR)
report(f"FE: baseline vs reduced {len(xa)}/{len(xb)}, duals {len(za)}/{len(zb)} (both 0/0)",
       not (xa or xb or za or zb))
report("verify.check accepts dual(reduced) against dual(baseline)", check(dB, dR) == (True, "ok"))
report("verify.check rejects dual(no flag) against dual(baseline)", not check(dB, dN)[0])

# 6
t0 = time.time()
r2 = ft_report(BASELINE, 2)
print(f"t=2 on the Steane baseline: fault sets {r2['fault_sets']}, unflagged {r2['undetected']}, "
      f"FT_2 {r2['FT']} ({len(r2['violations'])} violating sets), "
      f"{time.time() - t0:.1f} s")
report("t=1 part of fe_compare_t(t=2) agrees with t=1 (baseline vs reduced)",
       {k for k, v in fe_compare_t(BASELINE, REDUCED, 2)[0].items() if v == 1} == set())

big = Code("random 17-qubit CSS", 17, [random.sample(range(17), 8) for _ in range(8)],
           [random.sample(range(17), 8) for _ in range(8)])
xs = [(random.getrandbits(17), random.getrandbits(17)) for _ in range(2000)]
t0 = time.time()
for x, z in xs:
    big.weight_mod(x, z)
dt = (time.time() - t0) / len(xs)
brute = min(bin((xs[0][0] ^ a) | (xs[0][1] ^ b)).count("1") for a, b in big.group)
report(f"17-qubit weight_mod: {1e3 * dt:.2f} ms per uncached call, {len(big.group)} group elements, "
       f"matches brute force", big.weight_mod(*xs[0]) == brute)

# 7. vectorised checker (verify.check) = slow t-generic functions, on the E1.2 output circuits
import json
import os

from ftrs.ftlib import code_of, fe_signature_t
from ftrs.multi import GOTO, check_factory, derive, effects
from ftrs.search import composite_reductions

e12 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "e1_2_shor_search.json")
circs = [{"ops": [tuple(o) for o in r["ops"]], "syndrome": r["syndrome"], "flags": r["flags"]}
         for r in json.load(open(e12))["circuits"]]
for t in (1, 2):
    diff, t0 = 0, time.time()
    for g in circs + [REDUCED, NO_FLAG]:
        a, b = fe_compare_t(BASELINE, g, t)
        slow = (not ideal_ok(g)) or bool(a or b) or not ft_report(g, t)["FT"]
        diff += slow == check(BASELINE, g, t)[0]
    report(f"verify.check(t={t}) = slow fe_compare_t + ft_report on {len(circs) + 2} circuits "
           f"({diff} differ, {time.time() - t0:.1f} s)", diff == 0)
sig_slow = fe_signature_t(REDUCED, 2)
sig_fast = {(k[0], k[1][0]): v for k, v in __import__("ftrs.ftlib", fromlist=["x"]).effect_signature(
    REDUCED, 2, [REDUCED["syndrome"]])[0].items()}
report(f"effect_signature(t=2) = fe_signature_t(t=2) on Fig. 1b ({len(sig_fast)} effects)", sig_fast == sig_slow)

# 8. multi.check_factory(t=1) = the old set-based single-fault checker, on Goto children
rp, _ = derive(GOTO)
rsig_old, _ = effects(GOTO, rp)
proj = lambda keys: {(k[0][0], k[1]) for k in keys}
rsig_old = proj(rsig_old)
ident = (STEANE.canonical_mod(0, 0)[0], tuple(0 for _ in rp))


def old_multi(g):
    p, good = derive(g)
    if not good:
        return False
    sig, bad = effects(g, p)
    sig = proj(sig)
    return not ([k for k in rsig_old if k not in sig and k != ident]
                or [k for k in sig if k not in rsig_old and k != ident] or bad)


new1, new2 = check_factory(GOTO, 1), check_factory(GOTO, 2)
kids, frontier = [], [GOTO]
while frontier and len(kids) < 400:
    g = frontier.pop(0)
    for child, _ in composite_reductions(g, 2):
        kids.append(child)
        if new1(GOTO, child)[0]:
            frontier.append(child)
kids = kids[:400]
diff = sum(old_multi(g) != new1(GOTO, g)[0] for g in kids)
n1, n2 = sum(new1(GOTO, g)[0] for g in kids), sum(new2(GOTO, g)[0] for g in kids)
report(f"multi check_factory(t=1) = old checker on {len(kids)} Goto children ({diff} differ; "
       f"accepted t=1 {n1}, t=2 {n2})", diff == 0 and n2 <= n1)

print("\nALL PASS" if ok else "\nSOME CHECKS FAILED")
