"""E1.5b: t = 2 negative control on [[17,1,5]] (supplementary plan, E1.5b; AUTOLOOP C1).

Under test: the t = 2 checker, multi.check_factory(LIOU, 2) (ideal check, 2-fault effect
equivalence to the root, FT_2), and its FT_2 part alone (effect_signature(..., 2)[1] == 0).
Truth: oracle_ft (stim tableaus, brute-force weight; no ftlib.propagate) for FT_1/FT_2 and
multi.derive (stim sampling) for the noiseless check.

0. Oracle self-test on Steane: FT_1 of the strict-mode E1.2 circuits and the 47 E1.5 variants
   must equal ftlib.ft1_report. (The plan names "E1.2's 7 circuits"; the strict mode now gives 5.)
1. Variants of the Liou Fig. 20 circuit: P (FE-preserving edits, must pass), A1 (remove a flag),
   A2 (reorder g8's data CNOTs, seeds 0-49), A3 (move one flag CNOT), B (break FT_1), C (break
   the circuit). K = A1-A3 variants the oracle calls ideal-correct, FT_1 and not FT_2.
Writes results/e1_5b_t2_negative_control.json.
Run on Fir:  sbatch e1_5b.sh
"""
import json
import os
import random
import time
from collections import Counter

from ftrs.codes import COLOR_17
from ftrs.ftlib import STEANE_SUPPORTS, effect_signature, ft1_report
from ftrs.multi import check_factory, derive
from oracle_ft import Weight, ft_oracle

here = os.path.dirname(os.path.abspath(__file__))
res_dir = os.path.join(here, "..", "results")
t_start = time.time()

# --- 0. oracle self-test on Steane ---------------------------------------------------------
w7 = Weight(7, STEANE_SUPPORTS, STEANE_SUPPORTS)
src = open(os.path.join(here, "e1_5_negative_control.py"), encoding="utf-8").read()
ns = {"__file__": os.path.join(here, "e1_5_negative_control.py")}
exec(src[:src.index("counts, passed = {}, []")], ns)
steane = [g for _, _, g in ns["mutants"]] + [ns["root"]]
strict = json.load(open(os.path.join(res_dir, "e1_2_shor_search_t2_jeon.json")))["circuits"]
steane += [{"ops": [tuple(o) for o in r["ops"]], "syndrome": r["syndrome"], "flags": r["flags"]} for r in strict]
agree = sum(ft_oracle(g["ops"], g["flags"], w7, 1)[0] == ft1_report(g)["FT1"] for g in steane)
print(f"oracle self-test on Steane: {agree}/{len(steane)} FT_1 verdicts agree with ftlib", flush=True)

# --- 1. Liou variants ----------------------------------------------------------------------
LIOU = json.load(open(os.path.join(here, "..", "data", "circuits", "liou_17_1_5.json")))
LIOU = {**LIOU, "ops": [tuple(o) for o in LIOU["ops"]], "targets": [tuple(t) for t in LIOU["targets"]],
        "code": COLOR_17}
meas, flags = LIOU["syndrome"], LIOU["flags"]
w17 = Weight(17, COLOR_17.x_supports, COLOR_17.z_supports)
check = check_factory(LIOU, 2)


def variant(ops, **kw):
    return {**LIOU, "ops": ops, **kw}


def relabel(seed):
    rng = random.Random(seed)
    perm_m, perm_f = meas[:], flags[:]
    rng.shuffle(perm_m)
    rng.shuffle(perm_f)
    m = dict(zip(meas + flags, perm_m + perm_f))
    ops = [tuple([op[0]] + [m.get(q, q) for q in op[1:]]) for op in LIOU["ops"]]
    return variant(ops, syndrome=[m[q] for q in meas], flags=[m[q] for q in flags])


def disjoint_swaps(seed, n=3):
    rng = random.Random(seed)
    ops = LIOU["ops"][:]
    idx = [i for i in range(len(ops) - 1) if ops[i][0] == ops[i + 1][0] == "CX"
           and not set(ops[i][1:]) & set(ops[i + 1][1:])]
    for i in rng.sample(idx, n):
        ops[i], ops[i + 1] = ops[i + 1], ops[i]
    return variant(ops)


def remove_flag(f):
    return variant([op for op in LIOU["ops"] if f not in op[1:]], flags=[x for x in flags if x != f])


def reorder_g8(seed):
    rng = random.Random(seed)
    m8 = meas[7]
    pos = [i for i, op in enumerate(LIOU["ops"]) if op[0] == "CX" and op[2] == m8 and op[1] < 17]
    data = [LIOU["ops"][i][1] for i in pos]
    rng.shuffle(data)
    ops = LIOU["ops"][:]
    for i, d in zip(pos, data):
        ops[i] = ("CX", d, m8)
    return variant(ops)


def move_flag_cnot():
    out = []
    for k, op in enumerate(LIOU["ops"]):
        if op[0] != "CX" or op[1] not in flags:
            continue
        m = op[2]
        rest = LIOU["ops"][:k] + LIOU["ops"][k + 1:]
        slots = [i for i, o in enumerate(rest) if o[0] == "CX" and o[2] == m and o[1] < 17] + \
                [max(i for i, o in enumerate(rest) if m in o[1:] and o[0] == "CX") + 1]
        for s in slots:
            ops = rest[:s] + [op] + rest[s:]
            if ops != LIOU["ops"]:
                out.append((f"move {op} to slot {s}", variant(ops)))
    return out


def break_ft1(seed):
    rng = random.Random(seed)
    ops = LIOU["ops"][:]
    k = rng.choice([i for i, op in enumerate(ops) if op[0] == "CX" and op[1] < 17])
    c, t = ops[k][1:]
    ops[k] = ("CX", c, rng.choice([m for m in meas if m != t]))
    return variant(ops)


def break_circuit(seed):
    rng = random.Random(seed)
    ops = LIOU["ops"][:]
    kind = seed % 3
    if kind == 0:
        del ops[rng.choice([i for i, op in enumerate(ops) if op[0] == "CX"])]
    elif kind == 1:
        del ops[rng.choice([i for i, op in enumerate(ops) if op[0] == "H"])]
    else:
        k = rng.choice([i for i, op in enumerate(ops) if op[0] == "CX"])
        c, t = ops[k][1:]
        ops[k] = ("CX", c, rng.choice([q for q in range(17) if q != c]))
    return variant(ops)


variants = [("P", "original", LIOU)]
variants += [("P", f"relabel ancillas seed {s}", relabel(s)) for s in range(5)]
variants += [("P", f"3 swaps of adjacent disjoint CNOTs, seed {s}", disjoint_swaps(s)) for s in range(5)]
variants += [("A1", f"remove flag q{f}", remove_flag(f)) for f in flags]
variants += [("A2", f"reorder g8 data CNOTs, seed {s}", reorder_g8(s)) for s in range(50)]
variants += [("A3", d, g) for d, g in move_flag_cnot()]
variants += [("B", f"data CNOT moved to another measurement qubit, seed {s}", break_ft1(s)) for s in range(10)]
variants += [("C", f"random edit seed {s}", break_circuit(s)) for s in range(30)]
print(f"{len(variants)} variants: {dict(Counter(v[0] for v in variants))}", flush=True)

rows = []
for grp, desc, g in variants:
    t0 = time.time()
    parities, ideal = derive(g, shots=2000)
    o1 = o2 = None
    viol = None
    if ideal:
        o1 = ft_oracle(g["ops"], g["flags"], w17, 1)[0]
        o2, viol, _ = ft_oracle(g["ops"], g["flags"], w17, 2)
    ok, why = check(LIOU, g)
    ft2_part = (effect_signature(g, 2, parities)[1] == 0) if ideal else None
    rows.append({"group": grp, "change": desc, "oracle_ideal": ideal, "oracle_FT1": o1, "oracle_FT2": o2,
                 "oracle_violation": viol, "checker_pass": ok, "checker_reason": why,
                 "checker_FT2_part": ft2_part, "seconds": round(time.time() - t0, 1)})
    print(f"[{grp}] {desc}: ideal {ideal}, oracle FT1 {o1} FT2 {o2} | checker {ok} ({why}) | "
          f"FT2 part {ft2_part} | {rows[-1]['seconds']} s", flush=True)

# --- summary -------------------------------------------------------------------------------
K = [r for r in rows if r["group"] in ("A1", "A2", "A3") and r["oracle_ideal"] and r["oracle_FT1"]
     and r["oracle_FT2"] is False]
conf = Counter((r["oracle_FT2"], r["checker_FT2_part"]) for r in rows if r["oracle_ideal"])
summary = {
    "oracle_self_test": f"{agree}/{len(steane)}",
    "P_pass": f"{sum(r['checker_pass'] for r in rows if r['group'] == 'P')}/{sum(r['group'] == 'P' for r in rows)}",
    "K_size": len(K), "K_detected": sum(not r["checker_pass"] for r in K),
    "K_detected_by_FT2_part": sum(r["checker_FT2_part"] is False for r in K),
    "B_rejected": f"{sum(not r['checker_pass'] for r in rows if r['group'] == 'B')}/{sum(r['group'] == 'B' for r in rows)}",
    "C_rejected": f"{sum(not r['checker_pass'] for r in rows if r['group'] == 'C')}/{sum(r['group'] == 'C' for r in rows)}",
    "confusion_oracleFT2_vs_checkerFT2part": {f"oracle {a} / checker {b}": n for (a, b), n in conf.items()},
    "wall_seconds": round(time.time() - t_start, 1),
}
print("\nsummary:", json.dumps(summary, indent=1))
json.dump({"summary": summary, "variants": rows}, open(os.path.join(res_dir, "e1_5b_t2_negative_control.json"), "w"), indent=1)
print("saved")
