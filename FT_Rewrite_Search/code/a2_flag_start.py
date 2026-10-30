"""Direction A, step 2: a Liou-like start (one measurement qubit per stabilizer, own flags) and
flag merging on it.

Z-type gadget per generator (as in Liou Fig. 3, before any sharing): measurement qubit m in |0>,
data CNOTs d_i -> m in support order; each flag f (|+>, measured in X) is a pair of CNOTs f -> m
placed after data CNOT a and after data CNOT b (a < b), so it brackets data CNOTs a+1..b.
Flag positions per weight are chosen greedily against the distance-5 condition of d1 (every set
of s <= 3 faults: s + flags fired >= residual weight). The eight gadgets form the Z half, its
Hadamard dual the X half. Then greedy flag merging (a1_flag_merge: accept only if noiseless,
single-round FT_2, 5-round memory distance 5), and setting-B LER of the start and the end.
Writes data/circuits/a2_flag_start.json, results/a2_flag_start.json.
Run:  python a2_flag_start.py [max_errors]
"""
import itertools
import json
import os
import sys
import time

from ftrs.codes import COLOR_17, COLOR_17_SUPPORTS
from ftrs.decoders import DemLookup
from ftrs.ftlib import dual, effect_signature, faults_upto, mask, resources
from ftrs.mem import wilson
from ftrs.memgen import build
from ftrs.multi import derive

here = os.path.dirname(os.path.abspath(__file__))
max_errors = int(sys.argv[1]) if len(sys.argv) > 1 else 400
W8_ORDER = sys.argv[2] if len(sys.argv) > 2 else "sorted"  # "liou": Liou Fig. 20 data order for g8
TAG = sys.argv[3] if len(sys.argv) > 3 else ""
t_start = time.time()


def gadget(support, flags, base):
    m = base
    fl = list(range(base + 1, base + 1 + len(flags)))
    after = {}
    for f, (a, b) in zip(fl, flags):
        after.setdefault(a, []).append(f)
        after.setdefault(b, []).append(f)
    ops = [("R", m)] + [("R", f) for f in fl] + [("H", f) for f in fl]
    for i, d in enumerate(support):
        ops.append(("CX", d, m))
        ops += [("CX", f, m) for f in after.get(i, [])]
    ops += [("H", f) for f in fl] + [("M", q) for q in [m] + fl]
    return {"ops": ops, "syndrome": [m], "flags": fl, "kind": "Z", "targets": [tuple(support)], "code": COLOR_17}


def violations(g, smax=3):
    fl, dm = mask(g["flags"]), COLOR_17.data_mask
    bad = 0
    for s, _, x, z in faults_upto(g["ops"], smax):
        f = bin(x & fl).count("1")
        if s + f < 4 and COLOR_17.weight_mod(x & dm, z & dm) > s + f:
            bad += 1
    return bad


def choose(support):
    w = len(support)
    cands = [(a, b) for a in range(w - 1) for b in range(a + 1, w - 1)]
    v = lambda fs: violations(gadget(support, fs, 17))
    fs, cur = [], v([])
    while cur:
        scores = {c: v(fs + [c]) for c in cands if c not in fs}
        best = min(scores, key=scores.get)
        fs.append(best)
        cur = scores[best]
        print(f"  weight {w}: +flag{best} -> {cur} violations ({time.time() - t_start:.0f} s)", flush=True)
        if len(fs) > 8:
            raise RuntimeError("no flag set")
    for c in fs[:]:
        if not v([x for x in fs if x != c]):
            fs.remove(c)
    return fs


LIOU_G8 = (13, 9, 14, 5, 2, 6, 3, 10)
per = {}
for sup in sorted({tuple(s) for s in COLOR_17_SUPPORTS}, key=len):
    if len(sup) not in per:
        per[len(sup)] = choose(LIOU_G8 if (W8_ORDER == "liou" and len(sup) == 8) else sup)
        print(f"weight {len(sup)}: flags {per[len(sup)]}", flush=True)

ops, syn, flags, base = [], [], [], 17
LIOU_G8 = (13, 9, 14, 5, 2, 6, 3, 10)  # m8's data CNOT order in Liou Fig. 20 (0-indexed)
order = lambda sup: LIOU_G8 if (W8_ORDER == "liou" and len(sup) == 8) else sup
for sup in COLOR_17_SUPPORTS:
    g = gadget(order(sup), per[len(sup)], base)
    ops += g["ops"]
    syn += g["syndrome"]
    flags += g["flags"]
    base += 1 + len(g["flags"])
Zh = {"name": "Liou-like start", "kind": "Z", "targets": [tuple(s) for s in COLOR_17_SUPPORTS],
      "syndrome": syn, "flags": flags, "ops": ops, "code": COLOR_17}


def halves(zh):
    return zh, {**dual(zh), "targets": [tuple(s) for s in COLOR_17.x_supports]}


def accepted(zh):
    zh, xh = halves(zh)
    pz, okz = derive(zh, shots=2000)
    px, okx = derive(xh, shots=2000)
    if not (okz and okx):
        return False, "ideal"
    if effect_signature(zh, 2, pz)[1] or effect_signature(xh, 2, px)[1]:
        return False, "not FT_2"
    for mem in ("Z", "X"):
        c, _ = build(COLOR_17, [zh, xh], rounds=5, memory=mem, p=1e-3)
        err = c.search_for_undetectable_logical_errors(dont_explore_detection_event_sets_with_size_above=4,
              dont_explore_edges_with_degree_above=4, dont_explore_edges_increasing_symptom_degree=False,
              canonicalize_circuit_errors=True)
        if len(err) < 5:
            return False, f"distance {len(err)}"
    return True, "ok"


def merge(zh, fa, fb):
    """Redirect fb's CNOTs to fa and read fa out after its last CNOT. (Before 2026-10-29 22:30
    fa's H and M stayed at the end of its own gadget, so in a serial round fb's CNOTs hit fa after
    its readout and cancelled in pairs: those "merges" were flag deletions.)"""
    ops = [tuple([op[0]] + [fa if q == fb else q for q in op[1:]]) for op in zh["ops"]
           if op not in (("R", fb), ("M", fb), ("H", fb))]
    cx_idx = [i for i, op in enumerate(ops) if op[0] == "CX" and fa in op[1:]]
    mid = [i for i, op in enumerate(ops) if op[0] != "CX" and op[1] == fa and cx_idx[0] < i < cx_idx[-1]]
    tail = [ops[i] for i in mid]  # fa's readout (H, M) that now sits before some of its CNOTs
    ops = [op for i, op in enumerate(ops) if i not in set(mid)]
    last = max(i for i, op in enumerate(ops) if op[0] == "CX" and fa in op[1:])
    ops = ops[:last + 1] + tail + ops[last + 1:]
    return {**zh, "ops": ops, "flags": [f for f in zh["flags"] if f != fb]}


def ler(zh, p, mem):
    zh, xh = halves(zh)
    c, _ = build(COLOR_17, [zh, xh], rounds=5, memory=mem, p=p)
    dec = DemLookup(c)
    sampler = c.compile_detector_sampler(seed=11)
    e = n = 0
    while e < max_errors and n < 2 * 10 ** 8:
        d, o = sampler.sample(20000, separate_observables=True)
        e += int((dec.decode_batch(d) != o).any(axis=1).sum())
        n += len(d)
    return {"p": p, "memory": mem, "errors": e, "shots": n, "ler": e / n, "ci95": list(wilson(e, n))}


ok, why = accepted(Zh)
print(f"start: {resources(ops)}, {len(flags)} flags, accepted {ok} ({why})", flush=True)
out = {"flags_per_weight": {str(k): v for k, v in per.items()}, "start": resources(ops), "start_flags": len(flags),
       "start_accepted": ok, "start_reason": why}
steps = []
cur = Zh
if ok:
    while True:
        best, verdicts = None, {}
        for fa, fb in itertools.combinations(cur["flags"], 2):
            cand = merge(cur, fa, fb)
            good, why = accepted(cand)
            verdicts[why] = verdicts.get(why, 0) + 1
            if good and (best is None or resources(cand["ops"])["cnot_depth"] < resources(best[0]["ops"])["cnot_depth"]):
                best = (cand, (fa, fb))
        print(f"merge step {len(steps) + 1}: {verdicts} ({time.time() - t_start:.0f} s)", flush=True)
        if best is None:
            break
        cur = best[0]
        steps.append({"merge": best[1], "resources": resources(cur["ops"]), "flags": len(cur["flags"]), "verdicts": verdicts})
        print(f"  merged {best[1]} -> {resources(cur['ops'])}, {len(cur['flags'])} flags", flush=True)
    out["merge_steps"] = steps
    out["final"] = resources(cur["ops"])
    out["final_flags"] = len(cur["flags"])
    out["ler"] = {"start": [ler(Zh, p, m) for p in (1e-3, 2e-3) for m in ("Z", "X")],
                  "final": [ler(cur, p, m) for p in (1e-3, 2e-3) for m in ("Z", "X")] if steps else None}
out["seconds"] = round(time.time() - t_start, 1)
print(json.dumps(out))
json.dump({k: v for k, v in cur.items() if k != "code"} | {"code": COLOR_17.name},
          open(os.path.join(here, "..", "data", "circuits", f"a2_flag_final{TAG}.json"), "w"), indent=1)
json.dump({k: v for k, v in Zh.items() if k != "code"} | {"code": COLOR_17.name},
          open(os.path.join(here, "..", "data", "circuits", f"a2_flag_start{TAG}.json"), "w"), indent=1)
json.dump(out, open(os.path.join(here, "..", "results", f"a2_flag_start{TAG}.json"), "w"), indent=1)
