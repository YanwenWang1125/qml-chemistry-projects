"""Direction A, step 1: a flag-merge rewrite (Liou-style flag sharing) on top of Jeon's rules.

Flag merge (fa <- fb): every CNOT onto flag fb is redirected to flag fa, and fb's preparation and
measurement are dropped. The shared flag then measures the product of both parity checks, which
is still deterministic without noise. This is not a fault-equivalent rewrite in general (a fault
that fired both flags now cancels on the shared one), so each merge is accepted only if
  - the noiseless check passes for the X half and its dual Z half,
  - both halves stay single-round FT_2,
  - the 5-round memory circuit keeps circuit distance 5 in both bases (stim search).
Greedy: at each step try every flag pair and keep the accepted merge that leaves the lowest CNOT
depth; stop when no merge is accepted. Starts from a round of e2_7_round_ler.py.
Writes results/a1_flag_merge_<variant>.json (steps, final X-half ops).
Run:  python a1_flag_merge.py <variant>      e.g. fusion_d5-jeon_d5
"""
import itertools
import json
import os
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
variant = sys.argv[1]
src = open(os.path.join(here, "e2_7_round_ler.py"), encoding="utf-8").read()
ns = {"__file__": os.path.join(here, "e2_7_round_ler.py")}
sys.argv = ["x", variant, "0.001", "Z"]
exec(src[:src.index("t0 = time.time()")], ns)
COLOR_17, build, derive, dual = ns["COLOR_17"], ns["build"], ns["derive"], ns["dual"]
effect_signature, resources = ns["effect_signature"], ns["resources"]
t_start = time.time()


def halves(xh):
    zh = {**dual(xh), "targets": [tuple(s) for s in COLOR_17.z_supports]}
    return xh, zh


def accepted(xh):
    xh, zh = halves(xh)
    px, okx = derive(xh, shots=2000)
    pz, okz = derive(zh, shots=2000)
    if not (okx and okz):
        return False, "ideal"
    if effect_signature(xh, 2, px)[1] or effect_signature(zh, 2, pz)[1]:
        return False, "not FT_2"
    for mem in ("Z", "X"):
        c, _ = build(COLOR_17, [zh, xh], rounds=5, memory=mem, p=1e-3)
        err = c.search_for_undetectable_logical_errors(dont_explore_detection_event_sets_with_size_above=4,
              dont_explore_edges_with_degree_above=4, dont_explore_edges_increasing_symptom_degree=False,
              canonicalize_circuit_errors=True)
        if len(err) < 5:
            return False, f"distance {len(err)} ({mem} memory)"
    return True, "ok"


def merge(xh, fa, fb):
    ops = []
    for op in xh["ops"]:
        if op in (("R", fb), ("M", fb)):
            continue
        ops.append(tuple([op[0]] + [fa if q == fb else q for q in op[1:]]))
    return {**xh, "ops": ops, "flags": [f for f in xh["flags"] if f != fb],
            "syndrome": [q for q in xh["syndrome"] if q != fb]}


cur = ns["Xh"]
ok, why = accepted(cur)
steps = [{"step": 0, "merge": None, "resources": resources(cur["ops"]), "flags": len(cur["flags"]), "accepted": ok}]
print(f"start {variant}: {resources(cur['ops'])}, {len(cur['flags'])} flags, accepted {ok} ({why})", flush=True)
assert ok
while True:
    best, tried, reasons = None, 0, {}
    for fa, fb in itertools.permutations(cur["flags"], 2):
        if fa > fb:
            continue  # merging fb into fa or fa into fb gives the same circuit up to labels
        cand = merge(cur, fa, fb)
        tried += 1
        ok, why = accepted(cand)
        reasons[why] = reasons.get(why, 0) + 1
        if ok and (best is None or resources(cand["ops"])["cnot_depth"] < resources(best[0]["ops"])["cnot_depth"]):
            best = (cand, (fa, fb))
    print(f"step {len(steps)}: tried {tried}, verdicts {reasons}, {time.time() - t_start:.0f} s", flush=True)
    if best is None:
        break
    cur = best[0]
    steps.append({"step": len(steps), "merge": best[1], "resources": resources(cur["ops"]), "flags": len(cur["flags"]),
                  "verdicts": reasons})
    print(f"  merged {best[1]} -> {resources(cur['ops'])}, {len(cur['flags'])} flags", flush=True)

out = {"variant": variant, "steps": steps, "final_resources": resources(cur["ops"]), "final_flags": len(cur["flags"]),
       "final_ops": cur["ops"], "final_flag_qubits": cur["flags"], "seconds": round(time.time() - t_start, 1)}
print(json.dumps({k: v for k, v in out.items() if k != "final_ops"}))
json.dump(out, open(os.path.join(here, "..", "results", f"a1_flag_merge_{variant}.json"), "w"), indent=1)
