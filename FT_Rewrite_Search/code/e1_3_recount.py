"""E1.3 recount: try coarser counting conventions on results/e1_3_steane_search.json to see
which one lands near Jeon's 391 candidates (Sec. V C). Our search dedups only by ASAP-layer
canonical form, which keeps circuits that differ by commuting CNOTs or ancilla labels apart.
Writes results/e1_3_recount.json.
Run on Fir:  sbatch e1_3_recount.sh
Also works on another search output:  python e1_3_recount.py <search.json> <out.json>
"""
import json
import os
import sys
from collections import Counter

DATA = 7  # qubits 0-6 are data, >= 7 are ancillas
here = os.path.dirname(os.path.abspath(__file__))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "results", "e1_3_steane_search.json")
res = json.load(open(path))
rows = res["circuits"]
for r in rows:  # E1.2 rows lack H count and total depth
    r.setdefault("H", sum(op[0] == "H" for op in r["ops"]))
    r.setdefault("depth", None)
print(f"distinct (ASAP-layer canonical) {len(rows)} | BFS nodes {res['bfs_nodes_dedup']} (cap 100000) "
      f"| children generated {res.get('children_generated')}")


def ancilla_labels(ops):
    """Relabel ancillas by a label-free signature: their sorted gate roles with data qubits."""
    sig = {}
    for op in ops:
        for i, q in enumerate(op[1:]):
            if q >= DATA:
                other = [x if x < DATA else -1 for x in op[1:] if x != q]
                sig.setdefault(q, []).append((op[0], i, tuple(other)))
    order = sorted(sig, key=lambda q: (sorted(sig[q]), q))
    return {q: DATA + k for k, q in enumerate(order)}


def gate_multiset(ops):
    """Order-free key: the multiset of gates after ancilla relabelling. Circuits that differ only
    in the order of (commuting or not) gates collapse, so this is a lower bound on true classes."""
    m = ancilla_labels(ops)
    return tuple(sorted((op[0], *(m.get(q, q) for q in op[1:])) for op in ops))


def cx_commute_key(ops):
    """Per-qubit gate sequences with runs of CNOTs that share this qubit in the same role sorted,
    i.e. merge CNOTs that commute because they share a control or share a target."""
    m = ancilla_labels(ops)
    seq = {}
    for op in ops:
        g = (op[0], *(m.get(q, q) for q in op[1:]))
        for i, q in enumerate(g[1:]):
            seq.setdefault(q, []).append((g, i))
    key = []
    for q in sorted(seq):
        out, run = [], []
        for g, i in seq[q]:
            if g[0] == "CX" and run and run[-1][1] == i:
                run.append((g, i))
            else:
                out += sorted(run)
                run = [(g, i)] if g[0] == "CX" else []
                if g[0] != "CX":
                    out.append((g, i))
        out += sorted(run)
        key.append((q, tuple(out)))
    return tuple(key)


keys = {
    "resources(anc,CNOT,H,depth,CNOTdepth)": lambda r: (r["ancillas"], r["cnots"], r["H"], r["depth"], r["cnot_depth"]),
    "resources(anc,CNOT,H)": lambda r: (r["ancillas"], r["cnots"], r["H"]),
    "syndrome+flags layout": lambda r: (json.dumps(r["syndrome"]), json.dumps(r["flags"])),
    "CX-commutation canonical": lambda r: cx_commute_key(r["ops"]),
    "gate multiset (order-free)": lambda r: gate_multiset(r["ops"]),
}
out = {"distinct_asap": len(rows), "bfs_nodes": res["bfs_nodes_dedup"], "counts": {}, "by_steps": {}}
for name, f in keys.items():
    n = len({f(r) for r in rows})
    out["counts"][name] = n
    print(f"{name:42s} {n}")

steps = Counter(r["steps"] for r in rows)
print("\nby Bell reductions:", dict(sorted(steps.items())))
for k in sorted(steps):
    sub = [r for r in rows if r["steps"] <= k]
    n_ms = len({gate_multiset(r["ops"]) for r in sub})
    n_cx = len({cx_commute_key(r["ops"]) for r in sub})
    out["by_steps"][k] = {"asap": len(sub), "cx_commute": n_cx, "multiset": n_ms}
    print(f"  steps <= {k}: asap {len(sub)}, CX-commutation {n_cx}, multiset {n_ms}")

dst = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "results", "e1_3_recount.json")
json.dump(out, open(dst, "w"), indent=1)
print("saved", os.path.normpath(dst))
