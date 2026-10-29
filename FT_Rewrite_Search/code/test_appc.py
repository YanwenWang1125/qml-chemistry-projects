"""E1.3-fix: can our rules turn each Appendix C circuit into the next one?

For each step F_i -> F_{i+1}: (1) both pass the checker against the Goto baseline; (2) a free BFS
over single rule applications (any rule, up to MAX_DEPTH rewrites, exactly one Bell reduction)
reaches a circuit isomorphic to F_{i+1}; (3) the composite search of ftrs/search.py does too.
Run:  C:/Users/Administrator/venvs/qml/Scripts/python.exe test_appc.py
"""
from ftrs.appc import F0, F1, F2, F3a, F3b, F4, isomorphic
from ftrs.ftlib import resources
from ftrs.multi import GOTO, check_factory
from ftrs.rules import basis_swap, bell_reduction, restricted_commutation, target_swap
from ftrs.search import composite_reductions

MAX_DEPTH = 4
check = check_factory(GOTO)

print("F0 isomorphic to GOTO:", isomorphic(F0, GOTO))
for g in (F0, F1, F2, F3a, F3b, F4):
    r = resources(g["ops"])
    print(f"{g['name']}: checker {check(GOTO, g)}, ancillas {r['ancillas']}, CNOT {r['cnots']}")


def free_reach(src, dst):
    frontier = [(src, [])]
    for depth in range(MAX_DEPTH):
        nxt = []
        for g, path in frontier:
            for rule in (restricted_commutation, target_swap, basis_swap, bell_reduction):
                if rule is bell_reduction and any(p.startswith("Bell") for p in path):
                    continue
                for h, d in rule(g):
                    if d.startswith("Bell") and isomorphic(h, dst):
                        return path + [d]
                    nxt.append((h, path + [d]))
        frontier = nxt
    return None


steps = [(F0, F1), (F1, F2), (F2, F3a), (F2, F3b), (F3a, F4), (F3b, F4)]
for a, b in steps:
    comp = [s for h, s in composite_reductions(a) if isomorphic(h, b)]
    print(f"\n{a['name']} -> {b['name']}")
    print("  composite search:", comp[0] if comp else "NOT reached")
    if not comp:
        print("  free BFS:", free_reach(a, b) or f"NOT reached within {MAX_DEPTH} rewrites")
