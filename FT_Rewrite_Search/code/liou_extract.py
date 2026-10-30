"""Read the CNOT list of Liou & Lai 2024 (arXiv:2407.00607) App. A figures from the PDF vector data.

Figures 19 and 20 (page 17) are PDF form XObjects /Im19, /Im20. Each wire is a long horizontal
stroke; each CNOT is a vertical stroke whose ends sit on two wires, with a filled dot (control)
on one end and a stroked circle (target) on the other; dashed verticals separate time steps.
Wire order from the top: data 1..17, then m1..m8, then the flag qubits.
Prints the CNOTs in left-to-right order with their time step, for transcription checking.
Run:  python liou_extract.py <liou2024.pdf> [Im20|Im19]
"""
import re
import sys

from pypdf import PdfReader

pdf, name = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else "Im20")
form = PdfReader(pdf).pages[16]["/Resources"]["/XObject"]["/" + name].get_object()
stream = form.get_data().decode("latin1")

num = r"-?\d*\.?\d+(?:[eE]-?\d+)?"
tok = re.compile(rf"({num})|(\[[^\]]*\]\s*{num}\s+d)|([A-Za-z*'\"]+)")

wires, stems, dots, rings, dashes = [], [], [], [], []
stack, path, dash = [], [], False
cur = None
for m in tok.finditer(stream):
    if m.group(1):
        stack.append(float(m.group(1)))
        continue
    if m.group(2):
        dash = not re.match(r"\[\s*\]", m.group(2))
        continue
    op = m.group(3)
    if op == "q":
        dash = False
    if op == "m":
        cur = [(stack[-2], stack[-1])]
        path.append(cur)
    elif op == "l":
        cur.append((stack[-2], stack[-1]))
    elif op == "c":
        cur.append((stack[-2], stack[-1]))
        cur.append(("c",))
    elif op in ("S", "s", "f", "f*", "B", "b", "F"):
        for seg in path:
            curved = any(len(p) == 1 for p in seg)
            pts = [p for p in seg if len(p) == 2]
            xs, ys = [p[0] for p in pts], [p[1] for p in pts]
            if curved:
                c = ((max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2, (max(xs) - min(xs)) / 2)
                (dots if op.startswith("f") or op == "F" or op == "B" else rings).append(c)
            elif len(pts) == 2:
                (x0, y0), (x1, y1) = pts
                if abs(y0 - y1) < 1e-3 and abs(x1 - x0) > 500:
                    wires.append(y0)
                elif abs(x0 - x1) < 0.5 and abs(y1 - y0) > 1:
                    (dashes if dash else stems).append((x0, min(y0, y1), max(y0, y1)))
        path = []
    if op in ("re", "W", "W*", "n"):
        path = [] if op == "n" else path
    if op not in ("m", "l", "c"):
        stack = []

wires = sorted(set(round(y, 1) for y in wires), reverse=True)
labels = [str(i) for i in range(1, 18)] + [f"m{i}" for i in range(1, 9)] + \
         [f"f{i}" for i in range(1, len(wires) - 24)]
wire_of = lambda y: min(range(len(wires)), key=lambda k: abs(wires[k] - y))
cuts = sorted(set(round(x, 1) for x, _, _ in dashes))
print(f"{len(wires)} wires, {len(stems)} stems, {len(dots)} dots, {len(rings)} rings, {len(cuts)} step separators")

gates, seen = [], set()
for x, ylo, yhi in sorted(stems):
    ends = (wire_of(ylo), wire_of(yhi))
    if ends[0] == ends[1] or (round(x, 1), ends) in seen:  # the bar inside a target ring, or a re-stroke
        continue
    seen.add((round(x, 1), ends))
    # a dot is filled and also stroked, so it shows up in both lists; the target end has a ring only
    d = [w for w in ends if any(abs(dx - x) < 2 and abs(dy - wires[w]) < 2 for dx, dy, _ in dots)]
    r = [w for w in ends if w not in d and any(abs(rx - x) < 2 and abs(ry - wires[w]) < 3 for rx, ry, _ in rings)]
    if len(d) == 1 and len(r) == 1:
        step = sum(c < x for c in cuts)
        gates.append((step, x, labels[d[0]], labels[r[0]]))
    else:
        print(f"  unresolved stem at x={x:.1f} wires {[labels[w] for w in ends]} dots {d} rings {r}")
for step, x, c, t in gates:
    print(f"step {step:2d}  x={x:7.1f}  CX {c:>3} -> {t}")
print(f"{len(gates)} CNOTs")
