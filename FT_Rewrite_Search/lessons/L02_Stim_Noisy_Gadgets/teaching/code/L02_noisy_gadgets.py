"""L02 verification: noisy Monte Carlo of the three gadgets with Stim.

Data start in |0>^7 and are measured in Z at the end, so X-type data errors are visible.
A shot is "malignant" when the flag stays 0 but the X error left on the data has weight >= 2
modulo the X-type stabilizers. For an FT_1 gadget this needs two faults, so the rate should
scale like p^2; without the flag one fault is enough, so it should scale like p.
Run from this folder:  C:/Users/Administrator/venvs/qml/Scripts/python.exe L02_noisy_gadgets.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                "L01_FT_and_Shor_Gadget", "teaching", "code"))
from ftlib import (BASELINE, NO_FLAG, REDUCED, STEANE_SUPPORTS, measured_qubits,  # noqa: E402
                   to_stim)

SHOTS = 1_000_000
PS = [2e-3, 5e-3, 1e-2, 2e-2, 5e-2]
rng_seed = 2026

# The 8 elements of the X-type stabilizer group, as boolean patterns on the 7 data qubits.
x_gens = np.zeros((3, 7), dtype=bool)
for i, s in enumerate(STEANE_SUPPORTS):
    x_gens[i, list(s)] = True
x_group = np.array([(b0 * x_gens[0]) ^ (b1 * x_gens[1]) ^ (b2 * x_gens[2])
                    for b0 in (0, 1) for b1 in (0, 1) for b2 in (0, 1)], dtype=bool)


def run(gadget, p, shots=SHOTS):
    c = to_stim(gadget["ops"], data_init="zero", p=p)
    s = c.compile_sampler(seed=rng_seed).sample(shots)
    meas = measured_qubits(gadget["ops"])
    n_g = len(meas)
    flag = (s[:, [meas.index(q) for q in gadget["flags"]]].any(axis=1)
            if gadget["flags"] else np.zeros(shots, bool))
    data = s[:, n_g:n_g + 7]
    # weight of the X error modulo the X stabilizers: min over the 8 coset elements
    w = np.min([(data ^ g).sum(axis=1) for g in x_group], axis=0)
    malignant = (~flag) & (w >= 2)
    return flag.mean(), malignant.sum()


print("Example: the noiseless baseline gadget as a Stim circuit (data in |0>^7)")
print(to_stim(BASELINE["ops"], data_init="zero"))
print()

results = {}
for g in (NO_FLAG, BASELINE, REDUCED):
    print("=" * 72)
    print(g["name"], f"({SHOTS} shots per point)")
    rows = []
    for p in PS:
        flag_rate, n_bad = run(g, p)
        rows.append((p, flag_rate, n_bad))
        print(f"  p = {p:.0e}: flag rate {flag_rate:.4f}, malignant shots {n_bad:7d}, "
              f"rate {n_bad / SHOTS:.2e}")
    ps = np.array([r[0] for r in rows if r[2] > 0])
    rates = np.array([r[2] / SHOTS for r in rows if r[2] > 0])
    low = ps <= 1e-2
    slope_all = np.polyfit(np.log10(ps), np.log10(rates), 1)[0]
    slope_low = np.polyfit(np.log10(ps[low]), np.log10(rates[low]), 1)[0] if low.sum() >= 2 else float("nan")
    print(f"  log-log slope, all points: {slope_all:.2f}; points with p <= 1e-2: {slope_low:.2f}")
    results[g["name"]] = rows

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for name, rows in results.items():
        pts = [(p, n / SHOTS) for p, _, n in rows if n > 0]
        ax.loglog([a for a, _ in pts], [b for _, b in pts], "o-", label=name)
    ax.set_xlabel("physical error rate p")
    ax.set_ylabel("P(flag silent and X error weight >= 2 mod S)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    out = os.path.join(os.path.dirname(__file__), "..", "figures", "L02_malignant_vs_p.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("\nfigure saved:", os.path.normpath(out))
except ImportError:
    print("matplotlib not available; no figure")
