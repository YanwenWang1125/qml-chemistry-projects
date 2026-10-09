"""Step 4: ADAPT-VQE. The circuit grows by one excitation gate per iteration, chosen by its gradient.

The operator pool is the set of gates that UCCSD uses, so the comparison with UCCSD is
"all of them once" against "only the ones selected".

Usage (from code/):
    python step04_adapt_vqe.py --molecule H2             # check: one gate should be enough
    python step04_adapt_vqe.py --molecule LiH --quick    # 1 bond length, 3 gates, 3 steps: checks and timings
    python step04_adapt_vqe.py --molecule LiH            # full experiment; rerunning skips finished bond lengths
    python step04_adapt_vqe.py --molecule LiH --plot-only
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pennylane as qml
from pennylane import numpy as pnp

from vqe_common import (
    CHEMICAL_ACCURACY_HA,
    beh2_geometry,
    build_hamiltonian,
    exact_energy,
    h2_geometry,
    hf_energy,
    lih_geometry,
    minimise_energy,
    tee_stdout,
)

# __file__ is undefined in Jupyter; there the working directory is code/ (needed for the vqe_common import).
try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

# Inner loop: the same budget as steps 2 and 3.
STEP_SIZE = 0.05
MAX_STEPS = 300
CONV_TOL = 1e-6
PATIENCE = 5

# Outer loop: stop when no gate in the pool has a gradient above GRAD_TOL (in Ha per radian),
# or when the circuit has MAX_GATES gates. Neither rule uses the exact energy.
GRAD_TOL = 1e-3
MAX_GATES = 40

# Bond lengths in Angstrom; "reference" is the one used by --quick.
MOLECULES = {
    "H2": {"geometry": h2_geometry, "bond_lengths": [0.74], "reference": 0.74},
    "LiH": {"geometry": lih_geometry, "bond_lengths": [1.2, 1.6, 2.0, 2.6, 3.2], "reference": 1.6},
    "BeH2": {"geometry": beh2_geometry, "bond_lengths": [1.0, 1.3, 1.7, 2.2, 2.8], "reference": 1.3},
}

COLUMNS = [
    "molecule", "bond_length_A", "iteration", "gate", "excitation", "max_gradient", "n_params", "steps",
    "e_hf", "e_exact", "e_vqe", "error_mHa", "seconds", "stop_reason",
]


def build_pool(n_electrons, n_qubits):
    """The single and double excitations of UCCSD, as (kind, spin orbitals, wires) entries.

    "spin orbitals" names the excitation, e.g. [0, 1, 4, 5]; "wires" is what the gate acts on.
    """
    singles, doubles = qml.qchem.excitations(n_electrons, n_qubits)
    s_wires, d_wires = qml.qchem.excitations_to_wires(singles, doubles)
    pool = [("D", excitation, wires) for excitation, wires in zip(doubles, d_wires)]
    pool += [("S", excitation, wires) for excitation, wires in zip(singles, s_wires)]
    return pool


def apply_gate(gate, angle):
    """The same gates that qml.UCCSD is built from. At angle 0 they do nothing."""
    kind, _, wires = gate
    if kind == "D":
        qml.FermionicDoubleExcitation(angle, wires1=wires[0], wires2=wires[1])
    else:
        qml.FermionicSingleExcitation(angle, wires=wires)


def run_adapt(hamiltonian, n_qubits, n_electrons, e_hf, e_exact, max_gates, max_steps):
    """Grow the circuit one gate at a time. Returns one record per iteration and the stop reason.

    Record 0 is the empty circuit (the HF state). e_exact is used only to report errors.
    """
    pool = build_pool(n_electrons, n_qubits)
    hf_state = qml.qchem.hf_state(n_electrons, n_qubits)
    dev = qml.device("lightning.qubit", wires=n_qubits)
    selected = []  # indices into pool, in the order the gates were added

    @qml.qnode(dev, diff_method="adjoint")
    def energy(theta):
        qml.BasisState(hf_state, wires=range(n_qubits))
        for position, index in enumerate(selected):
            apply_gate(pool[index], theta[position])
        return qml.expval(hamiltonian)

    @qml.qnode(dev, diff_method="adjoint")
    def energy_with_candidates(theta, candidates):
        # The current circuit followed by every pool gate. With all candidate angles at 0 the
        # state is unchanged, and the derivative with respect to candidate k is the gradient
        # that gate k would have if it alone were appended at angle 0.
        qml.BasisState(hf_state, wires=range(n_qubits))
        for position, index in enumerate(selected):
            apply_gate(pool[index], theta[position])
        for k, gate in enumerate(pool):
            apply_gate(gate, candidates[k])
        return qml.expval(hamiltonian)

    score = qml.grad(energy_with_candidates, argnums=1)

    theta = np.zeros(0)
    records = [{"iteration": 0, "gate": "", "excitation": "", "max_gradient": "", "n_params": 0, "steps": 0,
                "e_vqe": e_hf, "seconds": 0.0, "history": [e_hf]}]
    print(f"  pool: {len(pool)} gates ({sum(g[0] == 'D' for g in pool)} doubles, {sum(g[0] == 'S' for g in pool)} singles)")
    print("  iter  gate  excitation        max |grad|  params  steps  E_VQE        error (mHa)  seconds", flush=True)
    print(f"  {0:4d}  {'HF':4s}  {'':16s}  {'':10s}  {0:6d}  {0:5d}  {e_hf:+.6f}  {(e_hf - e_exact) * 1000:11.4f}  {0.0:7.1f}",
          flush=True)

    stop_reason = "max_gates"
    for iteration in range(1, max_gates + 1):
        start = time.perf_counter()

        # 1. Score every gate in the pool. No parameter is updated here.
        gradients = np.abs(np.asarray(score(pnp.array(theta, requires_grad=False),
                                            pnp.zeros(len(pool), requires_grad=True)), dtype=float))
        best = int(np.argmax(gradients))

        # 2. Stop when even the best gate has almost no slope.
        if gradients[best] < GRAD_TOL:
            stop_reason = "gradient"
            print(f"  stop: largest gradient {gradients[best]:.2e} is below {GRAD_TOL:.0e}", flush=True)
            break

        # 3. Append the best gate with its angle at 0: the energy does not change at this moment.
        selected.append(best)
        theta = pnp.array(np.append(np.asarray(theta, dtype=float), 0.0), requires_grad=True)

        # 4. Train all angles, old and new.
        history, theta = minimise_energy(energy, theta, STEP_SIZE, max_steps, CONV_TOL, PATIENCE)

        # 5. Record.
        kind, excitation, _ = pool[best]
        elapsed = time.perf_counter() - start
        records.append({"iteration": iteration, "gate": kind, "excitation": " ".join(str(i) for i in excitation),
                        "max_gradient": float(gradients[best]), "n_params": len(selected), "steps": len(history) - 1,
                        "e_vqe": history[-1], "seconds": elapsed, "history": history})
        print(f"  {iteration:4d}  {kind:4s}  {str(list(excitation)):16s}  {gradients[best]:10.6f}  {len(selected):6d}  "
              f"{len(history) - 1:5d}  {history[-1]:+.6f}  {(history[-1] - e_exact) * 1000:11.4f}  {elapsed:7.1f}", flush=True)
    else:
        print(f"  stop: reached the limit of {max_gates} gates", flush=True)

    return records, stop_reason


def read_rows(csv_path):
    if not csv_path.exists():
        return []
    with open(csv_path, newline="") as f:
        return list(csv.DictReader(f))


def run_experiment(molecule, bond_lengths, max_gates, max_steps, csv_path, history_path):
    done = {f"{float(row['bond_length_A']):.2f}" for row in read_rows(csv_path)}
    if not csv_path.exists():
        with open(csv_path, "w", newline="") as f:
            csv.writer(f).writerow(COLUMNS)

    for r in bond_lengths:
        if f"{r:.2f}" in done:
            print(f"R = {r:.2f} A: already in {csv_path.name}, skipped")
            continue
        symbols, coordinates = MOLECULES[molecule]["geometry"](r)
        hamiltonian, n_qubits, n_electrons = build_hamiltonian(symbols, coordinates)
        e_hf = hf_energy(hamiltonian, n_qubits, n_electrons)
        e_exact = exact_energy(hamiltonian, n_qubits, n_electrons)
        print(f"\n{molecule} at R = {r:.2f} A  ({n_qubits} qubits, {n_electrons} electrons)")
        print(f"  E_HF = {e_hf:+.6f}   E_exact = {e_exact:+.6f} Ha")

        records, stop_reason = run_adapt(hamiltonian, n_qubits, n_electrons, e_hf, e_exact, max_gates, max_steps)

        # A bond length is written only when it is finished: a restarted job redoes an interrupted one.
        with open(csv_path, "a", newline="") as f:
            writer = csv.writer(f)
            for record in records:
                last = record is records[-1]
                writer.writerow([molecule, r, record["iteration"], record["gate"], record["excitation"],
                                 record["max_gradient"], record["n_params"], record["steps"], e_hf, e_exact,
                                 record["e_vqe"], (record["e_vqe"] - e_exact) * 1000, record["seconds"],
                                 stop_reason if last else ""])
        with open(history_path, "a") as f:
            for record in records:
                f.write(json.dumps({"bond_length_A": r, "iteration": record["iteration"], "e_exact": e_exact,
                                    "history": record["history"]}) + "\n")


def uccsd_rows(molecule):
    """UCCSD results of step 3 for the same molecule, if that experiment has been run."""
    return [row for row in read_rows(RESULTS / f"{molecule.lower()}_ansatz.csv") if row["ansatz"] == "UCCSD"]


def summarise(molecule, rows):
    uccsd = {f"{float(row['bond_length_A']):.2f}": row for row in uccsd_rows(molecule)}
    print("\nSummary: parameters needed for chemical accuracy (1.6 mHa)")
    print("  R (A)  first within 1.6 mHa  final params  final error (mHa)  stop       UCCSD params  UCCSD error (mHa)")
    for r in sorted({float(row["bond_length_A"]) for row in rows}):
        selected = sorted((row for row in rows if float(row["bond_length_A"]) == r), key=lambda row: int(row["iteration"]))
        within = [int(row["n_params"]) for row in selected if abs(float(row["error_mHa"])) < CHEMICAL_ACCURACY_HA * 1000]
        first = f"{within[0]} params" if within else "not reached"
        last = selected[-1]
        reference = uccsd.get(f"{r:.2f}")
        uccsd_text = f"{int(reference['n_params']):12d}  {float(reference['error_mHa']):17.4f}" if reference else "           -                  -"
        print(f"  {r:5.2f}  {first:20s}  {int(last['n_params']):12d}  {float(last['error_mHa']):17.4f}  "
              f"{last['stop_reason']:9s}  {uccsd_text}")


def plot_error_vs_params(molecule, rows, figure_path):
    """Error after each added gate. Stars are the UCCSD results of step 3 at the same bond lengths."""
    uccsd = {f"{float(row['bond_length_A']):.2f}": row for row in uccsd_rows(molecule)}
    plt.figure(figsize=(7, 3.8))
    for i, r in enumerate(sorted({float(row["bond_length_A"]) for row in rows})):
        selected = sorted((row for row in rows if float(row["bond_length_A"]) == r), key=lambda row: int(row["iteration"]))
        n_params = [int(row["n_params"]) for row in selected]
        errors = [abs(float(row["error_mHa"])) + 1e-9 for row in selected]
        plt.semilogy(n_params, errors, marker="o", markersize=3, color=f"C{i}", label=f"ADAPT, {r:.2f} A")
        reference = uccsd.get(f"{r:.2f}")
        if reference:
            plt.semilogy([int(reference["n_params"])], [abs(float(reference["error_mHa"])) + 1e-9], marker="*",
                         markersize=9, linestyle="none", color=f"C{i}", label=f"UCCSD, {r:.2f} A")
    plt.axhline(CHEMICAL_ACCURACY_HA * 1000, color="gray", linestyle="--", label="1.6 mHa")
    plt.xlabel("Number of parameters")
    plt.ylabel("|Error| vs exact (mHa)")
    plt.title(f"{molecule}, STO-3G: ADAPT-VQE, error after each added gate")
    plt.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    plt.tight_layout()
    plt.savefig(figure_path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--molecule", choices=list(MOLECULES), default="LiH")
    parser.add_argument("--quick", action="store_true", help="1 bond length, 3 gates, 3 steps; writes *_quick files")
    parser.add_argument("--plot-only", action="store_true", help="skip the runs and redraw from the CSV")
    parser.add_argument("--bond-lengths", type=float, nargs="+", help="in Angstrom; default: the molecule's list")
    parser.add_argument("--max-gates", type=int, default=MAX_GATES)
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    molecule = args.molecule
    settings = MOLECULES[molecule]
    name = molecule.lower() + ("_quick" if args.quick else "")
    csv_path = RESULTS / f"{name}_adapt.csv"
    history_path = RESULTS / f"{name}_adapt_histories.jsonl"

    if args.plot_only:
        summarise(molecule, read_rows(csv_path))
    else:
        # The console output of every run is appended to this file, so it can be committed with the results.
        with tee_stdout(RESULTS / f"step04_{name}_output.txt"):
            print(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')}  step04_adapt_vqe.py {' '.join(sys.argv[1:])}")
            print(f"PennyLane {qml.__version__}; Adam, stepsize {STEP_SIZE}, conv_tol {CONV_TOL}, patience {PATIENCE}; "
                  f"grad_tol {GRAD_TOL}, max_gates {args.max_gates}")
            if args.quick:
                run_experiment(molecule, [settings["reference"]], 3, 3, csv_path, history_path)
            else:
                bond_lengths = args.bond_lengths or settings["bond_lengths"]
                run_experiment(molecule, bond_lengths, args.max_gates, MAX_STEPS, csv_path, history_path)
            summarise(molecule, read_rows(csv_path))

    plot_error_vs_params(molecule, read_rows(csv_path), FIGURES / f"{name}_adapt_error_vs_params.png")


if __name__ == "__main__":
    main()
