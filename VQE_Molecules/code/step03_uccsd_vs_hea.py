"""Step 3: UCCSD against a hardware-efficient ansatz (HEA) on LiH or BeH2, under one optimisation budget.

Usage (from code/):
    python step03_uccsd_vs_hea.py --molecule BeH2 --quick      # 1 bond length, 3 steps, 1 seed: checks and timings
    python step03_uccsd_vs_hea.py --molecule BeH2              # full experiment; rerunning skips runs in the CSV
    python step03_uccsd_vs_hea.py --molecule LiH --plot-only   # redraw the figures from the CSV
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
    build_hamiltonian,
    beh2_geometry,
    exact_energy,
    hf_energy,
    lih_geometry,
    lowest_energy,
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

# Same budget as step02_h2_vqe.py, shared by every run.
STEP_SIZE = 0.05
MAX_STEPS = 300
CONV_TOL = 1e-6
PATIENCE = 5

# Half-width of the initial angle range for the HEA that starts from the HF state.
HF_START_ANGLE = 0.1

# Bond lengths in Angstrom. "reference" is close to the equilibrium bond length; it is the one
# used by --quick and by the convergence figure.
MOLECULES = {
    "LiH": {"geometry": lih_geometry, "bond_lengths": [1.2, 1.6, 2.0, 2.6, 3.2], "reference": 1.6, "bond": "Li-H"},
    "BeH2": {"geometry": beh2_geometry, "bond_lengths": [1.0, 1.3, 1.7, 2.2, 2.8], "reference": 1.3, "bond": "Be-H"},
}
HEA_LAYERS = [2, 4, 6]
SEEDS = [0, 1, 2, 3, 4]

COLUMNS = [
    "molecule", "bond_length_A", "ansatz", "layers", "seed", "n_params", "steps",
    "e_hf", "e_exact", "e_lowest_any_sector", "e_vqe", "error_mHa", "seconds",
]


def run_uccsd(hamiltonian, n_qubits, n_electrons, max_steps):
    """Minimise <H> over UCCSD parameters, starting from the HF state (all parameters 0)."""
    singles, doubles = qml.qchem.excitations(n_electrons, n_qubits)
    s_wires, d_wires = qml.qchem.excitations_to_wires(singles, doubles)
    hf_state = qml.qchem.hf_state(n_electrons, n_qubits)
    dev = qml.device("lightning.qubit", wires=n_qubits)

    @qml.qnode(dev, diff_method="adjoint")
    def energy(weights):
        qml.UCCSD(weights, wires=range(n_qubits), s_wires=s_wires, d_wires=d_wires, init_state=hf_state)
        return qml.expval(hamiltonian)

    weights = pnp.zeros(len(singles) + len(doubles), requires_grad=True)
    history, weights = minimise_energy(energy, weights, STEP_SIZE, max_steps, CONV_TOL, PATIENCE)
    return history, weights.size


def hf_preimage(n_electrons, n_qubits, n_layers):
    """Basis state that the CNOT rings of an n_layers HEA map onto the HF state.

    With every angle at 0 the rotations are identities and the circuit is only CNOTs, which send
    one basis state to another. Undoing those CNOTs on the HF bit string, last gate first, gives
    the state to start from so that the circuit with all angles at 0 outputs the HF state.
    """
    bits = [int(b) for b in qml.qchem.hf_state(n_electrons, n_qubits)]
    for _ in range(n_layers):
        for wire in reversed(range(n_qubits)):
            bits[(wire + 1) % n_qubits] ^= bits[wire]
    return np.array(bits)


def run_hea(hamiltonian, n_qubits, n_layers, seed, max_steps, start_bits=None, angle_range=np.pi):
    """Minimise <H> over a layered circuit: RY and RZ on every qubit, then a ring of CNOTs.

    Starts from |0...0>, or from the basis state start_bits, with angles drawn uniformly
    from [-angle_range, angle_range).
    """
    dev = qml.device("lightning.qubit", wires=n_qubits)

    @qml.qnode(dev, diff_method="adjoint")
    def energy(weights):
        if start_bits is not None:
            qml.BasisState(start_bits, wires=range(n_qubits))
        for layer in range(n_layers):
            for wire in range(n_qubits):
                qml.RY(weights[layer, wire, 0], wires=wire)
                qml.RZ(weights[layer, wire, 1], wires=wire)
            for wire in range(n_qubits):
                qml.CNOT(wires=[wire, (wire + 1) % n_qubits])
        return qml.expval(hamiltonian)

    rng = np.random.default_rng(seed)
    weights = pnp.array(rng.uniform(-angle_range, angle_range, size=(n_layers, n_qubits, 2)), requires_grad=True)
    history, weights = minimise_energy(energy, weights, STEP_SIZE, max_steps, CONV_TOL, PATIENCE)
    return history, weights.size


def run_key(row):
    """Identifies one run, so that a restarted job can skip what is already in the CSV."""
    return (f"{float(row['bond_length_A']):.2f}", row["ansatz"], int(row["layers"]), int(row["seed"]))


def read_rows(csv_path):
    if not csv_path.exists():
        return []
    with open(csv_path, newline="") as f:
        return list(csv.DictReader(f))


def run_experiment(molecule, bond_lengths, hea_layers, seeds, max_steps, csv_path, history_path):
    done = {run_key(row) for row in read_rows(csv_path)}
    if not csv_path.exists():
        with open(csv_path, "w", newline="") as f:
            csv.writer(f).writerow(COLUMNS)

    for r in bond_lengths:
        # HEA starts from random angles on |0...0>; HEA_HF is the same circuit started next to the HF state.
        jobs = ([("UCCSD", 0, 0)]
                + [("HEA", layers, seed) for layers in hea_layers for seed in seeds]
                + [("HEA_HF", layers, seed) for layers in hea_layers for seed in seeds])
        jobs = [job for job in jobs if (f"{r:.2f}", *job) not in done]
        if not jobs:
            print(f"R = {r:.2f} A: all runs already in {csv_path.name}, skipped")
            continue

        start = time.perf_counter()
        symbols, coordinates = MOLECULES[molecule]["geometry"](r)
        hamiltonian, n_qubits, n_electrons = build_hamiltonian(symbols, coordinates)
        e_hf = hf_energy(hamiltonian, n_qubits, n_electrons)
        e_exact = exact_energy(hamiltonian, n_qubits, n_electrons)
        e_lowest = lowest_energy(hamiltonian, n_qubits)
        print(f"\n{molecule} at R = {r:.2f} A  ({n_qubits} qubits, {n_electrons} electrons, "
              f"{len(hamiltonian.terms()[0])} Pauli terms, built in {time.perf_counter() - start:.0f} s)")
        print(f"  E_HF = {e_hf:+.6f}   E_exact = {e_exact:+.6f}   lowest over all sectors = {e_lowest:+.6f} Ha")
        if e_lowest < e_exact - 1e-8:
            print("  Note: a state with a different electron number lies below E_exact; HEA can go below it.")
        print("  ansatz  layers  seed  params  steps  start (mHa)  E_VQE        error (mHa)  seconds", flush=True)

        for ansatz, layers, seed in jobs:
            start = time.perf_counter()
            if ansatz == "UCCSD":
                history, n_params = run_uccsd(hamiltonian, n_qubits, n_electrons, max_steps)
            elif ansatz == "HEA_HF":
                start_bits = hf_preimage(n_electrons, n_qubits, layers)
                history, n_params = run_hea(hamiltonian, n_qubits, layers, seed, max_steps, start_bits, HF_START_ANGLE)
            else:
                history, n_params = run_hea(hamiltonian, n_qubits, layers, seed, max_steps)
            elapsed = time.perf_counter() - start
            error_mha = (history[-1] - e_exact) * 1000
            print(f"  {ansatz:6s}  {layers:6d}  {seed:4d}  {n_params:6d}  {len(history) - 1:5d}  "
                  f"{(history[0] - e_exact) * 1000:11.4f}  {history[-1]:+.6f}  {error_mha:11.4f}  {elapsed:7.1f}", flush=True)

            # Written after every run, so a job that hits its time limit keeps what it finished.
            with open(csv_path, "a", newline="") as f:
                csv.writer(f).writerow([molecule, r, ansatz, layers, seed, n_params, len(history) - 1,
                                        e_hf, e_exact, e_lowest, history[-1], error_mha, elapsed])
            with open(history_path, "a") as f:
                f.write(json.dumps({"bond_length_A": r, "ansatz": ansatz, "layers": layers, "seed": seed,
                                    "e_exact": e_exact, "history": history}) + "\n")


ANSATZ_ORDER = ["UCCSD", "HEA", "HEA_HF"]


def label_of(ansatz, layers):
    if ansatz == "UCCSD":
        return "UCCSD"
    return f"HEA from HF, {layers} layers" if ansatz == "HEA_HF" else f"HEA, {layers} layers"


def labels_in(rows):
    """Ansatz labels in plotting order, with the (ansatz, layers) pair each one stands for."""
    pairs = sorted({(row["ansatz"], int(row["layers"])) for row in rows},
                   key=lambda p: (ANSATZ_ORDER.index(p[0]), p[1]))
    return [(label_of(ansatz, layers), ansatz, layers) for ansatz, layers in pairs]


def summarise(rows):
    print("\nSummary: |error| in mHa over seeds, and runs within chemical accuracy")
    print("  R (A)  ansatz                  runs  mean      std       min       max       within 1.6 mHa")
    for r in sorted({float(row["bond_length_A"]) for row in rows}):
        for label, ansatz, layers in labels_in(rows):
            errors = np.array([abs(float(row["error_mHa"])) for row in rows
                               if float(row["bond_length_A"]) == r and row["ansatz"] == ansatz and int(row["layers"]) == layers])
            if errors.size == 0:
                continue
            n_ok = int(np.sum(errors < CHEMICAL_ACCURACY_HA * 1000))
            print(f"  {r:5.2f}  {label:22s}  {errors.size:4d}  {errors.mean():8.3f}  {errors.std():8.3f}  "
                  f"{errors.min():8.3f}  {errors.max():8.3f}  {n_ok} of {errors.size}")


def plot_error_vs_bond_length(molecule, rows, figure_path):
    plt.figure(figsize=(7, 3.8))
    for label, ansatz, layers in labels_in(rows):
        selected = [row for row in rows if row["ansatz"] == ansatz and int(row["layers"]) == layers]
        bond_lengths = sorted({float(row["bond_length_A"]) for row in selected})
        errors = [np.array([abs(float(row["error_mHa"])) for row in selected if float(row["bond_length_A"]) == r])
                  for r in bond_lengths]
        mean = np.array([e.mean() for e in errors])
        low = mean - np.array([e.min() for e in errors])
        high = np.array([e.max() for e in errors]) - mean
        # Bars span the smallest and largest error over seeds.
        plt.errorbar(bond_lengths, mean, yerr=[low, high], marker="o", markersize=3, capsize=2, label=label)
    plt.axhline(CHEMICAL_ACCURACY_HA * 1000, color="gray", linestyle="--", label="1.6 mHa")
    plt.yscale("log")
    plt.xlabel(f"{MOLECULES[molecule]['bond']} bond length (Angstrom)")
    plt.ylabel("|Error| vs exact (mHa)")
    plt.title(f"{molecule}, STO-3G: mean over seeds, bars min to max")
    # Outside the axes: with seven series a legend inside covers the curves.
    plt.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    plt.tight_layout()
    plt.savefig(figure_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_convergence(molecule, history_path, figure_path):
    """Training curves at the bond length closest to the molecule's reference one: one line per run."""
    bond_length = MOLECULES[molecule]["reference"]
    with open(history_path) as f:
        runs = [json.loads(line) for line in f]
    r = min({run["bond_length_A"] for run in runs}, key=lambda x: abs(x - bond_length))
    runs = [run for run in runs if run["bond_length_A"] == r]
    colours = {}
    plt.figure(figsize=(7, 3.8))
    for run in runs:
        label = label_of(run["ansatz"], run["layers"])
        first = label not in colours
        colours.setdefault(label, f"C{len(colours)}")
        error = np.abs(np.array(run["history"]) - run["e_exact"]) * 1000 + 1e-9
        plt.semilogy(error, color=colours[label], linewidth=0.8, label=label if first else None)
    plt.axhline(CHEMICAL_ACCURACY_HA * 1000, color="gray", linestyle="--", label="1.6 mHa")
    plt.xlabel("Optimiser step")
    plt.ylabel("|Error| vs exact (mHa)")
    plt.title(f"{molecule}, STO-3G, {r:.2f} Angstrom")
    # Outside the axes: with seven series a legend inside covers the curves.
    plt.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    plt.tight_layout()
    plt.savefig(figure_path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--molecule", choices=list(MOLECULES), default="LiH")
    parser.add_argument("--quick", action="store_true", help="1 bond length, 3 steps, 1 seed; writes *_quick files")
    parser.add_argument("--plot-only", action="store_true", help="skip the runs and redraw from the CSV")
    parser.add_argument("--bond-lengths", type=float, nargs="+", help="in Angstrom; default: the molecule's list")
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    molecule = args.molecule
    settings = MOLECULES[molecule]
    name = molecule.lower() + ("_quick" if args.quick else "")
    csv_path = RESULTS / f"{name}_ansatz.csv"
    history_path = RESULTS / f"{name}_histories.jsonl"

    if args.plot_only:
        summarise(read_rows(csv_path))
    else:
        # The console output of every run is appended to this file, so it can be committed with the results.
        with tee_stdout(RESULTS / f"step03_{name}_output.txt"):
            print(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')}  step03_uccsd_vs_hea.py {' '.join(sys.argv[1:])}")
            print(f"PennyLane {qml.__version__}; Adam, stepsize {STEP_SIZE}, conv_tol {CONV_TOL}, patience {PATIENCE}")
            if args.quick:
                run_experiment(molecule, [settings["reference"]], HEA_LAYERS, SEEDS[:1], 3, csv_path, history_path)
            else:
                bond_lengths = args.bond_lengths or settings["bond_lengths"]
                run_experiment(molecule, bond_lengths, HEA_LAYERS, SEEDS, MAX_STEPS, csv_path, history_path)
            summarise(read_rows(csv_path))

    rows = read_rows(csv_path)
    plot_error_vs_bond_length(molecule, rows, FIGURES / f"{name}_error_vs_bond_length.png")
    plot_convergence(molecule, history_path, FIGURES / f"{name}_convergence.png")


if __name__ == "__main__":
    main()
