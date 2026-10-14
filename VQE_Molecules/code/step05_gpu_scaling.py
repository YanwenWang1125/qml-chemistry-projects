"""Step 5: statevector simulation time on CPU and GPU against the number of qubits.

Linear hydrogen chains H2, H4, ..., H12 in STO-3G give 4 to 24 qubits. For each chain the script
times one energy evaluation and one adjoint gradient of a fixed hardware-efficient circuit (the
cost of one optimiser step), on lightning.qubit and on lightning.gpu. Nothing is trained.

The GPU device needs a GPU node and:  pip install pennylane-lightning-gpu==0.45.0
A device that cannot be created is reported and skipped, so the script also runs on a CPU-only machine.

Usage (from code/):
    python step05_gpu_scaling.py --quick                 # H2 and H4, 1 repeat: checks the script
    python step05_gpu_scaling.py                         # H2 to H12; rerunning skips rows already in the CSV
    python step05_gpu_scaling.py --atoms 2 4 6 8 10      # stop at 20 qubits
    python step05_gpu_scaling.py --atoms 2 4 6 8 --devices lightning.qubit --redo   # measure these rows again
    python step05_gpu_scaling.py --atoms 12 --devices lightning.gpu --repeats 1     # 24 qubits on GPU only
    python step05_gpu_scaling.py --plot-only
"""
import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pennylane as qml
from pennylane import numpy as pnp

from vqe_common import ANGSTROM_TO_BOHR, build_hamiltonian, tee_stdout

# __file__ is undefined in Jupyter; there the working directory is code/ (needed for the vqe_common import).
try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

ATOMS = [2, 4, 6, 8, 10, 12]
SPACING_A = 1.0
LAYERS = 4
REPEATS = 3
DEVICES = ["lightning.qubit", "lightning.gpu"]

COLUMNS = ["n_atoms", "n_qubits", "n_terms", "layers", "n_params", "device", "build_seconds",
           "energy_seconds", "gradient_seconds", "repeats"]


def chain_geometry(n_atoms, spacing_angstrom):
    """n_atoms hydrogen atoms on the z axis, equally spaced. Coordinates in bohr."""
    symbols = ["H"] * n_atoms
    coordinates = np.array([[0.0, 0.0, i * spacing_angstrom * ANGSTROM_TO_BOHR] for i in range(n_atoms)])
    return symbols, coordinates


def describe_hardware():
    cores = os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count()
    print(f"CPU cores available: {cores}; OMP_NUM_THREADS = {os.environ.get('OMP_NUM_THREADS', 'not set')}")
    try:
        gpus = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
                              capture_output=True, text=True, timeout=30).stdout.strip()
        print(f"GPU: {gpus or 'none reported by nvidia-smi'}")
    except (OSError, subprocess.SubprocessError):
        print("GPU: nvidia-smi not found")


def time_device(device_name, hamiltonian, n_qubits, n_electrons, repeats):
    """Smallest time over `repeats` of one energy evaluation and of one adjoint gradient, in seconds."""
    dev = qml.device(device_name, wires=n_qubits)
    hf_state = qml.qchem.hf_state(n_electrons, n_qubits)

    @qml.qnode(dev, diff_method="adjoint")
    def energy(weights):
        qml.BasisState(hf_state, wires=range(n_qubits))
        for layer in range(LAYERS):
            for wire in range(n_qubits):
                qml.RY(weights[layer, wire, 0], wires=wire)
                qml.RZ(weights[layer, wire, 1], wires=wire)
            for wire in range(n_qubits):
                qml.CNOT(wires=[wire, (wire + 1) % n_qubits])
        return qml.expval(hamiltonian)

    # The same angles on every device and every chain length that shares a shape.
    rng = np.random.default_rng(0)
    weights = pnp.array(rng.uniform(-np.pi, np.pi, size=(LAYERS, n_qubits, 2)), requires_grad=True)
    gradient = qml.grad(energy)

    # One untimed call of each, so that set-up and compilation are not counted.
    energy(weights)
    gradient(weights)

    energy_times, gradient_times = [], []
    for _ in range(repeats):
        start = time.perf_counter()
        energy(weights)
        energy_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        gradient(weights)
        gradient_times.append(time.perf_counter() - start)
    return min(energy_times), min(gradient_times), weights.size


def read_rows(csv_path):
    if not csv_path.exists():
        return []
    with open(csv_path, newline="") as f:
        return list(csv.DictReader(f))


def write_row(csv_path, row, replace):
    """Appends one measurement. With replace, an earlier row for the same chain and device is dropped first."""
    if replace:
        kept = [r for r in read_rows(csv_path) if (int(r["n_atoms"]), r["device"]) != (row[0], row[5])]
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(COLUMNS)
            writer.writerows([[r[c] for c in COLUMNS] for r in kept])
    with open(csv_path, "a", newline="") as f:
        csv.writer(f).writerow(row)


def run_experiment(atoms, devices, repeats, csv_path, redo=False):
    if not csv_path.exists():
        with open(csv_path, "w", newline="") as f:
            csv.writer(f).writerow(COLUMNS)

    unavailable = set()
    print("  atoms  qubits  terms   device           build (s)  energy (s)  gradient (s)")
    for n_atoms in atoms:
        # With redo, rows already in the CSV are measured again and replaced once the new value exists.
        done = set() if redo else {(int(row["n_atoms"]), row["device"]) for row in read_rows(csv_path)}
        todo = [d for d in devices if (n_atoms, d) not in done and d not in unavailable]
        if not todo:
            print(f"  {n_atoms:5d}  all devices already in {csv_path.name} or unavailable, skipped")
            continue

        start = time.perf_counter()
        hamiltonian, n_qubits, n_electrons = build_hamiltonian(*chain_geometry(n_atoms, SPACING_A))
        build_seconds = time.perf_counter() - start
        n_terms = len(hamiltonian.terms()[0])

        for device_name in todo:
            try:
                energy_seconds, gradient_seconds, n_params = time_device(
                    device_name, hamiltonian, n_qubits, n_electrons, repeats)
            except Exception as error:  # a missing plugin or no GPU: report it and go on with the other device
                unavailable.add(device_name)
                print(f"  {device_name} is not available, skipped from here on: {type(error).__name__}: {error}")
                continue
            print(f"  {n_atoms:5d}  {n_qubits:6d}  {n_terms:6d}  {device_name:15s}  {build_seconds:9.1f}  "
                  f"{energy_seconds:10.4f}  {gradient_seconds:12.4f}", flush=True)
            # Written after every measurement, so a job that hits its time limit keeps what it finished.
            write_row(csv_path, [n_atoms, n_qubits, n_terms, LAYERS, n_params, device_name, build_seconds,
                                 energy_seconds, gradient_seconds, repeats], replace=redo)


def summarise(rows):
    by_qubits = {}
    for row in rows:
        by_qubits.setdefault(int(row["n_qubits"]), {})[row["device"]] = float(row["gradient_seconds"])
    print("\nSummary: seconds for one adjoint gradient, and CPU time divided by GPU time")
    print("  qubits  lightning.qubit  lightning.gpu  CPU / GPU")
    for n_qubits in sorted(by_qubits):
        cpu = by_qubits[n_qubits].get("lightning.qubit")
        gpu = by_qubits[n_qubits].get("lightning.gpu")
        cpu_text = f"{cpu:15.4f}" if cpu is not None else f"{'-':>15s}"
        gpu_text = f"{gpu:13.4f}" if gpu is not None else f"{'-':>13s}"
        ratio = f"{cpu / gpu:9.2f}" if cpu is not None and gpu is not None else f"{'-':>9s}"
        print(f"  {n_qubits:6d}  {cpu_text}  {gpu_text}  {ratio}")


def plot_scaling(rows, figure_path):
    if not rows:
        return
    plt.figure(figsize=(6, 3.8))
    for device_name in sorted({row["device"] for row in rows}):
        selected = sorted((int(row["n_qubits"]), float(row["gradient_seconds"]))
                          for row in rows if row["device"] == device_name)
        plt.semilogy([q for q, _ in selected], [t for _, t in selected], marker="o", markersize=4, label=device_name)
    plt.xlabel("Qubits (hydrogen chain, STO-3G)")
    plt.ylabel("Seconds per adjoint gradient")
    plt.title(f"HEA, {LAYERS} layers: one gradient of <H>")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figure_path, dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quick", action="store_true", help="H2 and H4, 1 repeat; writes *_quick files")
    parser.add_argument("--plot-only", action="store_true", help="skip the timings and redraw from the CSV")
    parser.add_argument("--atoms", type=int, nargs="+", help="even chain lengths; default 2 4 6 8 10 12")
    parser.add_argument("--devices", nargs="+", default=DEVICES)
    parser.add_argument("--repeats", type=int, default=REPEATS, help="timed repeats per measurement; default 3")
    parser.add_argument("--redo", action="store_true", help="measure again rows already in the CSV and replace them")
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    name = "gpu_scaling" + ("_quick" if args.quick else "")
    csv_path = RESULTS / f"{name}.csv"

    if not args.plot_only:
        with tee_stdout(RESULTS / f"step05_{name}_output.txt"):
            print(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')}  step05_gpu_scaling.py {' '.join(sys.argv[1:])}")
            print(f"PennyLane {qml.__version__}; HEA with {LAYERS} layers; H-H spacing {SPACING_A} A; "
                  f"times are the smallest of the repeats and exclude building the Hamiltonian")
            describe_hardware()
            if args.quick:
                run_experiment([2, 4], args.devices, 1, csv_path)
            else:
                run_experiment(args.atoms or ATOMS, args.devices, args.repeats, csv_path, args.redo)
            summarise(read_rows(csv_path))
    else:
        summarise(read_rows(csv_path))

    plot_scaling(read_rows(csv_path), FIGURES / f"{name}.png")


if __name__ == "__main__":
    main()
