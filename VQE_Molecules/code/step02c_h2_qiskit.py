"""Step 2c: the H2 VQE again in Qiskit, as a cross-check of the PennyLane results of step 2.

The qubit Hamiltonian is the one PennyLane builds (same coefficients as step 1), converted to a
Qiskit SparsePauliOp. The ansatz is a hand-built circuit with one parameter and three CNOTs that
prepares cos(theta/2)|1100> + sin(theta/2)|0011>, the only two basis states in the H2 ground state.

Needs Qiskit 1.0 or later:  pip install qiskit

Usage (from code/):
    python step02c_h2_qiskit.py                       # 6 bond lengths
    python step02c_h2_qiskit.py --bond-lengths 0.74   # one geometry
"""
import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import pennylane as qml
import qiskit
from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import SparsePauliOp
from scipy.optimize import minimize_scalar

from vqe_common import CHEMICAL_ACCURACY_HA, build_hamiltonian, exact_energy, h2_geometry, hf_energy, tee_stdout

# __file__ is undefined in Jupyter; there the working directory is code/ (needed for the vqe_common import).
try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().resolve().parent
RESULTS = ROOT / "results"

BOND_LENGTHS = [0.5, 0.74, 1.0, 1.5, 2.0, 2.5]
COLUMNS = ["bond_length_A", "n_terms", "e_hf", "e_exact", "e_hf_qiskit", "theta", "e_vqe_qiskit", "error_mHa",
           "evaluations"]


def to_sparse_pauli_op(hamiltonian, n_qubits):
    """PennyLane Hamiltonian -> Qiskit SparsePauliOp, with PennyLane wire i as Qiskit qubit i."""
    terms = []
    for word, coeff in qml.pauli.pauli_sentence(hamiltonian).items():
        # A PauliWord maps wire -> "X", "Y" or "Z"; the empty word is the identity.
        terms.append(("".join(word.values()), list(word.keys()), complex(coeff)))
    return SparsePauliOp.from_sparse_list(terms, num_qubits=n_qubits).simplify()


def h2_ansatz():
    """cos(theta/2)|1100> + sin(theta/2)|0011> on qubits 0..3; theta = 0 is the Hartree-Fock state."""
    theta = Parameter("theta")
    circuit = QuantumCircuit(4)
    circuit.x(0)
    circuit.x(1)
    circuit.ry(theta, 2)
    circuit.cx(2, 3)
    circuit.cx(2, 0)
    circuit.cx(2, 1)
    return circuit


def run(bond_lengths, csv_path):
    estimator = StatevectorEstimator()
    circuit = h2_ansatz()
    rows = []
    print("  R (A)  terms  E_HF         E_HF (Qiskit)  E_exact      theta    E_VQE (Qiskit)  error (mHa)  evals")
    for r in bond_lengths:
        hamiltonian, n_qubits, n_electrons = build_hamiltonian(*h2_geometry(r))
        e_hf = hf_energy(hamiltonian, n_qubits, n_electrons)
        e_exact = exact_energy(hamiltonian, n_qubits, n_electrons)
        operator = to_sparse_pauli_op(hamiltonian, n_qubits)

        def energy(theta):
            result = estimator.run([(circuit, operator, [theta])]).result()
            return float(np.asarray(result[0].data.evs).ravel()[0])

        # The energy is a sinusoid in theta with period 2*pi, so one bounded 1-D search finds its minimum.
        best = minimize_scalar(energy, bounds=(-np.pi, np.pi), method="bounded", options={"xatol": 1e-8})
        e_hf_qiskit = energy(0.0)
        error_mha = (best.fun - e_exact) * 1000
        print(f"  {r:5.2f}  {len(operator):5d}  {e_hf:+.6f}  {e_hf_qiskit:+.6f}      {e_exact:+.6f}  "
              f"{best.x:+.4f}  {best.fun:+.6f}       {error_mha:11.4f}  {best.nfev:5d}", flush=True)
        rows.append([r, len(operator), e_hf, e_exact, e_hf_qiskit, best.x, best.fun, error_mha, best.nfev])

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(COLUMNS)
        writer.writerows(rows)

    errors = np.array([abs(row[7]) for row in rows])
    hf_gaps = np.array([abs(row[2] - row[4]) for row in rows])
    print(f"\nWithin chemical accuracy: {int(np.sum(errors < CHEMICAL_ACCURACY_HA * 1000))} of {len(rows)}; "
          f"largest error {errors.max():.4f} mHa")
    print(f"Largest |E_HF (PennyLane matrix) - E_HF (Qiskit circuit)|: {hf_gaps.max():.2e} Ha")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bond-lengths", type=float, nargs="+", default=BOND_LENGTHS, help="in Angstrom")
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    with tee_stdout(RESULTS / "step02c_output.txt"):
        print(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')}  step02c_h2_qiskit.py {' '.join(sys.argv[1:])}")
        print(f"Qiskit {qiskit.__version__}; PennyLane {qml.__version__} (Hamiltonian only)")
        print(h2_ansatz().draw(output="text"))
        run(args.bond_lengths, RESULTS / "h2_qiskit.csv")


if __name__ == "__main__":
    main()
