"""Shared helpers for the VQE_Molecules project: Hamiltonians and reference energies."""
import contextlib
import sys

import numpy as np
import pennylane as qml
from scipy.sparse.linalg import eigsh

ANGSTROM_TO_BOHR = 1.8897259886
CHEMICAL_ACCURACY_HA = 1.6e-3


def h2_geometry(bond_length_angstrom):
    """H2 along the z axis. Coordinates are returned in bohr, the unit qml.qchem expects."""
    half = 0.5 * bond_length_angstrom * ANGSTROM_TO_BOHR
    symbols = ["H", "H"]
    coordinates = np.array([[0.0, 0.0, -half], [0.0, 0.0, half]])
    return symbols, coordinates


def lih_geometry(bond_length_angstrom):
    """LiH along the z axis, Li at the origin. Coordinates are returned in bohr."""
    symbols = ["Li", "H"]
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, bond_length_angstrom * ANGSTROM_TO_BOHR]])
    return symbols, coordinates


def beh2_geometry(bond_length_angstrom):
    """Linear H-Be-H along the z axis, Be at the origin, both Be-H bonds of the given length. In bohr."""
    d = bond_length_angstrom * ANGSTROM_TO_BOHR
    symbols = ["H", "Be", "H"]
    coordinates = np.array([[0.0, 0.0, -d], [0.0, 0.0, 0.0], [0.0, 0.0, d]])
    return symbols, coordinates


def build_hamiltonian(symbols, coordinates, basis="sto-3g"):
    """Qubit Hamiltonian (Jordan-Wigner) in Hartree, with the number of qubits and electrons."""
    molecule = qml.qchem.Molecule(symbols, coordinates, basis_name=basis)
    hamiltonian, n_qubits = qml.qchem.molecular_hamiltonian(molecule)
    return hamiltonian, n_qubits, molecule.n_electrons


def hf_energy(hamiltonian, n_qubits, n_electrons):
    """<HF|H|HF>, read off the diagonal of the Hamiltonian matrix."""
    hf_state = qml.qchem.hf_state(n_electrons, n_qubits)
    # Wire 0 is the most significant bit of the basis-state index.
    index = int("".join(str(int(b)) for b in hf_state), 2)
    matrix = hamiltonian.sparse_matrix(wire_order=range(n_qubits)).tocsr()
    return float(np.real(matrix[index, index]))


def exact_energy(hamiltonian, n_qubits, n_electrons):
    """Lowest eigenvalue among basis states with n_electrons occupied spin orbitals.

    This equals the FCI energy in the same basis set.
    """
    matrix = hamiltonian.sparse_matrix(wire_order=range(n_qubits)).tocsr()
    sector = [i for i in range(2**n_qubits) if bin(i).count("1") == n_electrons]
    block = matrix[sector, :][:, sector]
    if block.shape[0] <= 512:
        return float(np.linalg.eigvalsh(block.toarray())[0])
    return float(eigsh(block, k=1, which="SA", return_eigenvectors=False)[0])


def lowest_energy(hamiltonian, n_qubits):
    """Lowest eigenvalue over all basis states, with no restriction on the number of electrons.

    A circuit that does not conserve the electron number (such as a hardware-efficient ansatz)
    is bounded below by this value, not by exact_energy.
    """
    matrix = hamiltonian.sparse_matrix(wire_order=range(n_qubits)).tocsr()
    return float(eigsh(matrix, k=1, which="SA", return_eigenvectors=False)[0])


class _Tee:
    """Writes to several streams at once."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


@contextlib.contextmanager
def tee_stdout(path):
    """Inside the block, everything printed also goes to the end of the file at `path`."""
    with open(path, "a", encoding="utf-8") as f:
        with contextlib.redirect_stdout(_Tee(sys.stdout, f)):
            yield


def minimise_energy(energy, weights, step_size, max_steps, conv_tol, patience=1):
    """Adam on a QNode. Stops after `patience` consecutive steps that change the energy by less than conv_tol.

    Returns the energy history (including the starting point) and the final weights.
    """
    optimizer = qml.AdamOptimizer(stepsize=step_size)
    history = [float(energy(weights))]
    small_steps = 0
    for _ in range(max_steps):
        weights, previous = optimizer.step_and_cost(energy, weights)
        current = float(energy(weights))
        history.append(current)
        small_steps = small_steps + 1 if abs(current - float(previous)) < conv_tol else 0
        if small_steps >= patience:
            break
    return history, weights
