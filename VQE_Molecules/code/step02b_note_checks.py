"""Step 2b: checks on the H2 Hamiltonian at 0.74 Angstrom.

Matrix elements between |1100> and |0011>, the ground-state amplitudes, and the energy
reached by a single DoubleExcitation gate.
"""
import numpy as np
import pennylane as qml
from scipy.optimize import minimize_scalar

from vqe_common import build_hamiltonian, exact_energy, h2_geometry

symbols, coordinates = h2_geometry(0.74)
hamiltonian, n_qubits, n_electrons = build_hamiltonian(symbols, coordinates)
matrix = hamiltonian.sparse_matrix(wire_order=range(n_qubits)).toarray().real

singles, doubles = qml.qchem.excitations(n_electrons, n_qubits)
print(f"singles: {singles}")
print(f"doubles: {doubles}")

print(f"<1100|H|1100> = {matrix[0b1100, 0b1100]:.6f} Ha")
print(f"<0011|H|0011> = {matrix[0b0011, 0b0011]:.6f} Ha")
print(f"<1100|H|0011> = {matrix[0b1100, 0b0011]:.6f} Ha")

e_sector = exact_energy(hamiltonian, n_qubits, n_electrons)
values, vectors = np.linalg.eigh(matrix)
print(f"lowest eigenvalue, 2-electron sector: {e_sector:.6f} Ha")
print(f"lowest eigenvalue, all 16 states    : {values[0]:.6f} Ha")

ground = vectors[:, 0]
print("ground-state amplitudes above 1e-6:")
for index in np.flatnonzero(np.abs(ground) > 1e-6):
    print(f"  |{index:04b}>  {ground[index]:+.6f}")

dev = qml.device("default.qubit", wires=n_qubits)


@qml.qnode(dev)
def energy(theta):
    qml.BasisState(np.array([1, 1, 0, 0]), wires=range(n_qubits))
    qml.DoubleExcitation(theta, wires=[0, 1, 2, 3])
    return qml.expval(hamiltonian)


best = minimize_scalar(lambda t: float(energy(t)), bounds=(-np.pi, np.pi), method="bounded")
print(f"single DoubleExcitation: theta = {best.x:.6f}, energy = {best.fun:.6f} Ha, "
      f"error = {(best.fun - e_sector) * 1000:.6f} mHa")
