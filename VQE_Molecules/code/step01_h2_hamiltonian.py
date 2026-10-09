"""Steps 0 and 1: build the H2 Hamiltonian at 0.74 Angstrom and compute the reference energies."""
import pennylane as qml

from vqe_common import build_hamiltonian, exact_energy, h2_geometry, hf_energy

BOND_LENGTH = 0.74

symbols, coordinates = h2_geometry(BOND_LENGTH)
hamiltonian, n_qubits, n_electrons = build_hamiltonian(symbols, coordinates)

print(f"PennyLane version : {qml.__version__}")
print(f"Bond length       : {BOND_LENGTH} Angstrom")
print(f"Qubits            : {n_qubits}")
print(f"Electrons         : {n_electrons}")
print(f"Pauli terms       : {len(hamiltonian.terms()[0])}")
print(f"HF state          : {qml.qchem.hf_state(n_electrons, n_qubits)}")

e_hf = hf_energy(hamiltonian, n_qubits, n_electrons)
e_exact = exact_energy(hamiltonian, n_qubits, n_electrons)
print(f"HF energy         : {e_hf:.6f} Ha")
print(f"Exact energy      : {e_exact:.6f} Ha")
print(f"Correlation energy: {(e_exact - e_hf) * 1000:.3f} mHa")

print("\nHamiltonian terms (coefficient, Pauli word):")
coeffs, ops = hamiltonian.terms()
for c, op in zip(coeffs, ops):
    print(f"  {float(c):+.6f}  {op}")
