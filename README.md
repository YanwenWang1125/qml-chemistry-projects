# qml-chemistry-projects

Small quantum-chemistry experiments on simulated quantum circuits, written in PennyLane. Each project has its own folder with code, raw results and figures.

## Projects

| Project | Question | Status |
|---|---|---|
| [VQE_Molecules](VQE_Molecules/) | How close do UCCSD, hardware-efficient and ADAPT-VQE ansätze get to the exact ground-state energy of small molecules, and at what circuit cost? | In progress: H2 done; LiH, BeH2, ADAPT-VQE and GPU scaling planned |

## Setup

```
python -m venv .venv
.venv\Scripts\activate        # Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

The results in this repository were produced with PennyLane 0.45.1, pennylane-lightning 0.45.0, NumPy 2.5.3 and SciPy 1.18.1 on Windows. Molecular Hamiltonians come from `qml.qchem`, which does not need PySCF.

Each script writes its CSV files and figures into the `results/` and `figures/` folders of its project. See the project README for the commands.
