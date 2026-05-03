# Quantum Tunnelling

Python project for simulating quantum tunneling through 1D potential barriers and connecting transmission behavior to an approximate tunnel-diode-style I-V curve.

## What This Project Covers

- Single rectangular barrier transmission
- Barrier height and width sweeps
- Wavefunction profile through the barrier
- Transmission vs energy
- Transmission vs barrier width
- Approximate current vs voltage using a Landauer-style integral
- Double-barrier resonant tunneling with a quantum well

## Project Layout

```text
quantum-tunnelling/
  requirements.txt
  run_demo.py
  quantum_tunnelling/
    __init__.py
    constants.py
    potentials.py
    solver.py
    simulation.py
```

## Setup

```bash
cd ~/quantum-tunnelling
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run

```bash
python run_demo.py
```

The script generates figures in `outputs/`:

- `single_barrier_transmission.png`
- `single_barrier_wavefunction.png`
- `width_sweep.png`
- `iv_curve.png`
- `resonant_tunnelling.png`

## Notes

- Energies are shown in `eV`, distances in `nm`.
- The transmission calculation uses a transfer-matrix solution of the 1D time-independent Schrodinger equation on a piecewise-constant potential.
- The current model is intentionally approximate and meant for coursework visualization rather than device-grade modeling.
