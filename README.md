# Quantum Tunnelling

Python project for simulating quantum tunneling through 1D potential barriers and connecting resonant-tunneling behavior to an approximate free-electron-mass RTD baseline I-V curve.

## What This Project Covers

- Single rectangular barrier transmission
- Barrier height and width sweeps
- Wavefunction profile through the barrier
- Transmission vs energy
- Transmission vs barrier width
- Temperature-dependent RTD baseline current density with negative differential resistance
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
- The RTD current-density baseline uses the free electron mass. For each applied bias, the double-barrier potential is tilted, the biased transmission spectrum `T(E, V)` is recomputed, and the current is integrated with a Tsu-Esaki-like finite-temperature supply term.
- The earlier resonance-shift I-V curve was only a qualitative NDR illustration; the current baseline is still approximate, but it is a better physical baseline for comparing against material-specific effective-mass models.
