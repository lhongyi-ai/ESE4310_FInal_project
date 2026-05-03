"""Transfer-matrix solver for 1D tunnelling problems."""

from __future__ import annotations

import cmath
from dataclasses import dataclass

import numpy as np

from .constants import ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, NM


@dataclass
class TransmissionResult:
    transmission: float
    reflection: float
    wavefunction: np.ndarray


def _wave_number(energy_ev: float, potential_ev: float, mass: float) -> complex:
    energy_joule = energy_ev * ELEMENTARY_CHARGE
    potential_joule = potential_ev * ELEMENTARY_CHARGE
    return cmath.sqrt(2.0 * mass * (energy_joule - potential_joule)) / HBAR


def solve_scattering(
    x_nm: np.ndarray,
    potential_ev: np.ndarray,
    energy_ev: float,
    mass: float = ELECTRON_MASS,
) -> TransmissionResult:
    if len(x_nm) != len(potential_ev):
        raise ValueError("x and potential arrays must have the same length")
    if len(x_nm) < 3:
        raise ValueError("Need at least 3 spatial points")

    dx_m = (x_nm[1] - x_nm[0]) * NM
    k_values = np.array([_wave_number(energy_ev, v, mass) for v in potential_ev], dtype=complex)

    total = np.identity(2, dtype=complex)
    for i in range(len(x_nm) - 1):
        k_left = k_values[i]
        k_right = k_values[i + 1]

        if abs(k_right) < 1e-18:
            k_right = 1e-18 + 0j
        if abs(k_left) < 1e-18:
            k_left = 1e-18 + 0j

        interface = 0.5 * np.array(
            [
                [1.0 + k_left / k_right, 1.0 - k_left / k_right],
                [1.0 - k_left / k_right, 1.0 + k_left / k_right],
            ],
            dtype=complex,
        )
        propagation = np.array(
            [
                [np.exp(1j * k_right * dx_m), 0.0],
                [0.0, np.exp(-1j * k_right * dx_m)],
            ],
            dtype=complex,
        )
        total = propagation @ interface @ total

    r = -total[1, 0] / total[1, 1]
    t = total[0, 0] + total[0, 1] * r
    k_left = k_values[0]
    k_right = k_values[-1]
    transmission = abs(t) ** 2 * (k_right.real / k_left.real)
    reflection = abs(r) ** 2

    states = [np.array([1.0 + 0j, r], dtype=complex)]
    state = states[0]
    for i in range(len(x_nm) - 1):
        k_left = k_values[i]
        k_right = k_values[i + 1]
        if abs(k_right) < 1e-18:
            k_right = 1e-18 + 0j
        if abs(k_left) < 1e-18:
            k_left = 1e-18 + 0j
        interface = 0.5 * np.array(
            [
                [1.0 + k_left / k_right, 1.0 - k_left / k_right],
                [1.0 - k_left / k_right, 1.0 + k_left / k_right],
            ],
            dtype=complex,
        )
        propagation = np.array(
            [
                [np.exp(1j * k_right * dx_m), 0.0],
                [0.0, np.exp(-1j * k_right * dx_m)],
            ],
            dtype=complex,
        )
        state = propagation @ interface @ state
        states.append(state)

    wavefunction = np.zeros_like(x_nm, dtype=complex)
    for i, coeffs in enumerate(states):
        wavefunction[i] = coeffs[0] + coeffs[1]

    transmission = float(np.clip(transmission.real, 0.0, 1.0))
    reflection = float(max(reflection.real, 0.0))
    return TransmissionResult(
        transmission=transmission,
        reflection=reflection,
        wavefunction=wavefunction,
    )
