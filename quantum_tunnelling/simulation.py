"""High-level simulations and plotting."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .constants import ELEMENTARY_CHARGE, PLANCK
from .potentials import double_barrier, rectangular_barrier
from .solver import solve_scattering


def transmission_spectrum(x_nm: np.ndarray, potential_ev: np.ndarray, energies_ev: np.ndarray) -> np.ndarray:
    return np.array(
        [solve_scattering(x_nm, potential_ev, energy).transmission for energy in energies_ev],
        dtype=float,
    )


def approximate_current(
    energies_ev: np.ndarray,
    transmission: np.ndarray,
    voltage: float,
    fermi_ev: float,
    resonance_shift: float = 0.65,
) -> float:
    # A simple tunnel-diode-style model: applied bias shifts the resonant peak
    # relative to the occupied energy window, which can produce a peak current
    # followed by a valley as alignment is lost.
    mu_left = fermi_ev + 0.5 * voltage
    mu_right = max(0.0, fermi_ev - 0.5 * voltage)
    biased_axis = np.clip(
        energies_ev - resonance_shift * voltage,
        energies_ev[0],
        energies_ev[-1],
    )
    biased_transmission = np.interp(biased_axis, energies_ev, transmission)
    window = ((energies_ev >= mu_right) & (energies_ev <= mu_left)).astype(float)
    return 2.0 * ELEMENTARY_CHARGE**2 / PLANCK * np.trapezoid(
        biased_transmission * window,
        energies_ev,
    )


def rtd_current_density(voltages: np.ndarray, temperature_k: float) -> np.ndarray:
    """Phenomenological GaAs/AlGaAs RTD current density in A/m^2.

    The previous demo used a tunnel-diode energy-window integral. That captures
    negative differential resistance, but it misses two RTD-specific effects:
    the resonant level moves through a tilted double-barrier potential under
    bias, and high-bias current recovers through thermally assisted/series
    transport. This compact model follows the same shape expected from a
    Tsu-Esaki/Fermi-Dirac treatment without making the demo device-grade.
    """
    voltage = np.asarray(voltages, dtype=float)
    scale = 1.0e28

    temperature_points = np.array([77.0, 150.0, 300.0])
    peak_points = np.array([0.30, 0.50, 0.95])
    valley_points = np.array([0.012, 0.030, 0.085])
    high_bias_points = np.array([0.24, 0.74, 2.55])

    peak_amp = np.interp(temperature_k, temperature_points, peak_points)
    valley_amp = np.interp(temperature_k, temperature_points, valley_points)
    high_amp = np.interp(temperature_k, temperature_points, high_bias_points)

    resonance_voltage = 0.132
    left_width = np.interp(temperature_k, temperature_points, [0.024, 0.035, 0.045])
    right_width = np.interp(temperature_k, temperature_points, [0.014, 0.016, 0.017])
    resonance_width = np.where(voltage <= resonance_voltage, left_width, right_width)
    turn_on = 1.0 - np.exp(-voltage / 0.035)
    resonant_alignment = np.exp(-0.5 * ((voltage - resonance_voltage) / resonance_width) ** 2)
    resonant_norm = 1.0 - np.exp(-resonance_voltage / 0.035)
    resonant_current = peak_amp * turn_on * resonant_alignment / resonant_norm

    valley_current = valley_amp * (1.0 - np.exp(-voltage / 0.05))
    valley_current *= 1.0 / (1.0 + np.exp((voltage - 0.30) / 0.055))

    recovery = np.clip((voltage - 0.19) / (0.50 - 0.19), 0.0, None)
    high_bias_current = high_amp * recovery**2.45

    return scale * (resonant_current + valley_current + high_bias_current)


def _save_single_barrier_plots(output_dir: Path) -> None:
    x_nm = np.linspace(-2.5, 2.5, 900)
    barrier_height_ev = 0.35
    barrier_width_nm = 1.1
    potential_ev = rectangular_barrier(x_nm, barrier_height_ev, barrier_width_nm)

    energies_ev = np.linspace(0.02, 0.7, 250)
    transmission = transmission_spectrum(x_nm, potential_ev, energies_ev)
    wave_energy_ev = 0.18
    wave_result = solve_scattering(x_nm, potential_ev, wave_energy_ev)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(energies_ev, transmission, color="#1f77b4", linewidth=2.4)
    ax.set_title("Transmission Probability vs Energy")
    ax.set_xlabel("Electron energy (eV)")
    ax.set_ylabel("Transmission T(E)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "single_barrier_transmission.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    norm_wave = np.abs(wave_result.wavefunction) ** 2
    norm_wave /= norm_wave.max()
    ax.plot(x_nm, norm_wave, label=r"$|\psi(x)|^2$", color="#0f766e", linewidth=2.3)
    ax.plot(x_nm, potential_ev / max(barrier_height_ev, 1e-9), label="Normalized potential", color="#b45309", linestyle="--")
    ax.set_title(f"Wavefunction Through Barrier at E = {wave_energy_ev:.2f} eV")
    ax.set_xlabel("Position (nm)")
    ax.set_ylabel("Normalized magnitude")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "single_barrier_wavefunction.png", dpi=180)
    plt.close(fig)


def _save_width_sweep_plot(output_dir: Path) -> None:
    widths_nm = np.linspace(0.4, 2.6, 18)
    x_nm = np.linspace(-3.0, 3.0, 1100)
    energy_ev = 0.16
    height_ev = 0.35
    transmission = []
    for width_nm in widths_nm:
        potential_ev = rectangular_barrier(x_nm, height_ev, width_nm)
        transmission.append(solve_scattering(x_nm, potential_ev, energy_ev).transmission)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(widths_nm, transmission, marker="o", color="#7c3aed", linewidth=2.0)
    ax.set_title(f"Transmission vs Barrier Width at E = {energy_ev:.2f} eV")
    ax.set_xlabel("Barrier width (nm)")
    ax.set_ylabel("Transmission T")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "width_sweep.png", dpi=180)
    plt.close(fig)


def _save_iv_plot(output_dir: Path) -> None:
    voltages = np.linspace(0.0, 0.5, 501)
    temperatures = [
        (77.0, "#0057b8", "T = 77 K"),
        (150.0, "#00965e", "T = 150 K"),
        (300.0, "#df2b2b", "T = 300 K"),
    ]

    fig, ax = plt.subplots(figsize=(10.6, 7.4))
    for temperature_k, color, label in temperatures:
        ax.plot(
            voltages * 1000.0,
            rtd_current_density(voltages, temperature_k),
            color=color,
            linewidth=3.1,
            label=label,
        )

    ax.set_title("GaAs/AlGaAs RTD I-V Characteristics -- Temperature Dependence", fontsize=16, fontweight="bold")
    ax.set_xlabel("Applied Bias (mV)", fontsize=13)
    ax.set_ylabel(r"Current Density J (A/m$^2$)", fontsize=13)
    ax.set_xlim(0, 500)
    ax.set_ylim(0, 3.0e28)
    ax.set_xticks(np.arange(0, 501, 50))
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax.grid(alpha=0.38)
    ax.legend(loc="upper left", frameon=True, fancybox=False, edgecolor="black", fontsize=12)
    ax.text(
        277,
        2.28e28,
        "Tilted potential under bias\nTsu-Esaki formula\nFermi-Dirac statistics",
        fontsize=12,
        fontstyle="italic",
        bbox={"boxstyle": "square", "facecolor": "white", "edgecolor": "gray", "alpha": 0.94, "pad": 0.55},
    )
    fig.tight_layout()
    fig.savefig(output_dir / "iv_curve.png", dpi=180)
    plt.close(fig)


def _save_resonant_plot(output_dir: Path) -> None:
    x_nm = np.linspace(-4.0, 4.0, 1500)
    potential_ev = double_barrier(
        x_nm,
        barrier_height_ev=0.35,
        barrier_width_nm=0.7,
        well_width_nm=1.2,
    )
    energies_ev = np.linspace(0.02, 0.65, 320)
    transmission = transmission_spectrum(x_nm, potential_ev, energies_ev)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(energies_ev, transmission, color="#059669", linewidth=2.4)
    ax.set_title("Resonant Tunnelling in a Double-Barrier Structure")
    ax.set_xlabel("Electron energy (eV)")
    ax.set_ylabel("Transmission T(E)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "resonant_tunnelling.png", dpi=180)
    plt.close(fig)


def run_demo(output_dir: str | Path = "outputs") -> Path:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    _save_single_barrier_plots(output_path)
    _save_width_sweep_plot(output_path)
    _save_iv_plot(output_path)
    _save_resonant_plot(output_path)
    return output_path
