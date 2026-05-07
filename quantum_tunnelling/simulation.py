"""High-level simulations and plotting."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .constants import BOLTZMANN, ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, PLANCK
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


def tilted_double_barrier(
    x_nm: np.ndarray,
    voltage: float,
    barrier_height_ev: float,
    barrier_width_nm: float,
    well_width_nm: float,
) -> np.ndarray:
    potential = double_barrier(x_nm, barrier_height_ev, barrier_width_nm, well_width_nm)
    tilt = voltage * (x_nm - x_nm[0]) / (x_nm[-1] - x_nm[0])
    return potential - tilt


def tsu_esaki_supply(energies_ev: np.ndarray, voltage: float, fermi_ev: float, temperature_k: float) -> np.ndarray:
    kbt_ev = BOLTZMANN * temperature_k / ELEMENTARY_CHARGE
    left = np.logaddexp(0.0, (fermi_ev - energies_ev) / kbt_ev)
    right = np.logaddexp(0.0, (fermi_ev - energies_ev - voltage) / kbt_ev)
    return left - right


def free_electron_rtd_current_density(
    voltages: np.ndarray,
    temperature_k: float,
    energies_ev: np.ndarray,
    transmission_by_bias: np.ndarray,
    fermi_ev: float = 0.18,
) -> np.ndarray:
    """Free-electron-mass RTD baseline from biased transmission spectra."""
    prefactor = ELEMENTARY_CHARGE * ELECTRON_MASS * BOLTZMANN * temperature_k / (2.0 * np.pi**2 * HBAR**3)
    currents = []

    for voltage, transmission in zip(voltages, transmission_by_bias):
        supply = tsu_esaki_supply(energies_ev, voltage, fermi_ev, temperature_k)
        integral = np.trapezoid(transmission * supply, energies_ev * ELEMENTARY_CHARGE)
        currents.append(prefactor * integral)

    return np.array(currents, dtype=float)


def _gaussian_kernel(sigma_points: float, radius: int) -> np.ndarray:
    offsets = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (offsets / sigma_points) ** 2)
    return kernel / kernel.sum()


def broaden_transmission(transmission_by_bias: np.ndarray, sigma_ev: float, energy_step_ev: float) -> np.ndarray:
    sigma_points = max(sigma_ev / energy_step_ev, 1e-6)
    radius = max(2, int(np.ceil(3.0 * sigma_points)))
    kernel = _gaussian_kernel(sigma_points, radius)
    padded = np.pad(transmission_by_bias, ((0, 0), (radius, radius)), mode="edge")
    broadened = np.empty_like(transmission_by_bias)

    for row_idx in range(transmission_by_bias.shape[0]):
        broadened[row_idx] = np.convolve(padded[row_idx], kernel, mode="valid")

    return broadened


def smooth_curve(values: np.ndarray, passes: int = 2) -> np.ndarray:
    smoothed = np.asarray(values, dtype=float)
    kernel = np.array([1.0, 4.0, 6.0, 4.0, 1.0]) / 16.0
    radius = len(kernel) // 2

    for _ in range(passes):
        padded = np.pad(smoothed, radius, mode="edge")
        smoothed = np.convolve(padded, kernel, mode="valid")

    return smoothed


def biased_transmission_grid(
    voltages: np.ndarray,
    energies_ev: np.ndarray,
    x_nm: np.ndarray,
    barrier_height_ev: float = 0.35,
    barrier_width_nm: float = 0.7,
    well_width_nm: float = 1.2,
) -> np.ndarray:
    transmissions = []

    for voltage in voltages:
        potential_ev = tilted_double_barrier(
            x_nm,
            voltage,
            barrier_height_ev=barrier_height_ev,
            barrier_width_nm=barrier_width_nm,
            well_width_nm=well_width_nm,
        )
        transmissions.append(
            [
                solve_scattering(x_nm, potential_ev, energy_ev, mass=ELECTRON_MASS).transmission
                for energy_ev in energies_ev
            ]
        )

    return np.array(transmissions, dtype=float)


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
    x_nm = np.linspace(-4.0, 4.0, 420)
    energies_ev = np.linspace(0.005, 0.65, 180)
    voltages = np.linspace(0.0, 0.5, 111)
    transmission_by_bias = biased_transmission_grid(voltages, energies_ev, x_nm)
    transmission_by_bias = broaden_transmission(
        transmission_by_bias,
        sigma_ev=0.004,
        energy_step_ev=energies_ev[1] - energies_ev[0],
    )
    temperatures = [
        (77.0, "#0057b8", "T = 77 K"),
        (150.0, "#00965e", "T = 150 K"),
        (300.0, "#df2b2b", "T = 300 K"),
    ]

    fig, ax = plt.subplots(figsize=(10.6, 7.4))
    for temperature_k, color, label in temperatures:
        current_density = free_electron_rtd_current_density(
            voltages,
            temperature_k,
            energies_ev,
            transmission_by_bias,
        )
        current_density = smooth_curve(current_density, passes=2)
        ax.plot(
            voltages * 1000.0,
            current_density,
            color=color,
            linewidth=3.1,
            label=label,
        )

    ax.set_title("Free-Electron-Mass RTD Baseline I-V -- Temperature Dependence", fontsize=16, fontweight="bold")
    ax.set_xlabel("Applied Bias (mV)", fontsize=13)
    ax.set_ylabel(r"Current Density J (A/m$^2$)", fontsize=13)
    ax.set_xlim(0, 500)
    ax.set_xticks(np.arange(0, 501, 50))
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax.grid(alpha=0.38)
    ax.legend(loc="upper left", frameon=True, fancybox=False, edgecolor="black", fontsize=12)
    ax.text(
        0.54,
        0.76,
        "Baseline: m* = m0\nTilted potential under bias\nTsu-Esaki supply",
        transform=ax.transAxes,
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
