"""Potential profile helpers."""

from __future__ import annotations

import numpy as np


def rectangular_barrier(x_nm: np.ndarray, height_ev: float, width_nm: float) -> np.ndarray:
    potential = np.zeros_like(x_nm, dtype=float)
    half_width = width_nm / 2.0
    mask = np.abs(x_nm) <= half_width
    potential[mask] = height_ev
    return potential


def double_barrier(
    x_nm: np.ndarray,
    barrier_height_ev: float,
    barrier_width_nm: float,
    well_width_nm: float,
) -> np.ndarray:
    potential = np.zeros_like(x_nm, dtype=float)
    half_well = well_width_nm / 2.0
    left_start = -half_well - barrier_width_nm
    left_end = -half_well
    right_start = half_well
    right_end = half_well + barrier_width_nm
    left_mask = (x_nm >= left_start) & (x_nm <= left_end)
    right_mask = (x_nm >= right_start) & (x_nm <= right_end)
    potential[left_mask] = barrier_height_ev
    potential[right_mask] = barrier_height_ev
    return potential
