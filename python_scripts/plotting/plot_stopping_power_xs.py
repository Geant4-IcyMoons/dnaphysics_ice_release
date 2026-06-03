#!/usr/bin/env python3
"""
Compute stopping power and W-value directly from Geant4-DNA cross-section tables (no MC).

Water model modes:
- dna_opt4 (default): emulate Geant4-DNA `dna_opt4` water stack for electrons
  (ionisation + electronic excitation only in stopping-power total).
- custom: legacy table mix (adds vib + attachment to stopping-power total).

Ice curves keep the project custom ice model composition.

Notes:
- The Emfietzoglou excitation/ionisation tables stop at 10 keV.
- Born excitation/ionisation tables extend to 1 MeV and are used above 10 keV
  when available.
- Vib and attachment are limited to lower energies.
  Outside the model ranges the stopping power is set to 0.
- W-value uses:
  w(E) = [(dE/dx)_ion + (dE/dx)_exc] / [dN_ion/dx].
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

from constants import (
    CROSS_SECTIONS_DIR,
    CUSTOM_DATA_ROOT_GEANT4,
    CUSTOM_DATA_ROOT_PROJECT,
    EMFI_EXCITATION_EEV,
    EMFI_ION_BINDING_EEV,
    EMFIETZOGLOU_SCALE_1E16,
    FONT_COURIER,
    FONTSIZE_24,
    ICE_AMORPHOUS_DENSITY_G_CM3,
    ICE_HEXAGONAL_DENSITY_G_CM3,
    N,
    OUTPUT_DIR,
    PROJECT_ROOT,
    RC_BASE_ELASTIC,
    MICHAUD_SIGMA_SCALE_CM2,
    MELTON_SIGMA_SCALE_CM2,
    rcparams_with_fontsize,
)

# --- Plot style (match plot_stopping_power_root.py) ---
plt.rcParams["font.family"] = FONT_COURIER
plt.rcParams["mathtext.rm"] = FONT_COURIER
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams.update(rcparams_with_fontsize(RC_BASE_ELASTIC, FONTSIZE_24))

# --- Unit conversion ---
EV_NM_TO_MEV_CM = 10.0
EV_NM_TO_EV_ANG = 0.1

CM_PER_NM = 1.0e7

WATER_DENSITY_G_CM3 = 1.0
ICE_DENSITY_BY_TYPE = {
    "amorphous": ICE_AMORPHOUS_DENSITY_G_CM3,
    "hexagonal": ICE_HEXAGONAL_DENSITY_G_CM3,
}

N_CM3_WATER = N / 1.0e6

VIB_ELOSS_EEV = np.array(
    [0.024, 0.061, 0.092, 0.205, 0.417, 0.460, 0.510, 0.834], dtype=float
)

HIGH_ENERGY_SWITCH_EEV = 1.0e4
WVALUE_MAX_WATER_EEV = 1.0e6
WVALUE_MAX_ICE_EEV = 1.0e7


def _data_root() -> Path:
    if CUSTOM_DATA_ROOT_GEANT4.exists():
        return CUSTOM_DATA_ROOT_GEANT4
    return CUSTOM_DATA_ROOT_PROJECT


def _ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required data file: {path}")


def _interp_loglog(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_new = np.asarray(x_new, dtype=float)

    if x.size == 0:
        return np.zeros_like(x_new)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    positive = y > 0
    if np.count_nonzero(positive) < 2:
        return np.zeros_like(x_new)

    x_pos = x[positive]
    y_pos = y[positive]

    out = np.zeros_like(x_new)
    in_range = (x_new >= x_pos[0]) & (x_new <= x_pos[-1])
    if np.any(in_range):
        out[in_range] = 10 ** np.interp(
            np.log10(x_new[in_range]),
            np.log10(x_pos),
            np.log10(y_pos),
        )
    out[~np.isfinite(out)] = 0.0
    return out


def _units_and_label(units: str, density: float, dedx: np.ndarray) -> tuple[np.ndarray, str]:
    if units == "ev_ang":
        return dedx * EV_NM_TO_EV_ANG, r"S (eV/$\AA$)"
    if units == "mev_cm":
        return dedx * EV_NM_TO_MEV_CM, "S (MeV/cm)"
    if units == "mev_cm2_g":
        return (dedx * EV_NM_TO_MEV_CM) / density, "Mass Stopping power (MeV cm^2/g)"
    return dedx, "Stopping power (eV/nm)"


def _load_table(path: Path) -> np.ndarray:
    _ensure_exists(path)
    return np.loadtxt(path, dtype=float)


def _resolve_wvalue_path(path_like: str | Path) -> Path | None:
    raw = Path(path_like)
    candidates: list[Path] = [raw]
    if raw.suffix == "":
        candidates.append(raw.with_suffix(".txt"))
    alias_map = {
        "wvalue_water_10MeV": "wvalue_water.txt",
        "wvalue_hex_10MeV": "wvalue_ice_hex.txt",
        "wvalue_am_10MeV": "wvalue_ice_am.txt",
    }
    raw_stem = raw.stem
    if raw_stem in alias_map:
        alias = Path(alias_map[raw_stem])
        candidates.append(alias)
    if not raw.is_absolute():
        prefixed = []
        for cand in candidates:
            prefixed.append(Path("build") / cand)
            prefixed.append(PROJECT_ROOT / cand)
            prefixed.append(PROJECT_ROOT / "build" / cand)
        candidates.extend(prefixed)

    seen = set()
    deduped: list[Path] = []
    for cand in candidates:
        key = str(cand)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(cand)

    for cand in deduped:
        if cand.exists():
            return cand
    return None


def _load_wvalue_curve(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    data = np.loadtxt(path, dtype=float)
    data = np.atleast_2d(data)
    if data.shape[1] < 4:
        raise ValueError(
            f"W-value file {path} must have at least 4 columns: E, Nion, rms, W."
        )

    # Keep the latest entry for each energy (files are append-only across runs).
    latest: dict[float, tuple[float, float | None]] = {}
    for row in data:
        e = float(row[0])
        w = float(row[3])
        w_err = float(row[-1]) if data.shape[1] >= 5 else None
        latest[e] = (w, w_err)

    energies = np.asarray(sorted(latest.keys()), dtype=float)
    wvals = np.asarray([latest[e][0] for e in energies], dtype=float)
    if data.shape[1] >= 5:
        werrs = np.asarray(
            [
                np.nan if latest[e][1] is None else float(latest[e][1])
                for e in energies
            ],
            dtype=float,
        )
    else:
        werrs = None

    valid = np.isfinite(energies) & np.isfinite(wvals) & (energies > 0) & (wvals > 0)
    if werrs is not None:
        valid &= np.isfinite(werrs) & (werrs >= 0)
        return energies[valid], wvals[valid], werrs[valid]
    return energies[valid], wvals[valid], None


def _max_table_energy(path: Path) -> float:
    _ensure_exists(path)
    data = np.loadtxt(path, dtype=float, usecols=[0])
    if data.size == 0:
        raise ValueError(f"No energy data in {path}")
    return float(np.max(data))


def _vib_energy_loss_xs(energy_grid: np.ndarray, path: Path) -> np.ndarray:
    data = _load_table(path)
    energies = data[:, 0]
    sigma_levels = data[:, 1:]
    if sigma_levels.shape[1] != VIB_ELOSS_EEV.size:
        raise ValueError(f"Unexpected vib table shape in {path}")

    sigma_cm2 = sigma_levels * MICHAUD_SIGMA_SCALE_CM2
    eloss_xs = np.zeros_like(energy_grid)
    for i, omega in enumerate(VIB_ELOSS_EEV):
        sigma_i = _interp_loglog(energies, sigma_cm2[:, i], energy_grid)
        eloss_xs += sigma_i * omega
    return eloss_xs


def _attachment_energy_loss_xs(
    energy_grid: np.ndarray, path: Path, sigma_scale_cm2: float
) -> np.ndarray:
    data = _load_table(path)
    energies = data[:, 0]
    sigma = data[:, 1] * sigma_scale_cm2
    sigma_interp = _interp_loglog(energies, sigma, energy_grid)
    return sigma_interp * energy_grid


def _vib_sigma_xs(energy_grid: np.ndarray, path: Path) -> np.ndarray:
    data = _load_table(path)
    energies = data[:, 0]
    sigma_levels = data[:, 1:]
    if sigma_levels.shape[1] != VIB_ELOSS_EEV.size:
        raise ValueError(f"Unexpected vib table shape in {path}")
    sigma_cm2 = sigma_levels * MICHAUD_SIGMA_SCALE_CM2
    sigma_total = np.sum(sigma_cm2, axis=1)
    return _interp_loglog(energies, sigma_total, energy_grid)


def _attachment_sigma_xs(
    energy_grid: np.ndarray, path: Path, sigma_scale_cm2: float
) -> np.ndarray:
    data = _load_table(path)
    energies = data[:, 0]
    sigma = data[:, 1] * sigma_scale_cm2
    return _interp_loglog(energies, sigma, energy_grid)


def _excitation_energy_loss_xs(energy_grid: np.ndarray, path: Path) -> np.ndarray:
    data = _load_table(path)
    energies = data[:, 0]
    sigma_levels = data[:, 1:1 + EMFI_EXCITATION_EEV.size]
    if sigma_levels.shape[1] != EMFI_EXCITATION_EEV.size:
        raise ValueError(f"Unexpected excitation table shape in {path}")

    scale_cm2 = EMFIETZOGLOU_SCALE_1E16 * 1.0e-16
    sigma_cm2 = sigma_levels * scale_cm2

    eloss_xs = np.zeros_like(energy_grid)
    for i, e_loss in enumerate(EMFI_EXCITATION_EEV):
        sigma_i = _interp_loglog(energies, sigma_cm2[:, i], energy_grid)
        eloss_xs += sigma_i * e_loss
    return eloss_xs


def _excitation_sigma_xs(energy_grid: np.ndarray, path: Path) -> np.ndarray:
    data = _load_table(path)
    energies = data[:, 0]
    sigma_levels = data[:, 1:1 + EMFI_EXCITATION_EEV.size]
    if sigma_levels.shape[1] != EMFI_EXCITATION_EEV.size:
        raise ValueError(f"Unexpected excitation table shape in {path}")

    scale_cm2 = EMFIETZOGLOU_SCALE_1E16 * 1.0e-16
    sigma_cm2 = sigma_levels * scale_cm2
    sigma_total = np.sum(sigma_cm2, axis=1)
    return _interp_loglog(energies, sigma_total, energy_grid)


def _trapz_nonuniform(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1])))


def _simpson_nonuniform(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 2:
        return 0.0

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if np.any(np.diff(x) == 0.0):
        unique_x, inv = np.unique(x, return_inverse=True)
        y_acc = np.zeros_like(unique_x)
        counts = np.zeros_like(unique_x)
        for i, idx in enumerate(inv):
            y_acc[idx] += y[i]
            counts[idx] += 1.0
        y = y_acc / np.where(counts == 0.0, 1.0, counts)
        x = unique_x

    n = x.size
    if n == 2:
        return _trapz_nonuniform(x, y)

    def _simpson_segment(x0, x1, x2, y0, y1, y2):
        h0 = x1 - x0
        h1 = x2 - x1
        if h0 <= 0.0 or h1 <= 0.0:
            return 0.0
        denom = h0 * h1
        if denom == 0.0:
            return 0.0
        w0 = 2.0 - (h1 / h0)
        w1 = ((h0 + h1) ** 2) / denom
        w2 = 2.0 - (h0 / h1)
        return (h0 + h1) / 6.0 * (w0 * y0 + w1 * y1 + w2 * y2)

    end = n if n % 2 == 1 else n - 1
    total = 0.0
    for i in range(0, end - 2, 2):
        total += _simpson_segment(
            x[i], x[i + 1], x[i + 2], y[i], y[i + 1], y[i + 2]
        )
    if n % 2 == 0:
        total += _trapz_nonuniform(x[-2:], y[-2:])
    return total


def _integrate_dcs_eloss(e_vals: list[float], sigma_vals: list[list[float]]) -> float:
    if len(e_vals) < 2:
        return 0.0
    e = np.asarray(e_vals, dtype=float)
    sigma = np.asarray(sigma_vals, dtype=float)
    scale_cm2 = EMFIETZOGLOU_SCALE_1E16 * 1.0e-16
    sigma *= scale_cm2

    eloss = 0.0
    for j in range(sigma.shape[0]):
        eloss += _simpson_nonuniform(e, e * sigma[j])
    return eloss


def _integrate_dcs_sigma(e_vals: list[float], sigma_vals: list[list[float]]) -> float:
    if len(e_vals) < 2:
        return 0.0
    e = np.asarray(e_vals, dtype=float)
    sigma = np.asarray(sigma_vals, dtype=float)
    scale_cm2 = EMFIETZOGLOU_SCALE_1E16 * 1.0e-16
    sigma *= scale_cm2

    sigma_total = 0.0
    for j in range(sigma.shape[0]):
        sigma_total += _simpson_nonuniform(e, sigma[j])
    return sigma_total


def _dcs_energy_loss_xs_table(path: Path) -> tuple[np.ndarray, np.ndarray]:
    _ensure_exists(path)
    energies = []
    eloss_xs = []

    current_e = None
    e_vals: list[float] = []
    sigma_vals: list[list[float]] = []
    n_channels = None

    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            e = float(parts[0])
            transfer = float(parts[1])
            if n_channels is None:
                n_channels = len(parts) - 2
                sigma_vals = [[] for _ in range(n_channels)]
            if len(parts) < 2 + n_channels:
                continue
            sigs = [float(val) for val in parts[2 : 2 + n_channels]]

            if current_e is None:
                current_e = e
            if e != current_e:
                energies.append(current_e)
                eloss_xs.append(_integrate_dcs_eloss(e_vals, sigma_vals))
                current_e = e
                e_vals = []
                sigma_vals = [[] for _ in range(n_channels)]

            e_vals.append(transfer)
            for j in range(n_channels):
                sigma_vals[j].append(sigs[j])

    if current_e is not None:
        energies.append(current_e)
        eloss_xs.append(_integrate_dcs_eloss(e_vals, sigma_vals))

    return np.asarray(energies, dtype=float), np.asarray(eloss_xs, dtype=float)


def _dcs_sigma_xs_table(path: Path) -> tuple[np.ndarray, np.ndarray]:
    _ensure_exists(path)
    energies = []
    sigma_xs = []

    current_e = None
    e_vals: list[float] = []
    sigma_vals: list[list[float]] = []
    n_channels = None

    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            e = float(parts[0])
            transfer = float(parts[1])
            if n_channels is None:
                n_channels = len(parts) - 2
                sigma_vals = [[] for _ in range(n_channels)]
            if len(parts) < 2 + n_channels:
                continue
            sigs = [float(val) for val in parts[2 : 2 + n_channels]]

            if current_e is None:
                current_e = e
            if e != current_e:
                energies.append(current_e)
                sigma_xs.append(_integrate_dcs_sigma(e_vals, sigma_vals))
                current_e = e
                e_vals = []
                sigma_vals = [[] for _ in range(n_channels)]

            e_vals.append(transfer)
            for j in range(n_channels):
                sigma_vals[j].append(sigs[j])

    if current_e is not None:
        energies.append(current_e)
        sigma_xs.append(_integrate_dcs_sigma(e_vals, sigma_vals))

    return np.asarray(energies, dtype=float), np.asarray(sigma_xs, dtype=float)


def _integrate_ion_dcs_moments(
    incident_energy_eV: float,
    transfer_vals: list[float],
    sigma_vals: list[list[float]],
    bindings_eV: np.ndarray,
) -> tuple[float, float]:
    if len(transfer_vals) < 2:
        return 0.0, 0.0

    omega = np.asarray(transfer_vals, dtype=float)
    sigma = np.asarray(sigma_vals, dtype=float)
    n_channels = sigma.shape[0]

    b = np.asarray(bindings_eV, dtype=float)
    if b.size < n_channels:
        b = np.pad(b, (0, n_channels - b.size), mode="edge")
    elif b.size > n_channels:
        b = b[:n_channels]

    scale_cm2 = EMFIETZOGLOU_SCALE_1E16 * 1.0e-16
    sigma *= scale_cm2

    sigma_ion = 0.0
    ion_eloss = 0.0
    for j in range(n_channels):
        s = sigma[j]
        # Electron-impact convention:
        # secondary kinetic W = omega - B_j with W in [0, (T - B_j)/2],
        # so omega in [B_j, (T + B_j)/2].
        omega_min = b[j]
        omega_max = min(incident_energy_eV, 0.5 * (incident_energy_eV + b[j]))
        if omega_max <= omega_min:
            continue

        mask = (omega >= omega_min) & (omega <= omega_max)
        if np.count_nonzero(mask) < 2:
            continue
        omega_eff = omega[mask]
        s_eff = s[mask]

        sigma_ion += _simpson_nonuniform(omega_eff, s_eff)
        ion_eloss += _simpson_nonuniform(omega_eff, omega_eff * s_eff)

    return sigma_ion, ion_eloss


def _ion_dcs_moments_table(
    path: Path, bindings_eV: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    _ensure_exists(path)
    energies = []
    sigma_ion_xs = []
    ion_eloss_xs = []

    current_e = None
    transfer_vals: list[float] = []
    sigma_vals: list[list[float]] = []
    n_channels = None

    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            e = float(parts[0])
            transfer = float(parts[1])
            if n_channels is None:
                n_channels = len(parts) - 2
                sigma_vals = [[] for _ in range(n_channels)]
            if len(parts) < 2 + n_channels:
                continue
            sigs = [float(val) for val in parts[2 : 2 + n_channels]]

            if current_e is None:
                current_e = e
            if e != current_e:
                sigma_i, eloss_i = _integrate_ion_dcs_moments(
                    current_e, transfer_vals, sigma_vals, bindings_eV
                )
                energies.append(current_e)
                sigma_ion_xs.append(sigma_i)
                ion_eloss_xs.append(eloss_i)
                current_e = e
                transfer_vals = []
                sigma_vals = [[] for _ in range(n_channels)]

            transfer_vals.append(transfer)
            for j in range(n_channels):
                sigma_vals[j].append(sigs[j])

    if current_e is not None:
        sigma_i, eloss_i = _integrate_ion_dcs_moments(
            current_e, transfer_vals, sigma_vals, bindings_eV
        )
        energies.append(current_e)
        sigma_ion_xs.append(sigma_i)
        ion_eloss_xs.append(eloss_i)

    return (
        np.asarray(energies, dtype=float),
        np.asarray(sigma_ion_xs, dtype=float),
        np.asarray(ion_eloss_xs, dtype=float),
    )


def _to_dedx_ev_nm(eloss_xs: np.ndarray, n_cm3: float) -> np.ndarray:
    return (eloss_xs * n_cm3) / CM_PER_NM


def _merge_low_high_grid(
    energy_grid: np.ndarray,
    low_vals: np.ndarray,
    high_vals: np.ndarray,
    switch_e: float,
) -> np.ndarray:
    out = np.asarray(high_vals, dtype=float).copy()
    low_vals = np.asarray(low_vals, dtype=float)
    mask = energy_grid <= switch_e
    out[mask] = low_vals[mask]
    return out


def _merge_low_high_tables(
    low_energy: np.ndarray,
    low_vals: np.ndarray,
    high_energy: np.ndarray,
    high_vals: np.ndarray,
    switch_e: float,
) -> tuple[np.ndarray, np.ndarray]:
    low_energy = np.asarray(low_energy, dtype=float)
    low_vals = np.asarray(low_vals, dtype=float)
    high_energy = np.asarray(high_energy, dtype=float)
    high_vals = np.asarray(high_vals, dtype=float)

    if low_energy.size == 0:
        return high_energy, high_vals
    if high_energy.size == 0:
        return low_energy, low_vals

    low_mask = low_energy <= switch_e
    high_mask = high_energy > switch_e
    merged_energy = np.concatenate([low_energy[low_mask], high_energy[high_mask]])
    merged_vals = np.concatenate([low_vals[low_mask], high_vals[high_mask]])

    order = np.argsort(merged_energy)
    return merged_energy[order], merged_vals[order]


def _compute_stopping_components(
    energy_grid: np.ndarray,
    n_cm3: float,
    vib_path: Path | None,
    attach_path: Path | None,
    attach_sigma_scale_cm2: float,
    exc_eloss_xs: np.ndarray,
    ion_energy: np.ndarray,
    ion_eloss_xs: np.ndarray,
) -> dict[str, np.ndarray]:
    if vib_path is not None:
        vib_eloss = _vib_energy_loss_xs(energy_grid, vib_path)
    else:
        vib_eloss = np.zeros_like(energy_grid)
    if attach_path is not None:
        attach_eloss = _attachment_energy_loss_xs(
            energy_grid, attach_path, attach_sigma_scale_cm2
        )
    else:
        attach_eloss = np.zeros_like(energy_grid)
    ion_eloss = _interp_loglog(ion_energy, ion_eloss_xs, energy_grid)

    vib_dedx = _to_dedx_ev_nm(vib_eloss, n_cm3)
    exc_dedx = _to_dedx_ev_nm(exc_eloss_xs, n_cm3)
    ion_dedx = _to_dedx_ev_nm(ion_eloss, n_cm3)
    attach_dedx = _to_dedx_ev_nm(attach_eloss, n_cm3)

    return {
        "vib": vib_dedx,
        "exc": exc_dedx,
        "ion": ion_dedx,
        "attach": attach_dedx,
        "total": vib_dedx + exc_dedx + ion_dedx + attach_dedx,
    }


def _compute_w_value(
    energy_grid: np.ndarray,
    n_cm3: float,
    ion_energy: np.ndarray,
    ion_sigma_xs: np.ndarray,
    ion_dedx: np.ndarray,
    exc_dedx: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    sigma_ion_interp = _interp_loglog(ion_energy, ion_sigma_xs, energy_grid)
    dNion_dx = (sigma_ion_interp * n_cm3) / CM_PER_NM

    w_value = np.full_like(energy_grid, np.nan, dtype=float)
    # Only compute W where ionization is significant (avoid low-energy artifacts)
    # Use threshold: energy > 12 eV (ionization threshold) AND dNion_dx meaningful
    min_energy_eV = 12.0  # First ionization threshold for water
    min_rate = 1.0e-8  # Minimum ionization rate (1/nm)
    valid = (energy_grid >= min_energy_eV) & (dNion_dx > min_rate)
    # W-value: (ionization + excitation) energy per ionization event
    # This accounts for energy spent on both ionization and excitation
    w_value[valid] = (ion_dedx[valid] + exc_dedx[valid]) / dNion_dx[valid]
    return w_value, dNion_dx


def _print_w_assumptions() -> None:
    print("\nW-value model assumptions:")
    print("  1) One ion pair per ionisation event (no explicit multiple-ionisation multiplicity).")
    print("  2) Ionisation DCS interpreted as dσ_j(T,ω)/dω; electron-impact limit ω in [B_j, (T+B_j)/2].")
    print("  3) W = [(dE/dx)_ion + (dE/dx)_exc] / (dN_ion/dx): includes excitation energy.")
    print("  4) Vibrational, elastic, and attachment channels excluded from W numerator.")
    print("  5) Completeness is limited by available tabulated cross sections and energy coverage.")


def _print_w_diagnostics(
    label: str,
    energy_grid: np.ndarray,
    ion_energy: np.ndarray,
    ion_sigma_xs: np.ndarray,
    dNion_dx: np.ndarray,
    ion_dedx: np.ndarray,
    exc_dedx: np.ndarray,
    w_value: np.ndarray,
    emax: float | None = None,
) -> None:
    sigma_interp = _interp_loglog(ion_energy, ion_sigma_xs, energy_grid)
    targets = np.array([20.0, 100.0, 1.0e3, 1.0e4], dtype=float)
    print(f"\n{label} diagnostics:")
    print("  E (eV) | sigma_ion (cm^2) | dNion/dx (1/nm) | dE/dx_ion (eV/nm) | dE/dx_exc (eV/nm) | W (eV/ion)")
    for target in targets:
        if target < energy_grid[0] or target > energy_grid[-1]:
            continue
        if emax is not None and target > emax:
            continue
        idx = int(np.argmin(np.abs(np.log10(energy_grid) - np.log10(target))))
        if not np.isfinite(w_value[idx]):
            continue
        print(
            f"  {energy_grid[idx]:8.2f} | "
            f"{sigma_interp[idx]:.3e} | "
            f"{dNion_dx[idx]:.3e} | "
            f"{ion_dedx[idx]:.3e} | "
            f"{exc_dedx[idx]:.3e} | "
            f"{w_value[idx]:.3f}"
        )


def _plot(
    energy_grid: np.ndarray,
    water_dedx: np.ndarray,
    water_w: np.ndarray,
    water_mfp_nm: np.ndarray,
    water_color,
    water_dom_sp: np.ndarray | None,
    water_dom_w: np.ndarray | None,
    water_dom_mfp: np.ndarray | None,
    am_dom_sp: np.ndarray | None,
    am_dom_w: np.ndarray | None,
    am_dom_mfp: np.ndarray | None,
    hex_dom_sp: np.ndarray | None,
    hex_dom_w: np.ndarray | None,
    hex_dom_mfp: np.ndarray | None,
    ice_series: list[dict],
    units: str,
    out_path: Path,
    rho_water: float,
    water_emax: float | None = None,
    water_w_emax: float | None = None,
    ice_w_emax_by_type: dict[str, float] | None = None,
    w_text_curves: dict[str, dict[str, np.ndarray]] | None = None,
) -> None:
    def _add_dominance_spans(ax, dom, ymin, ymax, colors, alpha=0.5):
        if dom is None:
            return
        dom = np.asarray(dom)
        valid = dom >= 0
        if not np.any(valid):
            return
        start = None
        current = None
        for i in range(dom.size):
            if not valid[i]:
                if start is not None:
                    end = i - 1
                    x0 = energy_grid[start]
                    x1 = energy_grid[end + 1] if end + 1 < energy_grid.size else energy_grid[end]
                    ax.axvspan(x0, x1, ymin=ymin, ymax=ymax, color=colors[int(current)], alpha=alpha, zorder=0)
                    start = None
                    current = None
                continue
            if current is None:
                start = i
                current = dom[i]
                continue
            if dom[i] != current:
                end = i - 1
                x0 = energy_grid[start]
                x1 = energy_grid[end + 1] if end + 1 < energy_grid.size else energy_grid[end]
                ax.axvspan(x0, x1, ymin=ymin, ymax=ymax, color=colors[int(current)], alpha=alpha, zorder=0)
                start = i
                current = dom[i]
        if start is not None:
            end = dom.size - 1
            x0 = energy_grid[start]
            x1 = energy_grid[end]
            ax.axvspan(x0, x1, ymin=ymin, ymax=ymax, color=colors[int(current)], alpha=alpha, zorder=0)

    water_plot, ylabel = _units_and_label(units, rho_water, water_dedx)
    fig, (ax_sp, ax_w, ax_mfp) = plt.subplots(
        3,
        1,
        figsize=(10, 13),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.0]},
    )
    # Grayscale dominance bands: lightest=vib/att, then excitation, then ionization.
    dom_colors = ["#f0f0f0", "#cfcfcf", "#a9a9a9"]
    # Top third: water, middle: amorphous, bottom: hexagonal
    _add_dominance_spans(ax_sp, water_dom_sp, 2.0 / 3.0, 1.0, dom_colors)
    _add_dominance_spans(ax_sp, am_dom_sp, 1.0 / 3.0, 2.0 / 3.0, dom_colors)
    _add_dominance_spans(ax_sp, hex_dom_sp, 0.0, 1.0 / 3.0, dom_colors)
    _add_dominance_spans(ax_w, water_dom_w, 2.0 / 3.0, 1.0, dom_colors)
    _add_dominance_spans(ax_w, am_dom_w, 1.0 / 3.0, 2.0 / 3.0, dom_colors)
    _add_dominance_spans(ax_w, hex_dom_w, 0.0, 1.0 / 3.0, dom_colors)
    _add_dominance_spans(ax_mfp, water_dom_mfp, 2.0 / 3.0, 1.0, dom_colors)
    _add_dominance_spans(ax_mfp, am_dom_mfp, 1.0 / 3.0, 2.0 / 3.0, dom_colors)
    _add_dominance_spans(ax_mfp, hex_dom_mfp, 0.0, 1.0 / 3.0, dom_colors)

    valid_w = np.isfinite(water_plot) & (water_plot > 0) & np.isfinite(energy_grid)
    if water_emax is not None:
        valid_w &= energy_grid <= water_emax
    ax_sp.loglog(
        energy_grid[valid_w],
        water_plot[valid_w],
        color=water_color,
        linewidth=3,
        label="Water",
        zorder=1,
    )

    for idx, series in enumerate(ice_series):
        ice_density = float(series.get("density_g_cm3", 1.0))
        ice_plot, _ = _units_and_label(units, ice_density, series["dedx"])
        valid_i = np.isfinite(ice_plot) & (ice_plot > 0) & np.isfinite(energy_grid)
        ice_emax = series.get("emax")
        if ice_emax is not None:
            valid_i &= energy_grid <= ice_emax
        ax_sp.loglog(
            energy_grid[valid_i],
            ice_plot[valid_i],
            color=series.get("color", f"0.{35 + idx * 20:02d}"),
            linewidth=3,
            ls=series.get("ls", "-"),
            label=series.get("label", "Ice"),
            zorder=3,
        )

    ax_sp.set_ylabel(ylabel)
    # Legend moved to figure-level below panels.

    water_w_energy = energy_grid
    water_w_vals = water_w
    water_w_err = None
    if w_text_curves and "water" in w_text_curves:
        water_w_energy = w_text_curves["water"]["energy"]
        water_w_vals = w_text_curves["water"]["w"]
        water_w_err = w_text_curves["water"].get("w_err")
    valid_wv = (
        np.isfinite(water_w_vals)
        & (water_w_vals > 0)
        & np.isfinite(water_w_energy)
        & (water_w_energy > 0)
    )
    if water_w_emax is not None:
        valid_wv &= water_w_energy <= water_w_emax
    if water_w_err is not None:
        valid_werr = (
            valid_wv
            & np.isfinite(water_w_err)
            & (water_w_vals - water_w_err > 0)
        )
        if np.any(valid_werr):
            ax_w.fill_between(
                water_w_energy[valid_werr],
                water_w_vals[valid_werr] - water_w_err[valid_werr],
                water_w_vals[valid_werr] + water_w_err[valid_werr],
                color=water_color,
                alpha=0.18,
                linewidth=0,
                zorder=0.8,
            )
    ax_w.loglog(
        water_w_energy[valid_wv],
        water_w_vals[valid_wv],
        color=water_color,
        linewidth=3,
        label="Water",
        zorder=1,
    )
    for idx, series in enumerate(ice_series):
        ice_w_energy = energy_grid
        ice_w_vals = series["w_value"]
        ice_w_err = None
        ice_type = series.get("ice_type")
        if w_text_curves and isinstance(ice_type, str) and ice_type in w_text_curves:
            ice_w_energy = w_text_curves[ice_type]["energy"]
            ice_w_vals = w_text_curves[ice_type]["w"]
            ice_w_err = w_text_curves[ice_type].get("w_err")
        valid_iw = (
            np.isfinite(ice_w_vals)
            & (ice_w_vals > 0)
            & np.isfinite(ice_w_energy)
            & (ice_w_energy > 0)
        )
        ice_emax = series.get("emax")
        ice_w_emax = None
        if isinstance(ice_type, str) and ice_w_emax_by_type is not None:
            ice_w_emax = ice_w_emax_by_type.get(ice_type)
        max_w_e = None
        if ice_emax is not None and ice_w_emax is not None:
            max_w_e = min(float(ice_emax), float(ice_w_emax))
        elif ice_emax is not None:
            max_w_e = float(ice_emax)
        elif ice_w_emax is not None:
            max_w_e = float(ice_w_emax)
        if max_w_e is not None:
            valid_iw &= ice_w_energy <= max_w_e
        if ice_w_err is not None:
            valid_iw_err = (
                valid_iw
                & np.isfinite(ice_w_err)
                & (ice_w_vals - ice_w_err > 0)
            )
            if np.any(valid_iw_err):
                ax_w.fill_between(
                    ice_w_energy[valid_iw_err],
                    ice_w_vals[valid_iw_err] - ice_w_err[valid_iw_err],
                    ice_w_vals[valid_iw_err] + ice_w_err[valid_iw_err],
                    color=series.get("color", f"0.{35 + idx * 20:02d}"),
                    alpha=0.18,
                    linewidth=0,
                    zorder=2.8,
                )
        ax_w.loglog(
            ice_w_energy[valid_iw],
            ice_w_vals[valid_iw],
            color=series.get("color", f"0.{35 + idx * 20:02d}"),
            linewidth=3,
            label=series.get("label", "Ice"),
            ls=series.get("ls", "-"),
            zorder=3,
        )
    ax_w.set_xlabel("")
    ax_w.set_ylabel("W (eV)")
    ax_w.set_xlim(energy_grid.min(), energy_grid.max())
    w_ticks = [25., 50., 100., 200.]
    ax_w.set_ylim(15.0, 250.0)
    ax_w.yaxis.set_major_locator(FixedLocator(w_ticks))
    ax_w.yaxis.set_major_formatter(
        FuncFormatter(lambda y, _: f"{int(round(y))}" if y >= 1 else f"{y:g}")
    )
    ax_w.yaxis.set_minor_locator(NullLocator())
    ax_sp.text(
        0.0,
        1.04,
        "(a)",
        transform=ax_sp.transAxes,
        ha="left",
        va="bottom",
        fontsize=float(FONTSIZE_24),
        clip_on=False,
    )
    ax_w.text(
        0.0,
        1.04,
        "(b)",
        transform=ax_w.transAxes,
        ha="left",
        va="bottom",
        fontsize=float(FONTSIZE_24),
        clip_on=False,
    )

    valid_mfp = np.isfinite(water_mfp_nm) & (water_mfp_nm > 0) & np.isfinite(energy_grid)
    if water_emax is not None:
        valid_mfp &= energy_grid <= water_emax
    ax_mfp.loglog(
        energy_grid[valid_mfp],
        water_mfp_nm[valid_mfp],
        color=water_color,
        linewidth=3,
        label="Water",
        zorder=1,
    )
    for idx, series in enumerate(ice_series):
        ice_mfp = series.get("inelastic_mfp_nm")
        if ice_mfp is None:
            continue
        valid_i = np.isfinite(ice_mfp) & (ice_mfp > 0) & np.isfinite(energy_grid)
        ice_emax = series.get("emax")
        if ice_emax is not None:
            valid_i &= energy_grid <= ice_emax
        ax_mfp.loglog(
            energy_grid[valid_i],
            ice_mfp[valid_i],
            color=series.get("color", f"0.{35 + idx * 20:02d}"),
            linewidth=3,
            label=series.get("label", "Ice"),
            ls=series.get("ls", "-"),
            zorder=3,
        )
    ax_mfp.set_xlabel("Electron energy ($T$; eV)")
    ax_mfp.set_ylabel("IMFP (nm)")
    ax_mfp.text(
        0.0,
        1.04,
        "(c)",
        transform=ax_mfp.transAxes,
        ha="left",
        va="bottom",
        fontsize=float(FONTSIZE_24),
        clip_on=False,
    )

    handles, labels = ax_sp.get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.02),
            bbox_transform=fig.transFigure,
            ncol=max(1, len(labels)),
            frameon=False,
            columnspacing=1.2,
            handlelength=2.4,
            handletextpad=0.6,
        )

    plt.tight_layout(rect=(0.0, 0.08, 1.0, 1.0))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, bbox_inches="tight")
    print(f"\nPlot saved to: {out_path}")
    plt.show()


def _write_table(
    out_path: Path,
    energy_grid: np.ndarray,
    water_dedx: np.ndarray,
    water_w: np.ndarray,
    water_mfp_nm: np.ndarray,
    ice_series: list[dict],
) -> None:
    lines = []
    header = ["E_eV", "water_dedx_eV_per_nm", "water_W_eV", "water_IMFP_nm"]
    for series in ice_series:
        label = series.get("ice_type", series.get("label", "ice")).replace(" ", "_")
        header.extend(
            [
                f"{label}_dedx_eV_per_nm",
                f"{label}_W_eV",
                f"{label}_IMFP_nm",
            ]
        )
    lines.append("# " + " ".join(header))

    for i, e in enumerate(energy_grid):
        row = [
            f"{e:.6e}",
            f"{water_dedx[i]:.6e}",
            f"{water_w[i]:.6e}",
            f"{water_mfp_nm[i]:.6e}",
        ]
        for series in ice_series:
            row.append(f"{series['dedx'][i]:.6e}")
            row.append(f"{series['w_value'][i]:.6e}")
            row.append(f"{series['inelastic_mfp_nm'][i]:.6e}")
        lines.append(" ".join(row))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nTable saved to: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute stopping power and W-value from cross sections for water vs ice."
    )
    parser.add_argument("--emin", type=float, default=1.0, help="Minimum energy (eV)")
    parser.add_argument("--emax", type=float, default=1.0e7, help="Maximum energy (eV)")
    parser.add_argument("--nbins", type=int, default=500, help="Number of log bins")
    parser.add_argument(
        "--units",
        type=str,
        default="ev_ang",
        choices=("ev_ang", "ev_nm", "mev_cm", "mev_cm2_g"),
        help="Output units for stopping power.",
    )
    parser.add_argument(
        "--rho-water",
        type=float,
        default=WATER_DENSITY_G_CM3,
        help="Water density (g/cm^3) for mass stopping power.",
    )
    parser.add_argument(
        "--rho-ice",
        type=float,
        default=None,
        help=(
            "Override density (g/cm^3) for all ice phases in computed curves. "
            "If omitted, uses 0.94 for amorphous and 0.917 for hexagonal."
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=OUTPUT_DIR / "stopping_power_xs_water_vs_ice.png",
        help="Output plot path.",
    )
    parser.add_argument(
        "--out-txt",
        type=Path,
        default=OUTPUT_DIR / "stopping_power_w_imfp_water_vs_ice.txt",
        help="Output text table path.",
    )
    parser.add_argument(
        "--ice-types",
        type=str,
        default="amorphous,hexagonal",
        help="Comma-separated ice types to plot (amorphous,hexagonal).",
    )
    parser.add_argument(
        "--wvalue-water",
        type=str,
        default="wvalue_water.txt",
        help="Water W-value text file (name or path).",
    )
    parser.add_argument(
        "--wvalue-hex",
        type=str,
        default="wvalue_ice_hex.txt",
        help="Hexagonal ice W-value text file (name or path).",
    )
    parser.add_argument(
        "--wvalue-am",
        type=str,
        default="wvalue_ice_am.txt",
        help="Amorphous ice W-value text file (name or path).",
    )
    parser.add_argument(
        "--water-model",
        type=str,
        default="dna_opt4",
        choices=("dna_opt4", "custom"),
        help=(
            "Water stopping-power model: "
            "dna_opt4 (ion+electronic excitation only) or custom (includes vib+attachment)."
        ),
    )
    args = parser.parse_args()

    if args.emin <= 0 or args.emax <= 0 or args.emin >= args.emax:
        raise ValueError("Invalid energy range.")

    data_root = _data_root()
    dna_dir = data_root / "G4EMLOW8.6.1" / "dna"

    water_vib = dna_dir / "sigma_excitationvib_e_michaud.dat"
    ice_vib = CROSS_SECTIONS_DIR / "sigma_excitationvib_e_michaud.dat"
    water_attach = dna_dir / "sigma_attachment_e_melton.dat"
    ice_attach = CROSS_SECTIONS_DIR / "sigma_attachment_e_michaud.dat"
    exc_path = dna_dir / "sigma_excitation_e_emfietzoglou.dat"
    ion_dcs_path = dna_dir / "sigmadiff_ionisation_e_emfietzoglou.dat"
    exc_born_path = dna_dir / "sigma_excitation_e_born.dat"
    ion_born_dcs_path = dna_dir / "sigmadiff_ionisation_e_born.dat"
    ice_types = [t.strip().lower() for t in args.ice_types.split(",") if t.strip()]
    ice_types = [t for t in ice_types if t in ("amorphous", "hexagonal")]
    if not ice_types:
        raise ValueError("No valid ice types specified (use amorphous and/or hexagonal).")

    max_low_water = min(
        _max_table_energy(exc_path),
        _max_table_energy(ion_dcs_path),
    )
    ice_info = []
    for ice_type in ice_types:
        ice_label = f"{ice_type}_ice"
        ice_exc_dcs_path = CROSS_SECTIONS_DIR / f"sigmadiff_excitation_e_{ice_label}_emfietzoglou_kyriakou.dat"
        ice_ion_dcs_path = CROSS_SECTIONS_DIR / f"sigmadiff_ionisation_e_{ice_label}_emfietzoglou_kyriakou.dat"
        if not ice_exc_dcs_path.exists():
            fallback = CROSS_SECTIONS_DIR / "sigmadiff_excitation_e_ice_emfietzoglou_kyriakou.dat"
            if fallback.exists():
                print(
                    f"Missing {ice_exc_dcs_path.name}; falling back to {fallback.name}."
                )
                ice_exc_dcs_path = fallback
        if not ice_ion_dcs_path.exists():
            fallback = CROSS_SECTIONS_DIR / "sigmadiff_ionisation_e_ice_emfietzoglou_kyriakou.dat"
            if fallback.exists():
                print(
                    f"Missing {ice_ion_dcs_path.name}; falling back to {fallback.name}."
                )
                ice_ion_dcs_path = fallback
        ice_info.append((ice_type, ice_label, ice_exc_dcs_path, ice_ion_dcs_path))
    switch_e = min(HIGH_ENERGY_SWITCH_EEV, max_low_water)

    born_available = exc_born_path.exists() and ion_born_dcs_path.exists()
    max_supported_water = max_low_water
    if born_available:
        max_supported_water = min(
            _max_table_energy(exc_born_path),
            _max_table_energy(ion_born_dcs_path),
        )
        if max_supported_water <= switch_e:
            born_available = False
            max_supported_water = max_low_water

    if args.emax > max_supported_water:
        print(
            f"Requested emax={args.emax:.3e} eV exceeds water table coverage; "
            f"water curve will be truncated at {max_supported_water:.3e} eV."
        )
    ice_max_by_type = {}
    for ice_type, ice_label, ice_exc_dcs_path, ice_ion_dcs_path in ice_info:
        max_ice = min(
            _max_table_energy(ice_exc_dcs_path),
            _max_table_energy(ice_ion_dcs_path),
        )
        ice_max_by_type[ice_type] = max_ice
        if args.emax > max_ice:
            print(
                f"Requested emax={args.emax:.3e} eV exceeds {ice_type} ice table coverage; "
                f"curve will be truncated at {max_ice:.3e} eV."
            )

    energy = np.logspace(np.log10(args.emin), np.log10(args.emax), args.nbins)

    exc_eloss_xs_water_low = _excitation_energy_loss_xs(energy, exc_path)
    exc_sigma_xs_water_low = _excitation_sigma_xs(energy, exc_path)
    (
        ion_energy_water_low,
        ion_sigma_xs_water_low,
        ion_eloss_xs_water_low,
    ) = _ion_dcs_moments_table(ion_dcs_path, EMFI_ION_BINDING_EEV)

    if born_available:
        exc_eloss_xs_born = _excitation_energy_loss_xs(energy, exc_born_path)
        exc_sigma_xs_born = _excitation_sigma_xs(energy, exc_born_path)
        (
            ion_energy_born,
            ion_sigma_xs_born,
            ion_eloss_xs_born,
        ) = _ion_dcs_moments_table(ion_born_dcs_path, EMFI_ION_BINDING_EEV)

        exc_eloss_xs_water = _merge_low_high_grid(
            energy, exc_eloss_xs_water_low, exc_eloss_xs_born, switch_e
        )
        exc_sigma_xs_water = _merge_low_high_grid(
            energy, exc_sigma_xs_water_low, exc_sigma_xs_born, switch_e
        )
        ion_energy_water, ion_sigma_xs_water = _merge_low_high_tables(
            ion_energy_water_low,
            ion_sigma_xs_water_low,
            ion_energy_born,
            ion_sigma_xs_born,
            switch_e,
        )
        ion_energy_water_eloss, ion_eloss_xs_water = _merge_low_high_tables(
            ion_energy_water_low,
            ion_eloss_xs_water_low,
            ion_energy_born,
            ion_eloss_xs_born,
            switch_e,
        )
        if not np.array_equal(ion_energy_water, ion_energy_water_eloss):
            ion_eloss_xs_water = _interp_loglog(
                ion_energy_water_eloss, ion_eloss_xs_water, ion_energy_water
            )
    else:
        exc_eloss_xs_water = exc_eloss_xs_water_low
        exc_sigma_xs_water = exc_sigma_xs_water_low
        ion_energy_water = ion_energy_water_low
        ion_sigma_xs_water = ion_sigma_xs_water_low
        ion_eloss_xs_water = ion_eloss_xs_water_low

    use_water_custom = args.water_model == "custom"
    if use_water_custom:
        print("Water stopping-power model: custom (ion + electronic excitation + vib + attachment).")
    else:
        print("Water stopping-power model: dna_opt4-like (ion + electronic excitation only).")

    water_components = _compute_stopping_components(
        energy_grid=energy,
        n_cm3=N_CM3_WATER,
        vib_path=water_vib if use_water_custom else None,
        attach_path=water_attach if use_water_custom else None,
        attach_sigma_scale_cm2=MELTON_SIGMA_SCALE_CM2,
        exc_eloss_xs=exc_eloss_xs_water,
        ion_energy=ion_energy_water,
        ion_eloss_xs=ion_eloss_xs_water,
    )
    water_dedx = water_components["total"]
    water_exc_for_w = water_components["exc"]  # Electronic excitation only (standard W-value)
    water_w, water_dNion_dx = _compute_w_value(
        energy_grid=energy,
        n_cm3=N_CM3_WATER,
        ion_energy=ion_energy_water,
        ion_sigma_xs=ion_sigma_xs_water,
        ion_dedx=water_components["ion"],
        exc_dedx=water_exc_for_w,
    )
    water_sigma_exc = exc_sigma_xs_water
    water_sigma_ion = _interp_loglog(ion_energy_water, ion_sigma_xs_water, energy)
    water_sigma_vib = _vib_sigma_xs(energy, water_vib) if use_water_custom else np.zeros_like(energy)
    water_sigma_att = (
        _attachment_sigma_xs(energy, water_attach, MELTON_SIGMA_SCALE_CM2)
        if use_water_custom
        else np.zeros_like(energy)
    )
    water_sigma_vib_att = water_sigma_vib + water_sigma_att
    water_dom_sp = None
    water_dom_w = None
    water_dom_mfp = None
    water_vib_att_dedx = water_components["vib"] + water_components["attach"]
    water_dom_stack_sp = np.vstack(
        [water_vib_att_dedx, water_components["exc"], water_components["ion"]]
    )
    water_dom_sp = np.full(energy.shape, -1, dtype=int)
    water_sum_sp = np.nansum(water_dom_stack_sp, axis=0)
    valid_dom = np.isfinite(water_sum_sp) & (water_sum_sp > 0)
    if np.any(valid_dom):
        water_dom_sp[valid_dom] = np.argmax(water_dom_stack_sp[:, valid_dom], axis=0)

    water_dom_stack_w = np.vstack([water_vib_att_dedx * 0.0, water_components["exc"], water_components["ion"]])
    water_dom_w = np.full(energy.shape, -1, dtype=int)
    water_sum_w = np.nansum(water_dom_stack_w, axis=0)
    valid_dom = np.isfinite(water_sum_w) & (water_sum_w > 0)
    if np.any(valid_dom):
        water_dom_w[valid_dom] = np.argmax(water_dom_stack_w[:, valid_dom], axis=0)

    water_dom_stack_mfp = np.vstack([water_sigma_vib_att, water_sigma_exc, water_sigma_ion])
    water_dom_mfp = np.full(energy.shape, -1, dtype=int)
    water_sum_mfp = np.nansum(water_dom_stack_mfp, axis=0)
    valid_dom = np.isfinite(water_sum_mfp) & (water_sum_mfp > 0)
    if np.any(valid_dom):
        water_dom_mfp[valid_dom] = np.argmax(water_dom_stack_mfp[:, valid_dom], axis=0)
    water_sigma_ion = _interp_loglog(ion_energy_water, ion_sigma_xs_water, energy)
    water_sigma_inel = water_sigma_ion + exc_sigma_xs_water
    if use_water_custom:
        water_sigma_inel += _vib_sigma_xs(energy, water_vib)
        water_sigma_inel += _attachment_sigma_xs(energy, water_attach, MELTON_SIGMA_SCALE_CM2)
    water_mfp_nm = np.full_like(energy, np.nan, dtype=float)
    valid_sigma = water_sigma_inel > 0
    water_mfp_nm[valid_sigma] = CM_PER_NM / (N_CM3_WATER * water_sigma_inel[valid_sigma])

    ice_series = []
    plasma = plt.get_cmap("plasma", 4)
    plasma_colors = plasma(np.linspace(0.0, 1.0, 4))
    water_color = plasma_colors[0]
    ice_styles = {
        "amorphous": {"color": plasma_colors[1], "ls": "-"},
        "hexagonal": {"color": plasma_colors[2], "ls": "--"},
    }
    for ice_type, ice_label, ice_exc_dcs_path, ice_ion_dcs_path in ice_info:
        phase_density = ICE_DENSITY_BY_TYPE.get(ice_type, ICE_HEXAGONAL_DENSITY_G_CM3)
        ice_density = float(args.rho_ice) if args.rho_ice is not None else float(phase_density)
        n_cm3_ice = N_CM3_WATER * (ice_density / WATER_DENSITY_G_CM3)

        exc_energy_ice, exc_eloss_xs_ice_table = _dcs_energy_loss_xs_table(ice_exc_dcs_path)
        exc_energy_sigma_ice, exc_sigma_xs_ice_table = _dcs_sigma_xs_table(ice_exc_dcs_path)
        ion_energy_ice, ion_sigma_xs_ice, ion_eloss_xs_ice = _ion_dcs_moments_table(
            ice_ion_dcs_path, EMFI_ION_BINDING_EEV
        )
        exc_eloss_xs_ice = _interp_loglog(exc_energy_ice, exc_eloss_xs_ice_table, energy)
        exc_sigma_xs_ice = _interp_loglog(
            exc_energy_sigma_ice, exc_sigma_xs_ice_table, energy
        )

        ice_components = _compute_stopping_components(
            energy_grid=energy,
            n_cm3=n_cm3_ice,
            vib_path=ice_vib,
            attach_path=ice_attach,
            attach_sigma_scale_cm2=MICHAUD_SIGMA_SCALE_CM2,
            exc_eloss_xs=exc_eloss_xs_ice,
            ion_energy=ion_energy_ice,
            ion_eloss_xs=ion_eloss_xs_ice,
        )
        ice_exc_for_w = ice_components["exc"]  # Electronic excitation only (standard W-value)
        ice_w, ice_dNion_dx = _compute_w_value(
            energy_grid=energy,
            n_cm3=n_cm3_ice,
            ion_energy=ion_energy_ice,
            ion_sigma_xs=ion_sigma_xs_ice,
            ion_dedx=ice_components["ion"],
            exc_dedx=ice_exc_for_w,
        )
        ice_sigma_ion = _interp_loglog(ion_energy_ice, ion_sigma_xs_ice, energy)
        ice_sigma_inel = ice_sigma_ion + exc_sigma_xs_ice
        ice_sigma_inel += _vib_sigma_xs(energy, ice_vib)
        ice_sigma_inel += _attachment_sigma_xs(energy, ice_attach, MICHAUD_SIGMA_SCALE_CM2)
        ice_mfp_nm = np.full_like(energy, np.nan, dtype=float)
        valid_sigma = ice_sigma_inel > 0
        ice_mfp_nm[valid_sigma] = CM_PER_NM / (n_cm3_ice * ice_sigma_inel[valid_sigma])

        style = ice_styles.get(ice_type, {"color": "0.6", "ls": "-"})
        ice_series.append(
            {
                "label": f"{ice_type.capitalize()} ice",
                "ice_type": ice_type,
                "dedx": ice_components["total"],
                "w_value": ice_w,
                "exc_sigma": exc_sigma_xs_ice,
                "ion_sigma": ice_sigma_ion,
                "vib_sigma": _vib_sigma_xs(energy, ice_vib),
                "attach_sigma": _attachment_sigma_xs(energy, ice_attach, MICHAUD_SIGMA_SCALE_CM2),
                "vib_dedx": ice_components["vib"],
                "exc_dedx": ice_components["exc"],
                "ion_dedx": ice_components["ion"],
                "attach_dedx": ice_components["attach"],
                "ion_energy": ion_energy_ice,
                "ion_sigma_xs": ion_sigma_xs_ice,
                "dNion_dx": ice_dNion_dx,
                "ion_dedx": ice_components["ion"],
                "exc_dedx_for_w": ice_exc_for_w,
                "inelastic_mfp_nm": ice_mfp_nm,
                "color": style["color"],
                "ls": style["ls"],
                "emax": ice_max_by_type.get(ice_type),
                "density_g_cm3": ice_density,
            }
        )

    _print_w_assumptions()
    water_w_plot_emax = min(max_supported_water, WVALUE_MAX_WATER_EEV)
    _print_w_diagnostics(
        "Water",
        energy,
        ion_energy_water,
        ion_sigma_xs_water,
        water_dNion_dx,
        water_components["ion"],
        water_exc_for_w,
        water_w,
        emax=water_w_plot_emax,
    )
    for series in ice_series:
        ice_w_plot_emax = min(float(series.get("emax", np.inf)), WVALUE_MAX_ICE_EEV)
        _print_w_diagnostics(
            series["label"],
            energy,
            series["ion_energy"],
            series["ion_sigma_xs"],
            series["dNion_dx"],
            series["ion_dedx"],
            series["exc_dedx_for_w"],
            series["w_value"],
            emax=ice_w_plot_emax,
        )

    w_text_curves: dict[str, dict[str, np.ndarray]] = {}
    w_sources = {
        "water": args.wvalue_water,
        "hexagonal": args.wvalue_hex,
        "amorphous": args.wvalue_am,
    }
    for phase, source in w_sources.items():
        resolved = _resolve_wvalue_path(source)
        if resolved is None:
            print(f"W-value file not found for {phase}: {source} (using computed curve).")
            continue
        try:
            w_e, w_v, w_err = _load_wvalue_curve(resolved)
        except Exception as exc:
            print(f"Failed to load W-value file for {phase} ({resolved}): {exc}")
            continue
        if w_e.size == 0:
            print(f"W-value file for {phase} has no valid rows: {resolved}")
            continue
        w_text_curves[phase] = {"energy": w_e, "w": w_v, "w_err": w_err}
        print(f"Using W-value text for {phase}: {resolved}")

    ice_w_emax_by_type = {}
    for ice_type in ice_types:
        max_ice = ice_max_by_type.get(ice_type)
        if max_ice is None:
            continue
        ice_w_emax_by_type[ice_type] = min(float(max_ice), WVALUE_MAX_ICE_EEV)

    am_dom_sp = None
    am_dom_w = None
    am_dom_mfp = None
    hex_dom_sp = None
    hex_dom_w = None
    hex_dom_mfp = None

    def _dominance_from_series(series):
        vib_att_dedx = series["vib_dedx"] + series["attach_dedx"]
        exc_dedx = series["exc_dedx"]
        ion_dedx = series["ion_dedx"]
        stack_sp = np.vstack([vib_att_dedx, exc_dedx, ion_dedx])
        dom_sp = np.full(energy.shape, -1, dtype=int)
        sum_sp = np.nansum(stack_sp, axis=0)
        valid = np.isfinite(sum_sp) & (sum_sp > 0)
        if np.any(valid):
            dom_sp[valid] = np.argmax(stack_sp[:, valid], axis=0)

        stack_w = np.vstack([vib_att_dedx * 0.0, exc_dedx, ion_dedx])
        dom_w = np.full(energy.shape, -1, dtype=int)
        sum_w = np.nansum(stack_w, axis=0)
        valid = np.isfinite(sum_w) & (sum_w > 0)
        if np.any(valid):
            dom_w[valid] = np.argmax(stack_w[:, valid], axis=0)

        vib_att_sigma = series["vib_sigma"] + series["attach_sigma"]
        exc_sigma = series["exc_sigma"]
        ion_sigma = series["ion_sigma"]
        stack_mfp = np.vstack([vib_att_sigma, exc_sigma, ion_sigma])
        dom_mfp = np.full(energy.shape, -1, dtype=int)
        sum_mfp = np.nansum(stack_mfp, axis=0)
        valid = np.isfinite(sum_mfp) & (sum_mfp > 0)
        if np.any(valid):
            dom_mfp[valid] = np.argmax(stack_mfp[:, valid], axis=0)
        return dom_sp, dom_w, dom_mfp

    for series in ice_series:
        ice_type = series.get("ice_type")
        if ice_type == "amorphous":
            am_dom_sp, am_dom_w, am_dom_mfp = _dominance_from_series(series)
        elif ice_type == "hexagonal":
            hex_dom_sp, hex_dom_w, hex_dom_mfp = _dominance_from_series(series)

    _plot(
        energy_grid=energy,
        water_dedx=water_dedx,
        water_w=water_w,
        water_mfp_nm=water_mfp_nm,
        water_color=water_color,
        water_dom_sp=water_dom_sp,
        water_dom_w=water_dom_w,
        water_dom_mfp=water_dom_mfp,
        am_dom_sp=am_dom_sp,
        am_dom_w=am_dom_w,
        am_dom_mfp=am_dom_mfp,
        hex_dom_sp=hex_dom_sp,
        hex_dom_w=hex_dom_w,
        hex_dom_mfp=hex_dom_mfp,
        ice_series=ice_series,
        units=args.units,
        out_path=args.out,
        rho_water=args.rho_water,
        water_emax=max_supported_water,
        water_w_emax=water_w_plot_emax,
        ice_w_emax_by_type=ice_w_emax_by_type,
        w_text_curves=w_text_curves,
    )
    _write_table(
        out_path=args.out_txt,
        energy_grid=energy,
        water_dedx=water_dedx,
        water_w=water_w,
        water_mfp_nm=water_mfp_nm,
        ice_series=ice_series,
    )


if __name__ == "__main__":
    main()
