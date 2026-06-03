#!/usr/bin/env python3
"""
Build flux-scaled Europa lat/lon metrics from per-energy ROOT simulations.

This script is designed around the 30x30 library generated in:
  geant4_projects/dnaphysics-ice/europa_energy_library

Workflow:
1) `prepare`
   - Read `latlon_cell_ranges.csv` + `latlon_energy_scaling.csv.gz`
   - Collapse cells into unique (e_min, e_max) ranges
   - Write a manifest JSON for downstream workers
2) `energy-worker` / `energy-merge` (recommended)
   - Compute per-energy observables independently (one NPZ per energy)
   - Merge per-energy NPZs into energy cache products
3) `build-energy-cache` (legacy one-shot)
   - Stream all per-energy ROOT files in one command
4) `range-worker`
   - Compute one unique-range result (PBS-array friendly)
5) `merge`
   - Map range results back to all lat/lon cells
   - Write database outputs and hemisphere plots
6) `range-count` / `energy-count`
   - Print number of unique ranges from manifest (for PBS submit scripts)
7) `plot-from-npz`
   - Regenerate hemisphere plots directly from saved map NPZ outputs

Unit convention:
- All depth / penetration quantities are in cm.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np

try:
    import uproot  # type: ignore
except Exception:  # pragma: no cover - optional import for non-analysis commands
    uproot = None

from constants import (
    FONT_COURIER,
    FONTSIZE_24,
    PROJECT_ROOT,
    RC_BASE_ELASTIC,
    rcparams_with_fontsize,
)
from root_utils import resolve_root_paths


SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
# Convert MeV/g/s -> MGy/yr.
# 1 MeV = 1.602176634e-13 J, and 1 g = 1e-3 kg.
# 1 MeV g^-1 = 1.602176634e-10 Gy.
MEV_TO_MGY_PER_YEAR = 1.602176634e-10 * SECONDS_PER_YEAR / 1.0e6

EVENT_FIELDS_REQUIRED = [
    "eventID",
    "primaryEnergy",
    "depositedEnergy",
    "escapedBackEnergy",
]
EVENT_FIELDS_OPTIONAL = [
    "escapedForwardEnergy",
    "escapedLateralEnergy",
]
ESCAPE_PARTICLE_ENERGY_FIELDS = [
    "escaped_electron_energy_eV_sum",
    "escaped_photon_energy_eV_sum",
    "backscatter_electron_energy_eV_sum",
    "backscatter_photon_energy_eV_sum",
    "escaped_forward_electron_energy_eV_sum",
    "escaped_forward_photon_energy_eV_sum",
    "escaped_lateral_electron_energy_eV_sum",
    "escaped_lateral_photon_energy_eV_sum",
]


def _safe_float(text: str, default: float = float("nan")) -> float:
    try:
        return float(text)
    except Exception:
        return default


def _safe_int(text: str, default: int = 0) -> int:
    try:
        return int(float(text))
    except Exception:
        return default


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="") as fh:
        return list(csv.DictReader(fh))


def _read_csv_gz_rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def _default_library_dir() -> Path:
    return PROJECT_ROOT / "europa_energy_library"


def _default_local_data_dir() -> Path:
    return Path(__file__).resolve().parent.parent / "data"


def _default_out_dir(library_dir: Path) -> Path:
    work = os.environ.get("WORK", "")
    if work:
        return Path(work) / "dnaphysics-ice" / "europa_energy_library" / "latlon_metrics"
    return library_dir / "latlon_metrics"


def _default_root_dir(library_dir: Path) -> Path:
    user = os.environ.get("USER", "")
    if user:
        candidate = Path("/work") / user / "dnaphysics-ice" / "europa_energy_library" / "root"
        if candidate.exists():
            return candidate
    work = os.environ.get("WORK", "")
    if work:
        candidate = Path(work) / "dnaphysics-ice" / "europa_energy_library" / "root"
        if candidate.exists():
            return candidate
    return library_dir / "root"


def _legacy_depth_edges_cm() -> np.ndarray:
    values_mm = np.round(
        np.array(
            list(np.arange(0.0, 0.05, 0.0001))
            + list(np.arange(0.05, 0.5, 0.001))
            + list(np.arange(0.5, 1.0, 0.01))
            + list(np.arange(1.0, 10.0, 0.1))
            + list(np.arange(10.0, 100.0, 0.5))
            + list(np.arange(100.0, 1000.0, 1.0)),
            dtype=float,
        ),
        8,
    )
    # Legacy grid was defined in mm. Convert to cm for this pipeline.
    values = values_mm * 0.1
    if values.size < 2:
        raise RuntimeError("Legacy depth grid unexpectedly short.")
    delta = values[-1] - values[-2]
    return np.append(values, values[-1] + delta)


def _depth_edges_from_file(path: Path) -> np.ndarray:
    values: list[float] = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        try:
            values.append(float(s))
        except ValueError:
            continue
    if len(values) < 2:
        raise ValueError(f"Need at least 2 depth values in {path}")
    vals = np.unique(np.asarray(values, dtype=float))
    delta = vals[-1] - vals[-2]
    return np.append(vals, vals[-1] + delta)


def _load_manifest(path: Path) -> dict[str, Any]:
    with path.open("r") as fh:
        return json.load(fh)


def _write_manifest(path: Path, obj: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


def _load_run_rows(run_table_csv: Path) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for row in _read_csv_rows(run_table_csv):
        out.append(
            {
                "energy_index": _safe_int(row.get("energy_index", "")),
                "E_center_MeV": _safe_float(row.get("E_center_MeV", "nan")),
                "sim_particles": _safe_int(row.get("sim_particles", "0")),
                "root_basename": str(row.get("root_basename", "")),
                "root_relpath": str(row.get("root_relpath", "")),
                "density_gcm3": _safe_float(row.get("density_gcm3", "nan")),
            }
        )
    out.sort(key=lambda r: int(r["energy_index"]))
    return out


def _resolve_root_files_for_row(
    run_row: dict[str, Any],
    library_dir: Path,
    root_dir_override: Path | None = None,
    work_root_override: Path | None = None,
) -> list[Path]:
    basename = str(run_row["root_basename"])
    rel = str(run_row["root_relpath"])
    candidates: list[Path] = []

    if root_dir_override is not None:
        candidates.append(root_dir_override / f"{basename}.root")

    if rel:
        candidates.append(library_dir / rel)
    candidates.append(library_dir / "root" / f"{basename}.root")

    if work_root_override is not None:
        candidates.append(work_root_override / f"{basename}.root")

    user = os.environ.get("USER", "")
    if user:
        candidates.append(
            Path("/work") / user / "dnaphysics-ice" / "europa_energy_library" / "root" / f"{basename}.root"
        )
    work = os.environ.get("WORK", "")
    if work:
        candidates.append(
            Path(work) / "dnaphysics-ice" / "europa_energy_library" / "root" / f"{basename}.root"
        )

    seen: set[str] = set()
    for cand in candidates:
        key = str(cand)
        if key in seen:
            continue
        seen.add(key)
        resolved = resolve_root_paths(str(cand))
        if resolved:
            return [Path(p) for p in resolved]
    return []


def _require_uproot() -> None:
    if uproot is None:
        raise RuntimeError(
            "uproot is not available in this Python environment. "
            "Install uproot or run with an interpreter that has it."
        )


def _update_event_max_depth(
    target: dict[int, float],
    event_ids: np.ndarray,
    depth_cm: np.ndarray,
) -> None:
    if event_ids.size == 0:
        return
    order = np.argsort(event_ids, kind="mergesort")
    ids_sorted = event_ids[order]
    dep_sorted = depth_cm[order]
    uniq, starts = np.unique(ids_sorted, return_index=True)
    chunk_max = np.maximum.reduceat(dep_sorted, starts)
    for eid, dmax in zip(uniq.tolist(), chunk_max.tolist()):
        eid_i = int(eid)
        dval = float(dmax)
        prev = target.get(eid_i)
        if prev is None or dval > prev:
            target[eid_i] = dval


def _update_event_minmax_depth(
    target_min: dict[int, float],
    target_max: dict[int, float],
    event_ids: np.ndarray,
    depth_cm: np.ndarray,
) -> None:
    if event_ids.size == 0:
        return
    order = np.argsort(event_ids, kind="mergesort")
    ids_sorted = event_ids[order]
    dep_sorted = depth_cm[order]
    uniq, starts = np.unique(ids_sorted, return_index=True)
    chunk_min = np.minimum.reduceat(dep_sorted, starts)
    chunk_max = np.maximum.reduceat(dep_sorted, starts)
    for eid, dmin, dmax in zip(uniq.tolist(), chunk_min.tolist(), chunk_max.tolist()):
        eid_i = int(eid)
        dmin_f = float(dmin)
        dmax_f = float(dmax)
        prev_min = target_min.get(eid_i)
        prev_max = target_max.get(eid_i)
        if prev_min is None or dmin_f < prev_min:
            target_min[eid_i] = dmin_f
        if prev_max is None or dmax_f > prev_max:
            target_max[eid_i] = dmax_f


def _accumulate_event_profile_histograms(
    target_profiles: dict[int, np.ndarray],
    event_ids: np.ndarray,
    depth_cm: np.ndarray,
    deposited_eV: np.ndarray,
    depth_edges_cm: np.ndarray,
) -> None:
    if event_ids.size == 0:
        return
    n_depth = int(depth_edges_cm.size - 1)
    bin_idx = np.searchsorted(depth_edges_cm, depth_cm, side="right") - 1
    valid = (
        np.isfinite(depth_cm)
        & np.isfinite(deposited_eV)
        & (bin_idx >= 0)
        & (bin_idx < n_depth)
    )
    if not np.any(valid):
        return
    event_ids_i64 = np.asarray(event_ids[valid], dtype=np.int64)
    bin_idx_i64 = np.asarray(bin_idx[valid], dtype=np.int64)
    deposited_f64 = np.asarray(deposited_eV[valid], dtype=np.float64)
    composite = event_ids_i64 * np.int64(n_depth) + bin_idx_i64
    order = np.argsort(composite, kind="mergesort")
    composite_sorted = composite[order]
    deposited_sorted = deposited_f64[order]
    uniq, starts = np.unique(composite_sorted, return_index=True)
    chunk_sum = np.add.reduceat(deposited_sorted, starts)
    event_sorted = (uniq // np.int64(n_depth)).astype(np.int64)
    bin_sorted = (uniq % np.int64(n_depth)).astype(np.int64)
    for eid, depth_bin, dep_sum in zip(
        event_sorted.tolist(),
        bin_sorted.tolist(),
        chunk_sum.tolist(),
    ):
        prof = target_profiles.get(int(eid))
        if prof is None:
            prof = np.zeros(n_depth, dtype=np.float64)
            target_profiles[int(eid)] = prof
        prof[int(depth_bin)] += float(dep_sum)


def _accumulate_profile_sumsq_from_events(
    target_sumsq: np.ndarray,
    event_profiles: dict[int, np.ndarray],
) -> None:
    if not event_profiles:
        return
    profiles_arr = np.stack(list(event_profiles.values()), axis=0)
    target_sumsq += np.einsum("ij,ij->j", profiles_arr, profiles_arr)


def _make_group_ids(event_ids: np.ndarray, track_ids: np.ndarray | None) -> np.ndarray:
    event_ids_i64 = np.asarray(event_ids, dtype=np.int64)
    if track_ids is None:
        return event_ids_i64
    track_ids_i64 = np.asarray(track_ids, dtype=np.int64)
    # Build a stable composite key so penetration is measured per track, not per event.
    return (
        (event_ids_i64.astype(np.int64) << 32)
        ^ (track_ids_i64.astype(np.int64) & np.int64(0xFFFFFFFF))
    )


def _sum_event_tree_budget(
    root_files: list[Path],
    step_size: str,
) -> dict[str, Any]:
    _require_uproot()
    n_events = 0
    primary_e_sum = 0.0
    deposited_e_sum = 0.0
    backscatter_e_sum = 0.0
    forward_e_sum = 0.0
    lateral_e_sum = 0.0
    saw_required = False
    saw_forward = False
    saw_lateral = False

    for root_path in root_files:
        with uproot.open(root_path) as rootf:
            if "event" not in rootf:
                continue
            tree = rootf["event"]
            keys = set(tree.keys())
            if not all(k in keys for k in EVENT_FIELDS_REQUIRED):
                continue
            saw_required = True
            fields = list(EVENT_FIELDS_REQUIRED)
            if "escapedForwardEnergy" in keys:
                fields.append("escapedForwardEnergy")
                saw_forward = True
            if "escapedLateralEnergy" in keys:
                fields.append("escapedLateralEnergy")
                saw_lateral = True
            for chunk in tree.iterate(fields, library="np", step_size=step_size):
                n_events += int(np.asarray(chunk["eventID"]).size)
                primary_e_sum += float(np.sum(np.asarray(chunk["primaryEnergy"], dtype=np.float64)))
                deposited_e_sum += float(np.sum(np.asarray(chunk["depositedEnergy"], dtype=np.float64)))
                backscatter_e_sum += float(np.sum(np.asarray(chunk["escapedBackEnergy"], dtype=np.float64)))
                if "escapedForwardEnergy" in chunk:
                    forward_e_sum += float(
                        np.sum(np.asarray(chunk["escapedForwardEnergy"], dtype=np.float64))
                    )
                if "escapedLateralEnergy" in chunk:
                    lateral_e_sum += float(
                        np.sum(np.asarray(chunk["escapedLateralEnergy"], dtype=np.float64))
                    )

    return {
        "ok": bool(saw_required and n_events > 0),
        "n_events": int(n_events),
        "primary_energy_eV_sum": float(primary_e_sum),
        "deposited_energy_eV_sum": float(deposited_e_sum),
        "backscatter_energy_eV_sum": float(backscatter_e_sum),
        "escaped_forward_energy_eV_sum": float(forward_e_sum) if saw_forward else float("nan"),
        "escaped_lateral_energy_eV_sum": float(lateral_e_sum) if saw_lateral else float("nan"),
    }


def _analyze_single_energy(
    run_row: dict[str, Any],
    root_files: list[Path],
    depth_edges_cm: np.ndarray,
    step_size: str,
    depth_origin_cm: float,
) -> tuple[
    dict[str, Any],
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    energy_index = int(run_row["energy_index"])
    e_center_mev = float(run_row["E_center_MeV"])
    sim_particles = int(run_row["sim_particles"])
    n_depth = depth_edges_cm.size - 1
    profile_eV = np.zeros(n_depth, dtype=np.float64)
    profile_sumsq_eV2 = np.zeros(n_depth, dtype=np.float64)
    primary_ke_sum_profile_eV = np.zeros(n_depth, dtype=np.float64)
    primary_ke_count_profile = np.zeros(n_depth, dtype=np.float64)
    secondary_ke_sum_profile_eV = np.zeros(n_depth, dtype=np.float64)
    secondary_ke_count_profile = np.zeros(n_depth, dtype=np.float64)

    metrics: dict[str, Any] = {
        "energy_index": energy_index,
        "E_center_MeV": e_center_mev,
        "sim_particles": sim_particles,
        "n_root_files": len(root_files),
        "status": "ok_event_tree",
        "budget_mode": "event_tree",
        "event_count": 0,
        "primary_energy_eV_sum": float("nan"),
        "deposited_energy_eV_sum": float("nan"),
        "backscatter_energy_eV_sum": float("nan"),
        "escaped_forward_energy_eV_sum": float("nan"),
        "escaped_lateral_energy_eV_sum": float("nan"),
        "escaped_electron_energy_eV_sum": float("nan"),
        "escaped_photon_energy_eV_sum": float("nan"),
        "backscatter_electron_energy_eV_sum": float("nan"),
        "backscatter_photon_energy_eV_sum": float("nan"),
        "escaped_forward_electron_energy_eV_sum": float("nan"),
        "escaped_forward_photon_energy_eV_sum": float("nan"),
        "escaped_lateral_electron_energy_eV_sum": float("nan"),
        "escaped_lateral_photon_energy_eV_sum": float("nan"),
        "primary_pen_sum_cm": 0.0,
        "primary_pen_sumsq_cm2": 0.0,
        "primary_pen_count": 0,
        "secondary_pen_sum_cm": 0.0,
        "secondary_pen_sumsq_cm2": 0.0,
        "secondary_pen_count": 0,
        "secondary_ke_sum_eV": 0.0,
        "secondary_ke_sumsq_eV2": float("nan"),
        "secondary_ke_count": 0,
        "error": "",
    }

    if not root_files:
        metrics["status"] = "missing_root"
        metrics["budget_mode"] = "none"
        metrics["error"] = "No ROOT file found"
        return (
            metrics,
            profile_eV,
            profile_sumsq_eV2,
            primary_ke_sum_profile_eV,
            primary_ke_count_profile,
            secondary_ke_sum_profile_eV,
            secondary_ke_count_profile,
        )

    # Event-tree budget (preferred)
    event_budget = _sum_event_tree_budget(root_files, step_size)

    # Step-tree quantities
    total_edep_step_eV = 0.0
    primary_pen_sum_cm = 0.0
    primary_pen_sumsq_cm2 = 0.0
    primary_pen_count = 0
    secondary_pen_sum_cm = 0.0
    secondary_pen_sumsq_cm2 = 0.0
    secondary_pen_count = 0
    secondary_ke_sum_eV = 0.0
    secondary_ke_count = 0
    saw_step = False
    saw_step_fields = False

    for root_path in root_files:
        with uproot.open(root_path) as rootf:
            if "step" not in rootf:
                continue
            saw_step = True
            step_tree = rootf["step"]
            keys = set(step_tree.keys())
            if "z" not in keys or "totalEnergyDeposit" not in keys:
                continue
            saw_step_fields = True
            primary_event_min: dict[int, float] = {}
            primary_event_max: dict[int, float] = {}
            secondary_event_min: dict[int, float] = {}
            secondary_event_max: dict[int, float] = {}
            event_profiles_eV: dict[int, np.ndarray] = {}
            step_fields = ["z", "totalEnergyDeposit"]
            if "eventID" in keys:
                step_fields.append("eventID")
            if "parentID" in keys:
                step_fields.append("parentID")
            if "flagParticle" in keys:
                step_fields.append("flagParticle")
            if "kineticEnergy" in keys:
                step_fields.append("kineticEnergy")
            if "trackID" in keys:
                step_fields.append("trackID")

            for chunk in step_tree.iterate(step_fields, library="np", step_size=step_size):
                z_cm = np.asarray(chunk["z"], dtype=np.float64) * 1.0e-7 - depth_origin_cm
                z_cm = np.maximum(z_cm, 0.0)
                edep = np.asarray(chunk["totalEnergyDeposit"], dtype=np.float64)
                total_edep_step_eV += float(np.sum(edep[np.isfinite(edep)]))

                if "flagParticle" in chunk:
                    flag = np.asarray(chunk["flagParticle"], dtype=np.int64)
                    m_e = flag == 1
                else:
                    m_e = np.ones_like(z_cm, dtype=bool)

                m_hist = m_e & np.isfinite(z_cm) & np.isfinite(edep)
                if np.any(m_hist):
                    hist, _ = np.histogram(z_cm[m_hist], bins=depth_edges_cm, weights=edep[m_hist])
                    profile_eV += hist
                    if "eventID" in chunk:
                        event_id_hist = np.asarray(chunk["eventID"], dtype=np.int64)
                        _accumulate_event_profile_histograms(
                            event_profiles_eV,
                            event_id_hist[m_hist],
                            z_cm[m_hist],
                            edep[m_hist],
                            depth_edges_cm,
                        )

                if "parentID" not in chunk:
                    continue

                parent_id = np.asarray(chunk["parentID"], dtype=np.int64)
                track_id = np.asarray(chunk["trackID"], dtype=np.int64) if "trackID" in chunk else None
                event_id = np.asarray(chunk["eventID"], dtype=np.int64) if "eventID" in chunk else None

                m_primary = m_e & (parent_id == 0)
                if track_id is not None:
                    m_primary &= track_id == 1

                m_secondary = m_e & (parent_id > 0)
                if event_id is not None:
                    if np.any(m_primary):
                        primary_group_id = _make_group_ids(
                            event_id[m_primary],
                            track_id[m_primary] if track_id is not None else None,
                        )
                        _update_event_minmax_depth(
                            primary_event_min,
                            primary_event_max,
                            primary_group_id,
                            z_cm[m_primary],
                        )
                    if np.any(m_secondary):
                        secondary_group_id = _make_group_ids(
                            event_id[m_secondary],
                            track_id[m_secondary] if track_id is not None else None,
                        )
                        _update_event_minmax_depth(
                            secondary_event_min,
                            secondary_event_max,
                            secondary_group_id,
                            z_cm[m_secondary],
                        )

                if "kineticEnergy" in chunk:
                    ke = np.asarray(chunk["kineticEnergy"], dtype=np.float64)
                    m_ke_primary = m_primary & np.isfinite(ke) & np.isfinite(z_cm)
                    if np.any(m_ke_primary):
                        ke_sum_hist, _ = np.histogram(
                            z_cm[m_ke_primary],
                            bins=depth_edges_cm,
                            weights=ke[m_ke_primary],
                        )
                        ke_cnt_hist, _ = np.histogram(
                            z_cm[m_ke_primary],
                            bins=depth_edges_cm,
                        )
                        primary_ke_sum_profile_eV += ke_sum_hist
                        primary_ke_count_profile += ke_cnt_hist

                    m_ke_secondary = m_secondary & np.isfinite(ke) & np.isfinite(z_cm)
                    if np.any(m_ke_secondary):
                        ke_sum_hist, _ = np.histogram(
                            z_cm[m_ke_secondary],
                            bins=depth_edges_cm,
                            weights=ke[m_ke_secondary],
                        )
                        ke_cnt_hist, _ = np.histogram(
                            z_cm[m_ke_secondary],
                            bins=depth_edges_cm,
                        )
                        secondary_ke_sum_profile_eV += ke_sum_hist
                        secondary_ke_count_profile += ke_cnt_hist
                        sec_vals = ke[m_ke_secondary]
                        secondary_ke_sum_eV += float(np.sum(sec_vals[np.isfinite(sec_vals)]))
                        secondary_ke_count += int(np.count_nonzero(np.isfinite(sec_vals)))

            _accumulate_profile_sumsq_from_events(profile_sumsq_eV2, event_profiles_eV)
            primary_dz_vals = []
            for eid, zmax in primary_event_max.items():
                zmin = primary_event_min.get(eid)
                if zmin is None:
                    continue
                dz = float(zmax) - float(zmin)
                if np.isfinite(dz) and dz >= 0.0:
                    primary_dz_vals.append(dz)
            if primary_dz_vals:
                primary_dz_arr = np.asarray(primary_dz_vals, dtype=np.float64)
                primary_pen_sum_cm += float(np.sum(primary_dz_arr))
                primary_pen_sumsq_cm2 += float(np.sum(primary_dz_arr * primary_dz_arr))
                primary_pen_count += int(primary_dz_arr.size)

            secondary_dz_vals = []
            for eid, zmax in secondary_event_max.items():
                zmin = secondary_event_min.get(eid)
                if zmin is None:
                    continue
                dz = float(zmax) - float(zmin)
                if np.isfinite(dz) and dz >= 0.0:
                    secondary_dz_vals.append(dz)
            if secondary_dz_vals:
                secondary_dz_arr = np.asarray(secondary_dz_vals, dtype=np.float64)
                secondary_pen_sum_cm += float(np.sum(secondary_dz_arr))
                secondary_pen_sumsq_cm2 += float(np.sum(secondary_dz_arr * secondary_dz_arr))
                secondary_pen_count += int(secondary_dz_arr.size)

    if not saw_step:
        metrics["status"] = "missing_step_tree"
        metrics["budget_mode"] = "none"
        metrics["error"] = "step tree missing in all ROOT files"
        return (
            metrics,
            profile_eV,
            profile_sumsq_eV2,
            primary_ke_sum_profile_eV,
            primary_ke_count_profile,
            secondary_ke_sum_profile_eV,
            secondary_ke_count_profile,
        )
    if not saw_step_fields:
        metrics["status"] = "missing_step_fields"
        metrics["budget_mode"] = "none"
        metrics["error"] = "Required step fields missing (need z,totalEnergyDeposit)"
        return (
            metrics,
            profile_eV,
            profile_sumsq_eV2,
            primary_ke_sum_profile_eV,
            primary_ke_count_profile,
            secondary_ke_sum_profile_eV,
            secondary_ke_count_profile,
        )

    if bool(event_budget["ok"]):
        metrics["status"] = "ok_event_tree"
        metrics["budget_mode"] = "event_tree"
        metrics["event_count"] = int(event_budget["n_events"])
        metrics["primary_energy_eV_sum"] = float(event_budget["primary_energy_eV_sum"])
        metrics["deposited_energy_eV_sum"] = float(event_budget["deposited_energy_eV_sum"])
        metrics["backscatter_energy_eV_sum"] = float(event_budget["backscatter_energy_eV_sum"])
        metrics["escaped_forward_energy_eV_sum"] = float(event_budget["escaped_forward_energy_eV_sum"])
        metrics["escaped_lateral_energy_eV_sum"] = float(event_budget["escaped_lateral_energy_eV_sum"])
        for field in ESCAPE_PARTICLE_ENERGY_FIELDS:
            if field in event_budget:
                metrics[field] = float(event_budget[field])
    else:
        # Fallback for older files that do not contain full event budget.
        metrics["status"] = "ok_fallback_step_tree"
        metrics["budget_mode"] = "step_only"
        metrics["event_count"] = int(sim_particles)
        metrics["primary_energy_eV_sum"] = e_center_mev * 1.0e6 * float(sim_particles)
        metrics["deposited_energy_eV_sum"] = total_edep_step_eV
        metrics["backscatter_energy_eV_sum"] = float("nan")
        metrics["escaped_forward_energy_eV_sum"] = float("nan")
        metrics["escaped_lateral_energy_eV_sum"] = float("nan")
        for field in ESCAPE_PARTICLE_ENERGY_FIELDS:
            metrics[field] = float("nan")

    metrics["primary_pen_sum_cm"] = float(primary_pen_sum_cm)
    metrics["primary_pen_sumsq_cm2"] = float(primary_pen_sumsq_cm2)
    metrics["primary_pen_count"] = int(primary_pen_count)
    metrics["secondary_pen_sum_cm"] = float(secondary_pen_sum_cm)
    metrics["secondary_pen_sumsq_cm2"] = float(secondary_pen_sumsq_cm2)
    metrics["secondary_pen_count"] = int(secondary_pen_count)
    metrics["secondary_ke_sum_eV"] = float(secondary_ke_sum_eV)
    metrics["secondary_ke_count"] = int(secondary_ke_count)

    return (
        metrics,
        profile_eV,
        profile_sumsq_eV2,
        primary_ke_sum_profile_eV,
        primary_ke_count_profile,
        secondary_ke_sum_profile_eV,
        secondary_ke_count_profile,
    )


def _save_energy_cache(
    metrics_rows: list[dict[str, Any]],
    deposited_profiles_eV: np.ndarray,
    deposited_profile_sumsq_eV2: np.ndarray,
    primary_ke_sum_profiles_eV: np.ndarray,
    primary_ke_count_profiles: np.ndarray,
    secondary_ke_sum_profiles_eV: np.ndarray,
    secondary_ke_count_profiles: np.ndarray,
    depth_edges_cm: np.ndarray,
    energy_metrics_csv: Path,
    energy_profiles_npz: Path,
    energy_feedback_npz: Path,
    energy_feedback_csv: Path,
    energy_observables_npz: Path,
) -> None:
    energy_metrics_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "energy_index",
        "E_center_MeV",
        "sim_particles",
        "n_root_files",
        "status",
        "budget_mode",
        "event_count",
        "primary_energy_eV_sum",
        "deposited_energy_eV_sum",
        "backscatter_energy_eV_sum",
        "escaped_forward_energy_eV_sum",
        "escaped_lateral_energy_eV_sum",
        *ESCAPE_PARTICLE_ENERGY_FIELDS,
        "primary_pen_sum_cm",
        "primary_pen_sumsq_cm2",
        "primary_pen_count",
        "secondary_pen_sum_cm",
        "secondary_pen_sumsq_cm2",
        "secondary_pen_count",
        "secondary_ke_sum_eV",
        "secondary_ke_sumsq_eV2",
        "secondary_ke_count",
        "error",
    ]
    with energy_metrics_csv.open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=fieldnames)
        wr.writeheader()
        wr.writerows(metrics_rows)

    energy_idx = np.asarray([int(r["energy_index"]) for r in metrics_rows], dtype=np.int64)
    np.savez_compressed(
        energy_profiles_npz,
        energy_index=energy_idx,
        depth_edges_cm=np.asarray(depth_edges_cm, dtype=np.float64),
        deposited_profile_eV=np.asarray(deposited_profiles_eV, dtype=np.float64),
        deposited_profile_sumsq_eV2=np.asarray(deposited_profile_sumsq_eV2, dtype=np.float64),
        primary_ke_sum_profile_eV=np.asarray(primary_ke_sum_profiles_eV, dtype=np.float64),
        primary_ke_count_profile=np.asarray(primary_ke_count_profiles, dtype=np.float64),
        secondary_ke_sum_profile_eV=np.asarray(secondary_ke_sum_profiles_eV, dtype=np.float64),
        secondary_ke_count_profile=np.asarray(secondary_ke_count_profiles, dtype=np.float64),
    )

    e_center = np.asarray([_safe_float(str(r.get("E_center_MeV", "nan"))) for r in metrics_rows], dtype=np.float64)
    sim_particles = np.asarray([_safe_int(str(r.get("sim_particles", "0"))) for r in metrics_rows], dtype=np.int64)
    status = np.asarray([str(r.get("status", "")) for r in metrics_rows], dtype="<U64")
    budget_mode = np.asarray([str(r.get("budget_mode", "")) for r in metrics_rows], dtype="<U64")
    event_count = np.asarray([_safe_int(str(r.get("event_count", "0"))) for r in metrics_rows], dtype=np.int64)

    primary_eV = np.asarray(
        [_safe_float(str(r.get("primary_energy_eV_sum", "nan"))) for r in metrics_rows],
        dtype=np.float64,
    )
    deposited_eV = np.asarray(
        [_safe_float(str(r.get("deposited_energy_eV_sum", "nan"))) for r in metrics_rows],
        dtype=np.float64,
    )
    backscatter_eV = np.asarray(
        [_safe_float(str(r.get("backscatter_energy_eV_sum", "nan"))) for r in metrics_rows],
        dtype=np.float64,
    )
    forward_eV = np.asarray(
        [_safe_float(str(r.get("escaped_forward_energy_eV_sum", "nan"))) for r in metrics_rows],
        dtype=np.float64,
    )
    lateral_eV = np.asarray(
        [_safe_float(str(r.get("escaped_lateral_energy_eV_sum", "nan"))) for r in metrics_rows],
        dtype=np.float64,
    )
    split_escape_eV = {
        field: np.asarray(
            [_safe_float(str(r.get(field, "nan"))) for r in metrics_rows],
            dtype=np.float64,
        )
        for field in ESCAPE_PARTICLE_ENERGY_FIELDS
    }
    deflected_eV = np.where(
        np.isfinite(primary_eV) & np.isfinite(deposited_eV),
        primary_eV - deposited_eV,
        np.nan,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        deposited_fraction = np.where(primary_eV > 0.0, deposited_eV / primary_eV, np.nan)
        backscatter_fraction = np.where(primary_eV > 0.0, backscatter_eV / primary_eV, np.nan)
        forward_fraction = np.where(primary_eV > 0.0, forward_eV / primary_eV, np.nan)
        lateral_fraction = np.where(primary_eV > 0.0, lateral_eV / primary_eV, np.nan)
        deflected_fraction = np.where(primary_eV > 0.0, deflected_eV / primary_eV, np.nan)

    primary_pen_sum_cm = np.asarray(
        [_safe_float(str(r.get("primary_pen_sum_cm", "0"))) for r in metrics_rows],
        dtype=np.float64,
    )
    primary_pen_sumsq_cm2 = np.asarray(
        [_safe_float(str(r.get("primary_pen_sumsq_cm2", "0"))) for r in metrics_rows],
        dtype=np.float64,
    )
    primary_pen_cnt = np.asarray(
        [_safe_float(str(r.get("primary_pen_count", "0"))) for r in metrics_rows],
        dtype=np.float64,
    )
    secondary_pen_sum_cm = np.asarray(
        [_safe_float(str(r.get("secondary_pen_sum_cm", "0"))) for r in metrics_rows],
        dtype=np.float64,
    )
    secondary_pen_sumsq_cm2 = np.asarray(
        [_safe_float(str(r.get("secondary_pen_sumsq_cm2", "0"))) for r in metrics_rows],
        dtype=np.float64,
    )
    secondary_pen_cnt = np.asarray(
        [_safe_float(str(r.get("secondary_pen_count", "0"))) for r in metrics_rows],
        dtype=np.float64,
    )
    secondary_ke_sum_eV = np.asarray(
        [_safe_float(str(r.get("secondary_ke_sum_eV", "0"))) for r in metrics_rows],
        dtype=np.float64,
    )
    secondary_ke_sumsq_eV2 = np.asarray(
        [_safe_float(str(r.get("secondary_ke_sumsq_eV2", "nan"))) for r in metrics_rows],
        dtype=np.float64,
    )
    secondary_ke_cnt = np.asarray(
        [_safe_float(str(r.get("secondary_ke_count", "0"))) for r in metrics_rows],
        dtype=np.float64,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_primary_pen_cm = np.where(primary_pen_cnt > 0.0, primary_pen_sum_cm / primary_pen_cnt, np.nan)
        mean_secondary_pen_cm = np.where(secondary_pen_cnt > 0.0, secondary_pen_sum_cm / secondary_pen_cnt, np.nan)
        mean_secondary_ke_mev = np.where(secondary_ke_cnt > 0.0, (secondary_ke_sum_eV / secondary_ke_cnt) * 1.0e-6, np.nan)
        primary_pen_var_cm2 = np.where(
            primary_pen_cnt > 0.0,
            primary_pen_sumsq_cm2 / primary_pen_cnt - np.square(mean_primary_pen_cm),
            np.nan,
        )
        secondary_pen_var_cm2 = np.where(
            secondary_pen_cnt > 0.0,
            secondary_pen_sumsq_cm2 / secondary_pen_cnt - np.square(mean_secondary_pen_cm),
            np.nan,
        )
    primary_pen_std_cm = np.sqrt(np.maximum(primary_pen_var_cm2, 0.0))
    secondary_pen_std_cm = np.sqrt(np.maximum(secondary_pen_var_cm2, 0.0))

    with np.errstate(divide="ignore", invalid="ignore"):
        mean_primary_ke_profile_eV = np.divide(
            primary_ke_sum_profiles_eV,
            primary_ke_count_profiles,
            out=np.full_like(primary_ke_sum_profiles_eV, np.nan),
            where=primary_ke_count_profiles > 0.0,
        )
        mean_secondary_ke_profile_eV = np.divide(
            secondary_ke_sum_profiles_eV,
            secondary_ke_count_profiles,
            out=np.full_like(secondary_ke_sum_profiles_eV, np.nan),
            where=secondary_ke_count_profiles > 0.0,
        )
    dep_cumsum = np.cumsum(deposited_profiles_eV, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cumulative_deposited_fraction_profile = np.where(
            primary_eV[:, np.newaxis] > 0.0,
            dep_cumsum / primary_eV[:, np.newaxis],
            np.nan,
        )

    np.savez_compressed(
        energy_feedback_npz,
        energy_index=energy_idx,
        E_center_MeV=e_center,
        sim_particles=sim_particles,
        event_count=event_count,
        status=status,
        budget_mode=budget_mode,
        primary_energy_MeV_sum=primary_eV * 1.0e-6,
        deposited_energy_MeV_sum=deposited_eV * 1.0e-6,
        backscatter_energy_MeV_sum=backscatter_eV * 1.0e-6,
        escaped_forward_energy_MeV_sum=forward_eV * 1.0e-6,
        escaped_lateral_energy_MeV_sum=lateral_eV * 1.0e-6,
        escaped_electron_energy_MeV_sum=split_escape_eV["escaped_electron_energy_eV_sum"] * 1.0e-6,
        escaped_photon_energy_MeV_sum=split_escape_eV["escaped_photon_energy_eV_sum"] * 1.0e-6,
        backscatter_electron_energy_MeV_sum=split_escape_eV["backscatter_electron_energy_eV_sum"] * 1.0e-6,
        backscatter_photon_energy_MeV_sum=split_escape_eV["backscatter_photon_energy_eV_sum"] * 1.0e-6,
        escaped_forward_electron_energy_MeV_sum=split_escape_eV[
            "escaped_forward_electron_energy_eV_sum"
        ]
        * 1.0e-6,
        escaped_forward_photon_energy_MeV_sum=split_escape_eV[
            "escaped_forward_photon_energy_eV_sum"
        ]
        * 1.0e-6,
        escaped_lateral_electron_energy_MeV_sum=split_escape_eV[
            "escaped_lateral_electron_energy_eV_sum"
        ]
        * 1.0e-6,
        escaped_lateral_photon_energy_MeV_sum=split_escape_eV[
            "escaped_lateral_photon_energy_eV_sum"
        ]
        * 1.0e-6,
        deflected_energy_MeV_sum=deflected_eV * 1.0e-6,
        deposited_fraction=deposited_fraction,
        backscatter_fraction=backscatter_fraction,
        escaped_forward_fraction=forward_fraction,
        escaped_lateral_fraction=lateral_fraction,
        deflected_fraction=deflected_fraction,
        mean_primary_penetration_cm=mean_primary_pen_cm,
        std_primary_penetration_cm=primary_pen_std_cm,
        mean_secondary_penetration_cm=mean_secondary_pen_cm,
        std_secondary_penetration_cm=secondary_pen_std_cm,
        mean_secondary_ke_mev=mean_secondary_ke_mev,
    )

    with energy_feedback_csv.open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(
            [
                "energy_index",
                "E_center_MeV",
                "sim_particles",
                "event_count",
                "status",
                "budget_mode",
                "primary_energy_MeV_sum",
                "deposited_energy_MeV_sum",
                "backscatter_energy_MeV_sum",
                "escaped_forward_energy_MeV_sum",
                "escaped_lateral_energy_MeV_sum",
                "escaped_electron_energy_MeV_sum",
                "escaped_photon_energy_MeV_sum",
                "backscatter_electron_energy_MeV_sum",
                "backscatter_photon_energy_MeV_sum",
                "escaped_forward_electron_energy_MeV_sum",
                "escaped_forward_photon_energy_MeV_sum",
                "escaped_lateral_electron_energy_MeV_sum",
                "escaped_lateral_photon_energy_MeV_sum",
                "deflected_energy_MeV_sum",
                "deposited_fraction",
                "backscatter_fraction",
                "escaped_forward_fraction",
                "escaped_lateral_fraction",
                "deflected_fraction",
                "mean_primary_penetration_cm",
                "std_primary_penetration_cm",
                "mean_secondary_penetration_cm",
                "std_secondary_penetration_cm",
                "mean_secondary_ke_mev",
            ]
        )
        for i in range(energy_idx.size):
            wr.writerow(
                [
                    int(energy_idx[i]),
                    f"{e_center[i]:.9g}" if np.isfinite(e_center[i]) else "nan",
                    int(sim_particles[i]),
                    int(event_count[i]),
                    str(status[i]),
                    str(budget_mode[i]),
                    f"{(primary_eV[i] * 1.0e-6):.9g}" if np.isfinite(primary_eV[i]) else "nan",
                    f"{(deposited_eV[i] * 1.0e-6):.9g}" if np.isfinite(deposited_eV[i]) else "nan",
                    f"{(backscatter_eV[i] * 1.0e-6):.9g}" if np.isfinite(backscatter_eV[i]) else "nan",
                    f"{(forward_eV[i] * 1.0e-6):.9g}" if np.isfinite(forward_eV[i]) else "nan",
                    f"{(lateral_eV[i] * 1.0e-6):.9g}" if np.isfinite(lateral_eV[i]) else "nan",
                    f"{(split_escape_eV['escaped_electron_energy_eV_sum'][i] * 1.0e-6):.9g}"
                    if np.isfinite(split_escape_eV["escaped_electron_energy_eV_sum"][i])
                    else "nan",
                    f"{(split_escape_eV['escaped_photon_energy_eV_sum'][i] * 1.0e-6):.9g}"
                    if np.isfinite(split_escape_eV["escaped_photon_energy_eV_sum"][i])
                    else "nan",
                    f"{(split_escape_eV['backscatter_electron_energy_eV_sum'][i] * 1.0e-6):.9g}"
                    if np.isfinite(split_escape_eV["backscatter_electron_energy_eV_sum"][i])
                    else "nan",
                    f"{(split_escape_eV['backscatter_photon_energy_eV_sum'][i] * 1.0e-6):.9g}"
                    if np.isfinite(split_escape_eV["backscatter_photon_energy_eV_sum"][i])
                    else "nan",
                    f"{(split_escape_eV['escaped_forward_electron_energy_eV_sum'][i] * 1.0e-6):.9g}"
                    if np.isfinite(split_escape_eV["escaped_forward_electron_energy_eV_sum"][i])
                    else "nan",
                    f"{(split_escape_eV['escaped_forward_photon_energy_eV_sum'][i] * 1.0e-6):.9g}"
                    if np.isfinite(split_escape_eV["escaped_forward_photon_energy_eV_sum"][i])
                    else "nan",
                    f"{(split_escape_eV['escaped_lateral_electron_energy_eV_sum'][i] * 1.0e-6):.9g}"
                    if np.isfinite(split_escape_eV["escaped_lateral_electron_energy_eV_sum"][i])
                    else "nan",
                    f"{(split_escape_eV['escaped_lateral_photon_energy_eV_sum'][i] * 1.0e-6):.9g}"
                    if np.isfinite(split_escape_eV["escaped_lateral_photon_energy_eV_sum"][i])
                    else "nan",
                    f"{(deflected_eV[i] * 1.0e-6):.9g}" if np.isfinite(deflected_eV[i]) else "nan",
                    f"{deposited_fraction[i]:.9g}" if np.isfinite(deposited_fraction[i]) else "nan",
                    f"{backscatter_fraction[i]:.9g}" if np.isfinite(backscatter_fraction[i]) else "nan",
                    f"{forward_fraction[i]:.9g}" if np.isfinite(forward_fraction[i]) else "nan",
                    f"{lateral_fraction[i]:.9g}" if np.isfinite(lateral_fraction[i]) else "nan",
                    f"{deflected_fraction[i]:.9g}" if np.isfinite(deflected_fraction[i]) else "nan",
                    f"{mean_primary_pen_cm[i]:.9g}" if np.isfinite(mean_primary_pen_cm[i]) else "nan",
                    f"{primary_pen_std_cm[i]:.9g}" if np.isfinite(primary_pen_std_cm[i]) else "nan",
                    f"{mean_secondary_pen_cm[i]:.9g}" if np.isfinite(mean_secondary_pen_cm[i]) else "nan",
                    f"{secondary_pen_std_cm[i]:.9g}" if np.isfinite(secondary_pen_std_cm[i]) else "nan",
                    f"{mean_secondary_ke_mev[i]:.9g}" if np.isfinite(mean_secondary_ke_mev[i]) else "nan",
                ]
            )

    np.savez_compressed(
        energy_observables_npz,
        energy_index=energy_idx,
        depth_edges_cm=np.asarray(depth_edges_cm, dtype=np.float64),
        deposited_profile_eV=np.asarray(deposited_profiles_eV, dtype=np.float64),
        deposited_profile_sumsq_eV2=np.asarray(deposited_profile_sumsq_eV2, dtype=np.float64),
        primary_ke_sum_profile_eV=np.asarray(primary_ke_sum_profiles_eV, dtype=np.float64),
        primary_ke_count_profile=np.asarray(primary_ke_count_profiles, dtype=np.float64),
        secondary_ke_sum_profile_eV=np.asarray(secondary_ke_sum_profiles_eV, dtype=np.float64),
        secondary_ke_count_profile=np.asarray(secondary_ke_count_profiles, dtype=np.float64),
        mean_primary_ke_profile_eV=np.asarray(mean_primary_ke_profile_eV, dtype=np.float64),
        mean_secondary_ke_profile_eV=np.asarray(mean_secondary_ke_profile_eV, dtype=np.float64),
        cumulative_deposited_fraction_profile=np.asarray(
            cumulative_deposited_fraction_profile, dtype=np.float64
        ),
        primary_energy_eV_sum=np.asarray(primary_eV, dtype=np.float64),
        deposited_energy_eV_sum=np.asarray(deposited_eV, dtype=np.float64),
        backscatter_energy_eV_sum=np.asarray(backscatter_eV, dtype=np.float64),
        escaped_forward_energy_eV_sum=np.asarray(forward_eV, dtype=np.float64),
        escaped_lateral_energy_eV_sum=np.asarray(lateral_eV, dtype=np.float64),
        **{field: np.asarray(values, dtype=np.float64) for field, values in split_escape_eV.items()},
        secondary_ke_sumsq_eV2=np.asarray(secondary_ke_sumsq_eV2, dtype=np.float64),
        deposited_fraction=np.asarray(deposited_fraction, dtype=np.float64),
        backscatter_fraction=np.asarray(backscatter_fraction, dtype=np.float64),
        escaped_forward_fraction=np.asarray(forward_fraction, dtype=np.float64),
        escaped_lateral_fraction=np.asarray(lateral_fraction, dtype=np.float64),
        event_count=np.asarray(event_count, dtype=np.int64),
    )


def _load_energy_cache(
    energy_metrics_csv: Path,
    energy_profiles_npz: Path,
) -> tuple[dict[int, dict[str, Any]], dict[int, np.ndarray], dict[int, np.ndarray], np.ndarray]:
    metric_rows = _read_csv_rows(energy_metrics_csv)
    metrics_by_idx: dict[int, dict[str, Any]] = {}
    for row in metric_rows:
        idx = _safe_int(row.get("energy_index", "-1"), default=-1)
        if idx < 0:
            continue
        metrics_by_idx[idx] = {
            "energy_index": idx,
            "status": str(row.get("status", "")),
            "budget_mode": str(row.get("budget_mode", "")),
            "sim_particles": _safe_float(row.get("sim_particles", "nan")),
            "event_count": _safe_float(row.get("event_count", row.get("sim_particles", "nan"))),
            "E_center_MeV": _safe_float(row.get("E_center_MeV", "nan")),
            "primary_energy_eV_sum": _safe_float(row.get("primary_energy_eV_sum", "nan")),
            "deposited_energy_eV_sum": _safe_float(row.get("deposited_energy_eV_sum", "nan")),
            "backscatter_energy_eV_sum": _safe_float(row.get("backscatter_energy_eV_sum", "nan")),
            "escaped_forward_energy_eV_sum": _safe_float(
                row.get("escaped_forward_energy_eV_sum", "nan")
            ),
            "escaped_lateral_energy_eV_sum": _safe_float(
                row.get("escaped_lateral_energy_eV_sum", "nan")
            ),
            **{
                field: _safe_float(row.get(field, "nan"))
                for field in ESCAPE_PARTICLE_ENERGY_FIELDS
            },
            "primary_pen_sum_cm": _safe_float(
                row.get("primary_pen_sum_cm", row.get("primary_pen_sum_mm", "0"))
            ),
            "primary_pen_sumsq_cm2": _safe_float(row.get("primary_pen_sumsq_cm2", "nan")),
            "primary_pen_count": _safe_float(row.get("primary_pen_count", "0")),
            "secondary_pen_sum_cm": _safe_float(
                row.get("secondary_pen_sum_cm", row.get("secondary_pen_sum_mm", "0"))
            ),
            "secondary_pen_sumsq_cm2": _safe_float(row.get("secondary_pen_sumsq_cm2", "nan")),
            "secondary_pen_count": _safe_float(row.get("secondary_pen_count", "0")),
            "secondary_ke_sum_eV": _safe_float(row.get("secondary_ke_sum_eV", "0")),
            "secondary_ke_sumsq_eV2": _safe_float(row.get("secondary_ke_sumsq_eV2", "nan")),
            "secondary_ke_count": _safe_float(row.get("secondary_ke_count", "0")),
        }

    npz = np.load(energy_profiles_npz)
    if "depth_edges_cm" in npz:
        depth_edges_cm = np.asarray(npz["depth_edges_cm"], dtype=np.float64)
    elif "depth_edges_mm" in npz:
        depth_edges_cm = np.asarray(npz["depth_edges_mm"], dtype=np.float64) * 0.1
    else:
        raise KeyError(f"Neither depth_edges_cm nor depth_edges_mm found in {energy_profiles_npz}")
    prof_idx = np.asarray(npz["energy_index"], dtype=np.int64)
    profiles = np.asarray(npz["deposited_profile_eV"], dtype=np.float64)
    profiles_sumsq = (
        np.asarray(npz["deposited_profile_sumsq_eV2"], dtype=np.float64)
        if "deposited_profile_sumsq_eV2" in npz
        else np.full_like(profiles, np.nan, dtype=np.float64)
    )
    profiles_by_idx: dict[int, np.ndarray] = {}
    profiles_sumsq_by_idx: dict[int, np.ndarray] = {}
    for i, eidx in enumerate(prof_idx.tolist()):
        profiles_by_idx[int(eidx)] = profiles[i]
        profiles_sumsq_by_idx[int(eidx)] = profiles_sumsq[i]

    return metrics_by_idx, profiles_by_idx, profiles_sumsq_by_idx, depth_edges_cm


def _compute_range_result(
    range_def: dict[str, Any],
    metrics_by_idx: dict[int, dict[str, Any]],
    profiles_by_idx: dict[int, np.ndarray],
    profile_sumsq_by_idx: dict[int, np.ndarray],
    depth_edges_cm: np.ndarray,
    density_gcm3: float,
    dose_depth_cm: float,
) -> dict[str, Any]:
    n_depth = depth_edges_cm.size - 1
    dose_profile = np.full(n_depth, np.nan, dtype=np.float64)
    dose_profile_std = np.full(n_depth, np.nan, dtype=np.float64)
    dep_profile_flux_eV = np.zeros(n_depth, dtype=np.float64)
    dep_profile_flux_var_eV2 = np.zeros(n_depth, dtype=np.float64)
    dz_cm = depth_edges_cm[1:] - depth_edges_cm[:-1]

    primary_flux_eV = 0.0
    deposited_flux_eV = 0.0
    back_flux_eV = 0.0
    forward_flux_eV = 0.0
    lateral_flux_eV = 0.0
    split_escape_flux_eV = {field: 0.0 for field in ESCAPE_PARTICLE_ENERGY_FIELDS}
    saw_split_escape = {field: False for field in ESCAPE_PARTICLE_ENERGY_FIELDS}
    missing_split_escape = {field: False for field in ESCAPE_PARTICLE_ENERGY_FIELDS}
    jde_flux_cm2_s = 0.0
    jde_energy_flux_model_mev_cm2_s = 0.0
    saw_forward_escape = False
    saw_lateral_escape = False
    missing_forward_escape = False
    missing_lateral_escape = False

    primary_pen_weighted_sum_cm = 0.0
    primary_pen_weight_particles = 0.0
    secondary_pen_weighted_sum_cm = 0.0
    secondary_pen_weight_particles = 0.0
    secondary_ke_sum_scaled = 0.0
    secondary_ke_cnt_scaled = 0.0

    used_energies = 0
    missing_energies = 0

    energy_rows = list(range_def.get("energy_rows", []))

    # Build G4beamline-style dE from neighboring sampled energies.
    # This mirrors:
    #   if first -> E[1]-E[0], if last -> E[-1]-E[-2], else -> E[i+1]-E[i]
    weighted_rows: list[dict[str, Any]] = []
    for erow in energy_rows:
        idx = _safe_int(str(erow.get("energy_index", "-1")), -1)
        if idx < 0:
            continue
        m = metrics_by_idx.get(idx)
        e_center = _safe_float(str(erow.get("E_center_MeV", "nan")))
        if (not np.isfinite(e_center)) and m is not None:
            e_center = _safe_float(str(m.get("E_center_MeV", "nan")))
        weighted_rows.append(
            {
                "idx": int(idx),
                "e_center_mev": float(e_center),
                "row": erow,
                "metrics": m,
            }
        )

    finite_energy_rows = [r for r in weighted_rows if np.isfinite(float(r["e_center_mev"]))]
    finite_energy_rows.sort(key=lambda r: float(r["e_center_mev"]))
    delta_e_by_idx: dict[int, float] = {}
    if len(finite_energy_rows) >= 2:
        e_sorted = np.asarray([float(r["e_center_mev"]) for r in finite_energy_rows], dtype=np.float64)
        for i, r in enumerate(finite_energy_rows):
            if i == 0:
                de = abs(float(e_sorted[1] - e_sorted[0]))
            elif i == len(finite_energy_rows) - 1:
                de = abs(float(e_sorted[-1] - e_sorted[-2]))
            else:
                de = abs(float(e_sorted[i + 1] - e_sorted[i]))
            delta_e_by_idx[int(r["idx"])] = de
    else:
        # Single-energy edge case: prefer overlap bin width if available.
        for r in finite_energy_rows:
            overlap_de = _safe_float(str(r["row"].get("overlap_dE_MeV", "nan")))
            if np.isfinite(overlap_de) and overlap_de > 0.0:
                delta_e_by_idx[int(r["idx"])] = float(overlap_de)

    for erow in energy_rows:
        idx = _safe_int(str(erow.get("energy_index", "-1")), -1)
        if idx < 0:
            missing_energies += 1
            continue
        m = metrics_by_idx.get(idx)
        sim_particles = _safe_float(str(erow.get("sim_particles", "nan")))
        if (not np.isfinite(sim_particles) or sim_particles <= 0.0) and m is not None:
            sim_particles = _safe_float(str(m.get("sim_particles", "nan")))

        # Legacy weight from generator table:
        #   scale_to_sim = (J(E_i) * overlap_dE_i) / N_sim(E_i)
        # G4beamline-style weight uses dE from neighboring sampled energies instead.
        scale_to_sim_legacy = _safe_float(str(erow.get("scale_to_sim", "nan")))
        overlap_de = _safe_float(str(erow.get("overlap_dE_MeV", "nan")))
        flux_model = _safe_float(str(erow.get("flux_model", "nan")))  # J(E_i)
        delta_e_g4bl = _safe_float(str(delta_e_by_idx.get(idx, float("nan"))))

        scale_to_sim = float("nan")
        if (
            np.isfinite(flux_model)
            and flux_model > 0.0
            and np.isfinite(delta_e_g4bl)
            and delta_e_g4bl > 0.0
            and np.isfinite(sim_particles)
            and sim_particles > 0.0
        ):
            scale_to_sim = (flux_model * delta_e_g4bl) / sim_particles
        elif (
            np.isfinite(scale_to_sim_legacy)
            and scale_to_sim_legacy > 0.0
            and np.isfinite(overlap_de)
            and overlap_de > 0.0
            and np.isfinite(delta_e_g4bl)
            and delta_e_g4bl > 0.0
        ):
            scale_to_sim = scale_to_sim_legacy * (delta_e_g4bl / overlap_de)
        else:
            scale_to_sim = scale_to_sim_legacy

        if not np.isfinite(scale_to_sim) or scale_to_sim <= 0.0:
            continue

        if np.isfinite(sim_particles) and sim_particles > 0.0:
            # Project directional intensity onto a surface over a hemisphere (x pi).
            jde_flux_cm2_s += sim_particles * scale_to_sim * math.pi
            e_center_mev = _safe_float(str(erow.get("E_center_MeV", "nan")))
            if (not np.isfinite(e_center_mev)) and m is not None:
                e_center_mev = _safe_float(str(m.get("E_center_MeV", "nan")))
            if np.isfinite(e_center_mev) and e_center_mev > 0.0:
                jde_energy_flux_model_mev_cm2_s += sim_particles * scale_to_sim * e_center_mev * math.pi

        prof = profiles_by_idx.get(idx)
        prof_sumsq = profile_sumsq_by_idx.get(idx)
        if m is None or prof is None:
            missing_energies += 1
            continue
        if not str(m.get("status", "")).startswith("ok"):
            missing_energies += 1
            continue

        used_energies += 1

        primary_e = float(m["primary_energy_eV_sum"])
        deposited_e = float(m["deposited_energy_eV_sum"])
        back_e = float(m["backscatter_energy_eV_sum"])
        forward_e = float(m.get("escaped_forward_energy_eV_sum", float("nan")))
        lateral_e = float(m.get("escaped_lateral_energy_eV_sum", float("nan")))
        split_escape_e = {
            field: float(m.get(field, float("nan")))
            for field in ESCAPE_PARTICLE_ENERGY_FIELDS
        }

        scale_to_surface = scale_to_sim * math.pi

        effective_particles = (
            sim_particles * scale_to_surface
            if np.isfinite(sim_particles) and sim_particles > 0.0
            else float("nan")
        )

        if np.isfinite(primary_e):
            primary_flux_eV += primary_e * scale_to_surface
        if np.isfinite(deposited_e):
            deposited_flux_eV += deposited_e * scale_to_surface
        if np.isfinite(back_e):
            back_flux_eV += back_e * scale_to_surface
        if np.isfinite(forward_e):
            forward_flux_eV += forward_e * scale_to_surface
            saw_forward_escape = True
        else:
            missing_forward_escape = True
        if np.isfinite(lateral_e):
            lateral_flux_eV += lateral_e * scale_to_surface
            saw_lateral_escape = True
        else:
            missing_lateral_escape = True
        for field, value in split_escape_e.items():
            if np.isfinite(value):
                split_escape_flux_eV[field] += value * scale_to_surface
                saw_split_escape[field] = True
            else:
                missing_split_escape[field] = True

        dep_profile_flux_eV += prof * scale_to_surface
        event_count = _safe_float(str(m.get("event_count", "nan")))
        if (
            prof_sumsq is not None
            and np.all(np.isfinite(prof_sumsq))
            and np.isfinite(event_count)
            and event_count > 0.0
        ):
            centered_sumsq = prof_sumsq - (prof * prof) / event_count
            centered_sumsq = np.maximum(centered_sumsq, 0.0)
            dep_profile_flux_var_eV2 += centered_sumsq * (scale_to_surface * scale_to_surface)

        primary_pen_count = float(m["primary_pen_count"])
        if np.isfinite(primary_pen_count) and primary_pen_count > 0.0:
            primary_pen_mean_cm = float(m["primary_pen_sum_cm"]) / primary_pen_count
            if (
                np.isfinite(primary_pen_mean_cm)
                and np.isfinite(effective_particles)
                and effective_particles > 0.0
            ):
                primary_pen_weighted_sum_cm += primary_pen_mean_cm * effective_particles
                primary_pen_weight_particles += effective_particles

        secondary_pen_count = float(m["secondary_pen_count"])
        if np.isfinite(secondary_pen_count) and secondary_pen_count > 0.0:
            secondary_pen_mean_cm = float(m["secondary_pen_sum_cm"]) / secondary_pen_count
            if (
                np.isfinite(secondary_pen_mean_cm)
                and np.isfinite(effective_particles)
                and effective_particles > 0.0
            ):
                secondary_pen_weighted_sum_cm += secondary_pen_mean_cm * effective_particles
                secondary_pen_weight_particles += effective_particles
        secondary_ke_sum_scaled += float(m["secondary_ke_sum_eV"]) * scale_to_surface
        secondary_ke_cnt_scaled += float(m["secondary_ke_count"]) * scale_to_surface

    if used_energies > 0 and density_gcm3 > 0.0:
        dep_profile_flux_mev = dep_profile_flux_eV * 1.0e-6
        with np.errstate(divide="ignore", invalid="ignore"):
            dose_profile = np.divide(
                dep_profile_flux_mev,
                density_gcm3 * dz_cm,
                out=np.zeros_like(dep_profile_flux_mev),
                where=dz_cm > 0,
            )
        dose_profile = dose_profile * MEV_TO_MGY_PER_YEAR
        dep_profile_flux_std_mev = np.sqrt(np.maximum(dep_profile_flux_var_eV2, 0.0)) * 1.0e-6
        with np.errstate(divide="ignore", invalid="ignore"):
            dose_profile_std = np.divide(
                dep_profile_flux_std_mev,
                density_gcm3 * dz_cm,
                out=np.full_like(dep_profile_flux_std_mev, np.nan),
                where=dz_cm > 0,
            )
        dose_profile_std = dose_profile_std * MEV_TO_MGY_PER_YEAR

    depth_centers = 0.5 * (depth_edges_cm[:-1] + depth_edges_cm[1:])
    dose_idx = int(np.argmin(np.abs(depth_centers - dose_depth_cm)))

    mean_primary_pen = (
        primary_pen_weighted_sum_cm / primary_pen_weight_particles
        if primary_pen_weight_particles > 0.0
        else float("nan")
    )
    mean_secondary_pen = (
        secondary_pen_weighted_sum_cm / secondary_pen_weight_particles
        if secondary_pen_weight_particles > 0.0
        else float("nan")
    )
    mean_secondary_ke_mev = (
        (secondary_ke_sum_scaled / secondary_ke_cnt_scaled) * 1.0e-6
        if secondary_ke_cnt_scaled > 0.0
        else float("nan")
    )
    deposited_fraction = (
        deposited_flux_eV / primary_flux_eV if primary_flux_eV > 0.0 else float("nan")
    )
    backscatter_fraction = (
        back_flux_eV / primary_flux_eV if primary_flux_eV > 0.0 else float("nan")
    )
    forward_escape_flux_mev_cm2_s = (
        float(forward_flux_eV * 1.0e-6)
        if saw_forward_escape and not missing_forward_escape
        else float("nan")
    )
    lateral_escape_flux_mev_cm2_s = (
        float(lateral_flux_eV * 1.0e-6)
        if saw_lateral_escape and not missing_lateral_escape
        else float("nan")
    )
    forward_escape_fraction = (
        forward_flux_eV / primary_flux_eV
        if saw_forward_escape and not missing_forward_escape and primary_flux_eV > 0.0
        else float("nan")
    )
    lateral_escape_fraction = (
        lateral_flux_eV / primary_flux_eV
        if saw_lateral_escape and not missing_lateral_escape and primary_flux_eV > 0.0
        else float("nan")
    )
    split_escape_flux_mev_cm2_s = {
        field.replace("_energy_eV_sum", "_flux_mev_cm2_s"): (
            float(value * 1.0e-6)
            if saw_split_escape[field] and not missing_split_escape[field]
            else float("nan")
        )
        for field, value in split_escape_flux_eV.items()
    }
    split_escape_fraction = {
        field.replace("_energy_eV_sum", "_fraction"): (
            value / primary_flux_eV
            if saw_split_escape[field] and not missing_split_escape[field] and primary_flux_eV > 0.0
            else float("nan")
        )
        for field, value in split_escape_flux_eV.items()
    }

    return {
        "range_id": int(range_def["range_id"]),
        "e_min_mev": float(range_def["e_min_mev"]),
        "e_max_mev": float(range_def["e_max_mev"]),
        "n_cells": int(len(range_def.get("cell_ids", []))),
        "used_energies": int(used_energies),
        "missing_energies": int(missing_energies),
        "dose_depth_cm": float(depth_centers[dose_idx]),
        "dose_at_depth_mgy_per_yr": float(dose_profile[dose_idx]) if used_energies > 0 else float("nan"),
        "mean_primary_penetration_cm": float(mean_primary_pen),
        "mean_secondary_ke_mev": float(mean_secondary_ke_mev),
        "mean_secondary_penetration_cm": float(mean_secondary_pen),
        "deposited_fraction": float(deposited_fraction),
        "backscatter_fraction": float(backscatter_fraction),
        "forward_escape_fraction": float(forward_escape_fraction),
        "lateral_escape_fraction": float(lateral_escape_fraction),
        **{k: float(v) for k, v in split_escape_fraction.items()},
        "backscatter_flux_mev_cm2_s": float(back_flux_eV * 1.0e-6),
        "forward_escape_flux_mev_cm2_s": forward_escape_flux_mev_cm2_s,
        "lateral_escape_flux_mev_cm2_s": lateral_escape_flux_mev_cm2_s,
        **{k: float(v) for k, v in split_escape_flux_mev_cm2_s.items()},
        "deposited_flux_mev_cm2_s": float(deposited_flux_eV * 1.0e-6),
        "primary_flux_mev_cm2_s": float(primary_flux_eV * 1.0e-6),
        "jde_flux_cm2_s": float(jde_flux_cm2_s),
        "jde_energy_flux_model_mev_cm2_s": float(jde_energy_flux_model_mev_cm2_s),
        "dose_profile_mgy_per_yr": dose_profile,
        "dose_profile_std_mgy_per_yr": dose_profile_std,
    }


def _range_result_path(range_results_dir: Path, range_id: int) -> Path:
    return range_results_dir / f"range_{range_id:04d}.npz"


def _save_range_result(path: Path, result: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        range_id=np.int64(result["range_id"]),
        e_min_mev=np.float64(result["e_min_mev"]),
        e_max_mev=np.float64(result["e_max_mev"]),
        n_cells=np.int64(result["n_cells"]),
        used_energies=np.int64(result["used_energies"]),
        missing_energies=np.int64(result["missing_energies"]),
        dose_depth_cm=np.float64(result["dose_depth_cm"]),
        dose_at_depth_mgy_per_yr=np.float64(result["dose_at_depth_mgy_per_yr"]),
        mean_primary_penetration_cm=np.float64(result["mean_primary_penetration_cm"]),
        mean_secondary_ke_mev=np.float64(result["mean_secondary_ke_mev"]),
        mean_secondary_penetration_cm=np.float64(result["mean_secondary_penetration_cm"]),
        deposited_fraction=np.float64(result["deposited_fraction"]),
        backscatter_fraction=np.float64(result["backscatter_fraction"]),
        forward_escape_fraction=np.float64(result["forward_escape_fraction"]),
        lateral_escape_fraction=np.float64(result["lateral_escape_fraction"]),
        escaped_electron_fraction=np.float64(result.get("escaped_electron_fraction", float("nan"))),
        escaped_photon_fraction=np.float64(result.get("escaped_photon_fraction", float("nan"))),
        backscatter_electron_fraction=np.float64(result.get("backscatter_electron_fraction", float("nan"))),
        backscatter_photon_fraction=np.float64(result.get("backscatter_photon_fraction", float("nan"))),
        escaped_forward_electron_fraction=np.float64(
            result.get("escaped_forward_electron_fraction", float("nan"))
        ),
        escaped_forward_photon_fraction=np.float64(
            result.get("escaped_forward_photon_fraction", float("nan"))
        ),
        escaped_lateral_electron_fraction=np.float64(
            result.get("escaped_lateral_electron_fraction", float("nan"))
        ),
        escaped_lateral_photon_fraction=np.float64(
            result.get("escaped_lateral_photon_fraction", float("nan"))
        ),
        backscatter_flux_mev_cm2_s=np.float64(result["backscatter_flux_mev_cm2_s"]),
        forward_escape_flux_mev_cm2_s=np.float64(result["forward_escape_flux_mev_cm2_s"]),
        lateral_escape_flux_mev_cm2_s=np.float64(result["lateral_escape_flux_mev_cm2_s"]),
        escaped_electron_flux_mev_cm2_s=np.float64(
            result.get("escaped_electron_flux_mev_cm2_s", float("nan"))
        ),
        escaped_photon_flux_mev_cm2_s=np.float64(
            result.get("escaped_photon_flux_mev_cm2_s", float("nan"))
        ),
        backscatter_electron_flux_mev_cm2_s=np.float64(
            result.get("backscatter_electron_flux_mev_cm2_s", float("nan"))
        ),
        backscatter_photon_flux_mev_cm2_s=np.float64(
            result.get("backscatter_photon_flux_mev_cm2_s", float("nan"))
        ),
        escaped_forward_electron_flux_mev_cm2_s=np.float64(
            result.get("escaped_forward_electron_flux_mev_cm2_s", float("nan"))
        ),
        escaped_forward_photon_flux_mev_cm2_s=np.float64(
            result.get("escaped_forward_photon_flux_mev_cm2_s", float("nan"))
        ),
        escaped_lateral_electron_flux_mev_cm2_s=np.float64(
            result.get("escaped_lateral_electron_flux_mev_cm2_s", float("nan"))
        ),
        escaped_lateral_photon_flux_mev_cm2_s=np.float64(
            result.get("escaped_lateral_photon_flux_mev_cm2_s", float("nan"))
        ),
        deposited_flux_mev_cm2_s=np.float64(result["deposited_flux_mev_cm2_s"]),
        primary_flux_mev_cm2_s=np.float64(result["primary_flux_mev_cm2_s"]),
        jde_flux_cm2_s=np.float64(result["jde_flux_cm2_s"]),
        jde_energy_flux_model_mev_cm2_s=np.float64(result["jde_energy_flux_model_mev_cm2_s"]),
        dose_profile_mgy_per_yr=np.asarray(result["dose_profile_mgy_per_yr"], dtype=np.float64),
        dose_profile_std_mgy_per_yr=np.asarray(result["dose_profile_std_mgy_per_yr"], dtype=np.float64),
    )


def _load_range_result(path: Path) -> dict[str, Any]:
    data = np.load(path)
    dose_depth_cm = (
        float(np.asarray(data["dose_depth_cm"]).item())
        if "dose_depth_cm" in data
        else float(np.asarray(data["dose_depth_mm"]).item()) * 0.1
    )
    primary_pen_cm = (
        float(np.asarray(data["mean_primary_penetration_cm"]).item())
        if "mean_primary_penetration_cm" in data
        else float(np.asarray(data["mean_primary_penetration_mm"]).item()) * 0.1
    )
    secondary_pen_cm = (
        float(np.asarray(data["mean_secondary_penetration_cm"]).item())
        if "mean_secondary_penetration_cm" in data
        else float(np.asarray(data["mean_secondary_penetration_mm"]).item()) * 0.1
    )
    return {
        "range_id": int(np.asarray(data["range_id"]).item()),
        "e_min_mev": float(np.asarray(data["e_min_mev"]).item()),
        "e_max_mev": float(np.asarray(data["e_max_mev"]).item()),
        "n_cells": int(np.asarray(data["n_cells"]).item()),
        "used_energies": int(np.asarray(data["used_energies"]).item()),
        "missing_energies": int(np.asarray(data["missing_energies"]).item()),
        "dose_depth_cm": dose_depth_cm,
        "dose_at_depth_mgy_per_yr": float(np.asarray(data["dose_at_depth_mgy_per_yr"]).item()),
        "mean_primary_penetration_cm": primary_pen_cm,
        "mean_secondary_ke_mev": float(np.asarray(data["mean_secondary_ke_mev"]).item()),
        "mean_secondary_penetration_cm": secondary_pen_cm,
        "deposited_fraction": float(np.asarray(data["deposited_fraction"]).item()),
        "backscatter_fraction": float(np.asarray(data["backscatter_fraction"]).item()),
        "forward_escape_fraction": float(np.asarray(data["forward_escape_fraction"]).item())
        if "forward_escape_fraction" in data
        else float("nan"),
        "lateral_escape_fraction": float(np.asarray(data["lateral_escape_fraction"]).item())
        if "lateral_escape_fraction" in data
        else float("nan"),
        "escaped_electron_fraction": float(np.asarray(data["escaped_electron_fraction"]).item())
        if "escaped_electron_fraction" in data
        else float("nan"),
        "escaped_photon_fraction": float(np.asarray(data["escaped_photon_fraction"]).item())
        if "escaped_photon_fraction" in data
        else float("nan"),
        "backscatter_electron_fraction": float(np.asarray(data["backscatter_electron_fraction"]).item())
        if "backscatter_electron_fraction" in data
        else float("nan"),
        "backscatter_photon_fraction": float(np.asarray(data["backscatter_photon_fraction"]).item())
        if "backscatter_photon_fraction" in data
        else float("nan"),
        "escaped_forward_electron_fraction": float(np.asarray(data["escaped_forward_electron_fraction"]).item())
        if "escaped_forward_electron_fraction" in data
        else float("nan"),
        "escaped_forward_photon_fraction": float(np.asarray(data["escaped_forward_photon_fraction"]).item())
        if "escaped_forward_photon_fraction" in data
        else float("nan"),
        "escaped_lateral_electron_fraction": float(np.asarray(data["escaped_lateral_electron_fraction"]).item())
        if "escaped_lateral_electron_fraction" in data
        else float("nan"),
        "escaped_lateral_photon_fraction": float(np.asarray(data["escaped_lateral_photon_fraction"]).item())
        if "escaped_lateral_photon_fraction" in data
        else float("nan"),
        "backscatter_flux_mev_cm2_s": float(np.asarray(data["backscatter_flux_mev_cm2_s"]).item()),
        "forward_escape_flux_mev_cm2_s": float(np.asarray(data["forward_escape_flux_mev_cm2_s"]).item())
        if "forward_escape_flux_mev_cm2_s" in data
        else float("nan"),
        "lateral_escape_flux_mev_cm2_s": float(np.asarray(data["lateral_escape_flux_mev_cm2_s"]).item())
        if "lateral_escape_flux_mev_cm2_s" in data
        else float("nan"),
        "escaped_electron_flux_mev_cm2_s": float(np.asarray(data["escaped_electron_flux_mev_cm2_s"]).item())
        if "escaped_electron_flux_mev_cm2_s" in data
        else float("nan"),
        "escaped_photon_flux_mev_cm2_s": float(np.asarray(data["escaped_photon_flux_mev_cm2_s"]).item())
        if "escaped_photon_flux_mev_cm2_s" in data
        else float("nan"),
        "backscatter_electron_flux_mev_cm2_s": float(
            np.asarray(data["backscatter_electron_flux_mev_cm2_s"]).item()
        )
        if "backscatter_electron_flux_mev_cm2_s" in data
        else float("nan"),
        "backscatter_photon_flux_mev_cm2_s": float(
            np.asarray(data["backscatter_photon_flux_mev_cm2_s"]).item()
        )
        if "backscatter_photon_flux_mev_cm2_s" in data
        else float("nan"),
        "escaped_forward_electron_flux_mev_cm2_s": float(
            np.asarray(data["escaped_forward_electron_flux_mev_cm2_s"]).item()
        )
        if "escaped_forward_electron_flux_mev_cm2_s" in data
        else float("nan"),
        "escaped_forward_photon_flux_mev_cm2_s": float(
            np.asarray(data["escaped_forward_photon_flux_mev_cm2_s"]).item()
        )
        if "escaped_forward_photon_flux_mev_cm2_s" in data
        else float("nan"),
        "escaped_lateral_electron_flux_mev_cm2_s": float(
            np.asarray(data["escaped_lateral_electron_flux_mev_cm2_s"]).item()
        )
        if "escaped_lateral_electron_flux_mev_cm2_s" in data
        else float("nan"),
        "escaped_lateral_photon_flux_mev_cm2_s": float(
            np.asarray(data["escaped_lateral_photon_flux_mev_cm2_s"]).item()
        )
        if "escaped_lateral_photon_flux_mev_cm2_s" in data
        else float("nan"),
        "deposited_flux_mev_cm2_s": float(np.asarray(data["deposited_flux_mev_cm2_s"]).item()),
        "primary_flux_mev_cm2_s": float(np.asarray(data["primary_flux_mev_cm2_s"]).item()),
        "jde_flux_cm2_s": float(np.asarray(data["jde_flux_cm2_s"]).item())
        if "jde_flux_cm2_s" in data
        else float("nan"),
        "jde_energy_flux_model_mev_cm2_s": float(np.asarray(data["jde_energy_flux_model_mev_cm2_s"]).item())
        if "jde_energy_flux_model_mev_cm2_s" in data
        else float("nan"),
        "dose_profile_mgy_per_yr": np.asarray(data["dose_profile_mgy_per_yr"], dtype=np.float64),
        "dose_profile_std_mgy_per_yr": (
            np.asarray(data["dose_profile_std_mgy_per_yr"], dtype=np.float64)
            if "dose_profile_std_mgy_per_yr" in data
            else np.full_like(np.asarray(data["dose_profile_mgy_per_yr"], dtype=np.float64), np.nan)
        ),
    }


def cmd_prepare(args: argparse.Namespace) -> None:
    library_dir = (args.library_dir or _default_library_dir()).resolve()
    out_dir = (args.out_dir or _default_out_dir(library_dir)).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    range_csv = library_dir / "latlon_cell_ranges.csv"
    scaling_csv = library_dir / "latlon_energy_scaling.csv.gz"
    run_table_csv = library_dir / "per_energy_runs.csv"

    if not range_csv.exists():
        raise FileNotFoundError(f"Missing {range_csv}")
    if not scaling_csv.exists():
        raise FileNotFoundError(f"Missing {scaling_csv}")
    if not run_table_csv.exists():
        raise FileNotFoundError(f"Missing {run_table_csv}")

    if args.z_range_file is not None:
        depth_edges_cm = _depth_edges_from_file(args.z_range_file)
        if str(args.z_range_unit).lower() == "mm":
            depth_edges_cm = depth_edges_cm * 0.1
    else:
        depth_edges_cm = _legacy_depth_edges_cm()
    depth_edges_file = out_dir / "depth_edges_cm.npy"
    np.save(depth_edges_file, depth_edges_cm)

    run_rows = _load_run_rows(run_table_csv)
    default_density = float(args.density_gcm3)
    if math.isnan(default_density):
        if run_rows:
            default_density = float(run_rows[0]["density_gcm3"])
        if not np.isfinite(default_density):
            default_density = 0.5

    raw_cells = _read_csv_rows(range_csv)
    raw_scaling = _read_csv_gz_rows(scaling_csv)

    scaling_by_cell: dict[int, list[dict[str, Any]]] = {}
    for row in raw_scaling:
        cid = _safe_int(row.get("cell_id", "-1"), -1)
        if cid < 0:
            continue
        scaling_by_cell.setdefault(cid, []).append(
            {
                "energy_index": _safe_int(row.get("energy_index", "-1"), -1),
                "E_center_MeV": _safe_float(row.get("E_center_MeV", "nan")),
                "overlap_dE_MeV": _safe_float(row.get("overlap_dE_MeV", "nan")),
                "flux_model": _safe_float(row.get("flux_model", "nan")),
                "scale_to_sim": _safe_float(row.get("scale_to_sim", "nan")),
                "sim_particles": _safe_float(row.get("sim_particles", "nan")),
            }
        )
    for cid in list(scaling_by_cell.keys()):
        scaling_by_cell[cid].sort(key=lambda r: int(r["energy_index"]))

    cells: list[dict[str, Any]] = []
    range_key_to_cells: dict[tuple[float, float], list[int]] = {}
    for row in raw_cells:
        cell_id = _safe_int(row.get("cell_id", "-1"), -1)
        lat = _safe_float(row.get("lat_deg", "nan"))
        lon = _safe_float(row.get("lon_w_deg", "nan"))
        hemisphere = str(row.get("hemisphere", ""))
        e_min_s = str(row.get("e_min_mev", "nan"))
        e_max_s = str(row.get("e_max_mev", "nan"))
        e_min = _safe_float(e_min_s)
        e_max = _safe_float(e_max_s)
        is_valid = _safe_int(row.get("is_valid", "0"), 0)
        valid = bool(is_valid == 1 and np.isfinite(e_min) and np.isfinite(e_max) and e_max > e_min)

        cell_obj = {
            "cell_id": cell_id,
            "lat_deg": lat,
            "lon_w_deg": lon,
            "hemisphere": hemisphere,
            "e_min_mev": e_min,
            "e_max_mev": e_max,
            "is_valid": int(valid),
            "range_id": -1,
        }
        cells.append(cell_obj)
        if valid:
            key = (round(e_min, 12), round(e_max, 12))
            range_key_to_cells.setdefault(key, []).append(cell_id)

    sorted_keys = sorted(range_key_to_cells.keys(), key=lambda t: (float(t[0]), float(t[1])))
    unique_ranges: list[dict[str, Any]] = []
    range_key_to_id: dict[tuple[float, float], int] = {}
    mismatch_count = 0
    for rid, key in enumerate(sorted_keys):
        range_key_to_id[key] = rid
        cell_ids = sorted(range_key_to_cells[key])
        representative = cell_ids[0] if cell_ids else -1
        ref_rows = scaling_by_cell.get(representative, [])
        ref_sig = [(int(r["energy_index"]), float(r["scale_to_sim"])) for r in ref_rows]
        for cid in cell_ids[1:]:
            cur_rows = scaling_by_cell.get(cid, [])
            cur_sig = [(int(r["energy_index"]), float(r["scale_to_sim"])) for r in cur_rows]
            if len(cur_sig) != len(ref_sig):
                mismatch_count += 1
                continue
            for (ei_ref, s_ref), (ei_cur, s_cur) in zip(ref_sig, cur_sig):
                if ei_ref != ei_cur or not np.isclose(s_ref, s_cur, rtol=1e-10, atol=1e-12):
                    mismatch_count += 1
                    break

        unique_ranges.append(
            {
                "range_id": rid,
                "e_min_mev": float(key[0]),
                "e_max_mev": float(key[1]),
                "cell_ids": cell_ids,
                "energy_rows": ref_rows,
            }
        )

    for cell in cells:
        if int(cell["is_valid"]) != 1:
            continue
        key = (round(float(cell["e_min_mev"]), 12), round(float(cell["e_max_mev"]), 12))
        if key in range_key_to_id:
            cell["range_id"] = int(range_key_to_id[key])
            continue
        matched = -1
        for raw_key, rid in range_key_to_id.items():
            if (
                np.isclose(float(raw_key[0]), float(cell["e_min_mev"]), rtol=0.0, atol=1e-12)
                and np.isclose(float(raw_key[1]), float(cell["e_max_mev"]), rtol=0.0, atol=1e-12)
            ):
                matched = int(rid)
                break
        cell["range_id"] = matched

    # Write helper tables for inspection.
    unique_ranges_csv = out_dir / "unique_ranges.csv"
    with unique_ranges_csv.open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["range_id", "e_min_mev", "e_max_mev", "n_cells", "n_energies", "cell_ids"])
        for r in unique_ranges:
            wr.writerow(
                [
                    int(r["range_id"]),
                    f"{float(r['e_min_mev']):.9g}",
                    f"{float(r['e_max_mev']):.9g}",
                    len(r["cell_ids"]),
                    len(r["energy_rows"]),
                    "|".join(str(c) for c in r["cell_ids"]),
                ]
            )

    cell_to_range_csv = out_dir / "cell_to_range.csv"
    with cell_to_range_csv.open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(
            [
                "cell_id",
                "lat_deg",
                "lon_w_deg",
                "hemisphere",
                "is_valid",
                "range_id",
                "e_min_mev",
                "e_max_mev",
            ]
        )
        for c in sorted(cells, key=lambda x: int(x["cell_id"])):
            wr.writerow(
                [
                    int(c["cell_id"]),
                    f"{float(c['lat_deg']):.9g}",
                    f"{float(c['lon_w_deg']):.9g}",
                    c["hemisphere"],
                    int(c["is_valid"]),
                    int(c["range_id"]),
                    f"{float(c['e_min_mev']):.9g}" if np.isfinite(float(c["e_min_mev"])) else "nan",
                    f"{float(c['e_max_mev']):.9g}" if np.isfinite(float(c["e_max_mev"])) else "nan",
                ]
            )

    manifest = {
        "version": 2,
        "library_dir": str(library_dir),
        "out_dir": str(out_dir),
        "range_csv": str(range_csv),
        "scaling_csv_gz": str(scaling_csv),
        "run_table_csv": str(run_table_csv),
        "energy_metrics_csv": str(out_dir / "energy_metrics.csv"),
        "energy_profiles_npz": str(out_dir / "energy_edep_profiles.npz"),
        "energy_feedback_npz": str(out_dir / "energy_feedback_summary.npz"),
        "energy_feedback_csv": str(out_dir / "energy_feedback_summary.csv"),
        "energy_observables_npz": str(out_dir / "energy_observables.npz"),
        "energy_results_dir": str(out_dir / "energy_results"),
        "range_results_dir": str(out_dir / "range_results"),
        "depth_edges_file": str(depth_edges_file),
        "dose_depth_cm_default": float(args.dose_depth_cm),
        "density_gcm3_default": float(default_density),
        "cells": sorted(cells, key=lambda x: int(x["cell_id"])),
        "unique_ranges": unique_ranges,
    }
    manifest_path = out_dir / "analysis_manifest.json"
    _write_manifest(manifest_path, manifest)

    print(f"Prepared manifest: {manifest_path}")
    print(f"Cells: {len(cells)}")
    print(f"Unique valid ranges: {len(unique_ranges)}")
    print(f"Depth bins: {depth_edges_cm.size - 1}")
    print(f"Default density: {default_density:.9g} g/cm3")
    print(f"Default dose depth: {float(args.dose_depth_cm):.9g} cm")
    if mismatch_count > 0:
        print(f"WARNING: found {mismatch_count} same-range cells with scaling-row mismatches.")


def _build_read_error_metrics(run_row: dict[str, Any], n_root_files: int, error: str) -> dict[str, Any]:
    return {
        "energy_index": int(run_row["energy_index"]),
        "E_center_MeV": float(run_row["E_center_MeV"]),
        "sim_particles": int(run_row["sim_particles"]),
        "n_root_files": int(n_root_files),
        "status": "read_error",
        "budget_mode": "none",
        "event_count": 0,
        "primary_energy_eV_sum": float("nan"),
        "deposited_energy_eV_sum": float("nan"),
        "backscatter_energy_eV_sum": float("nan"),
        "escaped_forward_energy_eV_sum": float("nan"),
        "escaped_lateral_energy_eV_sum": float("nan"),
        **{field: float("nan") for field in ESCAPE_PARTICLE_ENERGY_FIELDS},
        "primary_pen_sum_cm": 0.0,
        "primary_pen_sumsq_cm2": 0.0,
        "primary_pen_count": 0,
        "secondary_pen_sum_cm": 0.0,
        "secondary_pen_sumsq_cm2": 0.0,
        "secondary_pen_count": 0,
        "secondary_ke_sum_eV": 0.0,
        "secondary_ke_sumsq_eV2": float("nan"),
        "secondary_ke_count": 0,
        "error": str(error),
    }


def _analyze_cache_run_row(
    run_row: dict[str, Any],
    library_dir_str: str,
    root_dir_override_str: str | None,
    work_root_override_str: str | None,
    depth_edges_cm: np.ndarray,
    step_size: str,
    depth_origin_cm: float,
) -> tuple[
    dict[str, Any],
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    library_dir = Path(library_dir_str)
    root_dir_override = Path(root_dir_override_str) if root_dir_override_str else None
    work_root_override = Path(work_root_override_str) if work_root_override_str else None
    root_files = _resolve_root_files_for_row(
        run_row,
        library_dir=library_dir,
        root_dir_override=root_dir_override,
        work_root_override=work_root_override,
    )
    try:
        return _analyze_single_energy(
            run_row,
            root_files,
            depth_edges_cm=depth_edges_cm,
            step_size=step_size,
            depth_origin_cm=depth_origin_cm,
        )
    except Exception as exc:  # pragma: no cover - defensive
        return (
            _build_read_error_metrics(run_row, len(root_files), str(exc)),
            np.zeros(depth_edges_cm.size - 1, dtype=np.float64),
            np.zeros(depth_edges_cm.size - 1, dtype=np.float64),
            np.zeros(depth_edges_cm.size - 1, dtype=np.float64),
            np.zeros(depth_edges_cm.size - 1, dtype=np.float64),
            np.zeros(depth_edges_cm.size - 1, dtype=np.float64),
            np.zeros(depth_edges_cm.size - 1, dtype=np.float64),
        )


def cmd_build_energy_cache(args: argparse.Namespace) -> None:
    _require_uproot()

    manifest = _load_manifest(args.manifest)
    library_dir = Path(manifest["library_dir"]).resolve()
    out_dir = Path(manifest["out_dir"]).resolve()
    run_table_csv = Path(manifest["run_table_csv"]).resolve()
    depth_edges_cm = np.load(Path(manifest["depth_edges_file"]).resolve())

    energy_metrics_csv = Path(manifest["energy_metrics_csv"]).resolve()
    energy_profiles_npz = Path(manifest["energy_profiles_npz"]).resolve()
    energy_feedback_npz = Path(manifest.get("energy_feedback_npz", out_dir / "energy_feedback_summary.npz")).resolve()
    energy_feedback_csv = Path(manifest.get("energy_feedback_csv", out_dir / "energy_feedback_summary.csv")).resolve()
    energy_observables_npz = Path(
        manifest.get("energy_observables_npz", out_dir / "energy_observables.npz")
    ).resolve()

    if (
        energy_metrics_csv.exists()
        and energy_profiles_npz.exists()
        and energy_feedback_npz.exists()
        and energy_feedback_csv.exists()
        and energy_observables_npz.exists()
        and not args.force
    ):
        print("Energy cache already exists. Use --force to recompute.")
        print(f"  {energy_metrics_csv}")
        print(f"  {energy_profiles_npz}")
        print(f"  {energy_feedback_npz}")
        print(f"  {energy_feedback_csv}")
        print(f"  {energy_observables_npz}")
        return

    run_rows = _load_run_rows(run_table_csv)
    n_energy = len(run_rows)
    if n_energy == 0:
        raise RuntimeError(f"No rows in run table: {run_table_csv}")

    root_dir_override = args.root_dir.resolve() if args.root_dir is not None else None
    work_root_override = args.work_root_dir.resolve() if args.work_root_dir is not None else None
    if root_dir_override is None and work_root_override is None:
        work_root_override = _default_root_dir(library_dir)

    deposited_profiles = np.zeros((n_energy, depth_edges_cm.size - 1), dtype=np.float64)
    deposited_profile_sumsq = np.zeros((n_energy, depth_edges_cm.size - 1), dtype=np.float64)
    primary_ke_sum_profiles = np.zeros((n_energy, depth_edges_cm.size - 1), dtype=np.float64)
    primary_ke_count_profiles = np.zeros((n_energy, depth_edges_cm.size - 1), dtype=np.float64)
    secondary_ke_sum_profiles = np.zeros((n_energy, depth_edges_cm.size - 1), dtype=np.float64)
    secondary_ke_count_profiles = np.zeros((n_energy, depth_edges_cm.size - 1), dtype=np.float64)
    metrics_rows: list[dict[str, Any]] = [{} for _ in range(n_energy)]

    n_workers = max(1, min(int(args.workers), n_energy))
    depth_origin_cm = float(args.depth_origin_cm)
    library_dir_str = str(library_dir)
    root_dir_override_str = str(root_dir_override) if root_dir_override is not None else None
    work_root_override_str = str(work_root_override) if work_root_override is not None else None

    print(f"Energy-cache process workers: {n_workers}")

    if n_workers == 1:
        for i, run_row in enumerate(run_rows):
            (
                metrics,
                deposited_profile,
                deposited_profile_sumsq_eV2,
                primary_ke_sum_profile,
                primary_ke_count_profile,
                secondary_ke_sum_profile,
                secondary_ke_count_profile,
            ) = _analyze_cache_run_row(
                run_row=run_row,
                library_dir_str=library_dir_str,
                root_dir_override_str=root_dir_override_str,
                work_root_override_str=work_root_override_str,
                depth_edges_cm=depth_edges_cm,
                step_size=args.step_size,
                depth_origin_cm=depth_origin_cm,
            )
            deposited_profiles[i, :] = deposited_profile
            deposited_profile_sumsq[i, :] = deposited_profile_sumsq_eV2
            primary_ke_sum_profiles[i, :] = primary_ke_sum_profile
            primary_ke_count_profiles[i, :] = primary_ke_count_profile
            secondary_ke_sum_profiles[i, :] = secondary_ke_sum_profile
            secondary_ke_count_profiles[i, :] = secondary_ke_count_profile
            metrics_rows[i] = metrics

            if (i + 1) % 20 == 0 or i == 0 or (i + 1) == n_energy:
                print(
                    f"[{i + 1:4d}/{n_energy}] "
                    f"E{int(run_row['energy_index']):05d} "
                    f"status={metrics['status']} "
                    f"roots={metrics['n_root_files']}"
                )
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        import multiprocessing as mp

        completed = 0
        try:
            mp_ctx = mp.get_context("fork")
        except Exception:
            mp_ctx = None
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=mp_ctx) as ex:
            future_to_idx = {
                ex.submit(
                    _analyze_cache_run_row,
                    run_row,
                    library_dir_str,
                    root_dir_override_str,
                    work_root_override_str,
                    depth_edges_cm,
                    args.step_size,
                    depth_origin_cm,
                ): i
                for i, run_row in enumerate(run_rows)
            }
            for fut in as_completed(future_to_idx):
                i = int(future_to_idx[fut])
                run_row = run_rows[i]
                try:
                    (
                        metrics,
                        deposited_profile,
                        deposited_profile_sumsq_eV2,
                        primary_ke_sum_profile,
                        primary_ke_count_profile,
                        secondary_ke_sum_profile,
                        secondary_ke_count_profile,
                    ) = fut.result()
                except Exception as exc:  # pragma: no cover - defensive
                    metrics = _build_read_error_metrics(run_row, 0, f"worker_failed: {exc}")
                    deposited_profile = np.zeros(depth_edges_cm.size - 1, dtype=np.float64)
                    deposited_profile_sumsq_eV2 = np.full(depth_edges_cm.size - 1, np.nan, dtype=np.float64)
                    primary_ke_sum_profile = np.zeros(depth_edges_cm.size - 1, dtype=np.float64)
                    primary_ke_count_profile = np.zeros(depth_edges_cm.size - 1, dtype=np.float64)
                    secondary_ke_sum_profile = np.zeros(depth_edges_cm.size - 1, dtype=np.float64)
                    secondary_ke_count_profile = np.zeros(depth_edges_cm.size - 1, dtype=np.float64)
                deposited_profiles[i, :] = deposited_profile
                deposited_profile_sumsq[i, :] = deposited_profile_sumsq_eV2
                primary_ke_sum_profiles[i, :] = primary_ke_sum_profile
                primary_ke_count_profiles[i, :] = primary_ke_count_profile
                secondary_ke_sum_profiles[i, :] = secondary_ke_sum_profile
                secondary_ke_count_profiles[i, :] = secondary_ke_count_profile
                metrics_rows[i] = metrics
                completed += 1
                if completed % 20 == 0 or completed == 1 or completed == n_energy:
                    print(
                        f"[{completed:4d}/{n_energy}] "
                        f"E{int(run_row['energy_index']):05d} "
                        f"status={metrics['status']} "
                        f"roots={metrics['n_root_files']}"
                    )

    metrics_rows.sort(key=lambda r: int(r["energy_index"]))
    order = np.argsort(np.asarray([int(r["energy_index"]) for r in metrics_rows], dtype=np.int64))
    deposited_profiles_sorted = deposited_profiles[order, :]
    deposited_profile_sumsq_sorted = deposited_profile_sumsq[order, :]
    primary_ke_sum_profiles_sorted = primary_ke_sum_profiles[order, :]
    primary_ke_count_profiles_sorted = primary_ke_count_profiles[order, :]
    secondary_ke_sum_profiles_sorted = secondary_ke_sum_profiles[order, :]
    secondary_ke_count_profiles_sorted = secondary_ke_count_profiles[order, :]
    metrics_sorted = [metrics_rows[int(i)] for i in order.tolist()]

    _save_energy_cache(
        metrics_rows=metrics_sorted,
        deposited_profiles_eV=deposited_profiles_sorted,
        deposited_profile_sumsq_eV2=deposited_profile_sumsq_sorted,
        primary_ke_sum_profiles_eV=primary_ke_sum_profiles_sorted,
        primary_ke_count_profiles=primary_ke_count_profiles_sorted,
        secondary_ke_sum_profiles_eV=secondary_ke_sum_profiles_sorted,
        secondary_ke_count_profiles=secondary_ke_count_profiles_sorted,
        depth_edges_cm=depth_edges_cm,
        energy_metrics_csv=energy_metrics_csv,
        energy_profiles_npz=energy_profiles_npz,
        energy_feedback_npz=energy_feedback_npz,
        energy_feedback_csv=energy_feedback_csv,
        energy_observables_npz=energy_observables_npz,
    )
    print(f"Wrote energy metrics CSV: {energy_metrics_csv}")
    print(f"Wrote energy profile NPZ: {energy_profiles_npz}")
    print(f"Wrote per-energy feedback CSV: {energy_feedback_csv}")
    print(f"Wrote per-energy feedback NPZ: {energy_feedback_npz}")
    print(f"Wrote per-energy observables NPZ: {energy_observables_npz}")


def _energy_result_path(energy_results_dir: Path, energy_index: int) -> Path:
    return energy_results_dir / f"energy_{energy_index:05d}.npz"


def _save_energy_result(
    path: Path,
    metrics: dict[str, Any],
    depth_edges_cm: np.ndarray,
    deposited_profile_eV: np.ndarray,
    deposited_profile_sumsq_eV2: np.ndarray,
    primary_ke_sum_profile_eV: np.ndarray,
    primary_ke_count_profile: np.ndarray,
    secondary_ke_sum_profile_eV: np.ndarray,
    secondary_ke_count_profile: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        metrics_json=np.asarray(json.dumps(metrics), dtype="<U8192"),
        depth_edges_cm=np.asarray(depth_edges_cm, dtype=np.float64),
        deposited_profile_eV=np.asarray(deposited_profile_eV, dtype=np.float64),
        deposited_profile_sumsq_eV2=np.asarray(deposited_profile_sumsq_eV2, dtype=np.float64),
        primary_ke_sum_profile_eV=np.asarray(primary_ke_sum_profile_eV, dtype=np.float64),
        primary_ke_count_profile=np.asarray(primary_ke_count_profile, dtype=np.float64),
        secondary_ke_sum_profile_eV=np.asarray(secondary_ke_sum_profile_eV, dtype=np.float64),
        secondary_ke_count_profile=np.asarray(secondary_ke_count_profile, dtype=np.float64),
    )


def _load_energy_result(
    path: Path,
) -> tuple[dict[str, Any], np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    data = np.load(path, allow_pickle=False)
    metrics = json.loads(str(np.asarray(data["metrics_json"]).item()))
    if "depth_edges_cm" in data:
        depth_edges_cm = np.asarray(data["depth_edges_cm"], dtype=np.float64)
    elif "depth_edges_mm" in data:
        depth_edges_cm = np.asarray(data["depth_edges_mm"], dtype=np.float64) * 0.1
    else:
        raise KeyError(f"Missing depth edges in {path}")
    return (
        metrics,
        depth_edges_cm,
        np.asarray(data["deposited_profile_eV"], dtype=np.float64),
        (
            np.asarray(data["deposited_profile_sumsq_eV2"], dtype=np.float64)
            if "deposited_profile_sumsq_eV2" in data
            else np.full_like(np.asarray(data["deposited_profile_eV"], dtype=np.float64), np.nan)
        ),
        np.asarray(data["primary_ke_sum_profile_eV"], dtype=np.float64),
        np.asarray(data["primary_ke_count_profile"], dtype=np.float64),
        np.asarray(data["secondary_ke_sum_profile_eV"], dtype=np.float64),
        np.asarray(data["secondary_ke_count_profile"], dtype=np.float64),
    )


def cmd_energy_count(args: argparse.Namespace) -> None:
    manifest = _load_manifest(args.manifest)
    run_rows = _load_run_rows(Path(manifest["run_table_csv"]).resolve())
    print(len(run_rows))


def cmd_energy_worker(args: argparse.Namespace) -> None:
    _require_uproot()

    manifest = _load_manifest(args.manifest)
    library_dir = Path(manifest["library_dir"]).resolve()
    out_dir = Path(manifest["out_dir"]).resolve()
    run_rows = _load_run_rows(Path(manifest["run_table_csv"]).resolve())
    depth_edges_cm = np.load(Path(manifest["depth_edges_file"]).resolve())
    energy_results_dir = Path(manifest.get("energy_results_dir", out_dir / "energy_results")).resolve()

    pos = int(args.energy_pos)
    if pos < 0 or pos >= len(run_rows):
        raise IndexError(f"energy_pos {pos} out of bounds [0, {len(run_rows)-1}]")
    run_row = run_rows[pos]
    energy_index = int(run_row["energy_index"])
    out_file = _energy_result_path(energy_results_dir, energy_index)
    if out_file.exists() and not args.force:
        print(f"Energy result exists, skipping E{energy_index:05d}: {out_file}")
        return

    root_dir_override = args.root_dir.resolve() if args.root_dir is not None else None
    work_root_override = args.work_root_dir.resolve() if args.work_root_dir is not None else None
    if root_dir_override is None and work_root_override is None:
        work_root_override = _default_root_dir(library_dir)

    (
        metrics,
        deposited_profile,
        deposited_profile_sumsq_eV2,
        primary_ke_sum_profile,
        primary_ke_count_profile,
        secondary_ke_sum_profile,
        secondary_ke_count_profile,
    ) = _analyze_cache_run_row(
        run_row=run_row,
        library_dir_str=str(library_dir),
        root_dir_override_str=str(root_dir_override) if root_dir_override is not None else None,
        work_root_override_str=str(work_root_override) if work_root_override is not None else None,
        depth_edges_cm=np.asarray(depth_edges_cm, dtype=np.float64),
        step_size=str(args.step_size),
        depth_origin_cm=float(args.depth_origin_cm),
    )

    _save_energy_result(
        path=out_file,
        metrics=metrics,
        depth_edges_cm=np.asarray(depth_edges_cm, dtype=np.float64),
        deposited_profile_eV=deposited_profile,
        deposited_profile_sumsq_eV2=deposited_profile_sumsq_eV2,
        primary_ke_sum_profile_eV=primary_ke_sum_profile,
        primary_ke_count_profile=primary_ke_count_profile,
        secondary_ke_sum_profile_eV=secondary_ke_sum_profile,
        secondary_ke_count_profile=secondary_ke_count_profile,
    )
    print(
        f"E{energy_index:05d} status={metrics.get('status','')} roots={metrics.get('n_root_files',0)} -> {out_file}"
    )


def cmd_energy_merge(args: argparse.Namespace) -> None:
    manifest = _load_manifest(args.manifest)
    out_dir = Path(manifest["out_dir"]).resolve()
    run_rows = _load_run_rows(Path(manifest["run_table_csv"]).resolve())
    energy_results_dir = Path(manifest.get("energy_results_dir", out_dir / "energy_results")).resolve()

    energy_metrics_csv = Path(manifest["energy_metrics_csv"]).resolve()
    energy_profiles_npz = Path(manifest["energy_profiles_npz"]).resolve()
    energy_feedback_npz = Path(manifest.get("energy_feedback_npz", out_dir / "energy_feedback_summary.npz")).resolve()
    energy_feedback_csv = Path(manifest.get("energy_feedback_csv", out_dir / "energy_feedback_summary.csv")).resolve()
    energy_observables_npz = Path(
        manifest.get("energy_observables_npz", out_dir / "energy_observables.npz")
    ).resolve()

    if (
        energy_metrics_csv.exists()
        and energy_profiles_npz.exists()
        and energy_feedback_npz.exists()
        and energy_feedback_csv.exists()
        and energy_observables_npz.exists()
        and not args.force
    ):
        print("Merged energy cache already exists. Use --force to recompute.")
        print(f"  {energy_metrics_csv}")
        print(f"  {energy_profiles_npz}")
        print(f"  {energy_feedback_npz}")
        print(f"  {energy_feedback_csv}")
        print(f"  {energy_observables_npz}")
        return

    metrics_rows: list[dict[str, Any]] = []
    deposited_profiles: list[np.ndarray] = []
    deposited_profile_sumsq: list[np.ndarray] = []
    primary_ke_sum_profiles: list[np.ndarray] = []
    primary_ke_count_profiles: list[np.ndarray] = []
    secondary_ke_sum_profiles: list[np.ndarray] = []
    secondary_ke_count_profiles: list[np.ndarray] = []
    depth_edges_ref: np.ndarray | None = None
    missing: list[int] = []

    for run_row in run_rows:
        energy_index = int(run_row["energy_index"])
        path = _energy_result_path(energy_results_dir, energy_index)
        if not path.exists():
            missing.append(energy_index)
            continue
        (
            metrics,
            depth_edges_cm,
            deposited_profile,
            deposited_profile_sumsq_eV2,
            primary_ke_sum_profile,
            primary_ke_count_profile,
            secondary_ke_sum_profile,
            secondary_ke_count_profile,
        ) = _load_energy_result(path)
        if depth_edges_ref is None:
            depth_edges_ref = depth_edges_cm
        elif not np.allclose(depth_edges_ref, depth_edges_cm, rtol=0.0, atol=1e-12):
            raise ValueError(f"Depth edges mismatch in {path}")

        metrics_rows.append(metrics)
        deposited_profiles.append(deposited_profile)
        deposited_profile_sumsq.append(deposited_profile_sumsq_eV2)
        primary_ke_sum_profiles.append(primary_ke_sum_profile)
        primary_ke_count_profiles.append(primary_ke_count_profile)
        secondary_ke_sum_profiles.append(secondary_ke_sum_profile)
        secondary_ke_count_profiles.append(secondary_ke_count_profile)

    if missing:
        preview = ",".join(f"{x:05d}" for x in missing[:15])
        more = "" if len(missing) <= 15 else f" ... (+{len(missing)-15} more)"
        raise FileNotFoundError(
            f"Missing {len(missing)} energy-result NPZ files in {energy_results_dir}: {preview}{more}"
        )
    if depth_edges_ref is None:
        raise RuntimeError("No energy result files found to merge.")

    metrics_rows.sort(key=lambda r: int(r["energy_index"]))
    order = np.argsort(np.asarray([int(r["energy_index"]) for r in metrics_rows], dtype=np.int64))
    dep_arr = np.asarray(deposited_profiles, dtype=np.float64)[order, :]
    dep_sumsq_arr = np.asarray(deposited_profile_sumsq, dtype=np.float64)[order, :]
    pk_sum_arr = np.asarray(primary_ke_sum_profiles, dtype=np.float64)[order, :]
    pk_cnt_arr = np.asarray(primary_ke_count_profiles, dtype=np.float64)[order, :]
    sk_sum_arr = np.asarray(secondary_ke_sum_profiles, dtype=np.float64)[order, :]
    sk_cnt_arr = np.asarray(secondary_ke_count_profiles, dtype=np.float64)[order, :]
    metrics_sorted = [metrics_rows[int(i)] for i in order.tolist()]

    _save_energy_cache(
        metrics_rows=metrics_sorted,
        deposited_profiles_eV=dep_arr,
        deposited_profile_sumsq_eV2=dep_sumsq_arr,
        primary_ke_sum_profiles_eV=pk_sum_arr,
        primary_ke_count_profiles=pk_cnt_arr,
        secondary_ke_sum_profiles_eV=sk_sum_arr,
        secondary_ke_count_profiles=sk_cnt_arr,
        depth_edges_cm=np.asarray(depth_edges_ref, dtype=np.float64),
        energy_metrics_csv=energy_metrics_csv,
        energy_profiles_npz=energy_profiles_npz,
        energy_feedback_npz=energy_feedback_npz,
        energy_feedback_csv=energy_feedback_csv,
        energy_observables_npz=energy_observables_npz,
    )
    print(f"Merged {len(metrics_sorted)} energy NPZ results from: {energy_results_dir}")
    print(f"Wrote energy metrics CSV: {energy_metrics_csv}")
    print(f"Wrote energy profile NPZ: {energy_profiles_npz}")
    print(f"Wrote per-energy feedback CSV: {energy_feedback_csv}")
    print(f"Wrote per-energy feedback NPZ: {energy_feedback_npz}")
    print(f"Wrote per-energy observables NPZ: {energy_observables_npz}")


def cmd_range_count(args: argparse.Namespace) -> None:
    manifest = _load_manifest(args.manifest)
    n_ranges = len(manifest.get("unique_ranges", []))
    print(n_ranges)


def cmd_range_worker(args: argparse.Namespace) -> None:
    manifest = _load_manifest(args.manifest)
    ranges = manifest.get("unique_ranges", [])
    if not isinstance(ranges, list) or len(ranges) == 0:
        raise RuntimeError("Manifest has no unique_ranges.")

    idx = int(args.range_index)
    if idx < 0 or idx >= len(ranges):
        raise IndexError(f"range_index {idx} out of bounds [0, {len(ranges)-1}]")

    energy_metrics_csv = Path(manifest["energy_metrics_csv"]).resolve()
    energy_profiles_npz = Path(manifest["energy_profiles_npz"]).resolve()
    if not energy_metrics_csv.exists() or not energy_profiles_npz.exists():
        raise FileNotFoundError(
            "Energy cache missing. Run `energy-merge` (recommended) or `build-energy-cache` first.\n"
            f"Expected:\n  {energy_metrics_csv}\n  {energy_profiles_npz}"
        )

    metrics_by_idx, profiles_by_idx, profile_sumsq_by_idx, depth_edges_cm = _load_energy_cache(
        energy_metrics_csv=energy_metrics_csv,
        energy_profiles_npz=energy_profiles_npz,
    )

    density = float(args.density_gcm3)
    if not np.isfinite(density):
        density = _safe_float(str(manifest.get("density_gcm3_default", "nan")))
    if not np.isfinite(density) or density <= 0.0:
        raise ValueError(f"Invalid density_gcm3={density}")

    dose_depth_cm = float(args.dose_depth_cm)
    if not np.isfinite(dose_depth_cm):
        dose_depth_cm = _safe_float(str(manifest.get("dose_depth_cm_default", "nan")))
    if not np.isfinite(dose_depth_cm):
        dose_depth_cm = 0.01

    range_def = ranges[idx]
    result = _compute_range_result(
        range_def=range_def,
        metrics_by_idx=metrics_by_idx,
        profiles_by_idx=profiles_by_idx,
        profile_sumsq_by_idx=profile_sumsq_by_idx,
        depth_edges_cm=depth_edges_cm,
        density_gcm3=density,
        dose_depth_cm=dose_depth_cm,
    )

    range_results_dir = Path(manifest["range_results_dir"]).resolve()
    out_file = _range_result_path(range_results_dir, int(result["range_id"]))
    _save_range_result(out_file, result)

    print(
        f"range_id={result['range_id']} "
        f"cells={result['n_cells']} "
        f"used_energies={result['used_energies']} "
        f"missing_energies={result['missing_energies']} "
        f"-> {out_file}"
    )


def _degree_formatter(value: float, _pos: int) -> str:
    if np.isclose(value, 0.0):
        value = 0.0
    if float(value).is_integer():
        label = f"{int(value)}"
    else:
        label = f"{value:g}"
    return rf"${label}\!^\circ$"


def _add_adaptive_grid(
    ax: Any,
    data: np.ndarray,
    im: Any,
    extent: tuple[float, float, float, float],
    linewidth: float = 0.30,
    alpha: float = 0.50,
) -> None:
    from matplotlib.collections import LineCollection

    if data.ndim != 2 or data.size == 0:
        return

    n_rows, n_cols = data.shape
    rgba = im.to_rgba(data)
    luminance = (
        0.2126 * rgba[..., 0]
        + 0.7152 * rgba[..., 1]
        + 0.0722 * rgba[..., 2]
    )

    x_edges = np.linspace(extent[0], extent[1], n_cols + 1)
    y_edges = np.linspace(extent[3], extent[2], n_rows + 1)

    segments: list[list[tuple[float, float]]] = []
    colors: list[tuple[float, float, float, float]] = []

    for i in range(n_rows):
        y_top = float(y_edges[i])
        y_bottom = float(y_edges[i + 1])
        for j in range(n_cols):
            x_left = float(x_edges[j])
            x_right = float(x_edges[j + 1])
            is_bright = bool(np.isfinite(luminance[i, j]) and luminance[i, j] > 0.5)
            color = (0.0, 0.0, 0.0, alpha) if is_bright else (1.0, 1.0, 1.0, alpha)
            segments.extend(
                [
                    [(x_left, y_top), (x_right, y_top)],
                    [(x_left, y_bottom), (x_right, y_bottom)],
                    [(x_left, y_top), (x_left, y_bottom)],
                    [(x_right, y_top), (x_right, y_bottom)],
                ]
            )
            colors.extend([color, color, color, color])

    grid = LineCollection(segments, colors=colors, linewidths=linewidth, zorder=3)
    ax.add_collection(grid)


def _robust_limits(
    arr: np.ndarray,
    floor_zero: bool = True,
    q_low: float = 2.0,
    q_high: float = 98.0,
) -> tuple[float, float]:
    vals = arr[np.isfinite(arr)]
    if vals.size == 0:
        return (0.0, 1.0)
    lo = float(np.percentile(vals, q_low))
    hi = float(np.percentile(vals, q_high))
    if floor_zero:
        lo = max(0.0, lo)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo = float(np.min(vals))
        hi = float(np.max(vals))
        if not np.isfinite(lo) or not np.isfinite(hi):
            return (0.0, 1.0)
        if hi <= lo:
            hi = lo + 1.0
    return lo, hi


def _format_sigfig_label(value: float, sigfigs: int = 4) -> str:
    if not np.isfinite(value):
        return "nan"
    if value == 0.0:
        return "0." + "0" * (sigfigs - 1)
    exponent = int(np.floor(np.log10(abs(value))))
    decimals = max(sigfigs - exponent - 1, 0)
    return f"{value:.{decimals}f}"


def _format_decimal_label(value: float, decimals: int = 5) -> str:
    if not np.isfinite(value):
        return "nan"
    text = f"{value:.{decimals}f}".rstrip("0").rstrip(".")
    return text if text else "0"


def _format_compact_sigfig_label(value: float, sigfigs: int = 3) -> str:
    if not np.isfinite(value):
        return "nan"
    return f"{value:.{sigfigs}g}"


def _configure_plot_style(force_interactive: bool = False) -> None:
    import matplotlib

    if force_interactive:
        if str(matplotlib.get_backend()).lower() == "agg":
            for candidate in ("MacOSX", "QtAgg", "TkAgg"):
                try:
                    matplotlib.use(candidate, force=True)
                    break
                except Exception:
                    continue
    else:
        # Keep headless/batch runs stable when no display is available.
        if (
            os.environ.get("DISPLAY", "") == ""
            and os.environ.get("WAYLAND_DISPLAY", "") == ""
            and sys.platform != "darwin"
        ):
            matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    plt.rcParams["font.family"] = FONT_COURIER
    plt.rcParams["mathtext.rm"] = FONT_COURIER
    plt.rcParams["mathtext.fontset"] = "custom"
    plt.rcParams.update(rcparams_with_fontsize(RC_BASE_ELASTIC, FONTSIZE_24))


def _plot_panel(
    fig: Any,
    ax: Any,
    data: np.ndarray,
    title: str | None,
    cbar_label: str,
    cmap_name: str,
    extent: tuple[float, float, float, float],
    xticks: list[float],
    fixed_limits: tuple[float, float] | None = None,
    floor_zero: bool = True,
    discrete: bool = False,
    low_threshold_positive: bool = False,
    n_ranks: int = 8,
    adaptive_grid: bool = False,
    colorbar_power_scale: bool = False,
    log_discrete: bool = False,
    show_empty_group: bool = False,
    empty_group_symbol: str = r"$\emptyset$",
    empty_group_color: str = "#bdbdbd",
    image_alpha: float = 1.0,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.ticker import FuncFormatter
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    vals = data[np.isfinite(data)]
    positive = vals[vals > 0.0]
    if fixed_limits is None:
        if low_threshold_positive and positive.size > 0:
            vmin = float(np.percentile(positive, 5.0))
            vmax = float(np.percentile(positive, 98.0))
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
                vmin = float(np.min(positive))
                vmax = float(np.max(positive))
                if vmax <= vmin:
                    vmax = vmin + 1.0
        else:
            vmin, vmax = _robust_limits(data, floor_zero=floor_zero)
    else:
        vmin, vmax = fixed_limits

    plot_data = np.asarray(data, dtype=np.float64)

    imshow_kwargs = {
        "origin": "upper",
        "extent": extent,
        "aspect": "equal",
        "alpha": image_alpha,
    }
    cbar_kwargs: dict[str, Any] = {}
    cbar_ticklabels: list[str] | None = None
    cbar_scale_title: str | None = None

    if discrete:
        if log_discrete:
            finite_positive = positive[np.isfinite(positive) & (positive > 0.0)]
            if finite_positive.size == 0:
                raise ValueError("log_discrete=True requires positive finite data.")
            vmin = max(float(vmin), float(np.min(finite_positive)))
            vmax = max(float(vmax), float(np.max(finite_positive)))
            if vmax <= vmin:
                vmax = vmin * (1.0 + 1e-6)
            base_bounds = np.geomspace(vmin, vmax, n_ranks + 1)
        else:
            base_bounds = np.linspace(vmin, vmax, n_ranks + 1)
        main_cmap = plt.get_cmap(cmap_name, n_ranks)
        plot_data_raw = np.asarray(plot_data, dtype=np.float64)
        if show_empty_group:
            if base_bounds.size >= 2:
                if log_discrete:
                    ratio = float(base_bounds[1] / base_bounds[0])
                    empty_lower = float(base_bounds[0] / ratio)
                else:
                    bin_width = float(base_bounds[1] - base_bounds[0])
                    empty_lower = float(vmin - bin_width)
            else:
                empty_lower = float(vmin - max(1.0, abs(vmax - vmin)))
            bounds = np.concatenate(([empty_lower], base_bounds))
            empty_center = 0.5 * (bounds[0] + bounds[1])
            cmap = ListedColormap([empty_group_color] + [main_cmap(i) for i in range(n_ranks)])
            norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=True)
            empty_mask = ~np.isfinite(plot_data_raw)
            plot_data = np.array(plot_data_raw, copy=True)
            plot_data = np.where(empty_mask, empty_center, plot_data)
            plot_data = np.where(~empty_mask & (plot_data < vmin), vmin, plot_data)
            upper_clip = np.nextafter(vmax, vmin)
            plot_data = np.where(~empty_mask & (plot_data >= vmax), upper_clip, plot_data)
            tick_positions = np.concatenate(([empty_center], base_bounds))
            if colorbar_power_scale and np.isfinite(vmax) and abs(vmax) >= 100.0:
                scale_exp = int(np.floor(np.log10(abs(vmax))))
                scale = 10.0 ** scale_exp
                cbar_scale_title = f"1e{scale_exp}"
                scaled_tick_positions = base_bounds / scale
                decimals = 2
                cbar_ticklabels = [empty_group_symbol] + [f"{x:.{decimals}f}" for x in scaled_tick_positions]
            else:
                cbar_ticklabels = [empty_group_symbol] + [_format_compact_sigfig_label(x, 3) for x in base_bounds]
            cbar_kwargs = {
                "boundaries": bounds,
                "ticks": tick_positions,
                "spacing": "uniform" if log_discrete else "proportional",
            }
        else:
            cmap = main_cmap.copy()
            under_color = cmap(0)
            cmap.set_under(under_color)
            cmap.set_bad(under_color)
            bounds = base_bounds
            norm = BoundaryNorm(bounds, ncolors=n_ranks, clip=False)
            tick_positions = bounds
            if colorbar_power_scale and np.isfinite(vmax) and abs(vmax) >= 100.0:
                scale_exp = int(np.floor(np.log10(abs(vmax))))
                scale = 10.0 ** scale_exp
                cbar_scale_title = f"1e{scale_exp}"
                scaled_tick_positions = tick_positions / scale
                decimals = 2
                cbar_ticklabels = [f"{x:.{decimals}f}" for x in scaled_tick_positions]
            else:
                cbar_ticklabels = [_format_compact_sigfig_label(x, 3) for x in tick_positions]
            under_value = float(vmin - max(1e-12, 1e-9 * max(abs(vmin), abs(vmax), 1.0)))
            plot_data = np.where(np.isfinite(plot_data_raw), plot_data_raw, under_value)
            plot_data = np.where(plot_data < vmin, under_value, plot_data)
            imshow_kwargs["cmap"] = cmap
            imshow_kwargs["norm"] = norm
            cbar_kwargs = {
                "extend": "min",
                "boundaries": bounds,
                "ticks": tick_positions,
                "spacing": "uniform" if log_discrete else "proportional",
            }
        imshow_kwargs["cmap"] = cmap
        imshow_kwargs["norm"] = norm
    else:
        cmap = plt.get_cmap(cmap_name).copy()
        cmap.set_bad("#d9d9d9")
        imshow_kwargs["cmap"] = cmap
        imshow_kwargs["vmin"] = vmin
        imshow_kwargs["vmax"] = vmax

    im = ax.imshow(plot_data, **imshow_kwargs)
    ax.set_facecolor("white")
    if adaptive_grid:
        _add_adaptive_grid(ax, plot_data, im, extent)
    if title:
        ax.set_title(title, pad=14)
    ax.set_xlabel(r"Longitude ($^\circ$ W)")
    ax.set_ylabel(r"Latitude ($^\circ$)")
    ax.set_xticks(xticks)
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.xaxis.set_major_formatter(FuncFormatter(_degree_formatter))
    ax.yaxis.set_major_formatter(FuncFormatter(_degree_formatter))
    ax.invert_xaxis()
    ax.tick_params(axis="x", pad=10)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, **cbar_kwargs)
    cbar.set_label(cbar_label)
    if cbar_ticklabels is not None:
        cbar.set_ticklabels(cbar_ticklabels)
    if cbar_scale_title is not None:
        cbar.ax.set_title(cbar_scale_title, pad=6)
    cbar.ax.tick_params(length=0)


def _plot_split_particle_panel(
    fig: Any,
    ax: Any,
    electron_data: np.ndarray,
    photon_data: np.ndarray,
    extent: tuple[float, float, float, float],
    xticks: list[float],
    title: str,
    n_ranks: int = 6,
) -> None:
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.ticker import FuncFormatter
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt

    electron = np.asarray(electron_data, dtype=np.float64)
    photon = np.asarray(photon_data, dtype=np.float64)
    if electron.shape != photon.shape:
        raise ValueError("electron_data and photon_data must have the same shape.")

    def discrete_style(data: np.ndarray, cmap_name: str) -> tuple[Any, Any, np.ndarray, np.ndarray, list[str], float, float]:
        vals = np.asarray(data, dtype=np.float64)
        vmin, vmax = _robust_limits(vals, floor_zero=True)
        if vmax <= vmin:
            vmax = vmin + 1.0e-6
        base_bounds = np.linspace(vmin, vmax, n_ranks + 1)
        bin_width = float(base_bounds[1] - base_bounds[0]) if base_bounds.size >= 2 else 1.0
        empty_lower = float(vmin - bin_width)
        bounds = np.concatenate(([empty_lower], base_bounds))
        empty_center = 0.5 * (bounds[0] + bounds[1])
        main_cmap = plt.get_cmap(cmap_name, n_ranks)
        cmap = ListedColormap(["#bdbdbd"] + [main_cmap(i) for i in range(n_ranks)])
        norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=True)
        tick_positions = np.concatenate(([empty_center], base_bounds))
        ticklabels = [r"$\emptyset$"] + [_format_compact_sigfig_label(x, 3) for x in base_bounds]
        return cmap, norm, bounds, tick_positions, ticklabels, vmin, vmax

    electron_cmap, electron_norm, electron_bounds, electron_ticks, electron_ticklabels, e_vmin, e_vmax = (
        discrete_style(electron, "cividis")
    )
    photon_cmap, photon_norm, photon_bounds, photon_ticks, photon_ticklabels, p_vmin, p_vmax = (
        discrete_style(photon, "plasma")
    )

    n_rows, n_cols = electron.shape
    x_edges = np.linspace(extent[0], extent[1], n_cols + 1)
    y_edges = np.linspace(extent[3], extent[2], n_rows + 1)

    empty_polys: list[list[tuple[float, float]]] = []
    electron_polys: list[list[tuple[float, float]]] = []
    electron_colors: list[tuple[float, float, float, float]] = []
    photon_polys: list[list[tuple[float, float]]] = []
    photon_colors: list[tuple[float, float, float, float]] = []

    for i in range(n_rows):
        y_top = float(y_edges[i])
        y_bottom = float(y_edges[i + 1])
        for j in range(n_cols):
            x_left = float(x_edges[j])
            x_right = float(x_edges[j + 1])
            e_val = float(electron[i, j])
            p_val = float(photon[i, j])
            has_e = np.isfinite(e_val)
            has_p = np.isfinite(p_val)
            if not has_e and not has_p:
                empty_polys.append(
                    [(x_left, y_top), (x_right, y_top), (x_right, y_bottom), (x_left, y_bottom)]
                )
                continue
            if has_e:
                e_val = max(e_vmin, min(np.nextafter(e_vmax, e_vmin), e_val))
                electron_polys.append([(x_left, y_bottom), (x_left, y_top), (x_right, y_bottom)])
                electron_colors.append(electron_cmap(electron_norm(e_val)))
            if has_p:
                p_val = max(p_vmin, min(np.nextafter(p_vmax, p_vmin), p_val))
                photon_polys.append([(x_left, y_top), (x_right, y_top), (x_right, y_bottom)])
                photon_colors.append(photon_cmap(photon_norm(p_val)))

    if empty_polys:
        ax.add_collection(
            PolyCollection(empty_polys, facecolors="#d9d9d9", edgecolors="none", zorder=0)
        )
    if electron_polys:
        ax.add_collection(
            PolyCollection(electron_polys, facecolors=electron_colors, edgecolors="none", zorder=1)
        )
    if photon_polys:
        ax.add_collection(
            PolyCollection(photon_polys, facecolors=photon_colors, edgecolors="none", zorder=2)
        )

    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal")
    ax.set_facecolor("white")
    ax.set_title(title, pad=18)
    ax.set_xlabel(r"Longitude ($^\circ$ W)")
    ax.set_ylabel(r"Latitude ($^\circ$)")
    ax.set_xticks(xticks)
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.xaxis.set_major_formatter(FuncFormatter(_degree_formatter))
    ax.yaxis.set_major_formatter(FuncFormatter(_degree_formatter))
    ax.invert_xaxis()
    ax.tick_params(axis="x", pad=10)
    ax.grid(color="white", linewidth=0.55, alpha=0.9)

    divider = make_axes_locatable(ax)
    cax_e = divider.append_axes("right", size="3.8%", pad=0.05)
    cax_p = divider.append_axes("right", size="3.8%", pad=0.80)
    cbar_e = fig.colorbar(
        plt.cm.ScalarMappable(norm=electron_norm, cmap=electron_cmap),
        cax=cax_e,
        boundaries=electron_bounds,
        ticks=electron_ticks,
        spacing="proportional",
    )
    cbar_p = fig.colorbar(
        plt.cm.ScalarMappable(norm=photon_norm, cmap=photon_cmap),
        cax=cax_p,
        boundaries=photon_bounds,
        ticks=photon_ticks,
        spacing="proportional",
    )
    cbar_e.ax.set_title("e-", pad=6, fontsize=float(FONTSIZE_24) * 0.55)
    cbar_p.ax.set_title("photon", pad=6, fontsize=float(FONTSIZE_24) * 0.55)
    cbar_e.set_ticklabels(electron_ticklabels)
    cbar_p.set_ticklabels(photon_ticklabels)
    cbar_e.ax.tick_params(length=0, labelsize=float(FONTSIZE_24) * 0.50, pad=2)
    cbar_p.ax.tick_params(length=0, labelsize=float(FONTSIZE_24) * 0.50, pad=2)


def _finite_limits_for_maps(maps: list[np.ndarray], floor_zero: bool = True) -> tuple[float, float]:
    vals = np.concatenate([np.asarray(m, dtype=np.float64).ravel() for m in maps])
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return (0.0, 1.0)
    vmin, vmax = _robust_limits(vals, floor_zero=floor_zero)
    if vmax <= vmin:
        vmax = vmin + 1.0e-6
    return vmin, vmax


def _subset_hemisphere(
    full_map: np.ndarray,
    lon_values: np.ndarray,
    hemisphere: str,
) -> tuple[np.ndarray, tuple[float, float, float, float], list[float]]:
    if hemisphere == "leading":
        cols = np.where(lon_values < 180.0)[0]
        extent = (0.0, 180.0, -90.0, 90.0)
        xticks = [0, 30, 60, 90, 120, 150, 180]
    else:
        cols = np.where(lon_values >= 180.0)[0]
        extent = (180.0, 360.0, -90.0, 90.0)
        xticks = [180, 210, 240, 270, 300, 330, 360]
    return full_map[:, cols], extent, xticks

def _npz_first_key(npz: Any, names: list[str]) -> str:
    for name in names:
        if name in npz:
            return name
    raise KeyError(f"None of keys found in NPZ: {names}")


def _npz_get_array(npz: Any, names: list[str]) -> np.ndarray:
    key = _npz_first_key(npz, names)
    return np.asarray(npz[key], dtype=np.float64)


def _npz_get_optional_array(npz: Any, names: list[str]) -> np.ndarray | None:
    for name in names:
        if name in npz:
            return np.asarray(npz[name], dtype=np.float64)
    return None


def _npz_get_scalar(npz: Any, names: list[str], default: float = float("nan")) -> float:
    for name in names:
        if name in npz:
            return float(np.asarray(npz[name]).item())
    return float(default)


def _load_density_gcm3_from_manifest(out_dir: Path, fallback_dir: Path | None = None) -> float:
    candidate_dirs = [out_dir]
    if fallback_dir is not None and fallback_dir not in candidate_dirs:
        candidate_dirs.append(fallback_dir)
    for base_dir in candidate_dirs:
        manifest_path = base_dir / "analysis_manifest.json"
        if not manifest_path.exists():
            continue
        try:
            manifest = json.loads(manifest_path.read_text())
        except Exception:
            continue
        density = _safe_float(str(manifest.get("density_gcm3_default", "nan")))
        if np.isfinite(density) and density > 0.0:
            return float(density)
    return float("nan")


def _resolve_manifest_cached_path(
    manifest_file: Path,
    configured_path: str | Path | None,
    local_name: str,
    expect_dir: bool = False,
) -> Path:
    manifest_dir = manifest_file.resolve().parent
    candidates: list[Path] = []
    if configured_path is not None:
        candidates.append(Path(configured_path).resolve())
    candidates.append((manifest_dir / local_name).resolve())
    for candidate in candidates:
        if expect_dir:
            if candidate.is_dir():
                return candidate
        else:
            if candidate.exists():
                return candidate
    return candidates[-1]


def _nearest_latitude_index(lat_values: np.ndarray, nominal_lat_deg: float) -> int:
    lat_arr = np.asarray(lat_values, dtype=np.float64)
    diffs = np.abs(lat_arr - float(nominal_lat_deg))
    best = float(np.nanmin(diffs))
    candidates = np.where(np.isclose(diffs, best, rtol=0.0, atol=1.0e-12))[0]
    if candidates.size == 1:
        return int(candidates[0])

    # Keep 0 aligned to the north-side grid point and prefer same-sign matches elsewhere.
    if np.isclose(nominal_lat_deg, 0.0):
        non_negative = candidates[lat_arr[candidates] >= 0.0]
        if non_negative.size > 0:
            return int(non_negative[0])
    same_sign = candidates[np.sign(lat_arr[candidates]) == np.sign(nominal_lat_deg)]
    if same_sign.size > 0:
        return int(same_sign[0])
    return int(candidates[0])


def _resolve_target_grid(args: argparse.Namespace) -> tuple[int, int] | None:
    size = int(getattr(args, "jde_grid_size", 0) or 0)
    nlat = int(getattr(args, "jde_grid_nlat", 0) or 0)
    nlon = int(getattr(args, "jde_grid_nlon", 0) or 0)
    if nlat > 0 and nlon > 0:
        return nlat, nlon
    if size > 0:
        return size, size
    return None


def _resample_2d_map(
    src_map: np.ndarray,
    src_lat: np.ndarray,
    src_lon: np.ndarray,
    nlat: int,
    nlon: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if nlat < 2 or nlon < 2:
        raise ValueError("Requested grid must have at least 2x2 points.")
    arr = np.asarray(src_map, dtype=np.float64)
    lat = np.asarray(src_lat, dtype=np.float64)
    lon = np.asarray(src_lon, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError("Source map must be 2D.")
    if arr.shape != (lat.size, lon.size):
        raise ValueError(
            f"Map shape {arr.shape} does not match lat/lon sizes ({lat.size}, {lon.size})."
        )

    lat_descending = bool(lat[0] > lat[-1])
    if lat_descending:
        lat_work = lat[::-1]
        arr_work = arr[::-1, :]
    else:
        lat_work = lat
        arr_work = arr

    if lon[0] > lon[-1]:
        lon_work = lon[::-1]
        arr_work = arr_work[:, ::-1]
    else:
        lon_work = lon

    tgt_lat_work = np.linspace(float(lat_work[0]), float(lat_work[-1]), int(nlat))
    tgt_lon = np.linspace(float(lon_work[0]), float(lon_work[-1]), int(nlon))

    interp_lon = np.full((arr_work.shape[0], tgt_lon.size), np.nan, dtype=np.float64)
    for i in range(arr_work.shape[0]):
        row = arr_work[i, :]
        valid = np.isfinite(row)
        if np.count_nonzero(valid) >= 2:
            interp_lon[i, :] = np.interp(
                tgt_lon,
                lon_work[valid],
                row[valid],
                left=np.nan,
                right=np.nan,
            )
        elif np.count_nonzero(valid) == 1:
            interp_lon[i, :] = float(row[valid][0])

    out_work = np.full((tgt_lat_work.size, tgt_lon.size), np.nan, dtype=np.float64)
    for j in range(tgt_lon.size):
        col = interp_lon[:, j]
        valid = np.isfinite(col)
        if np.count_nonzero(valid) >= 2:
            out_work[:, j] = np.interp(
                tgt_lat_work,
                lat_work[valid],
                col[valid],
                left=np.nan,
                right=np.nan,
            )
        elif np.count_nonzero(valid) == 1:
            out_work[:, j] = float(col[valid][0])

    if lat_descending:
        return out_work[::-1, :], tgt_lat_work[::-1], tgt_lon
    return out_work, tgt_lat_work, tgt_lon


def _compute_cumulative_fraction_depth_map(
    depth_edges_cm: np.ndarray,
    dose_profile_map_mgy_per_yr: np.ndarray,
    target_fraction: float = 1.0 - 1.0 / math.e,
) -> np.ndarray:
    depth_edges = np.asarray(depth_edges_cm, dtype=np.float64)
    profiles = np.asarray(dose_profile_map_mgy_per_yr, dtype=np.float64)
    if profiles.ndim != 3:
        raise ValueError("dose_profile_map_mgy_per_yr must be a 3D lat-lon-depth array.")
    if depth_edges.size != profiles.shape[2] + 1:
        raise ValueError("depth_edges_cm must bracket the depth-profile bins.")
    if not np.isfinite(target_fraction) or not (0.0 < target_fraction < 1.0):
        raise ValueError("target_fraction must be finite and strictly between 0 and 1.")

    dz_cm = depth_edges[1:] - depth_edges[:-1]
    out = np.full(profiles.shape[:2], np.nan, dtype=np.float64)

    for ii in range(profiles.shape[0]):
        for jj in range(profiles.shape[1]):
            profile = np.asarray(profiles[ii, jj, :], dtype=np.float64)
            valid = (
                np.isfinite(profile)
                & np.isfinite(dz_cm)
                & (profile > 0.0)
                & (dz_cm > 0.0)
            )
            if not np.any(valid):
                continue

            shell_weight = np.where(valid, profile * dz_cm, 0.0)
            total_weight = float(np.sum(shell_weight))
            if not np.isfinite(total_weight) or total_weight <= 0.0:
                continue
            target_weight = total_weight * float(target_fraction)
            cumulative = np.cumsum(shell_weight)
            crossed = np.flatnonzero(cumulative >= target_weight)
            if crossed.size == 0:
                out[ii, jj] = float(depth_edges[-1])
                continue
            hit_idx = int(crossed[0])
            prev_weight = float(cumulative[hit_idx - 1]) if hit_idx > 0 else 0.0
            bin_weight = float(shell_weight[hit_idx])
            x0 = float(depth_edges[hit_idx])
            x1 = float(depth_edges[hit_idx + 1])
            if not np.isfinite(bin_weight) or bin_weight <= 0.0:
                out[ii, jj] = x1
                continue
            frac = (target_weight - prev_weight) / bin_weight
            frac = float(np.clip(frac, 0.0, 1.0))
            out[ii, jj] = x0 + frac * (x1 - x0)

    return out


def _save_hemisphere_plots(
    out_dir: Path,
    hemisphere: str,
    lon_values: np.ndarray,
    energy_panel_flux_map: np.ndarray,
    primary_flux_map: np.ndarray,
    primary_pen_map: np.ndarray,
    secondary_ke_map: np.ndarray,
    secondary_pen_map: np.ndarray,
    cumulative_depth_map: np.ndarray | None,
    dep_fraction_map: np.ndarray,
    backscatter_amount_map: np.ndarray,
    dose_depth_cm: float,
    backscatter_fraction_map: np.ndarray | None = None,
    lateral_escape_fraction_map: np.ndarray | None = None,
    forward_escape_fraction_map: np.ndarray | None = None,
    escaped_electron_fraction_map: np.ndarray | None = None,
    escaped_photon_fraction_map: np.ndarray | None = None,
    backscatter_electron_fraction_map: np.ndarray | None = None,
    backscatter_photon_fraction_map: np.ndarray | None = None,
    lateral_escape_electron_fraction_map: np.ndarray | None = None,
    lateral_escape_photon_fraction_map: np.ndarray | None = None,
    forward_escape_electron_fraction_map: np.ndarray | None = None,
    forward_escape_photon_fraction_map: np.ndarray | None = None,
    jde_energy_flux_map_mev_cm2_s: np.ndarray | None = None,
    show_plots: bool = False,
) -> None:
    import matplotlib.pyplot as plt

    energy_h, extent, xticks = _subset_hemisphere(energy_panel_flux_map, lon_values, hemisphere)
    ppen_h, _, _ = _subset_hemisphere(primary_pen_map, lon_values, hemisphere)
    ske_h, _, _ = _subset_hemisphere(secondary_ke_map, lon_values, hemisphere)
    ske_h_ev = ske_h * 1.0e6
    if cumulative_depth_map is not None:
        spen_h, _, _ = _subset_hemisphere(cumulative_depth_map, lon_values, hemisphere)
        penetration_label = "e-fold deposition depth (cm)"
    else:
        spen_h, _, _ = _subset_hemisphere(secondary_pen_map, lon_values, hemisphere)
        penetration_label = "Penetration depth (cm)"
    depf_h, _, _ = _subset_hemisphere(dep_fraction_map, lon_values, hemisphere)
    depf_4panel_h = np.asarray(depf_h, dtype=np.float64)

    # Cells with no relevant electrons in the modeled range should render as an explicit
    # empty category across the 4-panel map, even if one derived map happens to contain 0.0.
    empty_mask_4panel = ~np.isfinite(ppen_h) | ~np.isfinite(ske_h_ev) | ~np.isfinite(spen_h)
    energy_h_4panel = np.where(empty_mask_4panel, np.nan, energy_h)
    depf_h_4panel = np.where(empty_mask_4panel, np.nan, depf_4panel_h)
    ske_h_ev_4panel = np.where(empty_mask_4panel, np.nan, ske_h_ev)
    spen_h_4panel = np.where(empty_mask_4panel, np.nan, spen_h)
    depf_finite_4panel = depf_h_4panel[np.isfinite(depf_h_4panel)]
    depf_vmin_4panel = float(np.min(depf_finite_4panel)) if depf_finite_4panel.size > 0 else 0.0
    depf_vmax_4panel = float(np.max(depf_finite_4panel)) if depf_finite_4panel.size > 0 else 1.0
    if depf_vmax_4panel <= depf_vmin_4panel:
        depf_vmax_4panel = depf_vmin_4panel + 1.0e-6

    # Match the 2-panel europa electron energy map geometry so stacked figures
    # align like a 3x2 matrix in document layouts.
    fig1, axes1 = plt.subplots(2, 2, figsize=(18, 18), facecolor="white")
    fig1.subplots_adjust(left=0.05, right=0.98, bottom=0.10, top=0.92, wspace=0.45, hspace=0.02)

    n_ranks = 6

    _plot_panel(
        fig1,
        axes1[0, 0],
        energy_h_4panel,
        None,
        r"Deposited energy flux (MeV cm$^{-2}$ s$^{-1}$)",
        "plasma",
        extent,
        xticks,
        floor_zero=True,
        discrete=True,
        n_ranks=n_ranks,
        adaptive_grid=True,
        colorbar_power_scale=True,
        show_empty_group=True,
    )
    axes1[0, 0].set_xlabel("")
    axes1[0, 0].tick_params(labelbottom=False)
    _plot_panel(
        fig1,
        axes1[0, 1],
        depf_h_4panel,
        None,
        "Deposited energy fraction",
        "cividis",
        extent,
        xticks,
        fixed_limits=(depf_vmin_4panel, depf_vmax_4panel),
        floor_zero=True,
        discrete=True,
        n_ranks=n_ranks,
        adaptive_grid=True,
        show_empty_group=True,
    )
    axes1[0, 1].set_xlabel("")
    axes1[0, 1].tick_params(labelbottom=False)
    axes1[0, 1].set_ylabel("")
    _plot_panel(
        fig1,
        axes1[1, 0],
        ske_h_ev_4panel,
        None,
        "$\\bar{T}_{\\rm{secondary}}$ (eV)",
        "Purples",
        extent,
        xticks,
        floor_zero=False,
        discrete=True,
        low_threshold_positive=True,
        n_ranks=n_ranks,
        adaptive_grid=True,
        colorbar_power_scale=True,
        show_empty_group=True,
    )
    _plot_panel(
        fig1,
        axes1[1, 1],
        spen_h_4panel,
        None,
        penetration_label,
        "Greens",
        extent,
        xticks,
        floor_zero=False,
        discrete=True,
        low_threshold_positive=True,
        n_ranks=n_ranks,
        adaptive_grid=True,
        show_empty_group=True,
    )
    axes1[1, 1].set_ylabel("")
    panel_labels = (("(a)", axes1[0, 0]), ("(b)", axes1[0, 1]), ("(c)", axes1[1, 0]), ("(d)", axes1[1, 1]))
    for label, ax in panel_labels:
        ax.text(
            0.0,
            1.04,
            label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=float(FONTSIZE_24),
            clip_on=False,
        )
    four_panel_out = out_dir / f"europa_{hemisphere}_4panel_metrics.png"
    fig1.savefig(four_panel_out, bbox_inches="tight")
    if show_plots:
        plt.show()

    if (
        backscatter_fraction_map is not None
        and lateral_escape_fraction_map is not None
        and forward_escape_fraction_map is not None
    ):
        primary_h, _, _ = _subset_hemisphere(primary_flux_map, lon_values, hemisphere)
        backf_h, _, _ = _subset_hemisphere(backscatter_fraction_map, lon_values, hemisphere)
        latf_h, _, _ = _subset_hemisphere(lateral_escape_fraction_map, lon_values, hemisphere)
        forf_h, _, _ = _subset_hemisphere(forward_escape_fraction_map, lon_values, hemisphere)
        directional_escape_available = (
            np.any(np.isfinite(latf_h)) and np.any(np.isfinite(forf_h))
        )
        if directional_escape_available:
            empty_mask_escape = ~np.isfinite(primary_h)
            depf_escape_h = np.where(empty_mask_escape, np.nan, depf_h)
            backf_escape_h = np.where(empty_mask_escape, np.nan, backf_h)
            latf_escape_h = np.where(empty_mask_escape, np.nan, latf_h)
            forf_escape_h = np.where(empty_mask_escape, np.nan, forf_h)

            fig2, axes2 = plt.subplots(2, 2, figsize=(18, 18), facecolor="white")
            fig2.subplots_adjust(left=0.05, right=0.98, bottom=0.10, top=0.92, wspace=0.45, hspace=0.02)

            _plot_panel(
                fig2,
                axes2[0, 0],
                depf_escape_h,
                None,
                "Deposited energy fraction",
                "cividis",
                extent,
                xticks,
                floor_zero=True,
                discrete=True,
                n_ranks=n_ranks,
                adaptive_grid=True,
                show_empty_group=True,
            )
            axes2[0, 0].set_xlabel("")
            axes2[0, 0].tick_params(labelbottom=False)
            _plot_panel(
                fig2,
                axes2[0, 1],
                backf_escape_h,
                None,
                "Backscattered fraction",
                "cividis",
                extent,
                xticks,
                floor_zero=True,
                discrete=True,
                n_ranks=n_ranks,
                adaptive_grid=True,
                show_empty_group=True,
            )
            axes2[0, 1].set_xlabel("")
            axes2[0, 1].tick_params(labelbottom=False)
            axes2[0, 1].set_ylabel("")
            _plot_panel(
                fig2,
                axes2[1, 0],
                latf_escape_h,
                None,
                "Lateral escape fraction",
                "cividis",
                extent,
                xticks,
                floor_zero=True,
                discrete=True,
                n_ranks=n_ranks,
                adaptive_grid=True,
                show_empty_group=True,
            )
            _plot_panel(
                fig2,
                axes2[1, 1],
                forf_escape_h,
                None,
                "Forward escape fraction",
                "cividis",
                extent,
                xticks,
                floor_zero=True,
                discrete=True,
                n_ranks=n_ranks,
                adaptive_grid=True,
                show_empty_group=True,
            )
            axes2[1, 1].set_ylabel("")
            axes2[0, 0].set_title("Total deposited", pad=18)
            axes2[0, 1].set_title("Backscattered", pad=18)
            axes2[1, 0].set_title("Lateral escape", pad=18)
            axes2[1, 1].set_title("Forward escape", pad=18)
            escape_panel_out = out_dir / f"europa_{hemisphere}_4panel_escape_fractions.png"
            fig2.savefig(escape_panel_out, bbox_inches="tight")
            if show_plots:
                plt.show()

    def save_particle_escape_plot(
        particle: str,
        total_map: np.ndarray | None,
        back_map: np.ndarray | None,
        lateral_map: np.ndarray | None,
        forward_map: np.ndarray | None,
    ) -> None:
        if total_map is None or back_map is None or lateral_map is None or forward_map is None:
            return
        primary_h, _, _ = _subset_hemisphere(primary_flux_map, lon_values, hemisphere)
        total_h, _, _ = _subset_hemisphere(total_map, lon_values, hemisphere)
        back_h, _, _ = _subset_hemisphere(back_map, lon_values, hemisphere)
        lateral_h, _, _ = _subset_hemisphere(lateral_map, lon_values, hemisphere)
        forward_h, _, _ = _subset_hemisphere(forward_map, lon_values, hemisphere)
        if not (
            np.any(np.isfinite(total_h))
            or np.any(np.isfinite(back_h))
            or np.any(np.isfinite(lateral_h))
            or np.any(np.isfinite(forward_h))
        ):
            return

        empty_mask = ~np.isfinite(primary_h)
        total_h = np.where(empty_mask, np.nan, total_h)
        back_h = np.where(empty_mask, np.nan, back_h)
        lateral_h = np.where(empty_mask, np.nan, lateral_h)
        forward_h = np.where(empty_mask, np.nan, forward_h)

        figp, axesp = plt.subplots(2, 2, figsize=(18, 18), facecolor="white")
        figp.subplots_adjust(left=0.05, right=0.98, bottom=0.10, top=0.92, wspace=0.45, hspace=0.02)
        panel_specs = [
            (axesp[0, 0], total_h, f"Escaped {particle} fraction", f"Total escaped {particle}"),
            (axesp[0, 1], back_h, f"Backscattered {particle} fraction", f"Backscattered {particle}"),
            (axesp[1, 0], lateral_h, f"Lateral escaped {particle} fraction", f"Lateral escaped {particle}"),
            (axesp[1, 1], forward_h, f"Forward escaped {particle} fraction", f"Forward escaped {particle}"),
        ]
        for ax, panel_data, cbar_label, title in panel_specs:
            _plot_panel(
                figp,
                ax,
                panel_data,
                None,
                cbar_label,
                "cividis",
                extent,
                xticks,
                floor_zero=True,
                discrete=True,
                n_ranks=n_ranks,
                adaptive_grid=True,
                show_empty_group=True,
            )
            ax.set_title(title, pad=18)
        axesp[0, 0].set_xlabel("")
        axesp[0, 0].tick_params(labelbottom=False)
        axesp[0, 1].set_xlabel("")
        axesp[0, 1].tick_params(labelbottom=False)
        axesp[0, 1].set_ylabel("")
        axesp[1, 1].set_ylabel("")
        particle_out = out_dir / f"europa_{hemisphere}_{particle}_4panel_escape_fractions.png"
        figp.savefig(particle_out, bbox_inches="tight")
        if show_plots:
            plt.show()

    def save_electron_photon_escape_8panel_plot() -> None:
        electron_maps = [
            escaped_electron_fraction_map,
            backscatter_electron_fraction_map,
            lateral_escape_electron_fraction_map,
            forward_escape_electron_fraction_map,
        ]
        photon_maps = [
            escaped_photon_fraction_map,
            backscatter_photon_fraction_map,
            lateral_escape_photon_fraction_map,
            forward_escape_photon_fraction_map,
        ]
        if any(m is None for m in electron_maps + photon_maps):
            return

        primary_h, _, _ = _subset_hemisphere(primary_flux_map, lon_values, hemisphere)
        empty_mask = ~np.isfinite(primary_h)

        electron_h = [
            np.where(empty_mask, np.nan, _subset_hemisphere(m, lon_values, hemisphere)[0])
            for m in electron_maps
            if m is not None
        ]
        photon_h = [
            np.where(empty_mask, np.nan, _subset_hemisphere(m, lon_values, hemisphere)[0])
            for m in photon_maps
            if m is not None
        ]
        if not (any(np.any(np.isfinite(m)) for m in electron_h) or any(np.any(np.isfinite(m)) for m in photon_h)):
            return

        figep, axesep = plt.subplots(4, 2, figsize=(18, 36), facecolor="white")
        figep.subplots_adjust(left=0.05, right=0.98, bottom=0.04, top=0.97, wspace=0.45, hspace=0.08)
        panel_specs = [
            (0, 0, "Total escaped energy (electrons)", electron_h[0]),
            (0, 1, "Backscattered energy (electrons)", electron_h[1]),
            (1, 0, "Lateral escape energy (electrons)", electron_h[2]),
            (1, 1, "Forward escape energy (electrons)", electron_h[3]),
            (2, 0, "Total escaped energy (photons)", photon_h[0]),
            (2, 1, "Backscattered energy (photons)", photon_h[1]),
            (3, 0, "Lateral escape energy (photons)", photon_h[2]),
            (3, 1, "Forward escape energy (photons)", photon_h[3]),
        ]
        for row, col, title, panel_data in panel_specs:
            _plot_panel(
                figep,
                axesep[row, col],
                panel_data,
                None,
                "Fraction",
                "cividis",
                extent,
                xticks,
                floor_zero=True,
                discrete=True,
                n_ranks=n_ranks,
                adaptive_grid=True,
                show_empty_group=True,
            )
            axesep[row, col].set_title(title, pad=18)
            if col == 1:
                axesep[row, col].set_ylabel("")
            if row < 3:
                axesep[row, col].set_xlabel("")
                axesep[row, col].tick_params(labelbottom=False)

        panel_labels = (
            ("(a)", axesep[0, 0]),
            ("(b)", axesep[0, 1]),
            ("(c)", axesep[1, 0]),
            ("(d)", axesep[1, 1]),
            ("(e)", axesep[2, 0]),
            ("(f)", axesep[2, 1]),
            ("(g)", axesep[3, 0]),
            ("(h)", axesep[3, 1]),
        )
        for label, ax in panel_labels:
            ax.text(
                -0.14,
                1.02,
                label,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=float(FONTSIZE_24),
                clip_on=False,
            )

        eight_panel_out = out_dir / f"europa_{hemisphere}_electron_photon_8panel_escape_fractions.png"
        figep.savefig(eight_panel_out, bbox_inches="tight")
        if show_plots:
            plt.show()

    save_particle_escape_plot(
        "electron",
        escaped_electron_fraction_map,
        backscatter_electron_fraction_map,
        lateral_escape_electron_fraction_map,
        forward_escape_electron_fraction_map,
    )
    save_particle_escape_plot(
        "photon",
        escaped_photon_fraction_map,
        backscatter_photon_fraction_map,
        lateral_escape_photon_fraction_map,
        forward_escape_photon_fraction_map,
    )
    save_electron_photon_escape_8panel_plot()

    if jde_energy_flux_map_mev_cm2_s is not None:
        jde_h, _, _ = _subset_hemisphere(jde_energy_flux_map_mev_cm2_s, lon_values, hemisphere)
        fig3, ax3 = plt.subplots(1, 1, figsize=(12, 7))
        fig3.subplots_adjust(left=0.08, right=0.95, bottom=0.10, top=0.78)
        _plot_panel(
            fig3,
            ax3,
            jde_h,
            r"$\Sigma\,J(E)\,dE \cdot \langle E \rangle$",
            r"MeV cm$^{-2}$ s$^{-1}$",
            "turbo",
            extent,
            xticks,
            floor_zero=True,
        )
        jde_out = out_dir / f"europa_{hemisphere}_jde_diagnostic.pdf"
        fig3.savefig(jde_out, bbox_inches="tight")
        if show_plots:
            plt.show()


def _save_energy_deposition_depth_profiles_plot(
    out_dir: Path,
    lat_values: np.ndarray,
    lon_values: np.ndarray,
    depth_edges_cm: np.ndarray,
    dose_profile_map_mgy_per_yr: np.ndarray,
    density_gcm3: float,
    dose_profile_std_map_mgy_per_yr: np.ndarray | None = None,
    nominal_latitudes_deg: tuple[float, ...] = (0.0, 30.0, 60.0),
    show_plots: bool = False,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, LogFormatterMathtext, LogLocator, NullFormatter

    lat_arr = np.asarray(lat_values, dtype=np.float64)
    lon_arr = np.asarray(lon_values, dtype=np.float64)
    depth_edges = np.asarray(depth_edges_cm, dtype=np.float64)
    dose_profiles = np.asarray(dose_profile_map_mgy_per_yr, dtype=np.float64)
    if dose_profiles.ndim != 3:
        raise ValueError("dose_profile_map_mgy_per_yr must be a 3D lat-lon-depth array.")
    if depth_edges.size != dose_profiles.shape[2] + 1:
        raise ValueError("depth_edges_cm must bracket the depth-profile bins.")
    dose_profile_std = None
    if dose_profile_std_map_mgy_per_yr is not None:
        dose_profile_std = np.asarray(dose_profile_std_map_mgy_per_yr, dtype=np.float64)
        if dose_profile_std.shape != dose_profiles.shape:
            raise ValueError("dose_profile_std_map_mgy_per_yr must match dose_profile_map_mgy_per_yr.")

    depth_centers_cm = 0.5 * (depth_edges[:-1] + depth_edges[1:])
    positive_depth_centers = depth_centers_cm[np.isfinite(depth_centers_cm) & (depth_centers_cm > 0.0)]
    if positive_depth_centers.size == 0:
        raise ValueError("depth_edges_cm must include positive depth bins for log-scale plotting.")
    x_min_cm = float(np.min(positive_depth_centers))
    dose_rate_profiles = dose_profiles
    base_fontsize = float(FONTSIZE_24)
    plotted_depth_max_cm = float("nan")
    for nominal_lon_deg in (90.0, 270.0):
        lon_idx = int(np.argmin(np.abs(lon_arr - nominal_lon_deg)))
        for nominal_lat in nominal_latitudes_deg:
            lat_idx = _nearest_latitude_index(lat_arr, nominal_lat)
            mean_profile = dose_rate_profiles[lat_idx, lon_idx, :]
            valid = (
                np.isfinite(depth_centers_cm)
                & np.isfinite(mean_profile)
                & (depth_centers_cm > 0.0)
                & (mean_profile > 0.0)
            )
            if np.any(valid):
                last_depth = float(np.max(depth_centers_cm[valid]))
                if not np.isfinite(plotted_depth_max_cm) or last_depth > plotted_depth_max_cm:
                    plotted_depth_max_cm = last_depth
    if not np.isfinite(plotted_depth_max_cm):
        plotted_depth_max_cm = float(np.max(positive_depth_centers))
    x_tick_step_cm = 10.0 ** max(np.floor(np.log10(plotted_depth_max_cm)) - 1.0, -5.0)
    x_max_cm = float(np.ceil(plotted_depth_max_cm / x_tick_step_cm) * x_tick_step_cm)
    x_max_cm = min(x_max_cm, float(depth_edges[-1]))

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 10),
        sharex=True,
        facecolor="white",
        gridspec_kw={"height_ratios": [1.0, 1.0]},
    )
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.11, top=0.94, hspace=0.18)

    plasma_colors = plt.get_cmap("plasma")(np.linspace(0.0, 1.0, 4))
    colors = [plasma_colors[i] for i in range(3)]
    linestyles = ["-", "--", ":"]
    hemispheres = [
        ("(a)", "Leading hemisphere", 90.0, axes[0]),
        ("(b)", "Trailing hemisphere", 270.0, axes[1]),
    ]

    for idx, (panel_label, title, nominal_lon_deg, ax) in enumerate(hemispheres):
        lon_idx = int(np.argmin(np.abs(lon_arr - nominal_lon_deg)))
        for color, linestyle, nominal_lat in zip(colors, linestyles, nominal_latitudes_deg):
            lat_idx = _nearest_latitude_index(lat_arr, nominal_lat)
            mean_profile = dose_rate_profiles[lat_idx, lon_idx, :]
            std_profile = (
                dose_profile_std[lat_idx, lon_idx, :]
                if dose_profile_std is not None
                else None
            )
            valid = (
                np.isfinite(depth_centers_cm)
                & np.isfinite(mean_profile)
                & (depth_centers_cm > 0.0)
                & (mean_profile > 0.0)
            )
            if not np.any(valid):
                continue
            actual_lat = float(lat_arr[lat_idx])
            if std_profile is not None:
                lower_profile = np.maximum(
                    mean_profile - std_profile,
                    np.finfo(np.float64).tiny,
                )
                upper_profile = mean_profile + std_profile
                valid_band = (
                    valid
                    & np.isfinite(std_profile)
                    & (std_profile > 0.0)
                    & np.isfinite(lower_profile)
                    & np.isfinite(upper_profile)
                    & (upper_profile > lower_profile)
                )
                if np.any(valid_band):
                    ax.fill_between(
                        depth_centers_cm[valid_band],
                        lower_profile[valid_band],
                        upper_profile[valid_band],
                        color=color,
                        alpha=0.28,
                        linewidth=0.0,
                        zorder=1,
                    )
            ax.plot(
                depth_centers_cm[valid],
                mean_profile[valid],
                color=color,
                linewidth=3.0,
                linestyle=linestyle,
                label=rf"{int(round(actual_lat))}$^\circ$",
            )
        ax.text(
            0.0,
            1.04,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=base_fontsize,
            clip_on=False,
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(x_min_cm, x_max_cm)
        ax.xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=50))
        ax.xaxis.set_minor_locator(
            LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=200)
        )
        ax.xaxis.set_minor_formatter(NullFormatter())
        ymin, ymax = ax.get_ylim()
        ylo = max(min(ymin, ymax), np.finfo(np.float64).tiny)
        yhi = max(ymin, ymax)
        exp_min = int(np.floor(np.log10(ylo)))
        exp_max = int(np.ceil(np.log10(yhi)))
        all_exponents = np.arange(exp_min, exp_max + 1, dtype=int)
        major_exponents = all_exponents[(exp_max - all_exponents) % 2 == 0]
        minor_exponents = np.setdiff1d(all_exponents, major_exponents)
        ax.yaxis.set_major_locator(FixedLocator(np.power(10.0, major_exponents)))
        ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
        if minor_exponents.size > 0:
            ax.yaxis.set_minor_locator(FixedLocator(np.power(10.0, minor_exponents)))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(
            axis="both",
            which="major",
            labelsize=base_fontsize,
            length=8,
            width=1.2,
            direction="in",
            bottom=True,
            top=False,
            left=True,
            right=False,
        )
        ax.tick_params(
            axis="both",
            which="minor",
            labelsize=base_fontsize,
            length=4,
            width=1.0,
            direction="in",
            bottom=True,
            top=False,
            left=True,
            right=False,
        )
        if idx == len(hemispheres) - 1:
            ax.legend(
                loc="lower center",
                bbox_to_anchor=(0.5, 0.015),
                ncol=3,
                frameon=False,
                fontsize=base_fontsize,
                handlelength=2.0,
                columnspacing=1.4,
                borderaxespad=0.0,
            )

    axes[-1].set_xlabel("Depth (cm)", fontsize=base_fontsize)
    fig.supylabel(r"Deposited dose rate (MGy yr$^{-1}$)", x=0.04, fontsize=base_fontsize)
    out_pdf = out_dir / "europa_energy_deposition_profiles_3lat_vstack.pdf"
    out_png = out_dir / "europa_energy_deposition_profiles_3lat_vstack.png"
    png_tmp = out_dir / "europa_energy_deposition_profiles_3lat_vstack_fullscale_tmp.png"
    png_target_width_px = 855
    png_dpi = 300
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(
        png_tmp,
        bbox_inches="tight",
        dpi=png_dpi,
        pil_kwargs={"optimize": True},
    )
    from PIL import Image

    with Image.open(png_tmp) as img:
        scale = png_target_width_px / float(img.size[0])
        scaled_size = (
            png_target_width_px,
            max(1, int(round(img.size[1] * scale))),
        )
        resampling = getattr(Image, "Resampling", Image).LANCZOS
        img.resize(scaled_size, resampling).save(
            out_png,
            dpi=(png_dpi, png_dpi),
            optimize=True,
        )
    png_tmp.unlink(missing_ok=True)
    if show_plots:
        plt.show()


def _infer_legacy_dose_scale(
    depth_edges_cm: np.ndarray,
    dose_profile_map_mgy_per_yr: np.ndarray | None,
    deposited_flux_map_mev_cm2_s: np.ndarray | None,
    density_gcm3: float,
) -> tuple[float, float]:
    if (
        dose_profile_map_mgy_per_yr is None
        or deposited_flux_map_mev_cm2_s is None
        or not np.isfinite(density_gcm3)
        or density_gcm3 <= 0.0
    ):
        return 1.0, float("nan")

    depth_edges = np.asarray(depth_edges_cm, dtype=np.float64)
    profile = np.asarray(dose_profile_map_mgy_per_yr, dtype=np.float64)
    deposited_flux = np.asarray(deposited_flux_map_mev_cm2_s, dtype=np.float64)
    if profile.ndim != 3 or deposited_flux.shape != profile.shape[:2]:
        return 1.0, float("nan")
    if depth_edges.size != profile.shape[2] + 1:
        return 1.0, float("nan")

    dz_cm = depth_edges[1:] - depth_edges[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        reconstructed_flux = np.nansum(
            profile * (density_gcm3 * dz_cm[np.newaxis, np.newaxis, :]) / MEV_TO_MGY_PER_YEAR,
            axis=2,
        )
        ratio = np.divide(
            reconstructed_flux,
            deposited_flux,
            out=np.full_like(reconstructed_flux, np.nan),
            where=np.isfinite(deposited_flux) & (deposited_flux > 0.0),
        )

    valid_ratio = ratio[np.isfinite(ratio) & (ratio > 0.0)]
    if valid_ratio.size < 10:
        return 1.0, float("nan")

    median_ratio = float(np.nanmedian(valid_ratio))
    if 5.0e5 <= median_ratio <= 2.0e6:
        return 1.0e-6, median_ratio
    return 1.0, median_ratio


def cmd_plot_from_npz(args: argparse.Namespace) -> None:
    maps_npz = args.maps_npz.resolve()
    if not maps_npz.exists():
        raise FileNotFoundError(f"Missing NPZ map file: {maps_npz}")

    out_dir = args.out_dir.resolve() if args.out_dir is not None else maps_npz.parent.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    data = np.load(maps_npz)

    lon_values = _npz_get_array(data, ["lon_values"])
    lat_values = _npz_get_array(data, ["lat_values"])
    deposited_flux_map = _npz_get_array(
        data,
        ["deposited_flux_map_mev_cm2_s", "dose_map_depth_integrated_mgy_per_yr", "dose_map_at_depth_mgy_per_yr"],
    )

    primary_pen_key = _npz_first_key(
        data,
        ["primary_penetration_map_cm", "primary_penetration_map_mm"],
    )
    primary_pen_map = np.asarray(data[primary_pen_key], dtype=np.float64)
    if primary_pen_key.endswith("_mm"):
        primary_pen_map = primary_pen_map * 0.1

    secondary_pen_key = _npz_first_key(
        data,
        ["secondary_penetration_map_cm", "secondary_penetration_map_mm"],
    )
    secondary_pen_map = np.asarray(data[secondary_pen_key], dtype=np.float64)
    if secondary_pen_key.endswith("_mm"):
        secondary_pen_map = secondary_pen_map * 0.1

    secondary_ke_map = _npz_get_array(data, ["secondary_ke_map_mev"])
    dep_fraction_map = _npz_get_array(
        data,
        ["deposited_fraction_map", "energy_deposited_fraction_map"],
    )
    backscatter_fraction_map = (
        _npz_get_array(data, ["backscatter_fraction_map"])
        if "backscatter_fraction_map" in data
        else None
    )
    lateral_escape_fraction_map = (
        _npz_get_array(data, ["lateral_escape_fraction_map"])
        if "lateral_escape_fraction_map" in data
        else None
    )
    forward_escape_fraction_map = (
        _npz_get_array(data, ["forward_escape_fraction_map"])
        if "forward_escape_fraction_map" in data
        else None
    )
    escaped_electron_fraction_map = _npz_get_optional_array(
        data,
        ["escaped_electron_fraction_map", "energy_escaped_electron_fraction_map"],
    )
    escaped_photon_fraction_map = _npz_get_optional_array(
        data,
        ["escaped_photon_fraction_map", "energy_escaped_photon_fraction_map"],
    )
    backscatter_electron_fraction_map = _npz_get_optional_array(
        data,
        ["backscatter_electron_fraction_map"],
    )
    backscatter_photon_fraction_map = _npz_get_optional_array(
        data,
        ["backscatter_photon_fraction_map"],
    )
    lateral_escape_electron_fraction_map = _npz_get_optional_array(
        data,
        ["lateral_escape_electron_fraction_map"],
    )
    lateral_escape_photon_fraction_map = _npz_get_optional_array(
        data,
        ["lateral_escape_photon_fraction_map"],
    )
    forward_escape_electron_fraction_map = _npz_get_optional_array(
        data,
        ["forward_escape_electron_fraction_map"],
    )
    forward_escape_photon_fraction_map = _npz_get_optional_array(
        data,
        ["forward_escape_photon_fraction_map"],
    )
    depth_edges_cm = _npz_get_array(data, ["depth_edges_cm"])
    dose_profile_map_mgy_per_yr = _npz_get_array(data, ["dose_profile_map_mgy_per_yr"])
    dose_profile_std_map_mgy_per_yr = (
        _npz_get_array(data, ["dose_profile_std_map_mgy_per_yr"])
        if "dose_profile_std_map_mgy_per_yr" in data
        else None
    )
    cumulative_depth_map = _compute_cumulative_fraction_depth_map(
        depth_edges_cm,
        dose_profile_map_mgy_per_yr,
    )
    primary_flux_map = _npz_get_array(data, ["primary_flux_map_mev_cm2_s"])
    backscatter_amount_map = _npz_get_array(data, ["backscatter_amount_map_mev_cm2_s"])
    if backscatter_fraction_map is None:
        with np.errstate(divide="ignore", invalid="ignore"):
            backscatter_fraction_map = np.divide(
                backscatter_amount_map,
                primary_flux_map,
                out=np.full_like(backscatter_amount_map, np.nan),
                where=primary_flux_map > 0.0,
            )
    energy_panel_flux_map = np.asarray(deposited_flux_map, dtype=np.float64)
    jde_energy_flux_map_mev_cm2_s: np.ndarray | None = None
    if "jde_energy_flux_model_map_mev_cm2_s" in data:
        jde_energy_flux_map_mev_cm2_s = np.asarray(
            data["jde_energy_flux_model_map_mev_cm2_s"], dtype=np.float64
        )
    elif "jde_energy_flux_map_mev_cm2_s" in data:
        jde_energy_flux_map_mev_cm2_s = np.asarray(
            data["jde_energy_flux_map_mev_cm2_s"], dtype=np.float64
        )
    elif "jde_flux_map_cm2_s" in data:
        jde_flux_map = np.asarray(data["jde_flux_map_cm2_s"], dtype=np.float64)
        mean_incident_energy_mev = _npz_get_array(
            data,
            ["mean_incident_energy_map_mev", "incident_energy_map_mev", "primary_mean_energy_map_mev"],
        ) if any(
            k in data for k in ["mean_incident_energy_map_mev", "incident_energy_map_mev", "primary_mean_energy_map_mev"]
        ) else None
        if mean_incident_energy_mev is None and "primary_flux_map_mev_cm2_s" in data:
            primary_flux_map = np.asarray(data["primary_flux_map_mev_cm2_s"], dtype=np.float64)
            with np.errstate(divide="ignore", invalid="ignore"):
                mean_incident_energy_mev = np.divide(
                    primary_flux_map,
                    jde_flux_map,
                    out=np.full_like(primary_flux_map, np.nan),
                    where=jde_flux_map > 0.0,
                )
        if mean_incident_energy_mev is not None:
            jde_energy_flux_map_mev_cm2_s = jde_flux_map * mean_incident_energy_mev

    dose_depth_cm = float(args.dose_depth_cm)
    if not np.isfinite(dose_depth_cm):
        if "dose_depth_cm_selected" in data:
            dose_depth_cm = float(np.asarray(data["dose_depth_cm_selected"]).item())
        elif "dose_depth_cm" in data:
            dose_depth_cm = float(np.asarray(data["dose_depth_cm"]).item())
        elif "dose_depth_mm" in data:
            dose_depth_cm = float(np.asarray(data["dose_depth_mm"]).item()) * 0.1
        else:
            dose_depth_cm = 0.01

    density_gcm3 = _load_density_gcm3_from_manifest(out_dir=out_dir, fallback_dir=maps_npz.parent)
    dose_scale, dose_median_ratio = _infer_legacy_dose_scale(
        depth_edges_cm=depth_edges_cm,
        dose_profile_map_mgy_per_yr=dose_profile_map_mgy_per_yr,
        deposited_flux_map_mev_cm2_s=deposited_flux_map,
        density_gcm3=density_gcm3,
    )
    if dose_scale != 1.0:
        print(
            "WARNING: detected legacy dose units in cached NPZ "
            f"(median profile/flux ratio ~ {dose_median_ratio:.6g}); "
            "rescaling plotted dose quantities by 1e-6 from Gy/yr to MGy/yr."
        )
        dose_profile_map_mgy_per_yr = dose_profile_map_mgy_per_yr * dose_scale
        if dose_profile_std_map_mgy_per_yr is not None:
            dose_profile_std_map_mgy_per_yr = dose_profile_std_map_mgy_per_yr * dose_scale

    has_leading = bool(np.any(lon_values < 180.0))
    has_trailing = bool(np.any(lon_values >= 180.0))
    if not has_leading and not has_trailing:
        raise RuntimeError("NPZ lon_values does not contain valid longitude columns.")

    _configure_plot_style(force_interactive=True)
    hemispheres_written: list[str] = []
    if has_leading:
        _save_hemisphere_plots(
            out_dir=out_dir,
            hemisphere="leading",
            lon_values=lon_values,
            energy_panel_flux_map=energy_panel_flux_map,
            primary_flux_map=primary_flux_map,
            primary_pen_map=primary_pen_map,
            secondary_ke_map=secondary_ke_map,
            secondary_pen_map=secondary_pen_map,
            cumulative_depth_map=cumulative_depth_map,
            dep_fraction_map=dep_fraction_map,
            backscatter_amount_map=backscatter_amount_map,
            dose_depth_cm=dose_depth_cm,
            backscatter_fraction_map=backscatter_fraction_map,
            lateral_escape_fraction_map=lateral_escape_fraction_map,
            forward_escape_fraction_map=forward_escape_fraction_map,
            escaped_electron_fraction_map=escaped_electron_fraction_map,
            escaped_photon_fraction_map=escaped_photon_fraction_map,
            backscatter_electron_fraction_map=backscatter_electron_fraction_map,
            backscatter_photon_fraction_map=backscatter_photon_fraction_map,
            lateral_escape_electron_fraction_map=lateral_escape_electron_fraction_map,
            lateral_escape_photon_fraction_map=lateral_escape_photon_fraction_map,
            forward_escape_electron_fraction_map=forward_escape_electron_fraction_map,
            forward_escape_photon_fraction_map=forward_escape_photon_fraction_map,
            jde_energy_flux_map_mev_cm2_s=jde_energy_flux_map_mev_cm2_s,
            show_plots=True,
        )
        hemispheres_written.append("leading")
    if has_trailing:
        _save_hemisphere_plots(
            out_dir=out_dir,
            hemisphere="trailing",
            lon_values=lon_values,
            energy_panel_flux_map=energy_panel_flux_map,
            primary_flux_map=primary_flux_map,
            primary_pen_map=primary_pen_map,
            secondary_ke_map=secondary_ke_map,
            secondary_pen_map=secondary_pen_map,
            cumulative_depth_map=cumulative_depth_map,
            dep_fraction_map=dep_fraction_map,
            backscatter_amount_map=backscatter_amount_map,
            dose_depth_cm=dose_depth_cm,
            backscatter_fraction_map=backscatter_fraction_map,
            lateral_escape_fraction_map=lateral_escape_fraction_map,
            forward_escape_fraction_map=forward_escape_fraction_map,
            escaped_electron_fraction_map=escaped_electron_fraction_map,
            escaped_photon_fraction_map=escaped_photon_fraction_map,
            backscatter_electron_fraction_map=backscatter_electron_fraction_map,
            backscatter_photon_fraction_map=backscatter_photon_fraction_map,
            lateral_escape_electron_fraction_map=lateral_escape_electron_fraction_map,
            lateral_escape_photon_fraction_map=lateral_escape_photon_fraction_map,
            forward_escape_electron_fraction_map=forward_escape_electron_fraction_map,
            forward_escape_photon_fraction_map=forward_escape_photon_fraction_map,
            jde_energy_flux_map_mev_cm2_s=jde_energy_flux_map_mev_cm2_s,
            show_plots=True,
        )
        hemispheres_written.append("trailing")

    _save_energy_deposition_depth_profiles_plot(
        out_dir=out_dir,
        lat_values=lat_values,
        lon_values=lon_values,
        depth_edges_cm=depth_edges_cm,
        dose_profile_map_mgy_per_yr=dose_profile_map_mgy_per_yr,
        dose_profile_std_map_mgy_per_yr=dose_profile_std_map_mgy_per_yr,
        density_gcm3=density_gcm3,
        show_plots=True,
    )

    import matplotlib
    target_grid = _resolve_target_grid(args)
    if target_grid is not None:
        if jde_energy_flux_map_mev_cm2_s is None:
            raise RuntimeError(
                "Requested J(E)E dE grid resampling but J(E) diagnostic map is unavailable in NPZ."
            )
        nlat, nlon = target_grid
        jde_resampled, lat_resampled, lon_resampled = _resample_2d_map(
            jde_energy_flux_map_mev_cm2_s,
            lat_values,
            lon_values,
            nlat,
            nlon,
        )
        resampled_out = out_dir / f"jde_energy_flux_map_resampled_{nlat}x{nlon}.npz"
        np.savez(
            resampled_out,
            lat_values=np.asarray(lat_resampled, dtype=np.float64),
            lon_values=np.asarray(lon_resampled, dtype=np.float64),
            jde_energy_flux_map_mev_cm2_s=np.asarray(jde_resampled, dtype=np.float64),
            source_maps_npz=str(maps_npz),
        )
        print(f"Wrote resampled J(E)E dE map NPZ: {resampled_out}")

    print(f"Matplotlib backend: {matplotlib.get_backend()}")
    print(f"Regenerated hemisphere plots from NPZ: {maps_npz}")
    if jde_energy_flux_map_mev_cm2_s is None:
        print("NOTE: J(E) diagnostic map was unavailable; generated all non-J(E) plots.")
    expected: list[Path] = []
    has_escape_fraction_maps = (
        lateral_escape_fraction_map is not None
        and forward_escape_fraction_map is not None
        and np.any(np.isfinite(lateral_escape_fraction_map))
        and np.any(np.isfinite(forward_escape_fraction_map))
    )
    has_electron_escape_fraction_maps = all(
        m is not None
        for m in [
            escaped_electron_fraction_map,
            backscatter_electron_fraction_map,
            lateral_escape_electron_fraction_map,
            forward_escape_electron_fraction_map,
        ]
    )
    has_photon_escape_fraction_maps = all(
        m is not None
        for m in [
            escaped_photon_fraction_map,
            backscatter_photon_fraction_map,
            lateral_escape_photon_fraction_map,
            forward_escape_photon_fraction_map,
        ]
    )
    for hemi in hemispheres_written:
        expected.append(out_dir / f"europa_{hemi}_4panel_metrics.png")
        if has_escape_fraction_maps:
            expected.append(out_dir / f"europa_{hemi}_4panel_escape_fractions.png")
        if has_electron_escape_fraction_maps:
            expected.append(out_dir / f"europa_{hemi}_electron_4panel_escape_fractions.png")
        if has_photon_escape_fraction_maps:
            expected.append(out_dir / f"europa_{hemi}_photon_4panel_escape_fractions.png")
        if has_electron_escape_fraction_maps and has_photon_escape_fraction_maps:
            expected.append(out_dir / f"europa_{hemi}_electron_photon_8panel_escape_fractions.png")
        if jde_energy_flux_map_mev_cm2_s is not None:
            expected.append(out_dir / f"europa_{hemi}_jde_diagnostic.pdf")
    if not has_escape_fraction_maps:
        print("NOTE: forward/lateral escape fraction maps were unavailable; skipped escape-fraction 4-panel plots.")
    expected.append(out_dir / "europa_energy_deposition_profiles_3lat_vstack.pdf")
    expected.append(out_dir / "europa_energy_deposition_profiles_3lat_vstack.png")
    missing = [p for p in expected if not p.exists()]
    if missing:
        raise RuntimeError(
            "Plot generation completed but some expected files are missing: "
            + ", ".join(str(p) for p in missing)
        )
    for p in expected:
        print(f"Wrote plot: {p}")
    print(f"Wrote plots under: {out_dir}")


def cmd_merge(args: argparse.Namespace) -> None:
    manifest = _load_manifest(args.manifest)
    manifest_file = args.manifest.resolve()
    out_dir = _resolve_manifest_cached_path(
        manifest_file,
        manifest.get("out_dir"),
        ".",
        expect_dir=True,
    )
    range_results_dir = _resolve_manifest_cached_path(
        manifest_file,
        manifest.get("range_results_dir"),
        "range_results",
        expect_dir=True,
    )
    depth_edges_file = _resolve_manifest_cached_path(
        manifest_file,
        manifest.get("depth_edges_file"),
        "depth_edges_cm.npy",
        expect_dir=False,
    )
    depth_edges_cm = np.load(depth_edges_file)
    depth_centers_cm = 0.5 * (depth_edges_cm[:-1] + depth_edges_cm[1:])

    ranges = manifest.get("unique_ranges", [])
    cells = manifest.get("cells", [])
    if not ranges:
        raise RuntimeError("Manifest has no unique ranges.")
    if not cells:
        raise RuntimeError("Manifest has no cells.")

    missing = []
    range_results: dict[int, dict[str, Any]] = {}
    for r in ranges:
        rid = int(r["range_id"])
        f = _range_result_path(range_results_dir, rid)
        if not f.exists():
            missing.append(str(f))
            continue
        range_results[rid] = _load_range_result(f)
    if missing:
        preview = "\n".join(missing[:10])
        more = "" if len(missing) <= 10 else f"\n... and {len(missing)-10} more"
        raise FileNotFoundError(
            "Missing range result files. Run all range workers first.\n"
            + preview
            + more
        )

    max_cell_id = max(int(c["cell_id"]) for c in cells)
    n_cells = max_cell_id + 1
    n_depth = depth_centers_cm.size
    density_gcm3 = _safe_float(str(manifest.get("density_gcm3_default", "nan")))

    requested_depth_cm = float(args.dose_depth_cm)
    if not np.isfinite(requested_depth_cm):
        any_range = next(iter(range_results.values()))
        requested_depth_cm = float(any_range["dose_depth_cm"])
    dose_idx = int(np.argmin(np.abs(depth_centers_cm - requested_depth_cm)))
    dose_depth_cm = float(depth_centers_cm[dose_idx])

    stale_escape_ranges = [
        int(r["range_id"])
        for r in range_results.values()
        if not np.isfinite(float(r.get("forward_escape_fraction", float("nan"))))
        or not np.isfinite(float(r.get("lateral_escape_fraction", float("nan"))))
        or not np.isfinite(float(r.get("escaped_electron_fraction", float("nan"))))
        or not np.isfinite(float(r.get("escaped_photon_fraction", float("nan"))))
    ]
    if stale_escape_ranges:
        energy_metrics_csv = Path(
            _resolve_manifest_cached_path(
                manifest_file,
                manifest.get("energy_metrics_csv"),
                "energy_metrics.csv",
                expect_dir=False,
            )
        ).resolve()
        energy_profiles_npz = Path(
            _resolve_manifest_cached_path(
                manifest_file,
                manifest.get("energy_profiles_npz"),
                "energy_edep_profiles.npz",
                expect_dir=False,
            )
        ).resolve()
        metrics_by_idx, profiles_by_idx, profile_sumsq_by_idx, cache_depth_edges_cm = _load_energy_cache(
            energy_metrics_csv,
            energy_profiles_npz,
        )
        if not np.allclose(
            np.asarray(cache_depth_edges_cm, dtype=np.float64),
            np.asarray(depth_edges_cm, dtype=np.float64),
            rtol=0.0,
            atol=1.0e-12,
            equal_nan=True,
        ):
            raise RuntimeError(
                "Energy cache depth edges do not match merge manifest depth edges."
            )
        range_defs = {
            int(r["range_id"]): r
            for r in manifest.get("unique_ranges", [])
        }
        refreshed = 0
        for rid in stale_escape_ranges:
            range_def = range_defs.get(int(rid))
            if range_def is None:
                continue
            refreshed_result = _compute_range_result(
                range_def=range_def,
                metrics_by_idx=metrics_by_idx,
                profiles_by_idx=profiles_by_idx,
                profile_sumsq_by_idx=profile_sumsq_by_idx,
                depth_edges_cm=depth_edges_cm,
                density_gcm3=density_gcm3,
                dose_depth_cm=dose_depth_cm,
            )
            range_results[int(rid)] = refreshed_result
            _save_range_result(_range_result_path(range_results_dir, int(rid)), refreshed_result)
            refreshed += 1
        if refreshed > 0:
            print(
                "Refreshed directional and particle-split escape fractions for "
                f"{refreshed} stale range-result files from the existing energy cache."
            )

    cell_dose_profile = np.full((n_cells, n_depth), np.nan, dtype=np.float64)
    cell_dose_profile_std = np.full((n_cells, n_depth), np.nan, dtype=np.float64)
    cell_dose_at_depth = np.full(n_cells, np.nan, dtype=np.float64)
    cell_primary_pen = np.full(n_cells, np.nan, dtype=np.float64)
    cell_secondary_ke = np.full(n_cells, np.nan, dtype=np.float64)
    cell_secondary_pen = np.full(n_cells, np.nan, dtype=np.float64)
    cell_deposited_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_backscatter_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_forward_escape_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_lateral_escape_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_escaped_electron_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_escaped_photon_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_backscatter_electron_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_backscatter_photon_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_forward_escape_electron_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_forward_escape_photon_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_lateral_escape_electron_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_lateral_escape_photon_fraction = np.full(n_cells, np.nan, dtype=np.float64)
    cell_backscatter_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_forward_escape_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_lateral_escape_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_escaped_electron_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_escaped_photon_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_backscatter_electron_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_backscatter_photon_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_forward_escape_electron_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_forward_escape_photon_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_lateral_escape_electron_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_lateral_escape_photon_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_primary_flux = np.full(n_cells, np.nan, dtype=np.float64)
    cell_deposited_flux = np.full(n_cells, np.nan, dtype=np.float64)
    cell_jde_flux = np.full(n_cells, np.nan, dtype=np.float64)
    cell_jde_energy_flux_model = np.full(n_cells, np.nan, dtype=np.float64)
    cell_deflected_amount = np.full(n_cells, np.nan, dtype=np.float64)
    cell_deflected_fraction = np.full(n_cells, np.nan, dtype=np.float64)

    lat_values = sorted({float(c["lat_deg"]) for c in cells}, reverse=True)
    lon_values = sorted({float(c["lon_w_deg"]) for c in cells})
    lat_index = {lat: i for i, lat in enumerate(lat_values)}
    lon_index = {lon: j for j, lon in enumerate(lon_values)}

    shape = (len(lat_values), len(lon_values))
    dose_map = np.full(shape, np.nan, dtype=np.float64)
    dose_map_depth_integrated = np.full(shape, np.nan, dtype=np.float64)
    jde_dose_equiv_map = np.full(shape, np.nan, dtype=np.float64)
    jde_energy_flux_map = np.full(shape, np.nan, dtype=np.float64)
    primary_pen_map = np.full(shape, np.nan, dtype=np.float64)
    secondary_ke_map = np.full(shape, np.nan, dtype=np.float64)
    secondary_pen_map = np.full(shape, np.nan, dtype=np.float64)
    dep_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    backscatter_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    forward_escape_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    lateral_escape_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    escaped_electron_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    escaped_photon_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    backscatter_electron_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    backscatter_photon_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    forward_escape_electron_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    forward_escape_photon_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    lateral_escape_electron_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    lateral_escape_photon_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    backscatter_amount_map = np.full(shape, np.nan, dtype=np.float64)
    forward_escape_amount_map = np.full(shape, np.nan, dtype=np.float64)
    lateral_escape_amount_map = np.full(shape, np.nan, dtype=np.float64)
    escaped_electron_amount_map = np.full(shape, np.nan, dtype=np.float64)
    escaped_photon_amount_map = np.full(shape, np.nan, dtype=np.float64)
    backscatter_electron_amount_map = np.full(shape, np.nan, dtype=np.float64)
    backscatter_photon_amount_map = np.full(shape, np.nan, dtype=np.float64)
    forward_escape_electron_amount_map = np.full(shape, np.nan, dtype=np.float64)
    forward_escape_photon_amount_map = np.full(shape, np.nan, dtype=np.float64)
    lateral_escape_electron_amount_map = np.full(shape, np.nan, dtype=np.float64)
    lateral_escape_photon_amount_map = np.full(shape, np.nan, dtype=np.float64)
    primary_flux_map = np.full(shape, np.nan, dtype=np.float64)
    deposited_flux_map = np.full(shape, np.nan, dtype=np.float64)
    jde_flux_map = np.full(shape, np.nan, dtype=np.float64)
    jde_energy_flux_model_map = np.full(shape, np.nan, dtype=np.float64)
    deflected_fraction_map = np.full(shape, np.nan, dtype=np.float64)
    deflected_amount_map = np.full(shape, np.nan, dtype=np.float64)
    dose_profile_map = np.full((shape[0], shape[1], n_depth), np.nan, dtype=np.float64)
    dose_profile_std_map = np.full((shape[0], shape[1], n_depth), np.nan, dtype=np.float64)

    for c in cells:
        cid = int(c["cell_id"])
        rid = int(c["range_id"])
        lat = float(c["lat_deg"])
        lon = float(c["lon_w_deg"])
        ii = lat_index[lat]
        jj = lon_index[lon]
        if rid < 0:
            continue
        res = range_results.get(rid)
        if res is None:
            continue

        profile = np.asarray(res["dose_profile_mgy_per_yr"], dtype=np.float64)
        profile_std = np.asarray(res["dose_profile_std_mgy_per_yr"], dtype=np.float64)
        cell_dose_profile[cid, :] = profile
        cell_dose_profile_std[cid, :] = profile_std
        if profile.size > dose_idx and np.isfinite(profile[dose_idx]):
            cell_dose_at_depth[cid] = float(profile[dose_idx])
        else:
            cell_dose_at_depth[cid] = float("nan")
        cell_primary_pen[cid] = float(res["mean_primary_penetration_cm"])
        cell_secondary_ke[cid] = float(res["mean_secondary_ke_mev"])
        cell_secondary_pen[cid] = float(res["mean_secondary_penetration_cm"])
        cell_deposited_fraction[cid] = float(res["deposited_fraction"])
        cell_backscatter_fraction[cid] = float(res["backscatter_fraction"])
        cell_forward_escape_fraction[cid] = float(res["forward_escape_fraction"])
        cell_lateral_escape_fraction[cid] = float(res["lateral_escape_fraction"])
        cell_escaped_electron_fraction[cid] = float(res.get("escaped_electron_fraction", float("nan")))
        cell_escaped_photon_fraction[cid] = float(res.get("escaped_photon_fraction", float("nan")))
        cell_backscatter_electron_fraction[cid] = float(res.get("backscatter_electron_fraction", float("nan")))
        cell_backscatter_photon_fraction[cid] = float(res.get("backscatter_photon_fraction", float("nan")))
        cell_forward_escape_electron_fraction[cid] = float(
            res.get("escaped_forward_electron_fraction", float("nan"))
        )
        cell_forward_escape_photon_fraction[cid] = float(
            res.get("escaped_forward_photon_fraction", float("nan"))
        )
        cell_lateral_escape_electron_fraction[cid] = float(
            res.get("escaped_lateral_electron_fraction", float("nan"))
        )
        cell_lateral_escape_photon_fraction[cid] = float(
            res.get("escaped_lateral_photon_fraction", float("nan"))
        )
        cell_backscatter_amount[cid] = float(res["backscatter_flux_mev_cm2_s"])
        cell_forward_escape_amount[cid] = float(res["forward_escape_flux_mev_cm2_s"])
        cell_lateral_escape_amount[cid] = float(res["lateral_escape_flux_mev_cm2_s"])
        cell_escaped_electron_amount[cid] = float(res.get("escaped_electron_flux_mev_cm2_s", float("nan")))
        cell_escaped_photon_amount[cid] = float(res.get("escaped_photon_flux_mev_cm2_s", float("nan")))
        cell_backscatter_electron_amount[cid] = float(
            res.get("backscatter_electron_flux_mev_cm2_s", float("nan"))
        )
        cell_backscatter_photon_amount[cid] = float(
            res.get("backscatter_photon_flux_mev_cm2_s", float("nan"))
        )
        cell_forward_escape_electron_amount[cid] = float(
            res.get("escaped_forward_electron_flux_mev_cm2_s", float("nan"))
        )
        cell_forward_escape_photon_amount[cid] = float(
            res.get("escaped_forward_photon_flux_mev_cm2_s", float("nan"))
        )
        cell_lateral_escape_electron_amount[cid] = float(
            res.get("escaped_lateral_electron_flux_mev_cm2_s", float("nan"))
        )
        cell_lateral_escape_photon_amount[cid] = float(
            res.get("escaped_lateral_photon_flux_mev_cm2_s", float("nan"))
        )
        cell_primary_flux[cid] = float(res["primary_flux_mev_cm2_s"])
        cell_deposited_flux[cid] = float(res["deposited_flux_mev_cm2_s"])
        cell_jde_flux[cid] = float(res["jde_flux_cm2_s"])
        cell_jde_energy_flux_model[cid] = float(res["jde_energy_flux_model_mev_cm2_s"])
        if np.isfinite(cell_primary_flux[cid]) and np.isfinite(cell_deposited_flux[cid]):
            cell_deflected_amount[cid] = cell_primary_flux[cid] - cell_deposited_flux[cid]
        else:
            cell_deflected_amount[cid] = float("nan")
        if np.isfinite(cell_primary_flux[cid]) and cell_primary_flux[cid] > 0.0 and np.isfinite(cell_deflected_amount[cid]):
            cell_deflected_fraction[cid] = cell_deflected_amount[cid] / cell_primary_flux[cid]
        else:
            cell_deflected_fraction[cid] = float("nan")

        dose_map[ii, jj] = cell_dose_at_depth[cid]
        primary_pen_map[ii, jj] = cell_primary_pen[cid]
        secondary_ke_map[ii, jj] = cell_secondary_ke[cid]
        secondary_pen_map[ii, jj] = cell_secondary_pen[cid]
        dep_fraction_map[ii, jj] = cell_deposited_fraction[cid]
        backscatter_fraction_map[ii, jj] = cell_backscatter_fraction[cid]
        forward_escape_fraction_map[ii, jj] = cell_forward_escape_fraction[cid]
        lateral_escape_fraction_map[ii, jj] = cell_lateral_escape_fraction[cid]
        escaped_electron_fraction_map[ii, jj] = cell_escaped_electron_fraction[cid]
        escaped_photon_fraction_map[ii, jj] = cell_escaped_photon_fraction[cid]
        backscatter_electron_fraction_map[ii, jj] = cell_backscatter_electron_fraction[cid]
        backscatter_photon_fraction_map[ii, jj] = cell_backscatter_photon_fraction[cid]
        forward_escape_electron_fraction_map[ii, jj] = cell_forward_escape_electron_fraction[cid]
        forward_escape_photon_fraction_map[ii, jj] = cell_forward_escape_photon_fraction[cid]
        lateral_escape_electron_fraction_map[ii, jj] = cell_lateral_escape_electron_fraction[cid]
        lateral_escape_photon_fraction_map[ii, jj] = cell_lateral_escape_photon_fraction[cid]
        backscatter_amount_map[ii, jj] = cell_backscatter_amount[cid]
        forward_escape_amount_map[ii, jj] = cell_forward_escape_amount[cid]
        lateral_escape_amount_map[ii, jj] = cell_lateral_escape_amount[cid]
        escaped_electron_amount_map[ii, jj] = cell_escaped_electron_amount[cid]
        escaped_photon_amount_map[ii, jj] = cell_escaped_photon_amount[cid]
        backscatter_electron_amount_map[ii, jj] = cell_backscatter_electron_amount[cid]
        backscatter_photon_amount_map[ii, jj] = cell_backscatter_photon_amount[cid]
        forward_escape_electron_amount_map[ii, jj] = cell_forward_escape_electron_amount[cid]
        forward_escape_photon_amount_map[ii, jj] = cell_forward_escape_photon_amount[cid]
        lateral_escape_electron_amount_map[ii, jj] = cell_lateral_escape_electron_amount[cid]
        lateral_escape_photon_amount_map[ii, jj] = cell_lateral_escape_photon_amount[cid]
        primary_flux_map[ii, jj] = cell_primary_flux[cid]
        deposited_flux_map[ii, jj] = cell_deposited_flux[cid]
        jde_flux_map[ii, jj] = cell_jde_flux[cid]
        jde_energy_flux_model_map[ii, jj] = cell_jde_energy_flux_model[cid]
        deflected_fraction_map[ii, jj] = cell_deflected_fraction[cid]
        deflected_amount_map[ii, jj] = cell_deflected_amount[cid]
        dose_profile_map[ii, jj, :] = profile
        dose_profile_std_map[ii, jj, :] = profile_std

    with np.errstate(divide="ignore", invalid="ignore"):
        cell_nondeposited_escaped_electron_fraction = np.divide(
            cell_escaped_electron_amount,
            cell_deflected_amount,
            out=np.full_like(cell_escaped_electron_amount, np.nan),
            where=cell_deflected_amount > 0.0,
        )
        cell_nondeposited_escaped_photon_fraction = np.divide(
            cell_escaped_photon_amount,
            cell_deflected_amount,
            out=np.full_like(cell_escaped_photon_amount, np.nan),
            where=cell_deflected_amount > 0.0,
        )
        cell_nondeposited_unresolved_fraction = np.divide(
            cell_deflected_amount - cell_escaped_electron_amount - cell_escaped_photon_amount,
            cell_deflected_amount,
            out=np.full_like(cell_deflected_amount, np.nan),
            where=cell_deflected_amount > 0.0,
        )
        nondeposited_escaped_electron_fraction_map = np.divide(
            escaped_electron_amount_map,
            deflected_amount_map,
            out=np.full_like(escaped_electron_amount_map, np.nan),
            where=deflected_amount_map > 0.0,
        )
        nondeposited_escaped_photon_fraction_map = np.divide(
            escaped_photon_amount_map,
            deflected_amount_map,
            out=np.full_like(escaped_photon_amount_map, np.nan),
            where=deflected_amount_map > 0.0,
        )
        nondeposited_unresolved_fraction_map = np.divide(
            deflected_amount_map - escaped_electron_amount_map - escaped_photon_amount_map,
            deflected_amount_map,
            out=np.full_like(deflected_amount_map, np.nan),
            where=deflected_amount_map > 0.0,
        )

    dose_scale, dose_median_ratio = _infer_legacy_dose_scale(
        depth_edges_cm=depth_edges_cm,
        dose_profile_map_mgy_per_yr=dose_profile_map,
        deposited_flux_map_mev_cm2_s=deposited_flux_map,
        density_gcm3=density_gcm3,
    )
    if dose_scale != 1.0:
        print(
            "WARNING: detected legacy dose units in range results "
            f"(median profile/flux ratio ~ {dose_median_ratio:.6g}); "
            "rescaling dose products by 1e-6 from Gy/yr to MGy/yr."
        )
        cell_dose_profile *= dose_scale
        cell_dose_profile_std *= dose_scale
        cell_dose_at_depth *= dose_scale
        dose_map *= dose_scale
        dose_profile_map *= dose_scale
        dose_profile_std_map *= dose_scale
        for res in range_results.values():
            if "dose_at_depth_mgy_per_yr" in res and np.isfinite(res["dose_at_depth_mgy_per_yr"]):
                res["dose_at_depth_mgy_per_yr"] = float(res["dose_at_depth_mgy_per_yr"]) * dose_scale
            if "dose_profile_mgy_per_yr" in res:
                res["dose_profile_mgy_per_yr"] = (
                    np.asarray(res["dose_profile_mgy_per_yr"], dtype=np.float64) * dose_scale
                )
            if "dose_profile_std_mgy_per_yr" in res:
                res["dose_profile_std_mgy_per_yr"] = (
                    np.asarray(res["dose_profile_std_mgy_per_yr"], dtype=np.float64) * dose_scale
                )

    slab_thickness_cm = float(depth_edges_cm[-1] - depth_edges_cm[0]) if depth_edges_cm.size >= 2 else float("nan")
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_incident_energy_mev = np.divide(
            primary_flux_map,
            jde_flux_map,
            out=np.full_like(primary_flux_map, np.nan),
            where=jde_flux_map > 0.0,
        )
    jde_energy_flux_map = np.array(jde_energy_flux_model_map, copy=True)
    missing_model = ~np.isfinite(jde_energy_flux_map)
    if np.any(missing_model):
        jde_energy_flux_fallback = jde_flux_map * mean_incident_energy_mev
        jde_energy_flux_map[missing_model] = jde_energy_flux_fallback[missing_model]

    if (
        np.isfinite(density_gcm3)
        and density_gcm3 > 0.0
        and np.isfinite(slab_thickness_cm)
        and slab_thickness_cm > 0.0
    ):
        with np.errstate(divide="ignore", invalid="ignore"):
            dose_map_depth_integrated = np.divide(
                deposited_flux_map,
                density_gcm3 * slab_thickness_cm,
                out=np.full_like(deposited_flux_map, np.nan),
            )
        dose_map_depth_integrated = dose_map_depth_integrated * MEV_TO_MGY_PER_YEAR
        # Convert the diagnostic J(E)dE-derived energy flux to a dose-equivalent for optional downstream use.
        with np.errstate(divide="ignore", invalid="ignore"):
            jde_dose_equiv_map = np.divide(
                jde_energy_flux_map,
                density_gcm3 * slab_thickness_cm,
                out=np.full_like(jde_energy_flux_map, np.nan),
            )
        jde_dose_equiv_map = jde_dose_equiv_map * MEV_TO_MGY_PER_YEAR

    cumulative_depth_map = _compute_cumulative_fraction_depth_map(
        depth_edges_cm,
        dose_profile_map,
    )

    summary_csv = out_dir / "latlon_metrics_summary.csv"
    with summary_csv.open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(
            [
                "cell_id",
                "lat_deg",
                "lon_w_deg",
                "hemisphere",
                "range_id",
                "dose_at_depth_mgy_per_yr",
                "mean_primary_penetration_cm",
                "mean_secondary_ke_mev",
                "mean_secondary_penetration_cm",
                "deposited_fraction",
                "backscatter_fraction",
                "lateral_escape_fraction",
                "forward_escape_fraction",
                "escaped_electron_fraction",
                "escaped_photon_fraction",
                "nondeposited_escaped_electron_fraction",
                "nondeposited_escaped_photon_fraction",
                "nondeposited_unresolved_fraction",
                "backscatter_electron_fraction",
                "backscatter_photon_fraction",
                "lateral_escape_electron_fraction",
                "lateral_escape_photon_fraction",
                "backscatter_flux_mev_cm2_s",
                "lateral_escape_flux_mev_cm2_s",
                "forward_escape_flux_mev_cm2_s",
                "escaped_electron_flux_mev_cm2_s",
                "escaped_photon_flux_mev_cm2_s",
                "backscatter_electron_flux_mev_cm2_s",
                "backscatter_photon_flux_mev_cm2_s",
                "lateral_escape_electron_flux_mev_cm2_s",
                "lateral_escape_photon_flux_mev_cm2_s",
                "jde_flux_cm2_s",
                "jde_energy_flux_model_mev_cm2_s",
                "deflected_fraction",
                "deflected_flux_mev_cm2_s",
            ]
        )
        for c in sorted(cells, key=lambda x: int(x["cell_id"])):
            cid = int(c["cell_id"])
            wr.writerow(
                [
                    cid,
                    f"{float(c['lat_deg']):.9g}",
                    f"{float(c['lon_w_deg']):.9g}",
                    c["hemisphere"],
                    int(c["range_id"]),
                    f"{cell_dose_at_depth[cid]:.9g}" if np.isfinite(cell_dose_at_depth[cid]) else "nan",
                    f"{cell_primary_pen[cid]:.9g}" if np.isfinite(cell_primary_pen[cid]) else "nan",
                    f"{cell_secondary_ke[cid]:.9g}" if np.isfinite(cell_secondary_ke[cid]) else "nan",
                    f"{cell_secondary_pen[cid]:.9g}" if np.isfinite(cell_secondary_pen[cid]) else "nan",
                    f"{cell_deposited_fraction[cid]:.9g}" if np.isfinite(cell_deposited_fraction[cid]) else "nan",
                    f"{cell_backscatter_fraction[cid]:.9g}" if np.isfinite(cell_backscatter_fraction[cid]) else "nan",
                    f"{cell_lateral_escape_fraction[cid]:.9g}" if np.isfinite(cell_lateral_escape_fraction[cid]) else "nan",
                    f"{cell_forward_escape_fraction[cid]:.9g}" if np.isfinite(cell_forward_escape_fraction[cid]) else "nan",
                    f"{cell_escaped_electron_fraction[cid]:.9g}" if np.isfinite(cell_escaped_electron_fraction[cid]) else "nan",
                    f"{cell_escaped_photon_fraction[cid]:.9g}" if np.isfinite(cell_escaped_photon_fraction[cid]) else "nan",
                    f"{cell_nondeposited_escaped_electron_fraction[cid]:.9g}"
                    if np.isfinite(cell_nondeposited_escaped_electron_fraction[cid])
                    else "nan",
                    f"{cell_nondeposited_escaped_photon_fraction[cid]:.9g}"
                    if np.isfinite(cell_nondeposited_escaped_photon_fraction[cid])
                    else "nan",
                    f"{cell_nondeposited_unresolved_fraction[cid]:.9g}"
                    if np.isfinite(cell_nondeposited_unresolved_fraction[cid])
                    else "nan",
                    f"{cell_backscatter_electron_fraction[cid]:.9g}" if np.isfinite(cell_backscatter_electron_fraction[cid]) else "nan",
                    f"{cell_backscatter_photon_fraction[cid]:.9g}" if np.isfinite(cell_backscatter_photon_fraction[cid]) else "nan",
                    f"{cell_lateral_escape_electron_fraction[cid]:.9g}"
                    if np.isfinite(cell_lateral_escape_electron_fraction[cid])
                    else "nan",
                    f"{cell_lateral_escape_photon_fraction[cid]:.9g}"
                    if np.isfinite(cell_lateral_escape_photon_fraction[cid])
                    else "nan",
                    f"{cell_backscatter_amount[cid]:.9g}" if np.isfinite(cell_backscatter_amount[cid]) else "nan",
                    f"{cell_lateral_escape_amount[cid]:.9g}" if np.isfinite(cell_lateral_escape_amount[cid]) else "nan",
                    f"{cell_forward_escape_amount[cid]:.9g}" if np.isfinite(cell_forward_escape_amount[cid]) else "nan",
                    f"{cell_escaped_electron_amount[cid]:.9g}" if np.isfinite(cell_escaped_electron_amount[cid]) else "nan",
                    f"{cell_escaped_photon_amount[cid]:.9g}" if np.isfinite(cell_escaped_photon_amount[cid]) else "nan",
                    f"{cell_backscatter_electron_amount[cid]:.9g}"
                    if np.isfinite(cell_backscatter_electron_amount[cid])
                    else "nan",
                    f"{cell_backscatter_photon_amount[cid]:.9g}"
                    if np.isfinite(cell_backscatter_photon_amount[cid])
                    else "nan",
                    f"{cell_lateral_escape_electron_amount[cid]:.9g}"
                    if np.isfinite(cell_lateral_escape_electron_amount[cid])
                    else "nan",
                    f"{cell_lateral_escape_photon_amount[cid]:.9g}"
                    if np.isfinite(cell_lateral_escape_photon_amount[cid])
                    else "nan",
                    f"{cell_jde_flux[cid]:.9g}" if np.isfinite(cell_jde_flux[cid]) else "nan",
                    f"{cell_jde_energy_flux_model[cid]:.9g}" if np.isfinite(cell_jde_energy_flux_model[cid]) else "nan",
                    f"{cell_deflected_fraction[cid]:.9g}" if np.isfinite(cell_deflected_fraction[cid]) else "nan",
                    f"{cell_deflected_amount[cid]:.9g}" if np.isfinite(cell_deflected_amount[cid]) else "nan",
                ]
            )

    range_metrics_csv = out_dir / "unique_range_metrics.csv"
    with range_metrics_csv.open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(
            [
                "range_id",
                "e_min_mev",
                "e_max_mev",
                "n_cells",
                "used_energies",
                "missing_energies",
                "dose_depth_cm",
                "dose_at_depth_mgy_per_yr",
                "mean_primary_penetration_cm",
                "mean_secondary_ke_mev",
                "mean_secondary_penetration_cm",
                "deposited_fraction",
                "backscatter_fraction",
                "lateral_escape_fraction",
                "forward_escape_fraction",
                "escaped_electron_fraction",
                "escaped_photon_fraction",
                "nondeposited_escaped_electron_fraction",
                "nondeposited_escaped_photon_fraction",
                "nondeposited_unresolved_fraction",
                "backscatter_electron_fraction",
                "backscatter_photon_fraction",
                "lateral_escape_electron_fraction",
                "lateral_escape_photon_fraction",
                "backscatter_flux_mev_cm2_s",
                "lateral_escape_flux_mev_cm2_s",
                "forward_escape_flux_mev_cm2_s",
                "escaped_electron_flux_mev_cm2_s",
                "escaped_photon_flux_mev_cm2_s",
                "backscatter_electron_flux_mev_cm2_s",
                "backscatter_photon_flux_mev_cm2_s",
                "lateral_escape_electron_flux_mev_cm2_s",
                "lateral_escape_photon_flux_mev_cm2_s",
                "jde_flux_cm2_s",
                "jde_energy_flux_model_mev_cm2_s",
                "deflected_fraction",
                "deflected_flux_mev_cm2_s",
            ]
        )
        for rid in sorted(range_results.keys()):
            r = range_results[rid]
            primary_flux = float(r["primary_flux_mev_cm2_s"])
            deposited_flux = float(r["deposited_flux_mev_cm2_s"])
            deflected_flux = (
                primary_flux - deposited_flux
                if np.isfinite(primary_flux) and np.isfinite(deposited_flux)
                else float("nan")
            )
            deflected_frac = (
                deflected_flux / primary_flux
                if np.isfinite(deflected_flux) and np.isfinite(primary_flux) and primary_flux > 0.0
                else float("nan")
            )
            escaped_electron_flux = float(r["escaped_electron_flux_mev_cm2_s"])
            escaped_photon_flux = float(r["escaped_photon_flux_mev_cm2_s"])
            nondep_electron_frac = (
                escaped_electron_flux / deflected_flux
                if np.isfinite(escaped_electron_flux) and np.isfinite(deflected_flux) and deflected_flux > 0.0
                else float("nan")
            )
            nondep_photon_frac = (
                escaped_photon_flux / deflected_flux
                if np.isfinite(escaped_photon_flux) and np.isfinite(deflected_flux) and deflected_flux > 0.0
                else float("nan")
            )
            nondep_unresolved_frac = (
                (deflected_flux - escaped_electron_flux - escaped_photon_flux) / deflected_flux
                if np.isfinite(deflected_flux)
                and np.isfinite(escaped_electron_flux)
                and np.isfinite(escaped_photon_flux)
                and deflected_flux > 0.0
                else float("nan")
            )
            wr.writerow(
                [
                    rid,
                    f"{r['e_min_mev']:.9g}",
                    f"{r['e_max_mev']:.9g}",
                    int(r["n_cells"]),
                    int(r["used_energies"]),
                    int(r["missing_energies"]),
                    f"{r['dose_depth_cm']:.9g}",
                    f"{r['dose_at_depth_mgy_per_yr']:.9g}" if np.isfinite(r["dose_at_depth_mgy_per_yr"]) else "nan",
                    f"{r['mean_primary_penetration_cm']:.9g}" if np.isfinite(r["mean_primary_penetration_cm"]) else "nan",
                    f"{r['mean_secondary_ke_mev']:.9g}" if np.isfinite(r["mean_secondary_ke_mev"]) else "nan",
                    f"{r['mean_secondary_penetration_cm']:.9g}" if np.isfinite(r["mean_secondary_penetration_cm"]) else "nan",
                    f"{r['deposited_fraction']:.9g}" if np.isfinite(r["deposited_fraction"]) else "nan",
                    f"{r['backscatter_fraction']:.9g}" if np.isfinite(r["backscatter_fraction"]) else "nan",
                    f"{r['lateral_escape_fraction']:.9g}" if np.isfinite(r["lateral_escape_fraction"]) else "nan",
                    f"{r['forward_escape_fraction']:.9g}" if np.isfinite(r["forward_escape_fraction"]) else "nan",
                    f"{r['escaped_electron_fraction']:.9g}" if np.isfinite(r["escaped_electron_fraction"]) else "nan",
                    f"{r['escaped_photon_fraction']:.9g}" if np.isfinite(r["escaped_photon_fraction"]) else "nan",
                    f"{nondep_electron_frac:.9g}" if np.isfinite(nondep_electron_frac) else "nan",
                    f"{nondep_photon_frac:.9g}" if np.isfinite(nondep_photon_frac) else "nan",
                    f"{nondep_unresolved_frac:.9g}" if np.isfinite(nondep_unresolved_frac) else "nan",
                    f"{r['backscatter_electron_fraction']:.9g}" if np.isfinite(r["backscatter_electron_fraction"]) else "nan",
                    f"{r['backscatter_photon_fraction']:.9g}" if np.isfinite(r["backscatter_photon_fraction"]) else "nan",
                    f"{r['escaped_lateral_electron_fraction']:.9g}"
                    if np.isfinite(r["escaped_lateral_electron_fraction"])
                    else "nan",
                    f"{r['escaped_lateral_photon_fraction']:.9g}"
                    if np.isfinite(r["escaped_lateral_photon_fraction"])
                    else "nan",
                    f"{r['backscatter_flux_mev_cm2_s']:.9g}" if np.isfinite(r["backscatter_flux_mev_cm2_s"]) else "nan",
                    f"{r['lateral_escape_flux_mev_cm2_s']:.9g}" if np.isfinite(r["lateral_escape_flux_mev_cm2_s"]) else "nan",
                    f"{r['forward_escape_flux_mev_cm2_s']:.9g}" if np.isfinite(r["forward_escape_flux_mev_cm2_s"]) else "nan",
                    f"{r['escaped_electron_flux_mev_cm2_s']:.9g}"
                    if np.isfinite(r["escaped_electron_flux_mev_cm2_s"])
                    else "nan",
                    f"{r['escaped_photon_flux_mev_cm2_s']:.9g}"
                    if np.isfinite(r["escaped_photon_flux_mev_cm2_s"])
                    else "nan",
                    f"{r['backscatter_electron_flux_mev_cm2_s']:.9g}"
                    if np.isfinite(r["backscatter_electron_flux_mev_cm2_s"])
                    else "nan",
                    f"{r['backscatter_photon_flux_mev_cm2_s']:.9g}"
                    if np.isfinite(r["backscatter_photon_flux_mev_cm2_s"])
                    else "nan",
                    f"{r['escaped_lateral_electron_flux_mev_cm2_s']:.9g}"
                    if np.isfinite(r["escaped_lateral_electron_flux_mev_cm2_s"])
                    else "nan",
                    f"{r['escaped_lateral_photon_flux_mev_cm2_s']:.9g}"
                    if np.isfinite(r["escaped_lateral_photon_flux_mev_cm2_s"])
                    else "nan",
                    f"{r['jde_flux_cm2_s']:.9g}" if np.isfinite(r["jde_flux_cm2_s"]) else "nan",
                    f"{r['jde_energy_flux_model_mev_cm2_s']:.9g}" if np.isfinite(r["jde_energy_flux_model_mev_cm2_s"]) else "nan",
                    f"{deflected_frac:.9g}" if np.isfinite(deflected_frac) else "nan",
                    f"{deflected_flux:.9g}" if np.isfinite(deflected_flux) else "nan",
                ]
            )

    database_npz = out_dir / "latlon_metrics_database.npz"
    np.savez_compressed(
        database_npz,
        dose_depth_cm_selected=np.float64(dose_depth_cm),
        depth_edges_cm=np.asarray(depth_edges_cm, dtype=np.float64),
        depth_centers_cm=np.asarray(depth_centers_cm, dtype=np.float64),
        lat_values=np.asarray(lat_values, dtype=np.float64),
        lon_values=np.asarray(lon_values, dtype=np.float64),
        dose_map_at_depth_mgy_per_yr=np.asarray(dose_map, dtype=np.float64),
        dose_map_depth_integrated_mgy_per_yr=np.asarray(dose_map_depth_integrated, dtype=np.float64),
        dose_profile_map_mgy_per_yr=np.asarray(dose_profile_map, dtype=np.float64),
        primary_penetration_map_cm=np.asarray(primary_pen_map, dtype=np.float64),
        cumulative_63pct_depth_map_cm=np.asarray(cumulative_depth_map, dtype=np.float64),
        secondary_ke_map_mev=np.asarray(secondary_ke_map, dtype=np.float64),
        secondary_penetration_map_cm=np.asarray(secondary_pen_map, dtype=np.float64),
        deposited_fraction_map=np.asarray(dep_fraction_map, dtype=np.float64),
        backscatter_fraction_map=np.asarray(backscatter_fraction_map, dtype=np.float64),
        forward_escape_fraction_map=np.asarray(forward_escape_fraction_map, dtype=np.float64),
        lateral_escape_fraction_map=np.asarray(lateral_escape_fraction_map, dtype=np.float64),
        escaped_electron_fraction_map=np.asarray(escaped_electron_fraction_map, dtype=np.float64),
        escaped_photon_fraction_map=np.asarray(escaped_photon_fraction_map, dtype=np.float64),
        nondeposited_escaped_electron_fraction_map=np.asarray(
            nondeposited_escaped_electron_fraction_map, dtype=np.float64
        ),
        nondeposited_escaped_photon_fraction_map=np.asarray(
            nondeposited_escaped_photon_fraction_map, dtype=np.float64
        ),
        nondeposited_unresolved_fraction_map=np.asarray(nondeposited_unresolved_fraction_map, dtype=np.float64),
        backscatter_electron_fraction_map=np.asarray(backscatter_electron_fraction_map, dtype=np.float64),
        backscatter_photon_fraction_map=np.asarray(backscatter_photon_fraction_map, dtype=np.float64),
        forward_escape_electron_fraction_map=np.asarray(forward_escape_electron_fraction_map, dtype=np.float64),
        forward_escape_photon_fraction_map=np.asarray(forward_escape_photon_fraction_map, dtype=np.float64),
        lateral_escape_electron_fraction_map=np.asarray(lateral_escape_electron_fraction_map, dtype=np.float64),
        lateral_escape_photon_fraction_map=np.asarray(lateral_escape_photon_fraction_map, dtype=np.float64),
        deflected_fraction_map=np.asarray(deflected_fraction_map, dtype=np.float64),
        primary_flux_map_mev_cm2_s=np.asarray(primary_flux_map, dtype=np.float64),
        deposited_flux_map_mev_cm2_s=np.asarray(deposited_flux_map, dtype=np.float64),
        jde_flux_map_cm2_s=np.asarray(jde_flux_map, dtype=np.float64),
        jde_energy_flux_model_map_mev_cm2_s=np.asarray(jde_energy_flux_model_map, dtype=np.float64),
        jde_energy_flux_map_mev_cm2_s=np.asarray(jde_energy_flux_map, dtype=np.float64),
        jde_dose_equiv_map_mgy_per_yr=np.asarray(jde_dose_equiv_map, dtype=np.float64),
        backscatter_amount_map_mev_cm2_s=np.asarray(backscatter_amount_map, dtype=np.float64),
        forward_escape_amount_map_mev_cm2_s=np.asarray(forward_escape_amount_map, dtype=np.float64),
        lateral_escape_amount_map_mev_cm2_s=np.asarray(lateral_escape_amount_map, dtype=np.float64),
        escaped_electron_amount_map_mev_cm2_s=np.asarray(escaped_electron_amount_map, dtype=np.float64),
        escaped_photon_amount_map_mev_cm2_s=np.asarray(escaped_photon_amount_map, dtype=np.float64),
        backscatter_electron_amount_map_mev_cm2_s=np.asarray(backscatter_electron_amount_map, dtype=np.float64),
        backscatter_photon_amount_map_mev_cm2_s=np.asarray(backscatter_photon_amount_map, dtype=np.float64),
        forward_escape_electron_amount_map_mev_cm2_s=np.asarray(forward_escape_electron_amount_map, dtype=np.float64),
        forward_escape_photon_amount_map_mev_cm2_s=np.asarray(forward_escape_photon_amount_map, dtype=np.float64),
        lateral_escape_electron_amount_map_mev_cm2_s=np.asarray(lateral_escape_electron_amount_map, dtype=np.float64),
        lateral_escape_photon_amount_map_mev_cm2_s=np.asarray(lateral_escape_photon_amount_map, dtype=np.float64),
        deflected_amount_map_mev_cm2_s=np.asarray(deflected_amount_map, dtype=np.float64),
        energy_deposited_fraction_map=np.asarray(dep_fraction_map, dtype=np.float64),
        energy_forward_escape_fraction_map=np.asarray(forward_escape_fraction_map, dtype=np.float64),
        energy_lateral_escape_fraction_map=np.asarray(lateral_escape_fraction_map, dtype=np.float64),
        energy_escaped_electron_fraction_map=np.asarray(escaped_electron_fraction_map, dtype=np.float64),
        energy_escaped_photon_fraction_map=np.asarray(escaped_photon_fraction_map, dtype=np.float64),
        energy_nondeposited_escaped_electron_fraction_map=np.asarray(
            nondeposited_escaped_electron_fraction_map, dtype=np.float64
        ),
        energy_nondeposited_escaped_photon_fraction_map=np.asarray(
            nondeposited_escaped_photon_fraction_map, dtype=np.float64
        ),
        energy_nondeposited_unresolved_fraction_map=np.asarray(nondeposited_unresolved_fraction_map, dtype=np.float64),
        energy_deflected_fraction_map=np.asarray(deflected_fraction_map, dtype=np.float64),
        energy_deflected_amount_map_mev_cm2_s=np.asarray(deflected_amount_map, dtype=np.float64),
        cell_dose_profile_mgy_per_yr=np.asarray(cell_dose_profile, dtype=np.float64),
        cell_dose_profile_std_mgy_per_yr=np.asarray(cell_dose_profile_std, dtype=np.float64),
        cell_dose_at_depth_mgy_per_yr=np.asarray(cell_dose_at_depth, dtype=np.float64),
        cell_primary_penetration_cm=np.asarray(cell_primary_pen, dtype=np.float64),
        cell_secondary_ke_mev=np.asarray(cell_secondary_ke, dtype=np.float64),
        cell_secondary_penetration_cm=np.asarray(cell_secondary_pen, dtype=np.float64),
        cell_deposited_fraction=np.asarray(cell_deposited_fraction, dtype=np.float64),
        cell_backscatter_fraction=np.asarray(cell_backscatter_fraction, dtype=np.float64),
        cell_forward_escape_fraction=np.asarray(cell_forward_escape_fraction, dtype=np.float64),
        cell_lateral_escape_fraction=np.asarray(cell_lateral_escape_fraction, dtype=np.float64),
        cell_escaped_electron_fraction=np.asarray(cell_escaped_electron_fraction, dtype=np.float64),
        cell_escaped_photon_fraction=np.asarray(cell_escaped_photon_fraction, dtype=np.float64),
        cell_nondeposited_escaped_electron_fraction=np.asarray(
            cell_nondeposited_escaped_electron_fraction, dtype=np.float64
        ),
        cell_nondeposited_escaped_photon_fraction=np.asarray(
            cell_nondeposited_escaped_photon_fraction, dtype=np.float64
        ),
        cell_nondeposited_unresolved_fraction=np.asarray(cell_nondeposited_unresolved_fraction, dtype=np.float64),
        cell_backscatter_electron_fraction=np.asarray(cell_backscatter_electron_fraction, dtype=np.float64),
        cell_backscatter_photon_fraction=np.asarray(cell_backscatter_photon_fraction, dtype=np.float64),
        cell_forward_escape_electron_fraction=np.asarray(cell_forward_escape_electron_fraction, dtype=np.float64),
        cell_forward_escape_photon_fraction=np.asarray(cell_forward_escape_photon_fraction, dtype=np.float64),
        cell_lateral_escape_electron_fraction=np.asarray(cell_lateral_escape_electron_fraction, dtype=np.float64),
        cell_lateral_escape_photon_fraction=np.asarray(cell_lateral_escape_photon_fraction, dtype=np.float64),
        cell_backscatter_amount_mev_cm2_s=np.asarray(cell_backscatter_amount, dtype=np.float64),
        cell_forward_escape_amount_mev_cm2_s=np.asarray(cell_forward_escape_amount, dtype=np.float64),
        cell_lateral_escape_amount_mev_cm2_s=np.asarray(cell_lateral_escape_amount, dtype=np.float64),
        cell_escaped_electron_amount_mev_cm2_s=np.asarray(cell_escaped_electron_amount, dtype=np.float64),
        cell_escaped_photon_amount_mev_cm2_s=np.asarray(cell_escaped_photon_amount, dtype=np.float64),
        cell_backscatter_electron_amount_mev_cm2_s=np.asarray(cell_backscatter_electron_amount, dtype=np.float64),
        cell_backscatter_photon_amount_mev_cm2_s=np.asarray(cell_backscatter_photon_amount, dtype=np.float64),
        cell_forward_escape_electron_amount_mev_cm2_s=np.asarray(cell_forward_escape_electron_amount, dtype=np.float64),
        cell_forward_escape_photon_amount_mev_cm2_s=np.asarray(cell_forward_escape_photon_amount, dtype=np.float64),
        cell_lateral_escape_electron_amount_mev_cm2_s=np.asarray(cell_lateral_escape_electron_amount, dtype=np.float64),
        cell_lateral_escape_photon_amount_mev_cm2_s=np.asarray(cell_lateral_escape_photon_amount, dtype=np.float64),
        cell_primary_flux_mev_cm2_s=np.asarray(cell_primary_flux, dtype=np.float64),
        cell_deposited_flux_mev_cm2_s=np.asarray(cell_deposited_flux, dtype=np.float64),
        cell_jde_flux_cm2_s=np.asarray(cell_jde_flux, dtype=np.float64),
        cell_deflected_fraction=np.asarray(cell_deflected_fraction, dtype=np.float64),
        cell_deflected_amount_mev_cm2_s=np.asarray(cell_deflected_amount, dtype=np.float64),
    )

    maps_npz = out_dir / "latlon_metric_maps.npz"
    np.savez_compressed(
        maps_npz,
        dose_depth_cm_selected=np.float64(dose_depth_cm),
        depth_edges_cm=np.asarray(depth_edges_cm, dtype=np.float64),
        depth_centers_cm=np.asarray(depth_centers_cm, dtype=np.float64),
        lat_values=np.asarray(lat_values, dtype=np.float64),
        lon_values=np.asarray(lon_values, dtype=np.float64),
        dose_map_at_depth_mgy_per_yr=np.asarray(dose_map, dtype=np.float64),
        dose_map_depth_integrated_mgy_per_yr=np.asarray(dose_map_depth_integrated, dtype=np.float64),
        dose_profile_map_mgy_per_yr=np.asarray(dose_profile_map, dtype=np.float64),
        dose_profile_std_map_mgy_per_yr=np.asarray(dose_profile_std_map, dtype=np.float64),
        primary_penetration_map_cm=np.asarray(primary_pen_map, dtype=np.float64),
        cumulative_63pct_depth_map_cm=np.asarray(cumulative_depth_map, dtype=np.float64),
        secondary_ke_map_mev=np.asarray(secondary_ke_map, dtype=np.float64),
        secondary_penetration_map_cm=np.asarray(secondary_pen_map, dtype=np.float64),
        deposited_fraction_map=np.asarray(dep_fraction_map, dtype=np.float64),
        backscatter_fraction_map=np.asarray(backscatter_fraction_map, dtype=np.float64),
        forward_escape_fraction_map=np.asarray(forward_escape_fraction_map, dtype=np.float64),
        lateral_escape_fraction_map=np.asarray(lateral_escape_fraction_map, dtype=np.float64),
        escaped_electron_fraction_map=np.asarray(escaped_electron_fraction_map, dtype=np.float64),
        escaped_photon_fraction_map=np.asarray(escaped_photon_fraction_map, dtype=np.float64),
        nondeposited_escaped_electron_fraction_map=np.asarray(
            nondeposited_escaped_electron_fraction_map, dtype=np.float64
        ),
        nondeposited_escaped_photon_fraction_map=np.asarray(
            nondeposited_escaped_photon_fraction_map, dtype=np.float64
        ),
        nondeposited_unresolved_fraction_map=np.asarray(nondeposited_unresolved_fraction_map, dtype=np.float64),
        backscatter_electron_fraction_map=np.asarray(backscatter_electron_fraction_map, dtype=np.float64),
        backscatter_photon_fraction_map=np.asarray(backscatter_photon_fraction_map, dtype=np.float64),
        forward_escape_electron_fraction_map=np.asarray(forward_escape_electron_fraction_map, dtype=np.float64),
        forward_escape_photon_fraction_map=np.asarray(forward_escape_photon_fraction_map, dtype=np.float64),
        lateral_escape_electron_fraction_map=np.asarray(lateral_escape_electron_fraction_map, dtype=np.float64),
        lateral_escape_photon_fraction_map=np.asarray(lateral_escape_photon_fraction_map, dtype=np.float64),
        deflected_fraction_map=np.asarray(deflected_fraction_map, dtype=np.float64),
        primary_flux_map_mev_cm2_s=np.asarray(primary_flux_map, dtype=np.float64),
        deposited_flux_map_mev_cm2_s=np.asarray(deposited_flux_map, dtype=np.float64),
        jde_flux_map_cm2_s=np.asarray(jde_flux_map, dtype=np.float64),
        jde_energy_flux_model_map_mev_cm2_s=np.asarray(jde_energy_flux_model_map, dtype=np.float64),
        jde_energy_flux_map_mev_cm2_s=np.asarray(jde_energy_flux_map, dtype=np.float64),
        jde_dose_equiv_map_mgy_per_yr=np.asarray(jde_dose_equiv_map, dtype=np.float64),
        backscatter_amount_map_mev_cm2_s=np.asarray(backscatter_amount_map, dtype=np.float64),
        forward_escape_amount_map_mev_cm2_s=np.asarray(forward_escape_amount_map, dtype=np.float64),
        lateral_escape_amount_map_mev_cm2_s=np.asarray(lateral_escape_amount_map, dtype=np.float64),
        escaped_electron_amount_map_mev_cm2_s=np.asarray(escaped_electron_amount_map, dtype=np.float64),
        escaped_photon_amount_map_mev_cm2_s=np.asarray(escaped_photon_amount_map, dtype=np.float64),
        backscatter_electron_amount_map_mev_cm2_s=np.asarray(backscatter_electron_amount_map, dtype=np.float64),
        backscatter_photon_amount_map_mev_cm2_s=np.asarray(backscatter_photon_amount_map, dtype=np.float64),
        forward_escape_electron_amount_map_mev_cm2_s=np.asarray(forward_escape_electron_amount_map, dtype=np.float64),
        forward_escape_photon_amount_map_mev_cm2_s=np.asarray(forward_escape_photon_amount_map, dtype=np.float64),
        lateral_escape_electron_amount_map_mev_cm2_s=np.asarray(lateral_escape_electron_amount_map, dtype=np.float64),
        lateral_escape_photon_amount_map_mev_cm2_s=np.asarray(lateral_escape_photon_amount_map, dtype=np.float64),
        deflected_amount_map_mev_cm2_s=np.asarray(deflected_amount_map, dtype=np.float64),
        energy_deposited_fraction_map=np.asarray(dep_fraction_map, dtype=np.float64),
        energy_forward_escape_fraction_map=np.asarray(forward_escape_fraction_map, dtype=np.float64),
        energy_lateral_escape_fraction_map=np.asarray(lateral_escape_fraction_map, dtype=np.float64),
        energy_escaped_electron_fraction_map=np.asarray(escaped_electron_fraction_map, dtype=np.float64),
        energy_escaped_photon_fraction_map=np.asarray(escaped_photon_fraction_map, dtype=np.float64),
        energy_nondeposited_escaped_electron_fraction_map=np.asarray(
            nondeposited_escaped_electron_fraction_map, dtype=np.float64
        ),
        energy_nondeposited_escaped_photon_fraction_map=np.asarray(
            nondeposited_escaped_photon_fraction_map, dtype=np.float64
        ),
        energy_nondeposited_unresolved_fraction_map=np.asarray(nondeposited_unresolved_fraction_map, dtype=np.float64),
        energy_deflected_fraction_map=np.asarray(deflected_fraction_map, dtype=np.float64),
        energy_deflected_amount_map_mev_cm2_s=np.asarray(deflected_amount_map, dtype=np.float64),
    )

    lon_arr = np.asarray(lon_values, dtype=np.float64)
    leading_cols = np.where(lon_arr < 180.0)[0]
    trailing_cols = np.where(lon_arr >= 180.0)[0]

    leading_profiles_npz = out_dir / "europa_leading_pixel_dose_profiles_mgyyr.npz"
    np.savez_compressed(
        leading_profiles_npz,
        dose_depth_cm_selected=np.float64(dose_depth_cm),
        depth_edges_cm=np.asarray(depth_edges_cm, dtype=np.float64),
        depth_centers_cm=np.asarray(depth_centers_cm, dtype=np.float64),
        lat_values=np.asarray(lat_values, dtype=np.float64),
        lon_values=np.asarray(lon_arr[leading_cols], dtype=np.float64),
        dose_profile_map_mgy_per_yr=np.asarray(dose_profile_map[:, leading_cols, :], dtype=np.float64),
        dose_profile_std_map_mgy_per_yr=np.asarray(dose_profile_std_map[:, leading_cols, :], dtype=np.float64),
        dose_map_at_depth_mgy_per_yr=np.asarray(dose_map[:, leading_cols], dtype=np.float64),
        dose_map_depth_integrated_mgy_per_yr=np.asarray(dose_map_depth_integrated[:, leading_cols], dtype=np.float64),
        jde_flux_map_cm2_s=np.asarray(jde_flux_map[:, leading_cols], dtype=np.float64),
        jde_energy_flux_model_map_mev_cm2_s=np.asarray(
            jde_energy_flux_model_map[:, leading_cols], dtype=np.float64
        ),
        jde_energy_flux_map_mev_cm2_s=np.asarray(
            jde_energy_flux_map[:, leading_cols], dtype=np.float64
        ),
        jde_dose_equiv_map_mgy_per_yr=np.asarray(
            jde_dose_equiv_map[:, leading_cols], dtype=np.float64
        ),
    )

    trailing_profiles_npz = out_dir / "europa_trailing_pixel_dose_profiles_mgyyr.npz"
    np.savez_compressed(
        trailing_profiles_npz,
        dose_depth_cm_selected=np.float64(dose_depth_cm),
        depth_edges_cm=np.asarray(depth_edges_cm, dtype=np.float64),
        depth_centers_cm=np.asarray(depth_centers_cm, dtype=np.float64),
        lat_values=np.asarray(lat_values, dtype=np.float64),
        lon_values=np.asarray(lon_arr[trailing_cols], dtype=np.float64),
        dose_profile_map_mgy_per_yr=np.asarray(dose_profile_map[:, trailing_cols, :], dtype=np.float64),
        dose_profile_std_map_mgy_per_yr=np.asarray(dose_profile_std_map[:, trailing_cols, :], dtype=np.float64),
        dose_map_at_depth_mgy_per_yr=np.asarray(dose_map[:, trailing_cols], dtype=np.float64),
        dose_map_depth_integrated_mgy_per_yr=np.asarray(dose_map_depth_integrated[:, trailing_cols], dtype=np.float64),
        jde_flux_map_cm2_s=np.asarray(jde_flux_map[:, trailing_cols], dtype=np.float64),
        jde_energy_flux_model_map_mev_cm2_s=np.asarray(
            jde_energy_flux_model_map[:, trailing_cols], dtype=np.float64
        ),
        jde_energy_flux_map_mev_cm2_s=np.asarray(
            jde_energy_flux_map[:, trailing_cols], dtype=np.float64
        ),
        jde_dose_equiv_map_mgy_per_yr=np.asarray(
            jde_dose_equiv_map[:, trailing_cols], dtype=np.float64
        ),
    )

    if not args.no_plots:
        _configure_plot_style()
        _save_hemisphere_plots(
            out_dir=out_dir,
            hemisphere="leading",
            lon_values=lon_arr,
            energy_panel_flux_map=deposited_flux_map,
            primary_flux_map=primary_flux_map,
            primary_pen_map=primary_pen_map,
            secondary_ke_map=secondary_ke_map,
            secondary_pen_map=secondary_pen_map,
            cumulative_depth_map=cumulative_depth_map,
            dep_fraction_map=dep_fraction_map,
            backscatter_amount_map=backscatter_amount_map,
            dose_depth_cm=dose_depth_cm,
            backscatter_fraction_map=backscatter_fraction_map,
            lateral_escape_fraction_map=lateral_escape_fraction_map,
            forward_escape_fraction_map=forward_escape_fraction_map,
            escaped_electron_fraction_map=escaped_electron_fraction_map,
            escaped_photon_fraction_map=escaped_photon_fraction_map,
            backscatter_electron_fraction_map=backscatter_electron_fraction_map,
            backscatter_photon_fraction_map=backscatter_photon_fraction_map,
            lateral_escape_electron_fraction_map=lateral_escape_electron_fraction_map,
            lateral_escape_photon_fraction_map=lateral_escape_photon_fraction_map,
            forward_escape_electron_fraction_map=forward_escape_electron_fraction_map,
            forward_escape_photon_fraction_map=forward_escape_photon_fraction_map,
            jde_energy_flux_map_mev_cm2_s=jde_energy_flux_map,
        )
        _save_hemisphere_plots(
            out_dir=out_dir,
            hemisphere="trailing",
            lon_values=lon_arr,
            energy_panel_flux_map=deposited_flux_map,
            primary_flux_map=primary_flux_map,
            primary_pen_map=primary_pen_map,
            secondary_ke_map=secondary_ke_map,
            secondary_pen_map=secondary_pen_map,
            cumulative_depth_map=cumulative_depth_map,
            dep_fraction_map=dep_fraction_map,
            backscatter_amount_map=backscatter_amount_map,
            dose_depth_cm=dose_depth_cm,
            backscatter_fraction_map=backscatter_fraction_map,
            lateral_escape_fraction_map=lateral_escape_fraction_map,
            forward_escape_fraction_map=forward_escape_fraction_map,
            escaped_electron_fraction_map=escaped_electron_fraction_map,
            escaped_photon_fraction_map=escaped_photon_fraction_map,
            backscatter_electron_fraction_map=backscatter_electron_fraction_map,
            backscatter_photon_fraction_map=backscatter_photon_fraction_map,
            lateral_escape_electron_fraction_map=lateral_escape_electron_fraction_map,
            lateral_escape_photon_fraction_map=lateral_escape_photon_fraction_map,
            forward_escape_electron_fraction_map=forward_escape_electron_fraction_map,
            forward_escape_photon_fraction_map=forward_escape_photon_fraction_map,
            jde_energy_flux_map_mev_cm2_s=jde_energy_flux_map,
        )
        _save_energy_deposition_depth_profiles_plot(
            out_dir=out_dir,
            lat_values=np.asarray(lat_values, dtype=np.float64),
            lon_values=lon_arr,
            depth_edges_cm=np.asarray(depth_edges_cm, dtype=np.float64),
            dose_profile_map_mgy_per_yr=np.asarray(dose_profile_map, dtype=np.float64),
            dose_profile_std_map_mgy_per_yr=np.asarray(dose_profile_std_map, dtype=np.float64),
            density_gcm3=density_gcm3,
            show_plots=False,
        )

    print(f"Wrote summary CSV: {summary_csv}")
    print(f"Wrote range metrics CSV: {range_metrics_csv}")
    print(f"Wrote database NPZ: {database_npz}")
    print(f"Wrote map-matrix NPZ: {maps_npz}")
    print(f"Wrote hemisphere dose-profile NPZ: {leading_profiles_npz}")
    print(f"Wrote hemisphere dose-profile NPZ: {trailing_profiles_npz}")
    if not args.no_plots:
        print(f"Wrote plots under: {out_dir}")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Analyze Europa per-energy ROOT library into lat-lon metrics maps."
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_prepare = sub.add_parser("prepare", help="Build unique-range manifest from library CSVs.")
    p_prepare.add_argument(
        "--library-dir",
        type=Path,
        default=None,
        help="Path to europa_energy_library directory.",
    )
    p_prepare.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for manifest/cache/results.",
    )
    p_prepare.add_argument(
        "--z-range-file",
        type=Path,
        default=None,
        help="Optional depth grid file values (use --z-range-unit to declare units).",
    )
    p_prepare.add_argument(
        "--z-range-unit",
        choices=["cm", "mm"],
        default="cm",
        help="Unit of --z-range-file values (default: cm).",
    )
    p_prepare.add_argument(
        "--dose-depth-cm",
        type=float,
        default=0.01,
        help="Depth used for the map value in MGy/yr panel.",
    )
    p_prepare.add_argument(
        "--density-gcm3",
        type=float,
        default=float("nan"),
        help="Default density (g/cm3). If omitted, infer from run table.",
    )
    p_prepare.set_defaults(func=cmd_prepare)

    p_cache = sub.add_parser(
        "build-energy-cache",
        help="Legacy one-shot cache build (local multiprocess).",
    )
    p_cache.add_argument(
        "--manifest",
        type=Path,
        required=True,
        help="Manifest path from `prepare`.",
    )
    p_cache.add_argument(
        "--root-dir",
        type=Path,
        default=None,
        help="Override ROOT directory containing dna_*.root files.",
    )
    p_cache.add_argument(
        "--work-root-dir",
        type=Path,
        default=None,
        help="Secondary ROOT directory fallback.",
    )
    p_cache.add_argument(
        "--step-size",
        type=str,
        default="50 MB",
        help="uproot iterate step_size (default: 50 MB).",
    )
    p_cache.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel workers for per-energy ROOT reads (default: 1).",
    )
    p_cache.add_argument(
        "--depth-origin-cm",
        type=float,
        default=0.0,
        help="Depth origin offset in cm (depth = z_cm - origin).",
    )
    p_cache.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing cache files.",
    )
    p_cache.set_defaults(func=cmd_build_energy_cache)

    p_ecount = sub.add_parser("energy-count", help="Print number of energies in run table.")
    p_ecount.add_argument("--manifest", type=Path, required=True)
    p_ecount.set_defaults(func=cmd_energy_count)

    p_ework = sub.add_parser("energy-worker", help="Compute per-energy observables for one energy position.")
    p_ework.add_argument("--manifest", type=Path, required=True)
    p_ework.add_argument(
        "--energy-pos",
        type=int,
        required=True,
        help="Row position in sorted run table (0-based).",
    )
    p_ework.add_argument(
        "--root-dir",
        type=Path,
        default=None,
        help="Override ROOT directory containing dna_*.root files.",
    )
    p_ework.add_argument(
        "--work-root-dir",
        type=Path,
        default=None,
        help="Secondary ROOT directory fallback.",
    )
    p_ework.add_argument(
        "--step-size",
        type=str,
        default="50 MB",
        help="uproot iterate step_size (default: 50 MB).",
    )
    p_ework.add_argument(
        "--depth-origin-cm",
        type=float,
        default=0.0,
        help="Depth origin offset in cm (depth = z_cm - origin).",
    )
    p_ework.add_argument(
        "--force",
        action="store_true",
        help="Overwrite per-energy NPZ if it already exists.",
    )
    p_ework.set_defaults(func=cmd_energy_worker)

    p_emerge = sub.add_parser("energy-merge", help="Merge per-energy NPZs into cache CSV/NPZ products.")
    p_emerge.add_argument("--manifest", type=Path, required=True)
    p_emerge.add_argument(
        "--force",
        action="store_true",
        help="Overwrite merged cache outputs even if they exist.",
    )
    p_emerge.set_defaults(func=cmd_energy_merge)

    p_count = sub.add_parser("range-count", help="Print number of unique ranges in manifest.")
    p_count.add_argument("--manifest", type=Path, required=True)
    p_count.set_defaults(func=cmd_range_count)

    p_worker = sub.add_parser("range-worker", help="Compute one unique range result.")
    p_worker.add_argument("--manifest", type=Path, required=True)
    p_worker.add_argument(
        "--range-index",
        type=int,
        required=True,
        help="Index in manifest.unique_ranges (PBS_ARRAY_INDEX compatible).",
    )
    p_worker.add_argument(
        "--density-gcm3",
        type=float,
        default=float("nan"),
        help="Override density used for dose conversion.",
    )
    p_worker.add_argument(
        "--dose-depth-cm",
        type=float,
        default=float("nan"),
        help="Override map depth for dose panel.",
    )
    p_worker.set_defaults(func=cmd_range_worker)

    p_merge = sub.add_parser("merge", help="Merge all range results, write DB, and plot maps.")
    p_merge.add_argument("--manifest", type=Path, required=True)
    p_merge.add_argument(
        "--dose-depth-cm",
        type=float,
        default=float("nan"),
        help="Override plotted dose depth.",
    )
    p_merge.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot generation (write CSV/NPZ only).",
    )
    p_merge.set_defaults(func=cmd_merge)

    p_plot_npz = sub.add_parser(
        "plot-from-npz",
        help="Regenerate hemisphere plots from a saved map NPZ file.",
    )
    p_plot_npz.add_argument(
        "--maps-npz",
        type=Path,
        required=True,
        help="NPZ file with map matrices (typically latlon_metric_maps.npz).",
    )
    p_plot_npz.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Plot output directory (default: maps-npz parent directory).",
    )
    p_plot_npz.add_argument(
        "--dose-depth-cm",
        type=float,
        default=float("nan"),
        help="Override dose panel depth in cm; defaults to NPZ stored depth.",
    )
    p_plot_npz.add_argument(
        "--jde-grid-size",
        type=int,
        default=0,
        help="If >0, also resample J(E)E dE map to NxN and save NPZ (e.g., 90).",
    )
    p_plot_npz.add_argument(
        "--jde-grid-nlat",
        type=int,
        default=0,
        help="If >0 with --jde-grid-nlon, resample J(E)E dE map to nlat x nlon.",
    )
    p_plot_npz.add_argument(
        "--jde-grid-nlon",
        type=int,
        default=0,
        help="If >0 with --jde-grid-nlat, resample J(E)E dE map to nlat x nlon.",
    )
    p_plot_npz.set_defaults(func=cmd_plot_from_npz)

    return ap


def main() -> None:
    if len(sys.argv) == 1:
        data_dir = _default_local_data_dir()
        maps_npz = data_dir / "latlon_metric_maps.npz"
        if not maps_npz.exists():
            raise FileNotFoundError(
                "No CLI command provided and default maps NPZ is missing: "
                f"{maps_npz}. Run with an explicit subcommand, or place "
                "latlon_metric_maps.npz in the sibling data directory."
            )
        args = argparse.Namespace(
            maps_npz=maps_npz,
            out_dir=data_dir,
            dose_depth_cm=float("nan"),
        )
        cmd_plot_from_npz(args)
        return

    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
