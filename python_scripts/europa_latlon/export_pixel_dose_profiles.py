#!/usr/bin/env python3
"""
Export Europa latitude-longitude-depth dose profile NPZ products.

This is a lightweight release helper for users who already have
`latlon_metric_maps.npz` or `latlon_metrics_database.npz` and want explicit
pixel-dose profile products without rerunning the full `merge` workflow.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


OPTIONAL_2D_FIELDS = [
    "dose_map_at_depth_mgy_per_yr",
    "dose_map_depth_integrated_mgy_per_yr",
    "jde_flux_map_cm2_s",
    "jde_energy_flux_model_map_mev_cm2_s",
    "jde_energy_flux_map_mev_cm2_s",
    "jde_dose_equiv_map_mgy_per_yr",
]


def _required_array(data: np.lib.npyio.NpzFile, name: str) -> np.ndarray:
    if name not in data:
        raise KeyError(f"Required field missing from NPZ: {name}")
    return np.asarray(data[name], dtype=np.float64)


def _optional_2d(
    data: np.lib.npyio.NpzFile,
    name: str,
    shape: tuple[int, int],
) -> np.ndarray:
    if name not in data:
        return np.full(shape, np.nan, dtype=np.float64)
    arr = np.asarray(data[name], dtype=np.float64)
    if arr.shape != shape:
        raise ValueError(f"{name} has shape {arr.shape}; expected {shape}")
    return arr


def _write_profile_npz(
    out_path: Path,
    data: np.lib.npyio.NpzFile,
    cols: np.ndarray,
    lat_values: np.ndarray,
    lon_values: np.ndarray,
    depth_edges_cm: np.ndarray,
    depth_centers_cm: np.ndarray,
    dose_profile: np.ndarray,
    dose_profile_std: np.ndarray,
    dose_depth_cm: float,
) -> None:
    base_shape = (lat_values.size, lon_values.size)
    payload = {
        "dose_depth_cm_selected": np.float64(dose_depth_cm),
        "depth_edges_cm": np.asarray(depth_edges_cm, dtype=np.float64),
        "depth_centers_cm": np.asarray(depth_centers_cm, dtype=np.float64),
        "lat_values": np.asarray(lat_values, dtype=np.float64),
        "lon_values": np.asarray(lon_values[cols], dtype=np.float64),
        "dose_profile_map_mgy_per_yr": np.asarray(dose_profile[:, cols, :], dtype=np.float64),
        "dose_profile_std_map_mgy_per_yr": np.asarray(dose_profile_std[:, cols, :], dtype=np.float64),
    }
    for name in OPTIONAL_2D_FIELDS:
        payload[name] = np.asarray(_optional_2d(data, name, base_shape)[:, cols], dtype=np.float64)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_path, **payload)


def export_pixel_dose_profiles(
    maps_npz: Path,
    out_dir: Path,
    prefix: str,
    split_lon_deg: float,
    write_combined: bool,
) -> list[Path]:
    with np.load(maps_npz) as data:
        lat_values = _required_array(data, "lat_values")
        lon_values = _required_array(data, "lon_values")
        depth_edges_cm = _required_array(data, "depth_edges_cm")
        depth_centers_cm = _required_array(data, "depth_centers_cm")
        dose_profile = _required_array(data, "dose_profile_map_mgy_per_yr")
        if dose_profile.ndim != 3:
            raise ValueError("dose_profile_map_mgy_per_yr must have shape (lat, lon, depth)")

        expected_shape = (lat_values.size, lon_values.size, depth_centers_cm.size)
        if dose_profile.shape != expected_shape:
            raise ValueError(
                "dose_profile_map_mgy_per_yr has shape "
                f"{dose_profile.shape}; expected {expected_shape}"
            )

        if "dose_profile_std_map_mgy_per_yr" in data:
            dose_profile_std = np.asarray(data["dose_profile_std_map_mgy_per_yr"], dtype=np.float64)
            if dose_profile_std.shape != expected_shape:
                raise ValueError(
                    "dose_profile_std_map_mgy_per_yr has shape "
                    f"{dose_profile_std.shape}; expected {expected_shape}"
                )
        else:
            dose_profile_std = np.full(expected_shape, np.nan, dtype=np.float64)

        if "dose_depth_cm_selected" in data:
            dose_depth_cm = float(np.asarray(data["dose_depth_cm_selected"]).item())
        else:
            dose_depth_cm = float("nan")

        leading_cols = np.where(lon_values < split_lon_deg)[0]
        trailing_cols = np.where(lon_values >= split_lon_deg)[0]

        written: list[Path] = []
        leading_path = out_dir / f"{prefix}_leading_pixel_dose_profiles_mgyyr.npz"
        trailing_path = out_dir / f"{prefix}_trailing_pixel_dose_profiles_mgyyr.npz"
        _write_profile_npz(
            leading_path,
            data,
            leading_cols,
            lat_values,
            lon_values,
            depth_edges_cm,
            depth_centers_cm,
            dose_profile,
            dose_profile_std,
            dose_depth_cm,
        )
        _write_profile_npz(
            trailing_path,
            data,
            trailing_cols,
            lat_values,
            lon_values,
            depth_edges_cm,
            depth_centers_cm,
            dose_profile,
            dose_profile_std,
            dose_depth_cm,
        )
        written.extend([leading_path, trailing_path])

        if write_combined:
            combined_path = out_dir / f"{prefix}_latlon_pixel_dose_profiles_mgyyr.npz"
            all_cols = np.arange(lon_values.size)
            _write_profile_npz(
                combined_path,
                data,
                all_cols,
                lat_values,
                lon_values,
                depth_edges_cm,
                depth_centers_cm,
                dose_profile,
                dose_profile_std,
                dose_depth_cm,
            )
            written.append(combined_path)

    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export latitude-longitude-depth Europa dose profile NPZ files."
    )
    parser.add_argument(
        "--maps-npz",
        type=Path,
        required=True,
        help="Input latlon_metric_maps.npz or latlon_metrics_database.npz.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to the input NPZ directory.",
    )
    parser.add_argument(
        "--prefix",
        default="europa",
        help="Output filename prefix. Default: europa.",
    )
    parser.add_argument(
        "--split-lon-deg",
        type=float,
        default=180.0,
        help="Longitude split between leading and trailing hemispheres. Default: 180.",
    )
    parser.add_argument(
        "--combined",
        action="store_true",
        help="Also write one combined lat-lon-depth NPZ covering all longitudes.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    out_dir = args.out_dir if args.out_dir is not None else args.maps_npz.parent
    written = export_pixel_dose_profiles(
        maps_npz=args.maps_npz,
        out_dir=out_dir,
        prefix=args.prefix,
        split_lon_deg=args.split_lon_deg,
        write_combined=args.combined,
    )
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
