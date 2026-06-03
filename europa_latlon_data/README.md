# Europa Latitude-Longitude Derived Products

This directory contains compact derived products from the Europa latitude-longitude analysis workflow. These files are intended for plotting, inspection, and post-processing. They are not the raw per-energy or per-range worker outputs.

## Main Products

- `latlon_metric_maps.npz`
  Map matrices used to regenerate the leading and trailing hemisphere diagnostic figures.
  Includes `dose_profile_map_mgy_per_yr`, the full `(latitude, longitude, depth)` dose cube.
- `latlon_metrics_database.npz`
  Compact database of latitude-longitude metrics.
- `latlon_metrics_summary.csv`
  Tabular summary of latitude-longitude metrics.
- `cell_to_range.csv`
  Mapping from each latitude-longitude cell to a unique energy range.
- `unique_ranges.csv`
  Unique energy ranges used by the post-processing workflow.
- `unique_range_metrics.csv`
  Metrics collapsed by unique energy range.
- `energy_metrics.csv`
  Per-energy metrics used in the latitude-longitude analysis.
- `energy_feedback_summary.csv` and `energy_feedback_summary.npz`
  Energy-budget and feedback summaries.
- `energy_edep_profiles.npz`
  Energy-resolved deposited-energy profile cache.
- `energy_observables.npz`
  Energy-resolved observables used by the map construction.
- `depth_edges_cm.npy`
  Depth grid used by the profile products.

## Pixel Dose Profiles

The hemisphere-specific pixel dose-profile caches are:

- `europa_leading_pixel_dose_profiles_mgyyr.npz`
- `europa_trailing_pixel_dose_profiles_mgyyr.npz`

A combined all-longitudes file can be generated from `latlon_metric_maps.npz` with `python_scripts/europa_latlon/export_pixel_dose_profiles.py`.

Each file stores latitude-longitude maps and depth-resolved dose profiles for one hemisphere of Europa.

Both NPZ files contain:

- `dose_depth_cm_selected`: scalar depth in cm used for the 2D dose map at depth
- `depth_edges_cm`: depth-bin edges, shape `(N_depth + 1,)`
- `depth_centers_cm`: depth-bin centers, shape `(N_depth,)`
- `lat_values`: latitude grid, shape `(N_lat,)`
- `lon_values`: longitude grid for the selected hemisphere, shape `(N_lon,)`
- `dose_profile_map_mgy_per_yr`: depth-resolved dose map, shape `(N_lat, N_lon, N_depth)`
- `dose_profile_std_map_mgy_per_yr`: standard deviation for the depth-resolved dose map
- `dose_map_at_depth_mgy_per_yr`: 2D dose map at `dose_depth_cm_selected`
- `dose_map_depth_integrated_mgy_per_yr`: depth-integrated dose proxy
- `jde_flux_map_cm2_s`: JDE flux map
- `jde_energy_flux_model_map_mev_cm2_s`: model-based JDE energy-flux map
- `jde_energy_flux_map_mev_cm2_s`: JDE energy-flux map
- `jde_dose_equiv_map_mgy_per_yr`: JDE-equivalent dose map

For the current cached products:

- `N_lat = 30`
- `N_lon = 15`
- `N_depth = 4070`

## Plot Regeneration

From the repository root:

```bash
python3 python_scripts/europa_latlon/analyze_europa_latlon_metrics.py \
  plot-from-npz \
  --maps-npz europa_latlon_data/latlon_metric_maps.npz \
  --out-dir europa_latlon_data
```

This regenerates the leading and trailing hemisphere diagnostic figures from the packaged map NPZ.

## Notes

- Dose quantities are in `MGy / yr`.
- Depth coordinates are in cm.
- Longitude values are hemisphere-local where applicable.
- Raw worker outputs, logs, and archival backup files are intentionally not distributed in this release package.
