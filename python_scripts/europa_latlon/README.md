# Europa Latitude-Longitude Analysis

This directory contains the post-processing code used to construct Europa latitude-longitude maps from the monoenergetic transport library.

The main script is:

- `analyze_europa_latlon_metrics.py`

The release also includes:

- `export_pixel_dose_profiles.py`
  Export explicit latitude-longitude-depth pixel dose-profile NPZ files from a saved map/database NPZ.

It supports a staged workflow:

- `prepare`
  Build a manifest of unique latitude-longitude energy ranges from the Europa energy-library CSV files.
- `energy-worker` and `energy-merge`
  Compute and merge per-energy observables from ROOT outputs.
- `range-worker` and `merge`
  Build latitude-longitude dose, escape, and range metrics from the energy cache.
- `plot-from-npz`
  Regenerate diagnostic plots from saved map NPZ products.

For the release data bundled with `v1.0.1`, the most common local command is:

```bash
python3 python_scripts/europa_latlon/analyze_europa_latlon_metrics.py \
  plot-from-npz \
  --maps-npz europa_latlon_data/latlon_metric_maps.npz \
  --out-dir europa_latlon_data
```

The script uses `numpy` and `matplotlib` for loading and plotting the packaged products. Commands that read ROOT files also require `uproot`.

The compact derived data products are distributed separately in:

- `europa_latlon_data/`

## Pixel Dose-Profile Export

The `merge` command in `analyze_europa_latlon_metrics.py` writes the pixel dose-profile NPZ files during a full analysis run. If you already have `latlon_metric_maps.npz`, use the standalone exporter instead:

```bash
python3 python_scripts/europa_latlon/export_pixel_dose_profiles.py \
  --maps-npz europa_latlon_data/latlon_metric_maps.npz \
  --out-dir europa_latlon_data \
  --combined
```

This writes:

- `europa_leading_pixel_dose_profiles_mgyyr.npz`
- `europa_trailing_pixel_dose_profiles_mgyyr.npz`
- `europa_latlon_pixel_dose_profiles_mgyyr.npz` when `--combined` is used

The profile array field is:

- `dose_profile_map_mgy_per_yr`, with shape `(latitude, longitude, depth)`
