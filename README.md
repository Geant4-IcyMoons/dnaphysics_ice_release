# Geant4-IcyMoons

<p align="center">
  <img src="icymoons_logo.png" alt="Geant4-IcyMoons logo" />
</p>

## Overview

Geant4-IcyMoons is an extension of Geant4-DNA for Monte Carlo simulation of electron transport in water ice and related condensed phases. The present release extends the Geant4-DNA electron-transport framework with new water-ice interaction physics, most notably phase-specific inelastic electronic excitation and ionization cross sections for amorphous and hexagonal ice, together with the supporting models, data products, examples, and utilities required to use them in practice. At the level of the distributed transport physics, the release models the elastic and inelastic interactions of electrons in amorphous and hexagonal water ice, with the principal phase-resolved extension entering through the electronic inelastic sector.

The distribution provides:

- a standalone Geant4 application, `dnaphysics`
- new water-ice inelastic cross-section datasets for amorphous and hexagonal ice
- Geant4 model implementations and physics-list integration for those new datasets and the supporting transport channels used by the release
- validated example macros and run scripts
- Python utilities for diagnostics, plotting, and cross-section generation
- a Europa energy-library generation workflow

This release is intended for users who require a documented, reproducible framework for:

- transport simulations in amorphous ice, hexagonal ice, and liquid water
- comparison of standard and custom low-energy interaction channels
- W-value calculations
- inspection and regeneration of the distributed cross-section tables
- generation of monoenergetic source libraries for Europa surface studies

## Scope Of The Release

The release extends standard Geant4-DNA electron transport to water-ice media through new models and distributed cross-section tables. In particular, it distributes:

- phase-specific electronic excitation tables for amorphous and hexagonal water ice
- phase-specific ionization tables for amorphous and hexagonal water ice
- differential `sigmadiff_*` products used where angular or energy-loss sampling requires more than total cross sections alone
- Michaud-based vibrational-excitation data for low-energy condensed-ice transport
- Michaud-based low-energy attachment data
- multiple elastic datasets and angular-sampling tables spanning the low- and high-energy transport treatment used by the release

The executable supports three runtime physics modes:

- `ice_am`: amorphous ice
- `ice_hex`: hexagonal ice
- `water`: liquid water

If no mode is specified, the application defaults to `ice_hex`.

In the present release, the phase dependence between amorphous and hexagonal ice enters primarily through the electronic inelastic sector, namely excitation and ionization. The supporting elastic, vibrational-excitation, and attachment channels remain distributed as the transport inputs used by the current implementation, but are not yet fully phase-dependent in the same way as the electronic inelastic channels.

## Repository Contents

The principal directories in the release are:

- `cross_sections/`
  Distributed `.dat` files required by the custom ice transport models.
- `src/`
  Geant4 model implementations, physics lists, and application code.
- `include/`
  Headers corresponding to the distributed source code.
- `examples/`
  Maintained example macros, run scripts, and example-specific plotting helpers.
- `python_scripts/plotting/`
  Diagnostics and plotting utilities for ROOT outputs and tabulated cross sections.
- `python_scripts/physics_ice/`
  Utilities for generating, transforming, and visualizing the distributed ice cross sections.
- `python_scripts/europa/`
  Europa energy-library and macro-generation utilities.
- `tabular/`
  Tabular source data used in the construction of selected distributed cross sections.

## Software Requirements

The release requires:

- Geant4 11.x with Geant4-DNA available
- CMake 3.16 or later
- a C++17-capable compiler
- Geant4 low-energy data, including a valid `G4EMLOW...` installation

The custom models resolve data files through the standard Geant4 environment variable `G4LEDATA`.

## Installation Of Distributed Cross Sections

The custom ice models do not read the repository-local `cross_sections/` directory directly at runtime. Instead, they resolve their input files through `G4LEDATA/dna`. Accordingly, the distributed `.dat` files must be copied or linked into the Geant4 low-energy data tree before execution.

Typical installation procedure:

```bash
export G4LEDATA=/path/to/Geant4/data/G4EMLOW8.6.1
mkdir -p "${G4LEDATA}/dna"
cp cross_sections/* "${G4LEDATA}/dna/"
```

If symbolic links are preferred:

```bash
export G4LEDATA=/path/to/Geant4/data/G4EMLOW8.6.1
mkdir -p "${G4LEDATA}/dna"
for f in cross_sections/*; do
  ln -sf "$(pwd)/$f" "${G4LEDATA}/dna/$(basename "$f")"
done
```

This step is mandatory for the custom ice models. If the required files are not visible under `G4LEDATA/dna`, the corresponding models will not initialize correctly.

## Build Procedure

From the repository root:

```bash
cmake -S . -B build -DGeant4_DIR=/path/to/geant4-install/lib/Geant4-11.x
cmake --build build -j
```

This generates:

- `build/dnaphysics`

The build system also copies top-level `.mac`, `.in`, and `.C` files into the build directory to facilitate execution from within `build/`.

## Execution

### Batch Execution

```bash
DNA_PHYSICS=ice_am ./build/dnaphysics path/to/macro.mac 4
```

The positional arguments are:

1. macro file
2. number of Geant4 threads

### Interactive Execution

```bash
./build/dnaphysics
```

When invoked without a macro argument, the application starts an interactive Geant4 UI session and executes `vis.mac`.

## Runtime Environment Variables

The principal runtime variables supported by this release are listed below.

| Variable | Purpose |
|---|---|
| `DNA_PHYSICS` | Selects `ice_am`, `ice_hex`, or `water`. |
| `DNA_ROOT_BASENAME` | Sets the base name of the ROOT output file. Default: `dna`. |
| `DNA_WVALUE_FILE` | Appends W-value rows to the specified text file. |
| `DNA_NTUPLE_FILES` | Sets the number of reduced ROOT ntuple files in merged multi-threaded output. |
| `DNA_NTUPLE_MERGE` | Enables or disables ntuple merging. |
| `DNA_NTUPLE_STRINGS` | Enables string ntuple columns where supported. |
| `DNA_KEEP_OLD_ROOT` | Preserves pre-existing matching ROOT outputs. |
| `DNA_ROOT_SPLIT_EVENTS` | Rotates output after a fixed number of events. |
| `DNA_ROOT_MAX_MB` | Rotates output after a size threshold. |

For water-mode studies, the water physics list also exposes optional feature toggles such as:

- `DNA_WATER_ENABLE_HIGH_EM`
- `DNA_WATER_ENABLE_ELASTIC`
- `DNA_WATER_ENABLE_EXCITATION`
- `DNA_WATER_ENABLE_IONISATION`

## Application Macro Commands

The release defines several convenience commands under `/dna/test/`:

- `/dna/test/setMat NAME`
  Select the Geant4 material.
- `/dna/test/setMatDens NAME VALUE UNIT`
  Override the density of a named material.
- `/dna/test/setIceSize X Y Z UNIT`
  Set the slab dimensions as full lengths.
- `/dna/test/setLogMode minimal|full`
  Select compact or detailed ROOT logging.
- `/dna/test/setMaxTheta true|false`
  For cosine-distributed GPS sources, automatically constrain `maxtheta` so sampled trajectories intersect the slab footprint.

These commands are used extensively in the distributed examples and are part of the intended user interface of the release.

## Output Products

The executable writes Geant4 analysis ntuples to ROOT. Depending on the selected logging mode and run configuration, the output may contain:

- `step`
- `event`
- `track`
- `config`

In general:

- `minimal` logging is intended for compact depth-profile and deposition studies
- `full` logging is intended for detailed process-level and track-level diagnostics

If `DNA_WVALUE_FILE` is set, the application also appends W-value rows in the form:

```text
E0[eV]  <Nion>  rms(Nion)  W[eV]  W_err[eV]
```

## Distributed Examples

The release includes six curated examples under `examples/`.

### Example 1

- amorphous ice
- minimal logging
- one electron each at `0.1`, `1`, and `10 MeV`
- depth-diagnostic plotting helper

This example is recommended as the minimal validation run for new users.

### Example 2

- amorphous ice
- full logging
- one electron each at `0.1`, `1`, and `10 MeV`
- process-diagnostic plotting helper

This example is recommended when process-by-process inspection is required.

### Example 3

- amorphous ice
- `100` primaries at `0.5 MeV`
- `10` threads
- cosine angular distribution
- geometry-limited `maxtheta`

This example is intended as a compact analog of the Europa-style source configuration.

### Examples 4, 5, and 6

- Example 4: amorphous-ice W-value scan
- Example 5: hexagonal-ice W-value scan
- Example 6: water W-value scan

These examples reproduce the current W-value workflow distributed with the release and include plotting helpers for the resulting text output.

For detailed run instructions, consult:

- `examples/README.md`
- the individual `README.md` file in each example directory

## Quick Verification Workflow

After building the executable and installing the distributed cross sections into `G4LEDATA/dna`, a minimal verification sequence is:

```bash
./examples/example_1/run_example_1.sh --binary "./build/dnaphysics"
python3 examples/example_1/plot_example_1_depth_diagnostics.py
```

This should produce:

- `examples/example_1/outputs/example_1.root`
- `examples/example_1/outputs/example_1_depth_diagnostics.png`

A process-level verification run may then be performed with:

```bash
./examples/example_2/run_example_2.sh --binary "./build/dnaphysics"
python3 examples/example_2/plot_example_2_process_diagnostics.py
```

## Python Utilities

### Plotting And Diagnostics

Under `python_scripts/plotting/`, the release provides:

- `plot_depth_diagnostics.py`
  depth-profile and deposited-energy diagnostics
- `plot_diagnostics_all_processes.py`
  process-resolved diagnostics from ROOT outputs
- `plot_elastic_cross_sections.py`
  inspection plots for bundled elastic, vibrational, and attachment datasets
- `plot_stopping_power_xs.py`
  stopping-power and W-value style plots from tabulated cross sections

These plotting scripts may obtain reference files from either the repository-local `cross_sections/` directory or `G4LEDATA/dna`, depending on the script and invocation mode.

### Cross-Section Generation And Inspection

Under `python_scripts/physics_ice/`, the release provides:

- `generate_vibExc_cumulative_dat.py`
  regenerates the vibrational differential and cumulative sampling tables from Michaud tabular inputs
- `generate_ice_cross_sections.py`
  generates phase-specific inelastic excitation and ionization products for water ice
- `generate_blended_elastic_dat.py`
  prepares blended elastic datasets
- `plot_ice_cross_sections.py`
  visualizes generated ice cross sections

The tabular input files used by these scripts are located in `tabular/`.

## Distributed Cross-Section Families

The release already includes the principal `.dat` tables required by the present custom models. The central physics extension introduced here is the electronic inelastic treatment of water ice, distributed separately for amorphous and hexagonal phases. The remaining bundled families provide the supporting transport description used by the current release.

Representative files include:

### Vibrational Excitation

- `sigma_excitationvib_e_michaud.dat`
- `sigmadiff_cumulated_excitationvib_e_michaud_hp.dat`

### Attachment

- `sigma_attachment_e_michaud.dat`

### Elastic Scattering

- `sigma_elastic_e_michaud.dat`
- `sigma_elastic_e_elsepa_muffin.dat`
- `sigmadiff_cumulated_elastic_e_elsepa_muffin.dat`
- multiple low-energy, high-energy, and corrected Michaud-ELSEPA variants

### Amorphous-Ice Excitation And Ionization

- `sigma_excitation_e_amorphous_ice_emfietzoglou_kyriakou.dat`
- `sigma_ionisation_e_amorphous_ice_emfietzoglou_kyriakou.dat`

### Hexagonal-Ice Excitation And Ionization

- `sigma_excitation_e_hexagonal_ice_emfietzoglou_kyriakou.dat`
- `sigma_ionisation_e_hexagonal_ice_emfietzoglou_kyriakou.dat`

Differential `sigmadiff_*` files are included where angular or energy-loss sampling requires more information than a total cross section alone.

## Europa Energy-Library Workflow

The release includes a Europa energy-library generation workflow under `python_scripts/europa/`.

The primary entry point is:

- `generate_europa_energy_library.py`

This workflow:

- constructs a shared global energy grid
- generates one Geant4 macro per global energy
- writes offline weighting tables that map loc
