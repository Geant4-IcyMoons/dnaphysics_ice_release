# Geant4-IcyMoons

<p align="center">
  <img src="icymoons_logo.png" alt="Geant4-IcyMoons logo" />
</p>

Geant4-IcyMoons is a Geant4-DNA-based application for electron transport in icy media. This release packages:

- a standalone `dnaphysics` executable
- custom low-energy electron cross sections for ice
- example macros and run scripts
- Python diagnostics and cross-section generation utilities
- a Europa energy-library generator

The release is intended for users who want to:

- run single-event or batch electron transport in amorphous ice, hexagonal ice, or water
- inspect deposited-energy and process diagnostics from ROOT outputs
- run W-value scans
- regenerate or inspect the bundled custom cross-section tables
- prepare Europa-style monoenergetic macro libraries for offline reweighting

## What This Release Adds To Geant4-DNA

This project extends standard Geant4-DNA electron transport with ice-specific data and models.

At a high level:

- `water` uses an explicit Geant4-DNA water physics list
- `ice_am` and `ice_hex` use the custom ice physics list
- the ice models combine Geant4-DNA style transport with bundled custom tables under `cross_sections/`

The release includes:

- Michaud-based vibrational excitation data
- Michaud-based low-energy attachment data
- multiple elastic options and angular tables
- amorphous-ice and hexagonal-ice excitation and ionization tables generated from the Emfietzoglou-Kyriakou model chain

## Repository Layout

- `examples/`: six ready-to-run user examples
- `cross_sections/`: bundled `.dat` tables used by the ice models
- `src/`: Geant4 model and application sources
- `include/`: headers
- `python_scripts/plotting/`: ROOT diagnostics and plotting utilities
- `python_scripts/physics_ice/`: cross-section generation and inspection tools
- `python_scripts/europa/`: Europa macro-library generation tools
- `tabular/`: source tabular data used to build some of the custom cross sections

## Physics Modes

The executable chooses the transport setup from the `DNA_PHYSICS` environment variable:

- `DNA_PHYSICS=ice_am`: amorphous ice
- `DNA_PHYSICS=ice_hex`: hexagonal ice
- `DNA_PHYSICS=water`: liquid water

If `DNA_PHYSICS` is omitted, the executable defaults to `ice_hex`.

These modes not only change density. They select different physics-list branches and, for the ice modes, the phase-aware custom excitation and ionization tables.

## Before You Build

You need:

- Geant4 11.x with Geant4-DNA available
- CMake 3.16 or newer
- a C++17-capable compiler
- Geant4 low-energy data, especially `G4EMLOW...`

This release expects Geant4 data through the standard `G4LEDATA` environment variable.

## Install The Bundled Cross Sections

The C++ ice models resolve their `.dat` files through `G4LEDATA/dna`, not directly from the repository. Before running the executable, make sure the bundled files in `cross_sections/` are available under your Geant4 low-energy data tree.

Typical workflow:

```bash
export G4LEDATA=/path/to/Geant4/data/G4EMLOW8.6.1
mkdir -p "${G4LEDATA}/dna"
cp cross_sections/* "${G4LEDATA}/dna/"
```

If you prefer not to copy, symlinks also work on most systems:

```bash
export G4LEDATA=/path/to/Geant4/data/G4EMLOW8.6.1
mkdir -p "${G4LEDATA}/dna"
for f in cross_sections/*; do
  ln -sf "$(pwd)/$f" "${G4LEDATA}/dna/$(basename "$f")"
done
```

This step is important. If the custom files are not visible in `G4LEDATA/dna`, the corresponding ice processes cannot load their data.

## Build

From the repository root:

```bash
cmake -S . -B build -DGeant4_DIR=/path/to/geant4-install/lib/Geant4-11.x
cmake --build build -j
```

This creates:

- `build/dnaphysics`: the executable

The build also copies top-level `.mac`, `.in`, and `.C` files into `build/` so they can be run from the build directory if desired.

## How To Run

Batch mode:

```bash
DNA_PHYSICS=ice_am ./build/dnaphysics path/to/macro.mac 4
```

Arguments:

- first argument: macro file
- second argument: number of Geant4 threads

Interactive mode:

```bash
./build/dnaphysics
```

With no macro argument, the application starts an interactive UI session and executes `vis.mac`.

## Common Runtime Environment Variables

The most useful runtime variables are:

- `DNA_PHYSICS`
  Selects `ice_am`, `ice_hex`, or `water`.
- `DNA_ROOT_BASENAME`
  Base name of the output ROOT file. Default is `dna`.
- `DNA_WVALUE_FILE`
  If set, append W-value rows to the given text file.
- `DNA_NTUPLE_FILES`
  Number of reduced ROOT ntuple files when running with ntuple merging.
- `DNA_NTUPLE_MERGE`
  Enable or disable merged ntuple output in multi-threaded runs.
- `DNA_NTUPLE_STRINGS`
  Include string columns such as process and model names where supported.
- `DNA_KEEP_OLD_ROOT`
  Keep older ROOT outputs instead of deleting matching files at run start.
- `DNA_ROOT_SPLIT_EVENTS`
  Rotate output files after a fixed number of events.
- `DNA_ROOT_MAX_MB`
  Rotate output files after a size threshold.

For water-only studies, the water physics list also supports optional feature toggles such as `DNA_WATER_ENABLE_HIGH_EM`, `DNA_WATER_ENABLE_ELASTIC`, `DNA_WATER_ENABLE_EXCITATION`, and `DNA_WATER_ENABLE_IONISATION`.

## Useful Macro Commands

The release adds several convenience commands under `/dna/test/`:

- `/dna/test/setMat NAME`
  Set the Geant4 material.
- `/dna/test/setMatDens NAME VALUE UNIT`
  Override density for a named material.
- `/dna/test/setIceSize X Y Z UNIT`
  Set slab dimensions as full lengths.
- `/dna/test/setLogMode minimal|full`
  Choose reduced or detailed ROOT logging.
- `/dna/test/setMaxTheta true|false`
  For cosine GPS sources, automatically constrain `maxtheta` so rays still intersect the slab footprint.

These are the main commands used in the bundled examples.

## Output Files

The executable writes ROOT output through Geant4 analysis ntuples. Depending on logging mode, the run may contain:

- `step` tree
- `event` tree
- `track` tree
- `config` tree

In broad terms:

- `minimal` logging keeps the output compact and is suited for depth-profile work
- `full` logging keeps richer process and track metadata and is better for process-by-process diagnostics

If `DNA_WVALUE_FILE` is set, a text file is also written in append mode with rows like:

```text
E0[eV]  <Nion>  rms(Nion)  W[eV]  W_err[eV]
```

## Bundled Examples

The release includes six maintained examples under `examples/`.

### Example 1: Minimal Logging, Three Energies

- `examples/example_1/`
- amorphous ice
- one electron each at `0.1`, `1`, and `10 MeV`
- minimal logging
- depth-diagnostic plotting helper included

Use this when you want a lightweight first run and a compact ROOT file.

### Example 2: Full Logging, Three Energies

- `examples/example_2/`
- amorphous ice
- one electron each at `0.1`, `1`, and `10 MeV`
- full logging
- process-diagnostic plotting helper included

Use this when you want richer process accounting per step.

### Example 3: Europa-Style Batch Run

- `examples/example_3/`
- amorphous ice
- `100` primaries at `0.5 MeV`
- `10` threads
- cosine angular distribution
- geometry-limited `maxtheta`

Use this as the closest small example to the production Europa-style source setup.

### Examples 4-6: W-Value Scans

- `examples/example_4/`: amorphous ice
- `examples/example_5/`: hexagonal ice
- `examples/example_6/`: water

These examples mirror the current W-value workflow and include plotting helpers for the resulting text files.

For exact instructions, use the README in each example directory or the example index in `examples/README.md`.

## Quick Start

After building and installing the data files into `G4LEDATA/dna`, a good first test is:

```bash
./examples/example_1/run_example_1.sh --binary "./build/dnaphysics"
python3 examples/example_1/plot_example_1_depth_diagnostics.py
```

Then inspect:

- `examples/example_1/outputs/example_1.root`
- `examples/example_1/outputs/example_1_depth_diagnostics.png`

For a process-level view:

```bash
./examples/example_2/run_example_2.sh --binary "./build/dnaphysics"
python3 examples/example_2/plot_example_2_process_diagnostics.py
```

## Python Utilities

### Plotting ROOT Diagnostics

Under `python_scripts/plotting/`:

- `plot_depth_diagnostics.py`
  depth profiles and deposited-energy diagnostics
- `plot_diagnostics_all_processes.py`
  process-by-process cross-section and angle diagnostics from ROOT
- `plot_elastic_cross_sections.py`
  bundled elastic, vibrational, and attachment reference plots
- `plot_stopping_power_xs.py`
  stopping-power and W-value style plots from tables

These plotting scripts resolve reference files either from the repository `cross_sections/` directory or from `G4LEDATA/dna`, depending on the script.

### Regenerating Or Inspecting Cross Sections

Under `python_scripts/physics_ice/`:

- `generate_vibExc_cumulative_dat.py`
  rebuilds the vibrational differential and cumulative sampling tables from Michaud tabular inputs
- `generate_ice_cross_sections.py`
  builds the phase-specific ice excitation and ionisation products
- `generate_blended_elastic_dat.py`
  prepares blended elastic datasets
- `plot_ice_cross_sections.py`
  visualizes generated ice cross sections

The source data used by these scripts live in `tabular/`.

## Cross-Section Files Included In This Release

The repository already ships the main tables needed by the current release.

Examples:

- vibrational excitation
  - `sigma_excitationvib_e_michaud.dat`
  - `sigmadiff_cumulated_excitationvib_e_michaud_hp.dat`
- attachment
  - `sigma_attachment_e_michaud.dat`
- elastic
  - `sigma_elastic_e_michaud.dat`
  - `sigma_elastic_e_elsepa_muffin.dat`
  - `sigmadiff_cumulated_elastic_e_elsepa_muffin.dat`
  - several low/high and corrected Michaud-ELSEPA variants
- amorphous-ice excitation and ionisation
  - `sigma_excitation_e_amorphous_ice_emfietzoglou_kyriakou.dat`
  - `sigma_ionisation_e_amorphous_ice_emfietzoglou_kyriakou.dat`
- hexagonal-ice excitation and ionisation
  - `sigma_excitation_e_hexagonal_ice_emfietzoglou_kyriakou.dat`
  - `sigma_ionisation_e_hexagonal_ice_emfietzoglou_kyriakou.dat`

The differential `sigmadiff_*` files are used where angular or energy-loss sampling needs more than a total cross section.

## Europa Energy-Library Workflow

The release also includes the Europa macro-library generator under `python_scripts/europa/`.

Main entry point:

- `generate_europa_energy_library.py`

This workflow:

- builds a shared global energy grid
- writes one Geant4 macro per energy
- writes reweighting tables that map local Europa surface cells onto the shared energy library

Use `python_scripts/europa/README.md` for the detailed parameter guide for that workflow.

## Troubleshooting

### The executable cannot find a custom data file

Check:

- `echo $G4LEDATA`
- whether the required file exists under `"$G4LEDATA/dna/"`

The release C++ models use `G4FindDataDir("G4LEDATA")`, so the data must be visible through that environment variable.

### The run starts but uses the wrong phase

Check:

- `DNA_PHYSICS`

Valid values are:

- `ice_am`
- `ice_hex`
- `water`

### The output ROOT file gets overwritten

Set:

- `DNA_ROOT_BASENAME` to a unique name
- or `DNA_KEEP_OLD_ROOT=1` to preserve matching older outputs

### Multi-threaded output creates many files or rotates unexpectedly

Check:

- `DNA_NTUPLE_FILES`
- `DNA_NTUPLE_MERGE`
- `DNA_ROOT_SPLIT_EVENTS`
- `DNA_ROOT_MAX_MB`

These control merged output and file rotation.

## References

- Geant4-DNA collaboration papers and documentation
- Michaud et al. low-energy electron interaction data for water ice
- Emfietzoglou and Kyriakou model chain for ice excitation and ionisation tables

## License

This project extends and depends on Geant4. Use and redistribution remain subject to the Geant4 Software License and the licenses associated with any bundled data sources.
