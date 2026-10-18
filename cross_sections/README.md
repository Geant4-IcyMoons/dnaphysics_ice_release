# Cross-Section Tables

**Regeneration pending:** these bundled DATs predate the generator's
`threshold-kk-moment-budget-v2` repair. They remain historical release products,
not newly corrected tables. Do not mix them with regenerated energy patches.
See [the repair and validation requirements](../python_scripts/physics_ice/CONSISTENCY.md).

This directory contains the `.dat` files distributed with Geant4-IcyMoons.

The custom ice models read these files through the Geant4 low-energy data path, not directly from this repository directory. Before running the executable, copy or link the files into:

```bash
${G4LEDATA}/dna/
```

The distributed files include:

- phase-specific electronic excitation and ionization tables for amorphous and hexagonal ice
- differential `sigmadiff_*` products for sampling
- Michaud-based vibrational-excitation and attachment tables
- elastic and angular-sampling tables used by the release transport models

See the top-level `README.md` for installation commands and runtime notes.
