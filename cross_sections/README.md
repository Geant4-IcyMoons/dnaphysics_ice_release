# Cross-Section Tables

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
