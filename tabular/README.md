# Tabular Source Data

This directory contains source tables and intermediate tabular products used in constructing selected distributed cross-section files.

Examples include:

- Michaud low-energy electron interaction tables
- optical and finite-q ice dielectric model tables
- auxiliary tabular datasets used by the Python generation scripts

The release transport executable does not read these files directly at runtime. They are included for reproducibility of the generated `cross_sections/*.dat` products.
