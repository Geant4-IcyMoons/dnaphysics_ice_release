# Plotting And Diagnostics Scripts

This directory contains plotting and diagnostic helpers for Geant4-IcyMoons outputs.

The scripts operate on ROOT outputs, tabulated cross sections, or derived analysis files depending on the workflow. Common tasks include:

- process-level ROOT diagnostics
- depth-profile plotting
- elastic and stopping-power comparisons
- Europa energy-map visualization
- vibrational-excitation angular-distribution plots

Most scripts require `numpy` and `matplotlib`. Scripts that read ROOT files also require `uproot`.
