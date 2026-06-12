# Stage A capture heatmaps

This folder contains a small standalone plotting workflow for Stage A neutron
capture results.

Files:

- `generate_stageA_capture_heatmaps.py`
  reads `Input/stageA/*/neutron_capture_positions/*.csv` and writes two heatmaps
  plus the flattened CSV tables used to build them.
- `output/`
  stores the generated PNG and CSV files.

Generated figures:

1. `p_cap_vs_thickness_heatmap.png`
   uses the `neutron_capture_positions` CSVs directly. The script infers the
   sample incident count from the largest `eventID` seen for each ratio. In the
   current dataset that evaluates to `500000`, so
   `P_cap = n_capture_positions / 500000`.
2. `geant4_vs_exponential_model_delta_heatmap.png`
   compares the sampled Geant4 result against
   `P_cap = 1 - exp(-Sigma_eff t)`, where `Sigma_eff` is taken from
   `geant4_macroscopic_cross_sections/geant4_macroscopic_cross_sections.csv`
   at `thickness_um = 1000` and `energy_eV = 0.0253`, using the
   `total_removal_macroscopic_xs_per_cm` column.

Notes:

- The comparison heatmap keeps the full ratio axis, but ratios without a
  matching `geant4_macroscopic_cross_sections.csv` are shown as blank cells.
- The current repository only contains that macroscopic-cross-section input for
  `Input/stageA/1-2/`, so only the `1-2` row is populated in the second figure.

Run:

```bash
python3 plots/stageA_capture_heatmaps/generate_stageA_capture_heatmaps.py
```
