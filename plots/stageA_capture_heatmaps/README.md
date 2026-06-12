# Stage A absorption figures

This folder contains a standalone plotting workflow for Stage A neutron
absorption results.

The script reads simulation data from:

- `Input/stageA/*/neutron_capture_absorption/*.csv`
- `Input/stageA/*/neutron_capture_positions/*.csv`

Optional experimental absorption points are read from the local spreadsheet:

- `plots/stageA_capture_heatmaps/input/neutron_absorption.xlsx`

The `input/` and `output/` directories are intentionally not tracked by Git.

## Run from the project root

Run the script from the same repository that contains `plots/` and `Input/`:

```bash
python3 plots/stageA_capture_heatmaps/generate_stageA_capture_heatmaps.py
```

If you activate a virtual environment, still run the command from the project
root. For example, if the repository is `~/g4work/M1`, do not run the command
from `~/g4work/MC`.

## Recommended isolated environment

On shared servers, system Python packages often mix incompatible NumPy and
Matplotlib builds. Use a project-local virtual environment:

```bash
cd ~/g4work/M1
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r plots/stageA_capture_heatmaps/requirements.txt
python3 plots/stageA_capture_heatmaps/generate_stageA_capture_heatmaps.py
```

The requirements intentionally pin `numpy<2` because older system Matplotlib
builds are commonly compiled against NumPy 1.x. Installing both packages inside
the same virtual environment avoids the `_ARRAY_API` and
`numpy.core.multiarray failed to import` errors.

## Generated figures

The script writes figures and derived CSV tables to:

- `plots/stageA_capture_heatmaps/output/`

Current main figures:

- `neutron_absorption_rate_vs_thickness.png`
- `neutron_absorption_rate_heatmap.png`
- `mean_absorption_depth_vs_thickness.png`
- `normalized_capture_depth_distribution_panels.png`
- `absolute_capture_depth_distribution_panels.png`
