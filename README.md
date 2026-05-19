
# HIREXSR_py
A python port of THACO IDL analysis scripts for the Alcator C-Mod HIREX Doppler spectrometry data.
The original files can be found in: /usr/local/cmod/idl/HIREXSR/
Some additional functionality has been added in to estimate Z_eff, also ported from IDL, and pull electron density and temperature directly from the Thomson laser data.


## Installation
The project is set up to run with uv: use ```uv sync``` to drop into the correct environment, and then uv run for the codes of your choice
You can perminantly install it to the current virtual environment with ```uv pip install .```
To verify that the code is working properly, type ```uv run pytest```

## Usage
See the jupyter notebook for run examples on pulling and plotting the line integrated or inverted data for $T_i, T_e, n_e, n_i$ and $v_\phi$ or checking a list of shots to see if good data exists.
Check the ```CMOD_XICIS_database.xlsx``` provided with this package for a partial list of "good" $T_i, v_\phi$ profiles

Most users only need this core functionality (loading and quality-checking HIREX-SR toroidal rotation and profile data for specific impurity lines).

## Optional Drift Submodule
An additional submodule is available at `HIREXSR_py.drifts` for comparing HIREX-SR toroidal rotation against diamagnetic drift-based corrections.

Key entry points:
- `HIREXSR_py.drifts.compute.load_profiles_for_shot`
- `HIREXSR_py.drifts.compute.load_equilibrium_for_shot`
- `HIREXSR_py.drifts.compute.compute_diamagnetic_drift_frequencies`
- `HIREXSR_py.drifts.plotting.plot_diamagnetic_vs_q_times`
- `HIREXSR_py.drifts.calc_diamagnetic_drift` (CLI)

The drifts submodule includes a local `cmod_helpers` implementation (`openTree`, `currentShot`, `YAG`) so it does not depend on `Synthetic_Mirnov/C-Mod/get_Cmod_Data.py`.

Example:

```python
from HIREXSR_py.drifts.compute import (
	load_profiles_for_shot,
	load_equilibrium_for_shot,
	compute_diamagnetic_drift_frequencies,
)

profiles = load_profiles_for_shot(1120906030, line=2, tht=0)
equilibrium = load_equilibrium_for_shot(1120906030, tree="efit20")
result = compute_diamagnetic_drift_frequencies(
	profiles=profiles,
	equilibrium=equilibrium,
	ion_charge_state=1.0,
	selected_times_s=[0.6, 1.0, 1.3],
)
```

CLI usage:

```bash
python -m HIREXSR_py.drifts.calc_diamagnetic_drift --shot 1120906030 --plot-times 0.6 1.0 1.3
```

Headless-safe usage (save figures without opening interactive windows):

```bash
python -m HIREXSR_py.drifts.calc_diamagnetic_drift --doSave /path/to/output/ --no-show
```

Notes:
- `--no-show`: disables interactive display while still creating/saving figures.
- `--headless`: forces the non-interactive `Agg` backend (also auto-selected when no display is detected).

Backward compatibility:
- `Synthetic_Mirnov/C-Mod/calc_diamagnetic_drift.py` remains as a thin wrapper that forwards to the new CLI module.

## Features
The code features a number of switches to account for data being missing on some MDS nodes, or nodes potentially being renamed.
Users may need to manually check a number of lines to find good data

## Configuration
Project is set up to run within a uv envionrment. If the user does not wish to do this, the only potentially nonstanard package is mdsthin, the python port of the original MDS+ code for data access.
