# Shell RCS Cases

This repository contains a small set of Python validation and comparison scripts for monostatic backscatter from spherical targets. The current workspace includes:

- a 2 to 18 GHz dielectric shell comparison between PyMieScatt and scattnlay
- a 1 to 400 MHz PEC-core shell validation against DGTD data
- canonical PEC-sphere validation plots from both PyMieScatt and scattnlay

## Repository Layout

- `case_shell_normalized_rcs_pymiescatt_2GHz_18GHz.py`
  Computes normalized backscatter for the dielectric shell from 2 to 18 GHz with PyMieScatt.
- `case_shell_normalized_rcs_scattnlay_2GHz_18GHz.py`
  Computes the same dielectric shell case with scattnlay.
- `case_shell_normalized_rcs_compare_2GHz_18GHz.py`
  Reads the two 2 to 18 GHz CSV files and generates a combined comparison plot.
- `case_shell_pec_core_dgtd_validation_1MHz_400MHz.py`
  Parses the DGTD material and backscatter files, evaluates the PEC-core shell with both Mie solvers on a denser theory grid, and overlays those curves with DGTD `g1` and `g2`.
- `case_pec_sphere_pymiescatt.py`
  Runs a canonical PEC-sphere validation with PyMieScatt from 10 MHz to 5 GHz.
- `case_pec_sphere_scattnlay.py`
  Runs the same PEC-sphere validation with scattnlay from 10 MHz to 5 GHz.
- `dgtd_simulations/`
  Contains the DGTD shell material definition and the `rcs_g1.dat` and `rcs_g2.dat` validation datasets.
- `external/scattnlay/`
  Git submodule checkout of the upstream scattnlay project.

## Cases

### Dielectric Shell Comparison, 2 to 18 GHz

These three scripts compare PyMieScatt and scattnlay for the same two-layer spherical shell.

- Core radius: 1.5 m
- Shell thickness: 3e-4 m
- Outer radius: 1.5003 m
- Relative permittivity: 4.0
- Loss tangent: 0.005
- Frequency sweep: 2 to 18 GHz with 1001 linearly spaced samples

Outputs:

- `Normalized_RCS_shell_pymiescatt_2-18GHz.csv`
- `Normalized_RCS_shell_pymiescatt_2-18GHz_logY.png`
- `Normalized_RCS_shell_scattnlay_2-18GHz.csv`
- `Normalized_RCS_shell_scattnlay_2-18GHz_logY.png`
- `Normalized_RCS_shell_compare_2-18GHz_logY.png`

The 2 to 18 GHz CSV files contain:

```text
frequency_hz,frequency_ghz,wavelength_m,cback_m2,rcs_normalized
```

### PEC-Core Shell Validation Against DGTD, 1 to 400 MHz

This validation case reads [dgtd_simulations/dgtd_material_shell.txt](dgtd_simulations/dgtd_material_shell.txt) and the DGTD backscatter datasets [dgtd_simulations/rcs_g1.dat](dgtd_simulations/rcs_g1.dat) and [dgtd_simulations/rcs_g2.dat](dgtd_simulations/rcs_g2.dat).

Current model setup:

- Core radius `a`: 0.95 m
- Shell thickness `d`: 0.05 m
- Outer radius: 1.0 m
- Relative permittivity: 4
- Relative permeability: 1
- Conductivity: 3 S/m
- DGTD frequency range: 1 MHz to 400 MHz
- Theory frequency grid: 2001 linearly spaced points over the same limits

Outputs:

- `Normalized_RCS_shell_pec_core_pymiescatt_validation_1-400MHz.csv`
- `Normalized_RCS_shell_pec_core_scattnlay_validation_1-400MHz.csv`
- `Normalized_RCS_shell_pec_core_validation_1-400MHz_logY.png`

The DGTD validation CSV files contain:

```text
frequency_hz,frequency_mhz,wavelength_m,cback_m2,rcs_normalized
```

### Canonical PEC Sphere Validation

The PEC sphere scripts generate log-log reference plots for a sphere of radius 1 m.

- Frequency sweep: 10 MHz to 5 GHz with 800 logarithmically spaced samples
- PyMieScatt output: `PEC_Sphere_Validation_pymiescatt.png`
- scattnlay output: `PEC_Sphere_Validation_scattnlay.png`

## Normalization

For the shell cases, normalized RCS is reported as

$$
\sigma_{norm} = \frac{\sigma_{back}}{\pi r_{out}^2}
$$

where $r_{out}$ is the outer radius of the shell.

For the canonical PEC-sphere validation, the plotted quantity is the usual backscatter efficiency

$$
Q_{back} = \frac{\sigma_{back}}{\pi r^2}
$$

for a sphere of radius $r = 1$ m.

## Requirements

Install the Python packages used by the scripts:

```bash
pip install numpy matplotlib PyMieScatt
```

The scattnlay cases also require a working local `scattnlay` Python installation.

## scattnlay Submodule

At the repository root, `./external/scattnlay/` is kept as a Git submodule.

After cloning the repository, initialize submodules with:

```bash
git submodule update --init --recursive
```

For most users, that checkout is still not enough to run the scattnlay scripts. You will usually need to build and install the Python package from `./external/scattnlay/` first.

Use the upstream instructions in `./external/scattnlay/README.md` for the build steps.

## Usage

From the repository root, run:

```bash
python case_shell_normalized_rcs_pymiescatt_2GHz_18GHz.py
python case_shell_normalized_rcs_scattnlay_2GHz_18GHz.py
python case_shell_normalized_rcs_compare_2GHz_18GHz.py
python case_shell_pec_core_dgtd_validation_1MHz_400MHz.py
python case_pec_sphere_pymiescatt.py
python case_pec_sphere_scattnlay.py
```

If you are running without a display, set Matplotlib to a non-interactive backend first:

```powershell
$env:MPLBACKEND='Agg'
```

## Notes

- The 2 to 18 GHz comparison script expects both shell CSV files to already exist.
- The DGTD validation script reads its inputs directly from `./dgtd_simulations/`.
- The shell validation and comparison plots use linear frequency on the x-axis and logarithmic normalized RCS on the y-axis.
- The PEC sphere validation scripts use logarithmic scales on both axes and include a top axis for electrical size.
