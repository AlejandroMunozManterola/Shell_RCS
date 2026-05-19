# PyMieScattering Cases

This folder contains a small set of scripts for normalized monostatic backscatter from a thin dielectric spherical shell, plus one PEC sphere reference case.

This README documents only the files that exist in this folder.

## What Is Here

- `case_shell_normalized_rcs_pymiescatt_2GHz_18GHz.py`
  Computes normalized backscatter for the shell from 2 to 18 GHz with PyMieScatt.
- `case_shell_normalized_rcs_scattnlay_2GHz_18GHz.py`
  Computes the same shell case with scattnlay.
- `case_shell_normalized_rcs_compare_2GHz_18GHz.py`
  Reads the two CSV files and generates a comparison plot.
- `case_pec_sphere_pymiescatt.py`
  Runs a PEC sphere reference calculation with PyMieScatt.

## Shell Case

- Core radius: 1.5 m
- Shell thickness: 3e-4 m
- Outer radius: 1.5003 m
- Relative permittivity: 4.0
- Loss tangent: 0.005
- Frequency sweep: 2 to 18 GHz, 1001 linearly spaced points

Normalized RCS is reported as

$$
\sigma_{norm} = \frac{\sigma_{back}}{\pi a^2}
$$

where $a$ is the outer radius.

## Generated Files

Running the shell scripts produces:

- `Normalized_RCS_shell_pymiescatt_2-18GHz.csv`
- `Normalized_RCS_shell_pymiescatt_2-18GHz_logY.png`
- `Normalized_RCS_shell_scattnlay_2-18GHz.csv`
- `Normalized_RCS_shell_scattnlay_2-18GHz_logY.png`
- `Normalized_RCS_shell_compare_2-18GHz_logY.png`

Running the PEC reference script produces:

- `PEC_Sphere_Validation.png`

The shell CSV files contain these columns:

```text
frequency_hz,frequency_ghz,wavelength_m,cback_m2,rcs_normalized
```

## Requirements

Install the Python packages used by the PyMieScatt scripts:

```bash
pip install PyMieScatt numpy matplotlib
```

The scattnlay script also requires a working local `scattnlay` Python installation.

## scattnlay Submodule

At the repository root, `./external/scattnlay/` is kept as a Git submodule.

After cloning the repository, initialize submodules with:

```bash
git submodule update --init --recursive
```

For most users, that checkout is still not enough to run the scattnlay case. You will usually need to compile and install your own local scattnlay Python library from `./external/scattnlay/` first.

Use the upstream instructions in `./external/scattnlay/README.md` for the build steps.

## Usage

From this folder, run:

```bash
python case_shell_normalized_rcs_pymiescatt_2GHz_18GHz.py
python case_shell_normalized_rcs_scattnlay_2GHz_18GHz.py
python case_shell_normalized_rcs_compare_2GHz_18GHz.py
python case_pec_sphere_pymiescatt.py
```

If you are running without a display, set Matplotlib to a non-interactive backend first:

```powershell
$env:MPLBACKEND='Agg'
```

## Notes

- The comparison script expects both shell CSV files to already exist.
- The shell plots use linear frequency on the x-axis and logarithmic normalized RCS on the y-axis.
- The repository root is one level above this folder. The scattnlay source checkout lives there at `./external/scattnlay/`.
