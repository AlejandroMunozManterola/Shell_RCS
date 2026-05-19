from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    pymiescatt_csv = Path("Normalized_RCS_shell_pymiescatt_2-18GHz.csv")
    scattnlay_csv = Path("Normalized_RCS_shell_scattnlay_2-18GHz.csv")

    if not pymiescatt_csv.exists():
        raise FileNotFoundError(f"Missing input CSV: {pymiescatt_csv}")
    if not scattnlay_csv.exists():
        raise FileNotFoundError(f"Missing input CSV: {scattnlay_csv}")

    pymiescatt = np.genfromtxt(pymiescatt_csv, delimiter=",", names=True)
    scattnlay = np.genfromtxt(scattnlay_csv, delimiter=",", names=True)

    if pymiescatt.shape != scattnlay.shape:
        raise RuntimeError("CSV files do not have the same number of samples.")
    if not np.allclose(pymiescatt["frequency_hz"], scattnlay["frequency_hz"]):
        raise RuntimeError("Frequency grids do not match between the two CSV files.")

    pymiescatt_rcs = np.where(
        pymiescatt["rcs_normalized"] > 0.0, pymiescatt["rcs_normalized"], np.nan
    )
    scattnlay_rcs = np.where(
        scattnlay["rcs_normalized"] > 0.0, scattnlay["rcs_normalized"], np.nan
    )

    if not np.any(np.isfinite(pymiescatt_rcs)):
        raise RuntimeError("PyMieScatt comparison data contains no positive values for log plotting.")
    if not np.any(np.isfinite(scattnlay_rcs)):
        raise RuntimeError("scattnlay comparison data contains no positive values for log plotting.")

    frequency_ghz = pymiescatt["frequency_ghz"]

    plt.figure(figsize=(10, 6))
    plt.semilogy(frequency_ghz, pymiescatt_rcs, color="tab:blue", linewidth=1.5, label="PyMieScatt")
    plt.semilogy(frequency_ghz, scattnlay_rcs, color="tab:green", linewidth=1.5, label="scattnlay")
    plt.xlabel("Frequency (GHz)")
    plt.ylabel(r"Normalized RCS ($\sigma_{back}/(\pi a^2)$)")
    plt.title("Normalized Backscatter RCS of Hollow Dielectric Shell (2 to 18 GHz)")
    plt.grid(True, which="both", linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()

    plot_filename = "Normalized_RCS_shell_compare_2-18GHz_logY.png"
    plt.savefig(plot_filename, dpi=300)
    plt.close()

    print(f"Plot saved: {plot_filename}")


if __name__ == "__main__":
    main()