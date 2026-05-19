from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import PyMieScatt as ps
from scattnlay import scattnlay


C0 = 299_792_458.0
EPSILON0 = 8.854_187_812_8e-12
PEC_REFRACTIVE_INDEX = 1000.0 + 1000.0j
THEORY_FREQUENCY_POINTS = 2001


@dataclass(frozen=True)
class ShellMaterial:
    core_radius_m: float
    shell_thickness_m: float
    eps_r_shell: float
    mu_r_shell: float
    sigma_shell_s_per_m: float

    @property
    def outer_radius_m(self) -> float:
        return self.core_radius_m + self.shell_thickness_m


@dataclass(frozen=True)
class DgtdBackscatterData:
    theta_rad: np.ndarray
    phi_rad: np.ndarray
    frequency_hz: np.ndarray
    rcs_m2: np.ndarray
    normalization_term: np.ndarray


def parse_material_file(material_path: Path) -> ShellMaterial:
    key_map = {
        "a": "core_radius_m",
        "d": "shell_thickness_m",
        "\u03b5r": "eps_r_shell",
        "\u03bcr": "mu_r_shell",
        "\u03c3": "sigma_shell_s_per_m",
    }
    values: dict[str, float] = {}

    for raw_line in material_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or "=" not in line:
            continue

        raw_key, raw_value = line.split("=", maxsplit=1)
        key = raw_key.strip().rstrip(",")
        number_match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", raw_value)
        if number_match is None or key not in key_map:
            continue

        values[key_map[key]] = float(number_match.group(0))

    required = (
        "core_radius_m",
        "shell_thickness_m",
        "eps_r_shell",
        "mu_r_shell",
        "sigma_shell_s_per_m",
    )
    missing = [name for name in required if name not in values]
    if missing:
        raise ValueError(f"Missing material parameters in {material_path}: {', '.join(missing)}")

    material = ShellMaterial(**values)
    if material.core_radius_m <= 0.0:
        raise ValueError("Core radius must be positive.")
    if material.shell_thickness_m <= 0.0:
        raise ValueError("Shell thickness must be positive.")
    if not np.isclose(material.mu_r_shell, 1.0):
        raise ValueError("This case assumes a non-magnetic shell (mu_r = 1).")

    return material


def load_dgtd_backscatter(data_path: Path) -> DgtdBackscatterData:
    data = np.loadtxt(data_path, skiprows=1)
    if data.ndim != 2 or data.shape[1] != 5:
        raise ValueError(f"Unexpected DGTD data layout in {data_path}.")

    return DgtdBackscatterData(
        theta_rad=np.asarray(data[:, 0], dtype=np.float64),
        phi_rad=np.asarray(data[:, 1], dtype=np.float64),
        frequency_hz=np.asarray(data[:, 2], dtype=np.float64),
        rcs_m2=np.asarray(data[:, 3], dtype=np.float64),
        normalization_term=np.asarray(data[:, 4], dtype=np.float64),
    )


def validate_dgtd_inputs(dgtd_g1: DgtdBackscatterData, dgtd_g2: DgtdBackscatterData) -> None:
    if dgtd_g1.frequency_hz.shape != dgtd_g2.frequency_hz.shape:
        raise RuntimeError("DGTD frequency grids do not have the same number of samples.")
    if not np.allclose(dgtd_g1.frequency_hz, dgtd_g2.frequency_hz):
        raise RuntimeError("DGTD frequency grids do not match between g1 and g2.")
    if not np.allclose(dgtd_g1.theta_rad, np.pi, atol=1e-5):
        raise RuntimeError("DGTD g1 is not a monostatic backscatter cut (theta != pi).")
    if not np.allclose(dgtd_g2.theta_rad, np.pi, atol=1e-5):
        raise RuntimeError("DGTD g2 is not a monostatic backscatter cut (theta != pi).")
    if not np.allclose(dgtd_g1.phi_rad, 0.0, atol=1e-12):
        raise RuntimeError("DGTD g1 is not a phi = 0 cut.")
    if not np.allclose(dgtd_g2.phi_rad, 0.0, atol=1e-12):
        raise RuntimeError("DGTD g2 is not a phi = 0 cut.")

    wavelength_m = C0 / dgtd_g1.frequency_hz
    if not np.allclose(dgtd_g1.normalization_term, wavelength_m**2, rtol=2e-3, atol=0.0):
        raise RuntimeError("DGTD g1 normalization_term does not match the expected free-space wavelength squared.")
    if not np.allclose(dgtd_g2.normalization_term, wavelength_m**2, rtol=2e-3, atol=0.0):
        raise RuntimeError("DGTD g2 normalization_term does not match the expected free-space wavelength squared.")


def shell_refractive_index(
    frequency_hz: np.ndarray, eps_r_shell: float, sigma_shell_s_per_m: float
) -> np.ndarray:
    omega = 2.0 * np.pi * frequency_hz
    eps_shell_complex = eps_r_shell + 1j * sigma_shell_s_per_m / (omega * EPSILON0)
    m_shell = np.sqrt(eps_shell_complex)
    return np.where(np.imag(m_shell) < 0.0, np.conj(m_shell), m_shell)


def compute_pymiescatt_rcs(
    frequency_hz: np.ndarray,
    material: ShellMaterial,
    m_shell: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    wavelength_m = C0 / frequency_hz
    wavelength_nm = wavelength_m * 1e9
    d_core_nm = (2.0 * material.core_radius_m) * 1e9
    d_shell_nm = (2.0 * material.outer_radius_m) * 1e9

    cback_m2 = np.empty_like(frequency_hz)
    for index, (wavelength_point_nm, shell_index) in enumerate(zip(wavelength_nm, m_shell)):
        result = ps.MieQCoreShell(
            PEC_REFRACTIVE_INDEX,
            shell_index,
            wavelength_point_nm,
            d_core_nm,
            d_shell_nm,
            asDict=True,
            asCrossSection=True,
        )
        cback_m2[index] = result["Cback"] * 1e-18

    projected_area_m2 = np.pi * material.outer_radius_m**2
    rcs_normalized = cback_m2 / projected_area_m2
    return cback_m2, rcs_normalized


def compute_scattnlay_rcs(
    frequency_hz: np.ndarray,
    material: ShellMaterial,
    m_shell: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    wavelength_m = C0 / frequency_hz
    x_core = 2.0 * np.pi * material.core_radius_m / wavelength_m
    x_shell = 2.0 * np.pi * material.outer_radius_m / wavelength_m
    x = np.column_stack((x_core, x_shell)).astype(np.float64)
    m = np.column_stack(
        (
            np.ones_like(m_shell, dtype=np.complex128),
            np.asarray(m_shell, dtype=np.complex128),
        )
    )

    _, _, _, _, qbk, _, _, _, _, _ = scattnlay(x, m, pl=0, mp=False)

    projected_area_m2 = np.pi * material.outer_radius_m**2
    rcs_normalized = np.asarray(qbk, dtype=np.float64)
    cback_m2 = rcs_normalized * projected_area_m2
    return cback_m2, rcs_normalized


def save_rcs_csv(
    csv_path: Path,
    frequency_hz: np.ndarray,
    cback_m2: np.ndarray,
    rcs_normalized: np.ndarray,
) -> None:
    wavelength_m = C0 / frequency_hz
    frequency_mhz = frequency_hz / 1e6
    csv_data = np.column_stack(
        (frequency_hz, frequency_mhz, wavelength_m, cback_m2, rcs_normalized)
    )
    np.savetxt(
        csv_path,
        csv_data,
        delimiter=",",
        header="frequency_hz,frequency_mhz,wavelength_m,cback_m2,rcs_normalized",
        comments="",
    )


def validate_rcs_series(name: str, cback_m2: np.ndarray, rcs_normalized: np.ndarray) -> np.ndarray:
    if not np.all(np.isfinite(cback_m2)) or not np.all(np.isfinite(rcs_normalized)):
        raise RuntimeError(f"Non-finite values detected in {name} outputs.")
    if np.any(cback_m2 < 0.0):
        raise RuntimeError(f"Negative cross-sections detected in {name} outputs.")

    rcs_for_plot = np.where(rcs_normalized > 0.0, rcs_normalized, np.nan)
    if not np.any(np.isfinite(rcs_for_plot)):
        raise RuntimeError(f"Cannot render log plot because {name} contains no positive values.")
    return rcs_for_plot


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    dgtd_dir = base_dir / "dgtd_simulations"

    material = parse_material_file(dgtd_dir / "dgtd_material_shell.txt")
    dgtd_g1 = load_dgtd_backscatter(dgtd_dir / "rcs_g1.dat")
    dgtd_g2 = load_dgtd_backscatter(dgtd_dir / "rcs_g2.dat")
    validate_dgtd_inputs(dgtd_g1, dgtd_g2)

    dgtd_frequency_hz = dgtd_g1.frequency_hz
    dgtd_frequency_mhz = dgtd_frequency_hz / 1e6
    theory_frequency_hz = np.linspace(
        float(dgtd_frequency_hz[0]), float(dgtd_frequency_hz[-1]), THEORY_FREQUENCY_POINTS
    )
    theory_frequency_mhz = theory_frequency_hz / 1e6
    m_shell = shell_refractive_index(
        theory_frequency_hz, material.eps_r_shell, material.sigma_shell_s_per_m
    )

    pymiescatt_cback_m2, pymiescatt_rcs_normalized = compute_pymiescatt_rcs(
        theory_frequency_hz, material, m_shell
    )
    scattnlay_cback_m2, scattnlay_rcs_normalized = compute_scattnlay_rcs(
        theory_frequency_hz, material, m_shell
    )

    projected_area_m2 = np.pi * material.outer_radius_m**2
    dgtd_g1_rcs_normalized = dgtd_g1.rcs_m2 / projected_area_m2
    dgtd_g2_rcs_normalized = dgtd_g2.rcs_m2 / projected_area_m2

    pymiescatt_rcs_for_plot = validate_rcs_series(
        "PyMieScatt", pymiescatt_cback_m2, pymiescatt_rcs_normalized
    )
    scattnlay_rcs_for_plot = validate_rcs_series(
        "scattnlay", scattnlay_cback_m2, scattnlay_rcs_normalized
    )
    dgtd_g1_rcs_for_plot = validate_rcs_series(
        "DGTD g1", dgtd_g1.rcs_m2, dgtd_g1_rcs_normalized
    )
    dgtd_g2_rcs_for_plot = validate_rcs_series(
        "DGTD g2", dgtd_g2.rcs_m2, dgtd_g2_rcs_normalized
    )

    pymiescatt_csv = base_dir / "Normalized_RCS_shell_pec_core_pymiescatt_validation_1-400MHz.csv"
    scattnlay_csv = base_dir / "Normalized_RCS_shell_pec_core_scattnlay_validation_1-400MHz.csv"
    save_rcs_csv(
        pymiescatt_csv,
        theory_frequency_hz,
        pymiescatt_cback_m2,
        pymiescatt_rcs_normalized,
    )
    save_rcs_csv(
        scattnlay_csv,
        theory_frequency_hz,
        scattnlay_cback_m2,
        scattnlay_rcs_normalized,
    )

    plt.figure(figsize=(10, 6))
    plt.semilogy(
        theory_frequency_mhz,
        pymiescatt_rcs_for_plot,
        color="tab:blue",
        linewidth=1.5,
        label="PyMieScatt",
    )
    plt.semilogy(
        theory_frequency_mhz,
        scattnlay_rcs_for_plot,
        color="tab:green",
        linewidth=1.5,
        linestyle="--",
        label="scattnlay",
    )
    plt.semilogy(
        dgtd_frequency_mhz,
        dgtd_g1_rcs_for_plot,
        color="tab:orange",
        linewidth=1.0,
        marker="o",
        markersize=3,
        markevery=12,
        label="DGTD g1",
    )
    plt.semilogy(
        dgtd_frequency_mhz,
        dgtd_g2_rcs_for_plot,
        color="tab:red",
        linewidth=1.0,
        marker="s",
        markersize=3,
        markevery=12,
        label="DGTD g2",
    )
    plt.xlabel("Frequency (MHz)")
    plt.ylabel(r"Normalized RCS ($\sigma_{back}/(\pi r_{out}^2)$)")
    plt.title("PEC-Core Dielectric Shell Validation Against DGTD (1 to 400 MHz)")
    plt.xlim(float(theory_frequency_mhz[0]), float(theory_frequency_mhz[-1]))
    plt.grid(True, which="both", linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()

    plot_path = base_dir / "Normalized_RCS_shell_pec_core_validation_1-400MHz_logY.png"
    plt.savefig(plot_path, dpi=300)
    plt.close()

    print(f"PyMieScatt CSV saved: {pymiescatt_csv.name}")
    print(f"scattnlay CSV saved: {scattnlay_csv.name}")
    print(f"Plot saved: {plot_path.name}")


if __name__ == "__main__":
    main()