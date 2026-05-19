import numpy as np
import matplotlib.pyplot as plt
from scattnlay import scattnlay


def main() -> None:
    core_radius_m = 1.5
    shell_thickness_m = 3e-4
    outer_radius_m = core_radius_m + shell_thickness_m
    eps_r_shell = 4.0
    loss_tangent_shell = 0.005

    frequency_start_hz = 2e9
    frequency_stop_hz = 18e9
    frequency_points = 1001
    frequency_hz = np.linspace(frequency_start_hz, frequency_stop_hz, frequency_points)

    use_multiprecision = False

    if shell_thickness_m <= 0:
        raise ValueError("Shell thickness must be positive.")

    projected_area_m2 = np.pi * outer_radius_m**2

    m_core = 1.0 + 0.0j
    m_shell = np.sqrt(eps_r_shell * (1.0 + 1j * loss_tangent_shell))
    if m_shell.imag < 0:
        m_shell = np.conj(m_shell)

    c0 = 299_792_458.0
    wavelength_m = c0 / frequency_hz
    x_core = 2.0 * np.pi * core_radius_m / wavelength_m
    x_shell = 2.0 * np.pi * outer_radius_m / wavelength_m
    x = np.column_stack((x_core, x_shell)).astype(np.float64)
    m = np.array((m_core, m_shell), dtype=np.complex128)

    _, _, _, _, qbk, _, _, _, _, _ = scattnlay(x, m, mp=use_multiprecision)

    rcs_normalized = np.asarray(qbk, dtype=np.float64)
    cback_m2 = rcs_normalized * projected_area_m2

    if frequency_hz.size != frequency_points:
        raise RuntimeError("Unexpected number of frequency points.")
    if not np.isclose(frequency_hz[0], frequency_start_hz) or not np.isclose(
        frequency_hz[-1], frequency_stop_hz
    ):
        raise RuntimeError("Frequency endpoints do not match requested range.")
    if not np.all(np.isfinite(cback_m2)) or not np.all(np.isfinite(rcs_normalized)):
        raise RuntimeError("Non-finite values detected in computed outputs.")
    if np.any(cback_m2 < 0):
        raise RuntimeError("Negative cross-sections detected, which is non-physical.")

    frequency_ghz = frequency_hz / 1e9
    csv_data = np.column_stack(
        (frequency_hz, frequency_ghz, wavelength_m, cback_m2, rcs_normalized)
    )
    csv_filename = "Normalized_RCS_shell_scattnlay_2-18GHz.csv"
    np.savetxt(
        csv_filename,
        csv_data,
        delimiter=",",
        header="frequency_hz,frequency_ghz,wavelength_m,cback_m2,rcs_normalized",
        comments="",
    )

    rcs_for_plot = np.where(rcs_normalized > 0.0, rcs_normalized, np.nan)
    if not np.any(np.isfinite(rcs_for_plot)):
        raise RuntimeError("Cannot render log plot because all normalized RCS values are non-positive.")

    plt.figure(figsize=(10, 6))
    plt.semilogy(frequency_ghz, rcs_for_plot, color="tab:green", linewidth=1.5)
    plt.xlabel("Frequency (GHz)")
    plt.ylabel(r"Normalized RCS ($\sigma_{back}/(\pi a^2)$)")
    plt.title("Normalized Backscatter RCS of Hollow Dielectric Shell (scattnlay, 2 to 18 GHz)")
    plt.grid(True, which="both", linestyle="--", alpha=0.6)
    plt.tight_layout()

    plot_filename = "Normalized_RCS_shell_scattnlay_2-18GHz_logY.png"
    plt.savefig(plot_filename, dpi=300)
    plt.close()

    print(f"CSV saved: {csv_filename}")
    print(f"Plot saved: {plot_filename}")


if __name__ == "__main__":
    main()