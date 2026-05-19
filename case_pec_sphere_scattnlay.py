import matplotlib.pyplot as plt
import numpy as np
from scattnlay import scattnlay


C0 = 299_792_458.0
RADIUS_M = 1.0
FREQUENCY_START_HZ = 1e7
FREQUENCY_STOP_HZ = 5e9
FREQUENCY_POINTS = 800


def main() -> None:
    freq_hz = np.logspace(np.log10(FREQUENCY_START_HZ), np.log10(FREQUENCY_STOP_HZ), FREQUENCY_POINTS)
    wavelengths_m = C0 / freq_hz

    x = (2.0 * np.pi * RADIUS_M / wavelengths_m).reshape(-1, 1)
    m = np.array([1.0 + 0.0j], dtype=np.complex128)

    print("Running canonical PEC sphere validation with scattnlay...")
    _, _, _, _, q_back_array, _, _, _, _, _ = scattnlay(x, m, pl=0, mp=False)
    q_back_array = np.asarray(q_back_array, dtype=np.float64)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.loglog(freq_hz, q_back_array, color="tab:green", linewidth=1.5)
    ax1.set_xlabel("Frequency (Hz)")
    ax1.set_ylabel(r"Normalized RCS ($\sigma / \pi r^2$)")
    ax1.set_title("Canonical PEC Sphere Validation (scattnlay, $r = 1$ m)")
    ax1.grid(True, which="both", ls="--", alpha=0.6)

    ax2 = ax1.twiny()
    ax2.set_xscale("log")
    ka_min = (2.0 * np.pi * RADIUS_M) / (C0 / FREQUENCY_START_HZ)
    ka_max = (2.0 * np.pi * RADIUS_M) / (C0 / FREQUENCY_STOP_HZ)
    ax2.set_xlim([ka_min, ka_max])
    ax2.set_xlabel(r"Electrical Size ($ka = 2\pi r / \lambda$)")

    plt.tight_layout()
    plt.savefig("PEC_Sphere_Validation_scattnlay.png", dpi=300)
    plt.close()
    print("Plot saved: PEC_Sphere_Validation_scattnlay.png")


if __name__ == "__main__":
    main()