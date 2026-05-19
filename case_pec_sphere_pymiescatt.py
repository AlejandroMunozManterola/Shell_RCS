import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps

# 1. Target Geometry
radius_m = 1.0
diameter_nm = 2.0 * radius_m * 1e9  # PyMieScatt takes diameter!

# 2. Emulate PEC in standard Mie Theory
# A massive complex refractive index ensures perfect reflection
m_pec = 1000 + 1000j 

# 3. Frequency Sweep: 10 MHz to 5 GHz
# Using logspace to evenly distribute points across multiple decades
freq_hz = np.logspace(7, np.log10(100e9), 800)
wavelengths_m = 3e8 / freq_hz
wavelengths_nm = wavelengths_m * 1e9

# 4. Calculate Normalized RCS (Q_back)
q_back_array = np.zeros(len(freq_hz))

print("Running canonical PEC sphere validation...")
for i, wl_nm in enumerate(wavelengths_nm):
    # ps.MieQ handles homogeneous spheres (no shell)
    # Returns: [qext, qsca, qabs, g, qpr, qback, qratio]
    mie_results = ps.MieQ(m_pec, wl_nm, diameter_nm)
    q_back_array[i] = mie_results[5]  # Q_back

# 5. Plotting setup
fig, ax1 = plt.subplots(figsize=(10, 6))

# Main plot: Q_back vs Frequency
ax1.loglog(freq_hz, q_back_array, color='red', linewidth=1.5)
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel(r"Normalized RCS ($\sigma / \pi r^2$)")
ax1.set_title("Canonical PEC Sphere Validation ($r = 1$ m)")
ax1.grid(True, which="both", ls="--", alpha=0.6)

# Create a secondary X-axis at the top to show the electrical size (ka)
# ka = 2 * pi * r / wavelength
ax2 = ax1.twiny()
ax2.set_xscale('log')
ka_min = (2 * np.pi * radius_m) / (3e8 / 1e7)
ka_max = (2 * np.pi * radius_m) / (3e8 / 5e9)
ax2.set_xlim([ka_min, ka_max])
ax2.set_xlabel(r"Electrical Size ($ka = 2\pi r / \lambda$)")

plt.tight_layout()
plt.savefig("PEC_Sphere_Validation.png", dpi=300)
plt.show()