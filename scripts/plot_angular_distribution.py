import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter

# File paths
all_angles_file = "plots/angular_output.txt"
detected_angles_file = "plots/det_angule_output.txt"

# Load angle data
all_angles = np.loadtxt(all_angles_file)
detected_angles = np.loadtxt(detected_angles_file)

# Histogram parameters
num_bins = 100
min_angle = 0
max_angle = np.pi
bins = np.linspace(min_angle, max_angle, num_bins + 1)
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Histogram counts (not density yet)
all_counts, _ = np.histogram(all_angles, bins=bins)
det_counts, _ = np.histogram(detected_angles, bins=bins)

# Normalize to get probability densities
all_density = all_counts / np.sum(all_counts)
det_density = det_counts / np.sum(all_counts)

# Custom formatter for x-axis ticks
def pi_formatter(x, _):
    frac = x / np.pi
    if np.isclose(frac, 0):
        return "0"
    elif np.isclose(frac, 1):
        return r"$\pi$"
    elif np.isclose(frac, 0.5):
        return r"$\frac{\pi}{2}$"
    else:
        return r"${:.2f}\pi$".format(frac)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(bin_centers, all_density, label="All First-Scatter Angles", marker='o', linestyle='-')
plt.plot(bin_centers, det_density, label="Detected First-Scatter Angles", marker='s', linestyle='-')

# Format x-axis
plt.gca().xaxis.set_major_locator(MultipleLocator(np.pi / 4))
plt.gca().xaxis.set_major_formatter(FuncFormatter(pi_formatter))
plt.xticks(fontsize=14)

plt.xlim(left=0)
plt.ylim(bottom=0)
plt.margins(x=0, y=0)

plt.tick_params(axis='both', direction='in')

# Labels and legend
plt.xlabel("Primary Electron Scattering Angle(radians)", fontsize = 14)
plt.ylabel("Fraction of All Primary Electrons", fontsize = 14)
plt.title("Angular Distributions of First Scatters", fontsize = 14)
plt.tight_layout()
plt.savefig("plots/angular_dist.png")
plt.show()


