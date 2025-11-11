import numpy as np
import matplotlib.pyplot as plt

# Load energy values from text files
all_energies = np.loadtxt("data/energy_dist.txt")
detected_energies = np.loadtxt("data/energy_det.txt")

"""print("detected energies > 250: " + str(np.sum(detected_energies > 250)/len(detected_energies)))
print("all energies > 250: " + str(np.sum(all_energies > 250)/len(all_energies)))""" 

# Define histogram bin edges and centers
min_energy = 0
max_energy = 520
num_bins = 50
bins = np.linspace(min_energy, max_energy, num_bins + 1)
bin_centers = 0.5 * (bins[1:] + bins[:-1])

# Compute histograms
all_counts, _ = np.histogram(all_energies, bins=bins)
det_counts, _ = np.histogram(detected_energies, bins=bins)

# Normalize counts to represent efficiency (fraction of total)
all_efficiency = all_counts / len(all_energies)
det_efficiency = det_counts / len(all_energies)

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(bin_centers, all_efficiency, marker='o', label='All Primary Electrons')
plt.plot(bin_centers, det_efficiency, marker='s', label='Detected Primary Electrons')

plt.xlim(left=0)
plt.ylim(bottom=0)
plt.margins(x=0, y=0)

plt.tick_params(axis='both', direction='in')

plt.xlabel('Primary Electron Energy (keV)', fontsize = 14)
plt.ylabel('Fraction of all Primary Electrons', fontsize = 14)
plt.title('Energy Distribution of First Scattered Electrons', fontsize = 14)
plt.tight_layout()
plt.savefig("plots/energy_distribution.png")
plt.show()