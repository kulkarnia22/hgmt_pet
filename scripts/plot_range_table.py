import numpy as np
import matplotlib.pyplot as plt

path = "simulation_materials/kapton_ranges.csv"
rho = 1.42  # g/cm^3

data = np.loadtxt(path, skiprows=8)
energy = 1000 * data[:, 0]                 # MeV → keV
range_gcm2 = data[:, 1]
range_mm = (range_gcm2 / rho) * 10.0       # g/cm² → mm

mask = energy <= 511                       # restrict ≤ 511 keV
energy = energy[mask]
range_mm = range_mm[mask]

plt.figure(figsize=(7,5))
plt.plot(energy, range_mm, marker="o")
plt.xlabel("Electron kinetic energy (keV)", fontsize=14)
plt.ylabel("CSDA range (mm)", fontsize=14)
plt.title("Electron CSDA Range in Kapton", fontsize=14)
plt.tick_params(axis='both', direction='in', top=True, right=True, length=8, width=1.4)
plt.minorticks_on()
plt.tick_params(which='minor', direction='in', length=4, width=1.0)
plt.tight_layout()
plt.savefig("plots/range_table_plot.png")
plt.show()

