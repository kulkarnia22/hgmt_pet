import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter
from matplotlib.ticker import LogLocator, NullFormatter

tau  = np.array([5, 4, 3, 2, 1, 0])
lors_per_annihilation = np.array([0.51, 0.57, 0.63, 0.69, 0.69 + (0.69 - 0.63), 0.69 + (0.69 - 0.63)*2])
lors11_per_annihilation = np.array([0.17, 0.19, 0.21, 0.23, 0.23 + (0.23 - 0.21), 0.23 + (0.23 - 0.21)*2])
lors12_per_annihilation = np.array([0.13, 0.14, 0.15, 0.16, 0.16 + (0.16 - 0.15), 0.16 + (0.16 - 0.15)*2])
lors22_per_annihilation = np.array([0.024, 0.026, 0.028, 0.029, 0.029 + (0.029 - 0.028), 0.029 + (0.029 - 0.028)*2])

fig, ax = plt.subplots()
ax.plot(tau, lors_per_annihilation, 'o-', ms = 10, linewidth=2)
ax.plot(tau, lors11_per_annihilation, '^-', ms = 10, linewidth=2)
ax.plot(tau, lors12_per_annihilation, 's-', ms = 10, linewidth=2)
ax.plot(tau, lors22_per_annihilation, 'x-', ms = 10, linewidth=2)

ax.set_xlabel("Tau (mil)", fontsize=16)
ax.set_ylabel(r"LORs Per Annhilation $\frac{LORs}{Annihilaitons}$", fontsize=16)
ax.set_title("LORs Per Annihilation v.s Tau", fontsize=18)

# --- Axes limits ---
#ax.set_xlim(0, 5.1)
#ax.set_ylim(np.log10(0.02), np.log10(1.00))   # accommodates full data range
ax.set_yscale("log")
ax.set_ylim(1e-2, 1)   # now you should see 10^-2, 10^-1, 10^0
ax.tick_params(axis="y", which="minor", labelsize=10, direction="in")
ax.tick_params(axis="y", which="major", labelsize=10, direction="in")
ax.tick_params(axis="x", which="major", labelsize=10, direction="in")

fig.tight_layout()
plt.savefig("plots/lors_per_annihil_total.png", dpi=300)
plt.show()
