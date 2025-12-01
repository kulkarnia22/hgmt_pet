import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter

tau  = np.array([5, 4, 3, 2])
lors_per_annihilation = np.array([0.51, 0.57, 0.63, 0.69])
lors11_per_annihilation = np.array([0.17, 0.19, 0.21, 0.23])
lors12_per_annihilation = np.array([0.13, 0.14, 0.15, 0.16])
lors22_per_annihilation = np.array([0.024, 0.026, 0.028, 0.029])

plt.rcParams.update({
    "figure.figsize": (8, 5),
    "axes.linewidth": 1.5,
    "font.size": 14,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

fig, ax = plt.subplots()
ax.plot(tau, lors_per_annihilation, 'o-', ms = 10, linewidth=2)
ax.plot(tau, lors11_per_annihilation, '^-', ms = 10, linewidth=2)
ax.plot(tau, lors12_per_annihilation, 's-', ms = 10, linewidth=2)
ax.plot(tau, lors22_per_annihilation, 'x-', ms = 10, linewidth=2)

ax.set_xlabel("Tau (mil)", fontsize=16)
ax.set_ylabel(r"LORs Per Annhilation $\frac{LORs}{Annihilaitons}$", fontsize=16)
ax.set_title("LORs Per Annihilation v.s Tau", fontsize=18)

# --- Axes limits ---
ax.set_xlim(2, 5)
ax.set_ylim(0.0, 0.70)   # accommodates full data range

# --- X ticks: integers, no decimal formatting ---
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))

# --- Y ticks ---
ax.yaxis.set_major_locator(MultipleLocator(0.02))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

ax.tick_params(which='major', length=8, width=1.5, labelsize=12)
ax.tick_params(which='minor', length=4, width=1.0)

fig.tight_layout()
plt.savefig("plots/lors_per_annihil.png", dpi=300)
plt.show()
