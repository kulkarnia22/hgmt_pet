import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter

tau  = np.array([5, 4, 3, 2])
chains_per_gamma = np.array([0.68, 0.72, 0.76, 0.79])

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
ax.plot(tau, chains_per_gamma, marker='o', ms = 10, linewidth=2)

ax.set_xlabel("Tau (mil)", fontsize=16)
ax.set_ylabel(r"Compton Chains per Gamma", fontsize=16)
ax.set_title("Compton Chains per 511 keV Gamma", fontsize=18)

# --- Axes limits ---
ax.set_xlim(2, 5)
ax.set_ylim(0.65, 0.79)   # accommodates full data range

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
plt.savefig("plots/chains_per_gamma.png", dpi=300)
plt.show()
