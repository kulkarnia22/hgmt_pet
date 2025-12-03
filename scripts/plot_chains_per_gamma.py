import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter

tau  = np.array([5, 4, 3, 2, 1, 0])
chains_per_gamma = np.array([0.68, 0.72, 0.76, 0.79, 0.79 + (0.79 - 0.76), 0.79 + (0.79 - 0.76)*2])

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
ax.set_ylabel(r"Compton Chains per Gamma($\frac{Compton \ Chains}{Gammas}$)", fontsize=16)
ax.set_title("Compton Chains per Gamma", fontsize=18)

# --- Axes limits ---
#ax.set_xlim(0.0, 5.1)
ax.set_ylim(0.0, 1.00)   # accommodates full data range


ax.tick_params(which='major', length=8, width=1.5, labelsize=12)
ax.tick_params(which='minor', length=4, width=1.0)

fig.tight_layout()
plt.savefig("plots/chains_per_gamma.png", dpi=300)
plt.show()
