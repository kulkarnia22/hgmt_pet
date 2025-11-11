import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter

hits = [0.67, 0.69, 0.71, 0.73]
tau  = [5, 4, 3, 2]

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
ax.plot(tau, hits, marker='o', linewidth=2)

ax.set_xlabel("Tau (Mil)", fontsize=16)
ax.set_ylabel("Pore Hits Per First Scatters", fontsize=16)
ax.set_title("Detector Efficiency vs Tau for Kapton LMCP", fontsize=18)

# Keep your zoom but give ticks that actually land inside it
ax.set_xlim(2, 5)
ax.set_ylim(0.67, 0.73)

# x: ticks at each sample
ax.xaxis.set_major_locator(MultipleLocator(1))      # 0.67, 0.68, ..., 0.73
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# y: nice dense ticks across the small range
ax.yaxis.set_major_locator(MultipleLocator(0.01))      # 0.67, 0.68, ..., 0.73
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax.tick_params(which='major', length=8, width=1.5, labelsize=12)
ax.tick_params(which='minor', length=4, width=1.0)

fig.tight_layout()
plt.savefig("plots/eff_vs_tau.png")
plt.show()





