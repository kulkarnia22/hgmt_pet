import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Example lists (replace with your actual data)
thickness_mm = [2, 3, 4, 5]
efficiency   = [0.16, 0.14, 0.12, 0.11]

hit_first_scatter_eff = [0.058, 0.054, 0.049, 0.045]
hit_second_scatter_eff = [0.038, 0.034, 0.030, 0.027]
hit_third_scatter_eff = [0.023, 0.019, 0.016, 0.014]

# APS-like styling
plt.rcParams.update({
    "figure.figsize": (8, 5),
    "axes.linewidth": 1.5,
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 6,
    "xtick.minor.size": 3,
    "ytick.major.size": 6,
    "ytick.minor.size": 3,
})

fig, ax = plt.subplots()

# Point plot
ax.plot(thickness_mm, efficiency, "o-", color="black", markersize=7, linewidth=1.5)

# Labels
ax.set_xlabel(r"Thickness Tau (mils)", fontsize = 14)
ax.set_ylabel(r"Hit Efficieny$(\frac{hits}{scatters})$", fontsize = 14)
ax.set_title("Hits per Scatter vs Tau for T = 1 inch")

# Minor ticks
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))

plt.tight_layout()
plt.savefig("plots/efficiency_vs_thickness.png", dpi=300)
plt.show()
