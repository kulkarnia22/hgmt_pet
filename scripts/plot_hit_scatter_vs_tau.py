import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

thickness_mm = [2, 3, 4, 5]

hit_first_scatter_eff = [0.51, 0.56, 0.61, 0.65]
hit_second_scatter_eff = [0.33, 0.36, 0.41, 0.46]
hit_third_scatter_eff = [0.18, 0.20, 0.24, 0.29]

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
ax.plot(thickness_mm, hit_first_scatter_eff, "o-", color="red", markersize=7, linewidth=1.5)
ax.plot(thickness_mm, hit_second_scatter_eff, "^-", color="green", markersize=7, linewidth=1.5)
ax.plot(thickness_mm, hit_third_scatter_eff, "s-", color="blue", markersize=7, linewidth=1.5)

# Labels
ax.set_xlabel(r"Thickness Tau (mils)", fontsize = 14)
ax.set_ylabel(r"Hit Efficieny$(\frac{hits}{scatters})$", fontsize = 14)
ax.set_title("Hits per Scatter vs Tau for T = 1 inch")

# Minor ticks
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))

plt.tight_layout()
plt.savefig("plots/hit_scatter_vs_thickness.png", dpi=300)
plt.show()
