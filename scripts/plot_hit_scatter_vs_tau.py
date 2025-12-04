import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

thickness_mm = np.array([5, 4, 3, 2, 1])

hit_first_scatter_eff = [0.51, 0.56, 0.61, 0.65, 0.65 + (0.65-0.61)]
hit_second_scatter_eff = [0.33, 0.36, 0.41, 0.46, 0.46 + (0.46-0.41)]
hit_third_scatter_eff = [0.18, 0.20, 0.24, 0.29, 0.29 + (0.29-0.24)]

first_hit_second_scatter_eff = [0.083, 0.085, 0.086, 0.087, 0.087 + (0.087-0.086)]
first_hit_third_scatter_eff = [0.032, 0.031, 0.030, 0.028, 0.028 + (0.028-0.030)]

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

m_first_scatter, b_first_scatter = np.polyfit(thickness_mm, hit_first_scatter_eff, 1)   # slope m, intercept b
x_line = np.linspace(0.5, thickness_mm.max(), 200)   # includes x = 0
first_scatter_line = m_first_scatter * x_line + b_first_scatter

m_second_scatter, b_second_scatter = np.polyfit(thickness_mm, hit_second_scatter_eff, 1)   # slope m, intercept b
x_line = np.linspace(0.5, thickness_mm.max(), 200)   # includes x = 0
second_scatter_line = m_second_scatter * x_line + b_second_scatter

m_first_hit_second_scatter, b_first_hit_second_scatter = np.polyfit(thickness_mm, first_hit_second_scatter_eff, 1)   # slope m, intercept b
x_line = np.linspace(0.5, thickness_mm.max(), 200)   # includes x = 0
first_hit_second_scatter_line = m_first_hit_second_scatter * x_line + b_first_hit_second_scatter

m_third_scatter, b_third_scatter = np.polyfit(thickness_mm, hit_third_scatter_eff, 1)   # slope m, intercept b
x_line = np.linspace(0.5, thickness_mm.max(), 200)   # includes x = 0
third_scatter_line = m_third_scatter * x_line + b_third_scatter

m_first_hit_third_scatter, b_first_hit_third_scatter = np.polyfit(thickness_mm, first_hit_third_scatter_eff, 1)   # slope m, intercept b
x_line = np.linspace(0.5, thickness_mm.max(), 200)   # includes x = 0
first_hit_third_scatter_line = m_first_hit_third_scatter * x_line + b_first_hit_third_scatter

# Point plot
ax.plot(thickness_mm, hit_first_scatter_eff, "o", color="red", markersize=7)
ax.plot(x_line, first_scatter_line, linestyle='-', color='red')

ax.plot(thickness_mm, hit_second_scatter_eff, "^", color="green", markersize=7)
ax.plot(x_line, second_scatter_line, linestyle='-', color='green')

#ax.plot(thickness_mm, first_hit_second_scatter_eff, "^", color="green", markerfacecolor="none", markersize="7")
#ax.plot(x_line, first_hit_second_scatter_line, linestyle='--', color='green')

ax.plot(thickness_mm, hit_third_scatter_eff, "s", color="blue", markersize=7)
ax.plot(x_line, third_scatter_line, linestyle='-', color='blue')

#ax.plot(thickness_mm, first_hit_third_scatter_eff, "s", color="blue", markerfacecolor="none", markersize="7")
#ax.plot(x_line, first_hit_third_scatter_line, linestyle='--', color='blue')

# Labels
ax.set_xlabel(r"Thickness Tau (mils)", fontsize = 14)
ax.set_ylabel(r"Hit Efficieny$(\frac{hits}{scatters})$", fontsize = 14)
ax.set_title("Hits per Scatter vs Tau")

# Minor ticks
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax.set_ylim(1e-2, 1)
ax.set_ylim(0, 1)
#ax.set_yscale("log")
ax.set_xlim(0, 5.1)

plt.tight_layout()
plt.savefig("plots/hit_scatter_vs_thickness.png", dpi=300)
plt.show()
