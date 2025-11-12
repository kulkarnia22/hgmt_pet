import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# --- read binary integers ---
data = np.fromfile("data/num_pores_crossed.data", dtype=np.int32)
max_val = data.max()
#data = data[data > 1]

# --- histogram with one bin per integer value ---
bins = np.arange(data.min(), data.max() + 2) - 0.5
counts, edges = np.histogram(data, bins=bins)
centers = (edges[:-1] + edges[1:]) / 2

print(f"Max value: {max_val}")

# --- APS-style defaults ---
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

# --- plot ---
fig, ax = plt.subplots()

ax.plot(centers, counts, 'o', markersize=6, color='red')
ax.set_xlabel("Number of pores crossed")
ax.set_ylabel("Count")
ax.set_title("Distribution of pores crossed by first hits")

# --- inward ticks on all sides ---
ax.tick_params(
    which='both',
    direction='in',     # all ticks go in
    top=True, right=True
)

# --- tick locators ---
ax.xaxis.set_major_locator(MultipleLocator(1))           # one tick per integer
ax.xaxis.set_minor_locator(AutoMinorLocator(2))          # 2 minor ticks per interval
ax.yaxis.set_minor_locator(AutoMinorLocator(2))


plt.tight_layout()
plt.savefig("plots/pores_crossed_dist.png", dpi=300)
plt.show()


