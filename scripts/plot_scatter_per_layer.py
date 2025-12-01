import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

filename = "data/scatter_layers.data"
normalize_const = 19869792    # <-- total number of scatters

layers_raw = np.fromfile(filename, dtype=np.float64)
print("length check = " + str(len(layers_raw)) + "\n")
print(f"Loaded {len(layers_raw)} entries before filtering")

valid_mask = layers_raw >= 0
layers_valid = layers_raw[valid_mask]

print(f"Kept {len(layers_valid)} valid entries, discarded {len(layers_raw) - len(layers_valid)}")

layers = layers_valid.astype(int) + 1

num_layers = 12
counts = np.zeros(num_layers, dtype=int)

for L in layers:
    counts[L - 1] += 1

# Normalize
norm_counts = counts / normalize_const

print("Counts per layer:", counts)
print("Normalized:", norm_counts)

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


layers_x = np.arange(1, num_layers + 1)
ax.plot(layers_x, norm_counts, "o", markersize=7, color="black")

ax.set_xlabel("Layer number")
ax.set_ylabel("Fraction of Total Scatters Per Layer")
ax.set_title("Scatter distribution by layer")

ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))

plt.tight_layout()
plt.savefig("plots/layer_scatter_distribution.png", dpi=300)
plt.show()

