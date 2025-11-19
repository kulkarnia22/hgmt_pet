import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ------------------------------
# Parameters
# ------------------------------

filename = "data/scatter_layers.data"
normalize_const = 19869792    # <-- total number of scatters

# ------------------------------
# 1. Load binary doubles
# ------------------------------

layers_raw = np.fromfile(filename, dtype=np.float64)
print("length check = " + str(len(layers_raw)) + "\n")
print(f"Loaded {len(layers_raw)} entries before filtering")

# ------------------------------
# 2. Discard -1 (not in detector)
# ------------------------------

valid_mask = layers_raw >= 0
layers_valid = layers_raw[valid_mask]

print(f"Kept {len(layers_valid)} valid entries, discarded {len(layers_raw) - len(layers_valid)}")

# ------------------------------
# 3. Convert layer values:
#    raw_value + 1 → layer index
# ------------------------------

# Example: 0.0→1, 1.0→2, ..., 11.0→12
layers = layers_valid.astype(int) + 1

# Only keep 1–12
"""final_mask = (layers >= 1) & (layers <= 12)
layers = layers[final_mask]"""

print(f"Final valid layer entries: {len(layers)}")

# ------------------------------
# 4. Count occurrences per layer
# ------------------------------

num_layers = 12
counts = np.zeros(num_layers, dtype=int)

for L in layers:
    counts[L - 1] += 1

# Normalize
norm_counts = counts / normalize_const

print("Counts per layer:", counts)
print("Normalized:", norm_counts)

# ------------------------------
# 5. Make point-style plot
# ------------------------------

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
ax.set_ylabel("Normalized counts")
ax.set_title("Scatter distribution by layer")

# Minor ticks
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))

plt.tight_layout()
plt.savefig("plots/layer_scatter_distribution.png", dpi=300)
plt.show()

