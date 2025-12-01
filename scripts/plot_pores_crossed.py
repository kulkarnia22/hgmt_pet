import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ------------------------------
# Parameters
# ------------------------------

filename = "data/num_pores_crossed.data"
normalize_const = 2000000      # you choose this

max_pores_in_file = 13        # data may go up to 13
max_k_to_plot = 6             # we only care about >=1..>=6 pores

# ------------------------------
# 1. Load binary integers
# ------------------------------

pores = np.fromfile(filename, dtype=np.int32)
print(f"Loaded {len(pores)} entries")

if len(pores) == 0:
    raise ValueError("File appears empty or incorrect format.")

# Optional: discard 0-pore events if you don't consider them "in detector"
pores = pores[pores > 0]
print(f"Kept {len(pores)} events with >=1 pore crossed")

# ------------------------------
# 2. Exact counts by pore number
#    counts_exact[k] = # hits that crossed exactly k pores
# ------------------------------

max_pores = min(pores.max(), max_pores_in_file)
counts_exact = np.bincount(pores, minlength=max_pores + 1)
# counts_exact[0] is for value 0, which we discarded; ignore index 0

print("Exact counts by pores crossed (k=1..max_pores):")
for k in range(1, max_pores + 1):
    print(f"k={k}: {counts_exact[k]}")

# ------------------------------
# 3. Cumulative counts: >=k pores
#    counts_at_least[k] = sum_{j>=k} counts_exact[j]
# ------------------------------

# Make a copy up to max_pores and compute reverse cumsum
exact = counts_exact[:max_pores + 1]  # indices 0..max_pores
cumulative_from_end = exact[::-1].cumsum()[::-1]
# cumulative_from_end[k] = sum_{j>=k} exact[j]

counts_at_least = cumulative_from_end

print("\nCounts for 'at least k pores crossed':")
for k in range(1, max_k_to_plot + 1):
    print(f">= {k}: {counts_at_least[k]}")

# ------------------------------
# 4. Prepare data for plotting (k = 1..6)
# ------------------------------

ks = np.arange(1, max_k_to_plot + 1)
y_vals = counts_at_least[1:max_k_to_plot + 1] / normalize_const

# ------------------------------
# 5. Plot: point-style, APS-ish formatting
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

ax.plot(ks, y_vals, "o", markersize=7, color="black")
ax.set_xlabel("Minimum Number of Pores Crossed")
ax.set_ylabel(r"$\frac{hits}{gamma}$")
ax.set_title("Hits Per Gamma vs The Minimum Number of Pores Crossed")

ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))

plt.tight_layout()
plt.savefig("plots/num_pores_crossed_cumulative.png", dpi=300)
plt.show()



