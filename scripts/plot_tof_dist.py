import numpy as np
import matplotlib.pyplot as plt

# ---- file containing angles ----
scattered_tof_file = "data/scattered_tof.data"
non_scattered_tof_file = "data/non_scattered_tof.data"

# ---- load tof diff vals ----
scattered_tof = np.fromfile(scattered_tof_file, dtype=np.float64)
non_scattered_tof = np.fromfile(non_scattered_tof_file, dtype=np.float64)

#getting max and min
max_scattered = np.max(scattered_tof)
min_scattered = np.min(scattered_tof)

max_non_scattered = np.max(non_scattered_tof)
min_non_scattered = np.min(non_scattered_tof)

# ---- clip to range of interest ----
"""MAX_TOF = 4.0  # ns
scattered_tof     = scattered_tof[scattered_tof <= MAX_TOF]
non_scattered_tof = non_scattered_tof[non_scattered_tof <= MAX_TOF]"""

print(f"Scattered events after clip:     {len(scattered_tof):,}")
print(f"Non-scattered events after clip: {len(non_scattered_tof):,}")

# ---- number of bins ----
n_bins = 400

# ---- histogram on shared range ----
scattered_counts, scattered_edges = np.histogram(
    scattered_tof, bins=n_bins, range=(-4, 4))
non_scattered_counts, non_scattered_edges = np.histogram(
    non_scattered_tof, bins=n_bins, range=(-4, 4))

# ---- bin centers for line plot ----
scattered_centers     = 0.5 * (scattered_edges[:-1]     + scattered_edges[1:])
non_scattered_centers = 0.5 * (non_scattered_edges[:-1] + non_scattered_edges[1:])

# ---- plot ----
plt.figure()
plt.plot(scattered_centers, scattered_counts,
         color="red", label="scattered")
plt.plot(non_scattered_centers, non_scattered_counts,
         color="black", label="non scattered")
plt.xlabel("tof diff (ns)")
plt.ylabel("Counts")
plt.title("TOF Distribution In Patient v.s Non In Patient")
plt.legend()
plt.savefig("plots/scattered_tof_plot.png")
plt.show()