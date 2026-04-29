import numpy as np
import matplotlib.pyplot as plt

# ---- file containing angles ----
scattered_impact_file = "data/scattered_impact_param.data"
non_scattered_impact_file = "data/non_scattered_impact_param.data"

# ---- load tof diff vals ----
scattered_impact = np.fromfile(scattered_impact_file, dtype=np.float64)
non_scattered_impact = np.fromfile(non_scattered_impact_file, dtype=np.float64)

#getting max and min
max_scattered = np.max(scattered_impact)
min_scattered = np.min(scattered_impact)

max_non_scattered = np.max(non_scattered_impact)
min_non_scattered = np.min(non_scattered_impact)

# ---- number of bins ----
n_bins = 400

# ---- histogram on shared range ----
scattered_counts, scattered_edges = np.histogram(
    scattered_impact, bins=n_bins, range=(min_scattered, max_scattered))
non_scattered_counts, non_scattered_edges = np.histogram(
    non_scattered_impact, bins=n_bins, range=(min_non_scattered, max_non_scattered))

# ---- bin centers for line plot ----
scattered_centers     = 0.5 * (scattered_edges[:-1]     + scattered_edges[1:])
non_scattered_centers = 0.5 * (non_scattered_edges[:-1] + non_scattered_edges[1:])

# ---- plot ----
plt.figure()
plt.yscale("log")
plt.plot(scattered_centers, scattered_counts,
         color="red", label="scattered")
plt.plot(non_scattered_centers, non_scattered_counts,
         color="black", label="non scattered")
plt.xlabel("impact parameter(cm)")
plt.ylabel("Counts")
plt.title("Impact Parameter Distribution In Patient v.s Non In Patient")
plt.legend()
plt.savefig("plots/scattered_impact_logplot.png")
plt.show()