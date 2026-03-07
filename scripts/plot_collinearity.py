import numpy as np
import matplotlib.pyplot as plt

# ---- file containing angles ----
filename = "data/collinearity_angles.data"

# ---- load angles ----
angles = np.fromfile(filename, dtype=np.float64)

# ---- determine range ----
min_angle = np.min(angles)
max_angle = np.max(angles)

mask = (angles >= 178.0) & (angles <= 180.0)
angles_zoom = angles[mask]
print("Number of angles in range:", len(angles_zoom))

print(f"Min angle: {min_angle}")
print(f"Max angle: {max_angle}")

# ---- number of bins ----
n_bins = 200

# ---- histogram ----
#counts, bin_edges = np.histogram(angles, bins=n_bins, range=(min_angle, max_angle))
counts, bin_edges = np.histogram(angles_zoom, bins=400, range=(178.0, 180.0))
# bin centers for line plot
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
# ---- plot ----
plt.figure()
#plt.yscale("log")
plt.plot(bin_centers, counts)
plt.xlabel("Angle")
plt.ylabel("Counts")
plt.title("Angle Distribution")
plt.savefig("plots/collinearity_plot.png")
plt.show()