import numpy as np

# ------------------------------
# Parameters
# ------------------------------

filename = "data/angle_inside_pore.data"   # binary file of C doubles
num_bins = 100             # choose however many you want

# ------------------------------
# Read binary file
# ------------------------------

angles = np.fromfile(filename, dtype=np.float64)

# Optional: safety filter
#angles = angles[(angles >= 0.0) & (angles <= np.pi)]

# ------------------------------
# Bin between 0 and pi
# ------------------------------

bins = np.linspace(0.0, np.pi, num_bins + 1)
counts, edges = np.histogram(angles, bins=bins)

# Bin centers (useful for point-style plots)
centers = 0.5 * (edges[:-1] + edges[1:])

# ------------------------------
# Done — these are your results
# ------------------------------

print("Counts:", counts)
print("Bin centers:", centers)
print("Bin edges:", edges)
