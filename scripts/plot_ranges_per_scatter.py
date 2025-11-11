import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_doubles(path):
    # Flat binary file of IEEE-754 float64 (C/C++ 'double')
    return np.fromfile(path, dtype=np.float64)

def make_bins(data, *, bins=None, bin_width=None):
    if bin_width is not None:
        dmin, dmax = np.min(data), np.max(data)
        # pad edges so max lands inside the last bin
        left  = np.floor(dmin / bin_width) * bin_width
        right = np.ceil(dmax / bin_width)  * bin_width
        edges = np.arange(left, right + bin_width, bin_width)
    else:
        # e.g., bins=50 by default (Freedman–Diaconis is another option)
        edges = np.histogram_bin_edges(data, bins=bins if bins else 50)
    centers = 0.5 * (edges[:-1] + edges[1:])
    counts, _ = np.histogram(data, bins=edges)
    return centers, counts, edges

def plot_points(centers, counts, *, label=None, marker='o'):
    plt.plot(centers, counts, linestyle='none', marker=marker, label=label)

# ---- Usage ----
first_path  = Path("data/first_scatter_ranges.data")
second_path = Path("data/second_scatter_ranges.data")  # if you have it
third_path  = Path("data/third_scatter_ranges.data")   # if you have it

first = load_doubles(first_path)
second = load_doubles(second_path)
third = load_doubles(third_path)

# Choose ONE of:
centers, counts, edges = make_bins(third, bin_width=0.05)   # fixed width (e.g., 0.05 mm)
# centers, counts, edges = make_bins(first, bins=60)        # or fixed number of bins

plt.figure(figsize=(8,5))
plot_points(centers, counts, marker='o')

plt.xlabel("Range (mm)", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.title("Electron Range Distribution — Third Scatter", fontsize=15)

# APS-like ticks/look
plt.tick_params(axis='both', direction='in', top=True, right=True, length=8, width=1.4)
plt.minorticks_on()
plt.tick_params(which='minor', direction = "in", length=4, width=1.0)

plt.tight_layout()
plt.savefig("plots/third_scatter_ranges.png")
plt.show()
