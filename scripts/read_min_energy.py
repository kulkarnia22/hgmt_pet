import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os

# ------------------------------
# 1. Read binary (int, double) pairs
# ------------------------------

filename = "data/min_energy.data"

# Each record: int32 (num pores), float64 (min energy)
dtype = np.dtype([("n_pores", np.int32), ("E_min", np.float64)])

if not os.path.exists(filename):
    raise FileNotFoundError(f"Could not find {filename}")

records = np.fromfile(filename, dtype=dtype)
print(f"Loaded {len(records)} records from {filename}")

if len(records) == 0:
    raise ValueError("File appears to be empty or format is incorrect.")

# ------------------------------
# 2. Split energies by number of pores crossed (1, 2, 3)
# ------------------------------

energies_by_pores = {}
for k in (1, 2, 3):
    mask = (records["n_pores"] == k)
    energies = records["E_min"][mask]
    energies_by_pores[k] = energies
    print(f"Pores crossed = {k}: {len(energies)} entries")

# ------------------------------
# 3. Plot helper: point-style histogram
# ------------------------------

def point_histogram(energies, title, outfile, nbins=50):
    if len(energies) == 0:
        print(f"No data for {title}, skipping plot.")
        return

    # Compute histogram
    emin, emax = energies.min(), energies.max()
    counts, edges = np.histogram(energies, bins=nbins, range=(emin, emax))
    centers = 0.5 * (edges[:-1] + edges[1:])

    # APS-ish styling
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

    # Point-style histogram: counts vs energy bin centers
    ax.plot(centers, counts, "o", markersize=5, color="black")

    ax.set_xlabel(r"Minimum energy (keV)")
    ax.set_ylabel("Count")
    ax.set_title(title)

    # Minor ticks
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    # Optional subtle grid
    #ax.grid(which="major", linestyle="--", alpha=0.4)
    #ax.grid(which="minor", linestyle=":", alpha=0.3)

    plt.tight_layout()
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    plt.savefig(outfile, dpi=300)
    plt.close(fig)
    print(f"Saved {outfile}")

# ------------------------------
# 4. Make plots for 1, 2, 3 pores crossed
# ------------------------------

point_histogram(
    energies_by_pores[1],
    title="Min energy for events crossing 1 pore",
    outfile="plots/min_energy_pores_1.png",
)

point_histogram(
    energies_by_pores[2],
    title="Min energy for events crossing 2 pores",
    outfile="plots/min_energy_pores_2.png",
)

point_histogram(
    energies_by_pores[3],
    title="Min energy for events crossing 3 pores",
    outfile="plots/min_energy_pores_3.png",
)
