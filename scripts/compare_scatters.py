import matplotlib.pyplot as plt
import numpy as np

def read_histogram(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    lines = lines[3:]  # Skip the first three metadata lines

    bin_centers = []
    counts = []

    for line in lines:
        try:
            range_part, count_part = line.split(":")
            low, high = map(float, range_part.split("-"))
            count = float(count_part.strip())
            center = high
            bin_centers.append(center)
            counts.append(count)
        except ValueError:
            continue  # Ignore malformed lines

    return np.array(bin_centers), np.array(counts)

# Set your paths
file_detected = "plots/first_scatter_detected_pointvac.txt"
file_all = "plots/first_scatter_detector_pointvac.txt"

# Read histogram data
centers_detected, counts_detected = read_histogram(file_detected)
centers_all, counts_all = read_histogram(file_all)

# Compute errors
errors_detected = np.sqrt(counts_detected)/np.sum(counts_all)
errors_all = np.sqrt(counts_all)/np.sum(counts_all)

#Scatter Plot
plt.errorbar(centers_all, counts_all/2000000, yerr=errors_all, fmt='-o', label='All First Scatters', capsize=3)
plt.errorbar(centers_detected, counts_detected/2000000, yerr=errors_detected, fmt='-o', label='Detected First Scatters', capsize=3)

plt.xlim(left=0)
plt.ylim(bottom=0)
plt.margins(x=0, y=0)

plt.tick_params(axis='both', direction='in')

plt.xlabel("Detector Layer", fontsize = 12)
plt.ylabel(r"Fraction of First Scatters($\frac{scatters}{events}$)", fontsize = 12)
plt.xticks(ticks=np.arange(1, 13), labels=[str(i) for i in range(1, 13)])
plt.title("Distribution of First Scatters by Detector Layer", fontsize = 12)
plt.tight_layout()

#last minute details 
# Make y-axis go to the top
plt.ylim(0, np.array(counts_all/2000000).max() * 1.05)  # add a bit of headroom (5%)
# Force ticks to span 0 → max
yticks = np.linspace(0, np.array(counts_all/2000000).max() * 1.05, 6)  # choose number of ticks
# Round each tick to 2 significant figures
def round_sig(x, sig=2):
    return float(f"{x:.{sig}g}")

ytick_labels = [round_sig(y, 2) for y in yticks]

plt.yticks(yticks, ytick_labels)
plt.savefig("plots/first_scatter_comparison.png")
plt.show()

"""#Integral Plots
integral_detected = [np.sum(counts_detected[:i+1]/2000000) for i in range(len(counts_detected))]
integral_all = [np.sum(counts_all[:i+1]/2000000) for i in range(len(counts_all))]

integral_detected_err = []
sum = 0
for error in errors_detected:
    sum += error**2
    integral_detected_err.append(np.sqrt(sum))


integral_all_err = []
sum = 0
for error in errors_all:
    sum += error**2
    integral_all_err.append(np.sqrt(sum))

plt.errorbar(centers_all, integral_all, integral_all_err, fmt = "-o", capsize = 2)
plt.errorbar(centers_detected, integral_detected, integral_detected_err, fmt = "-o", capsize = 2)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.margins(x=0, y=0)
plt.tick_params(axis='both', direction='in')
plt.xlabel("Cumulative Detector Layer", fontsize = 12)
plt.ylabel(r"Summed Fraction of First Scatters($\frac{scatters}{events}$)", fontsize = 12)
plt.title("Integrated Distribution of First Scatters by Cumulative Detector Layer", fontsize = 12)
plt.xticks(ticks=np.arange(1, 13), labels=[str(i) for i in range(1, 13)])

#last minute details 
# Make y-axis go to the top
plt.ylim(0, np.array(integral_all).max() * 1.05)  # add a bit of headroom (5%)
# Force ticks to span 0 → max
yticks = np.linspace(0, np.array(integral_all).max() * 1.05, 6)  # choose number of ticks
# Round each tick to 2 significant figures
def round_sig(x, sig=2):
    return float(f"{x:.{sig}g}")

ytick_labels = [round_sig(y, 2) for y in yticks]

plt.yticks(yticks, ytick_labels)

plt.savefig("plots/integral_scatter_comparison.png")
plt.show()"""
