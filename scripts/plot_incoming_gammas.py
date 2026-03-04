"""import numpy as np
import matplotlib.pyplot as plt

# Read the binary file of doubles
file_path = "data/Incoming_Photons.data" 
file_path_first = "data/Incoming_Photons_First.data"
file_path_first_firstid = "data/Incoming_Photons_Firstid.data"
file_path_first_secondid = "data/Incoming_Photons_Secondid.data"
file_path_first_thirdid = "data/Incoming_Photons_Thirdid.data"
file_path_second = "data/Incoming_Photons_Second.data"
file_path_third = "data/Incoming_Photons_Third.data"
file_path_fourth = "data/Incoming_Photons_Fourth.data"
file_path_fifth = "data/Incoming_Photons_Fifth.data"
file_path_sixth = "data/Incoming_Photons_Sixth.data"
energies = np.fromfile(file_path, dtype=np.float64)
energies_first = np.fromfile(file_path_first, dtype=np.float64)
energies_first_firstid = np.fromfile(file_path_first_firstid, dtype=np.float64)
energies_first_secondid = np.fromfile(file_path_first_secondid, dtype=np.float64)
energies_first_thirdid = np.fromfile(file_path_first_thirdid, dtype=np.float64)
energies_second = np.fromfile(file_path_second, dtype=np.float64)
energies_third = np.fromfile(file_path_third, dtype=np.float64)
energies_fourth = np.fromfile(file_path_fourth, dtype=np.float64)
energies_fifth = np.fromfile(file_path_fifth, dtype=np.float64)
energies_sixth = np.fromfile(file_path_sixth, dtype=np.float64)

# Define bin edges: 300 bins between 0 and 520 keV
num_bins = 500
bin_edges = np.linspace(0, 520, num_bins + 1)

# Histogram: counts in each bin
counts, _ = np.histogram(energies, bins=bin_edges)
counts_first, _ = np.histogram(energies_first, bins=bin_edges)
counts_first_firstid, _ = np.histogram(energies_first_firstid, bins=bin_edges)
counts_first_secondid, _ = np.histogram(energies_first_secondid, bins=bin_edges)
counts_first_thirdid, _ = np.histogram(energies_first_thirdid, bins=bin_edges)
counts_second, _ = np.histogram(energies_second, bins=bin_edges)
counts_third, _ = np.histogram(energies_third, bins=bin_edges)
counts_fourth, _ = np.histogram(energies_fourth, bins=bin_edges)
counts_fifth, _ = np.histogram(energies_fifth, bins=bin_edges)
counts_sixth, _ = np.histogram(energies_sixth, bins=bin_edges)

# Bin centers for bar positions
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
bin_width = bin_edges[1] - bin_edges[0]

#Plot
fig, ax = plt.subplots(figsize=(8,6))
eps = 1e-10
#ax.errorbar(bin_centers, counts/10000000, np.sqrt(counts)/10000000, capsize = 2)
ax.plot(bin_centers, counts/20000000 + eps, label = 'all first hits')
#ax.errorbar(bin_centers, counts_first/10000000, np.sqrt(counts_first)/10000000, capsize = 2)
ax.plot(bin_centers, counts_first/20000000 + eps, label = 'i = 1')
#ax.errorbar(bin_centers, counts_first_firstid/10000000, np.sqrt(counts_first_firstid)/10000000, capsize = 2)
#ax.plot(bin_centers, np.log(counts_first_firstid/10000000 + eps))
#ax.errorbar(bin_centers, counts_first_secondid/10000000, np.sqrt(counts_first_secondid)/10000000, capsize = 2)
#ax.plot(bin_centers, np.log(counts_first_secondid/10000000 + eps))
#ax.errorbar(bin_centers, counts_first_thirdid/10000000, np.sqrt(counts_first_thirdid)/10000000, capsize = 2)
#ax.plot(bin_centers, np.log(counts_first_thirdid/10000000 + eps))
ax.plot(bin_centers, counts_second/20000000 + eps, label = 'i = 2')
ax.plot(bin_centers, counts_third/20000000 + eps, label = 'i = 3')
ax.plot(bin_centers, counts_fourth/20000000 + eps, label = 'i = 4')
ax.plot(bin_centers, counts_fifth/20000000 + eps, label = 'i = 5')
ax.plot(bin_centers, counts_sixth/20000000 + eps, label = 'i = 6')
ax.set_xlabel("Energy (keV)", fontsize = 14)
ax.set_ylabel(r"Log(Fraction of Incoming Photons($\frac{photons}{all gammas}$))", fontsize = 14)
ax.set_title(r"Energy Distribution of Photons Resulting in $PC^{i,1}_{N}$", fontsize = 14)
# Inward ticks
ax.tick_params(which='both', direction='in')
ax.set_yscale("log")
# Subticks on the x-axis
plt.minorticks_on()  # turns on subticks
plt.tick_params(axis='x', which='minor', bottom=True)  # show subticks

# Make y-axis go to the top
#plt.ylim(-10, 0)  # add a bit of headroom (5%)
#(np.log(counts/10000000 + eps)).max() * 1.05
# Force ticks to span 0 → max
#yticks = np.linspace(-10, 0, 6)  # choose number of ticks
# Round each tick to 2 significant figures
plt.xlim(0, 520)
def round_sig(x, sig=2):
    return float(f"{x:.{sig}g}")

#ytick_labels = [round_sig(y, 2) for y in yticks]

#plt.yticks(yticks, ytick_labels)
#plt.legend()
ax.legend(
    loc='upper center',
    bbox_to_anchor=(0, 0.8, 0.5, 0.2),   # (x0, y0, width, height) in axes coords
    mode='expand',                      # stretch to fill the bbox width
    ncol=3,                             # put entries in columns to use the width
    frameon=True
)
plt.show()
plt.savefig("plots/Incomming_Photons")"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


def plot_prettier(dpi=150, fontsize=15):
	plt.rcParams['figure.dpi'] = dpi
	plt.rc("savefig", dpi=dpi)
	plt.rc('font', size=fontsize)
	plt.rc('xtick', direction='in')
	plt.rc('ytick', direction='in')
	plt.rc('xtick.major', pad=5)
	plt.rc('xtick.minor', pad=5)
	plt.rc('ytick.major', pad=5)
	plt.rc('ytick.minor', pad=5)
	plt.rc('lines', dotted_pattern = [2., 2.])
	# plt.rc('text', usetex=True)

plot_prettier()

def add_minor_ticks(plot, ticks, bot=True, tp=True, lft=True, rght=True, xticks=None, yticks=None):
	plot.tick_params(which = 'minor', bottom=bot, top=tp, left=lft, right=rght)
	plot.tick_params(bottom=True, top=True, left=True, right=True)
	if (bot or tp):
		if (xticks != None):
			plot.xaxis.set_minor_locator(AutoMinorLocator(int(xticks)))
		else:
			plot.xaxis.set_minor_locator(AutoMinorLocator(int(ticks)))
	if (lft or rght):
		if (yticks != None):
			plot.yaxis.set_minor_locator(AutoMinorLocator(int(yticks)))
		else:
			plot.yaxis.set_minor_locator(AutoMinorLocator(int(ticks)))	

ex_fig = plt.figure()
ex_plt = ex_fig.add_subplot()
ex_fig.subplots_adjust(top=0.88,bottom=0.13,left=0.17,right=0.9,hspace=0.2,wspace=0.2)

file_path = "data/Incoming_Photons.data" 
file_path_first = "data/Incoming_Photons_First.data"
file_path_second = "data/Incoming_Photons_Second.data"
file_path_third = "data/Incoming_Photons_Third.data"
file_path_fourth = "data/Incoming_Photons_Fourth.data"
file_path_fifth = "data/Incoming_Photons_Fifth.data"
file_path_sixth = "data/Incoming_Photons_Sixth.data"

energies = np.fromfile(file_path, dtype=np.float64)
energies_first = np.fromfile(file_path_first, dtype=np.float64)
energies_second = np.fromfile(file_path_second, dtype=np.float64)
energies_third = np.fromfile(file_path_third, dtype=np.float64)
energies_fourth = np.fromfile(file_path_fourth, dtype=np.float64)
energies_fifth = np.fromfile(file_path_fifth, dtype=np.float64)
energies_sixth = np.fromfile(file_path_sixth, dtype=np.float64)

# Define bin edges: 300 bins between 0 and 520 keV
num_bins = 500
bin_edges = np.linspace(0, 520, num_bins + 1)


# Histogram: counts in each bin
counts, _ = np.histogram(energies, bins=bin_edges)
counts_first, _ = np.histogram(energies_first, bins=bin_edges)
counts_second, _ = np.histogram(energies_second, bins=bin_edges)
counts_third, _ = np.histogram(energies_third, bins=bin_edges)
counts_fourth, _ = np.histogram(energies_fourth, bins=bin_edges)
counts_fifth, _ = np.histogram(energies_fifth, bins=bin_edges)
counts_sixth, _ = np.histogram(energies_sixth, bins=bin_edges)

# Bin centers for bar positions
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
bin_width = bin_edges[1] - bin_edges[0]
num_hits = 2*753758

ex_plt.plot(bin_centers, counts/num_hits, '-', label="all hits", ms=10, color="black")
ex_plt.plot(bin_centers, counts_first/num_hits, '-', label="first scattered hits", ms=10, color="red")
ex_plt.plot(bin_centers, counts_second/num_hits, '-', label="second scattered hits", ms=10, color="blue")
ex_plt.plot(bin_centers, counts_third/num_hits, '-', label="third scattered hits", ms=10, color="green")
ex_plt.plot(bin_centers, counts_fourth/num_hits, '-', label="fourth scattered hits", ms=10, color="orange")
ex_plt.plot(bin_centers, counts_fifth/num_hits, '-', label="fifth scattered hits", ms=10, color="pink")
ex_plt.plot(bin_centers, counts_sixth/num_hits, '-', label="sixth scattered hits", ms=10, color="purple")

"""ex_plt.set_xlim(left=0.0)
ex_plt.set_ylim(bottom=0.0)
ex_plt.set_ylim(top=1.0)"""

ex_plt.set_xlabel("Incoming Energy of Scattered Gammas(keV)")
ex_plt.set_ylabel(r"Fraction of Hits Resulting in an LOR")
ex_plt.set_title("Log Plot Effective Energy Resolution")

add_minor_ticks(ex_plt, 5)
ex_plt.set_yscale("log")
plt.legend(fontsize = '8', loc = 'upper left')

plt.savefig("plots/log_incomming_gammas.png")


plt.show()
