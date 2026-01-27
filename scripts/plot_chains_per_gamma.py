"""import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter

tau  = np.array([5, 4, 3, 2, 1])
chains_per_gamma = np.array([0.68, 0.72, 0.76, 0.79, 0.79 + (0.79 - 0.76)])

m, b = np.polyfit(tau, chains_per_gamma, 1)   # slope m, intercept b
x_line = np.linspace(0.5, tau.max(), 200)   # includes x = 0
y_line = m * x_line + b

plt.rcParams.update({
    "figure.figsize": (8, 5),
    "axes.linewidth": 1.5,
    "font.size": 14,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

fig, ax = plt.subplots()
ax.plot(tau, chains_per_gamma, marker='o', ms = 10, color="black")
ax.plot(x_line, y_line, linestyle="-", color="black")


ax.set_xlabel("Tau (mil)", fontsize=16)
ax.set_ylabel(r"Compton Chains per Gamma($\frac{Compton \ Chains}{Gammas}$)", fontsize=16)
ax.set_title("Compton Chains per Gamma", fontsize=18)

# --- Axes limits ---
ax.set_xlim(0.0, 5.1)
ax.set_ylim(0.0, 1.00)   # accommodates full data range


ax.tick_params(which='major', length=8, width=1.5, labelsize=12)
ax.tick_params(which='minor', length=4, width=1.0)

fig.tight_layout()
plt.savefig("plots/chains_per_gamma.png", dpi=300)
plt.show()"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


def plot_prettier(dpi=150, fontsize=17):
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

tau  = np.array([5, 4, 3, 2, 1])

chains_per_gamma = np.array([0.68, 0.72, 0.76, 0.79, 0.79 + (0.79 - 0.76)])

ex_fig = plt.figure()
ex_plt = ex_fig.add_subplot()
ex_fig.subplots_adjust(top=0.88,bottom=0.13,left=0.17,right=0.9,hspace=0.2,wspace=0.2)

ex_plt.plot(tau, chains_per_gamma, 'o-', label=r"$ e^{-x} $", ms=10, color="black")

ex_plt.set_xlim(left=0)
ex_plt.set_ylim(bottom=0)
ex_plt.set_ylim(top=1.0)

ex_plt.set_xlabel("Tau(mils)")
ex_plt.set_ylabel(r"Compton Chains per Gamma")
ex_plt.set_title("Chains Per Gamma")

add_minor_ticks(ex_plt, 5)
"""ex_plt.set_yscale("log")
ex_plt.set_ylim(1e-2, 1)"""

plt.savefig("plots/chains_per_gamma.png")


plt.show()
