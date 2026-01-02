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



example_x = np.linspace(0, 2 * np.e, 100)
example_y = np.exp(-example_x)

tau  = np.array([5, 4, 3, 2, 1])

chains_per_gamma = np.array([0.68, 0.72, 0.76, 0.79, 0.79 + (0.79 - 0.76)])

hit_first_scatter_eff = [0.51, 0.56, 0.61, 0.65, 0.65 + (0.65-0.61)]
hit_second_scatter_eff = [0.33, 0.36, 0.41, 0.46, 0.46 + (0.46-0.41)]
hit_third_scatter_eff = [0.18, 0.20, 0.24, 0.29, 0.29 + (0.29-0.24)]
first_hit_second_scatter_eff = [0.083, 0.085, 0.086, 0.087, 0.087 + (0.087-0.086)]
first_hit_third_scatter_eff = [0.032, 0.031, 0.030, 0.028, 0.028 + (0.028-0.030)]

lors_per_annihilation = np.array([0.51, 0.57, 0.63, 0.69, 0.69 + (0.69 - 0.63)])
lors11_per_annihilation = np.array([0.17, 0.19, 0.21, 0.23, 0.23 + (0.23 - 0.21)])
lors12_per_annihilation = np.array([0.13, 0.14, 0.15, 0.16, 0.16 + (0.16 - 0.15)])
lors22_per_annihilation = np.array([0.024, 0.026, 0.028, 0.029, 0.029 + (0.029 - 0.028)])

ex_fig = plt.figure()
ex_plt = ex_fig.add_subplot()
ex_fig.subplots_adjust(top=0.88,bottom=0.13,left=0.17,right=0.9,hspace=0.2,wspace=0.2)

"""ex_plt.plot(tau, hit_first_scatter_eff, 'o-', label=r"$ e^{-x} $", ms=10, color="black")
ex_plt.plot(tau, hit_second_scatter_eff, '^-', label=r"$ e^{-x} $", ms=10, color="green")
ex_plt.plot(tau, hit_third_scatter_eff, 's-', label=r"$ e^{-x} $", ms=10, color="blue")
ex_plt.plot(tau, first_hit_second_scatter_eff, '^--', label=r"$ e^{-x} $", ms=10, markerfacecolor="none",color="green")
ex_plt.plot(tau, first_hit_third_scatter_eff, 's--', label=r"$ e^{-x} $", ms=10, markerfacecolor="none",color="blue")"""

ex_plt.plot(tau, lors_per_annihilation, 'o-', label=r"$ e^{-x} $", ms=10, color="black")
ex_plt.plot(tau, lors11_per_annihilation, '^-', label=r"$ e^{-x} $", ms=10, color="red")
ex_plt.plot(tau, lors12_per_annihilation, 's-', label=r"$ e^{-x} $", ms=10, color="blue")
ex_plt.plot(tau, lors22_per_annihilation, 'x-', label=r"$ e^{-x} $", ms=10, color="green")

ex_plt.set_xlim(left=0)
ex_plt.set_ylim(bottom=0)
ex_plt.set_ylim(top=1.0)

ex_plt.set_xlabel("Tau(mils)")
ex_plt.set_ylabel(r"LORs per Annihilation($\frac{LORs}{Annihilations}$)")
ex_plt.set_title("LORs per Annihilation vs Tau")

add_minor_ticks(ex_plt, 5)
"""ex_plt.set_yscale("log")
ex_plt.set_ylim(1e-2, 1)"""

plt.savefig("plots/lors_per_annihil.png")


plt.show()


