import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter
from matplotlib.ticker import LogLocator, NullFormatter

"""tau  = np.array([5, 4, 3, 2, 1])
lors_per_annihilation = np.array([0.51, 0.57, 0.63, 0.69, 0.69 + (0.69 - 0.63)])
lors11_per_annihilation = np.array([0.17, 0.19, 0.21, 0.23, 0.23 + (0.23 - 0.21)])
lors12_per_annihilation = np.array([0.13, 0.14, 0.15, 0.16, 0.16 + (0.16 - 0.15)])
lors22_per_annihilation = np.array([0.024, 0.026, 0.028, 0.029, 0.029 + (0.029 - 0.028)])

m_lor, b_lor = np.polyfit(tau, lors_per_annihilation, 1)   # slope m, intercept b
x_line = np.linspace(0.5, tau.max(), 200)   # includes x = 0
lors_line = m_lor * x_line + b_lor

m_lor11, b_lor11 = np.polyfit(tau, lors11_per_annihilation, 1)
lors11_line = m_lor11 * x_line + b_lor11

m_lor12, b_lor12 = np.polyfit(tau, lors12_per_annihilation, 1)
lors12_line = m_lor12 * x_line + b_lor12

m_lor22, b_lor22 = np.polyfit(tau, lors22_per_annihilation, 1)
lors22_line = m_lor22 * x_line + b_lor22

fig, ax = plt.subplots()
ax.plot(tau, lors_per_annihilation, 'o', ms = 10, color = 'black')
ax.plot(x_line, lors_line, linestyle="-", color="black")

ax.plot(tau, lors11_per_annihilation, '^', ms = 10, color='red')
ax.plot(x_line, lors11_line, linestyle="-", color="red")

ax.plot(tau, lors12_per_annihilation, 's', ms = 10, color='blue')
ax.plot(x_line, lors12_line, linestyle='-', color='blue')

ax.plot(tau, lors22_per_annihilation, 'x', ms = 10, color='green')
ax.plot(x_line, lors22_line, linestyle='-', color='green')

ax.set_xlabel("Tau (mil)", fontsize=16)
ax.set_ylabel(r"LORs Per Annhilation $\frac{LORs}{Annihilaitons}$", fontsize=16)
ax.set_title("LORs Per Annihilation v.s Tau", fontsize=18)

# --- Axes limits ---
ax.set_yscale("log")
ax.set_xlim(0, 5.1)
ax.set_ylim(1e-2, 1)   # now you should see 10^-2, 10^-1, 10^0
ax.tick_params(axis="y", which="minor", labelsize=10, direction="in")
ax.tick_params(axis="y", which="major", labelsize=10, direction="in")
ax.tick_params(axis="x", which="major", labelsize=10, direction="in")

fig.tight_layout()
plt.savefig("plots/lors_per_annihil_log.png", dpi=300)
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
lors_per_annihilation = np.array([0.54, 0.60, 0.66, 0.71, 0.75])
lors11_per_annihilation = np.array([0.17, 0.19, 0.21, 0.22, 0.20])
lors12_per_annihilation = np.array([0.13, 0.14, 0.15, 0.16, 0.14])
lors22_per_annihilation = np.array([0.025, 0.026, 0.027, 0.027, 0.026])

ex_fig = plt.figure()
ex_plt = ex_fig.add_subplot()
ex_fig.subplots_adjust(top=0.88,bottom=0.13,left=0.17,right=0.9,hspace=0.2,wspace=0.2)

ex_plt.plot(tau, lors_per_annihilation, 'o-', label=r"$ e^{-x} $", ms=10, color="black")
ex_plt.plot(tau, lors11_per_annihilation, '^-', label=r"$ e^{-x} $", ms=10, color="red")
ex_plt.plot(tau, lors12_per_annihilation, 's-', label=r"$ e^{-x} $", ms=10, color="blue")
ex_plt.plot(tau, lors22_per_annihilation, 'x-', label=r"$ e^{-x} $", ms=10, color="green")

ex_plt.set_xlim(left=0)
ex_plt.set_ylim(bottom=0)
ex_plt.set_ylim(top=1.0)

ex_plt.set_xlabel("Tau(mils)")
ex_plt.set_ylabel(r"Log LORs per Annihilation($\frac{LORs}{Annihilations}$)")
ex_plt.set_title("LORs per Annihilation vs Tau")

add_minor_ticks(ex_plt, 5)
ex_plt.set_yscale("log")
ex_plt.set_ylim(1e-2, 1)

plt.savefig("plots/lors_per_annihil_log.png")


plt.show()