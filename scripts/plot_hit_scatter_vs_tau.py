'''import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

thickness_mm = np.array([5, 4, 3, 2, 1])

hit_first_scatter_eff = [0.51, 0.56, 0.61, 0.65, 0.65 + (0.65-0.61)]
hit_second_scatter_eff = [0.33, 0.36, 0.41, 0.46, 0.46 + (0.46-0.41)]
hit_third_scatter_eff = [0.18, 0.20, 0.24, 0.29, 0.29 + (0.29-0.24)]

first_hit_second_scatter_eff = [0.083, 0.085, 0.086, 0.087, 0.087 + (0.087-0.086)]
first_hit_third_scatter_eff = [0.032, 0.031, 0.030, 0.028, 0.028 + (0.028-0.030)]

# APS-like styling
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

m_first_scatter, b_first_scatter = np.polyfit(thickness_mm, hit_first_scatter_eff, 1)   # slope m, intercept b
x_line = np.linspace(0.5, thickness_mm.max(), 200)   # includes x = 0
first_scatter_line = m_first_scatter * x_line + b_first_scatter

m_second_scatter, b_second_scatter = np.polyfit(thickness_mm, hit_second_scatter_eff, 1)   # slope m, intercept b
x_line = np.linspace(0.5, thickness_mm.max(), 200)   # includes x = 0
second_scatter_line = m_second_scatter * x_line + b_second_scatter

m_first_hit_second_scatter, b_first_hit_second_scatter = np.polyfit(thickness_mm, first_hit_second_scatter_eff, 1)   # slope m, intercept b
x_line = np.linspace(0.5, thickness_mm.max(), 200)   # includes x = 0
first_hit_second_scatter_line = m_first_hit_second_scatter * x_line + b_first_hit_second_scatter

m_third_scatter, b_third_scatter = np.polyfit(thickness_mm, hit_third_scatter_eff, 1)   # slope m, intercept b
x_line = np.linspace(0.5, thickness_mm.max(), 200)   # includes x = 0
third_scatter_line = m_third_scatter * x_line + b_third_scatter

m_first_hit_third_scatter, b_first_hit_third_scatter = np.polyfit(thickness_mm, first_hit_third_scatter_eff, 1)   # slope m, intercept b
x_line = np.linspace(0.5, thickness_mm.max(), 200)   # includes x = 0
first_hit_third_scatter_line = m_first_hit_third_scatter * x_line + b_first_hit_third_scatter

# Point plot
ax.plot(thickness_mm, hit_first_scatter_eff, "o", color="red", markersize=7)
ax.plot(x_line, first_scatter_line, linestyle='-', color='red')

ax.plot(thickness_mm, hit_second_scatter_eff, "^", color="green", markersize=7)
ax.plot(x_line, second_scatter_line, linestyle='-', color='green')

#ax.plot(thickness_mm, first_hit_second_scatter_eff, "^", color="green", markerfacecolor="none", markersize="7")
#ax.plot(x_line, first_hit_second_scatter_line, linestyle='--', color='green')

ax.plot(thickness_mm, hit_third_scatter_eff, "s", color="blue", markersize=7)
ax.plot(x_line, third_scatter_line, linestyle='-', color='blue')

#ax.plot(thickness_mm, first_hit_third_scatter_eff, "s", color="blue", markerfacecolor="none", markersize="7")
#ax.plot(x_line, first_hit_third_scatter_line, linestyle='--', color='blue')

# Labels
ax.set_xlabel(r"Thickness Tau (mils)", fontsize = 14)
ax.set_ylabel(r"Hit Efficieny$(\frac{hits}{scatters})$", fontsize = 14)
ax.set_title("Hits per Scatter vs Tau")

# Minor ticks
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax.set_ylim(1e-2, 1)
ax.set_ylim(0, 1)
#ax.set_yscale("log")
ax.set_xlim(0, 5.1)

plt.tight_layout()
plt.savefig("plots/hit_scatter_vs_thickness.png", dpi=300)
plt.show()'''

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
hit_first_scatter_eff = [0.51, 0.56, 0.61, 0.65, 0.686545]
hit_second_scatter_eff = [0.33, 0.36, 0.41, 0.46, 0.540820]
hit_third_scatter_eff = [0.18, 0.20, 0.24, 0.29, 0.399484]
first_hit_second_scatter_eff = [0.083, 0.085, 0.086, 0.087, 0.089053]
first_hit_third_scatter_eff = [0.032, 0.031, 0.030, 0.028, 0.028398]

ex_fig = plt.figure()
ex_plt = ex_fig.add_subplot()
ex_fig.subplots_adjust(top=0.88,bottom=0.13,left=0.17,right=0.9,hspace=0.2,wspace=0.2)

ex_plt.plot(tau, hit_first_scatter_eff, 'o-', label=r"$ e^{-x} $", ms=10, color="black")
ex_plt.plot(tau, hit_second_scatter_eff, '^-', label=r"$ e^{-x} $", ms=10, color="green")
ex_plt.plot(tau, hit_third_scatter_eff, 's-', label=r"$ e^{-x} $", ms=10, color="blue")
ex_plt.plot(tau, first_hit_second_scatter_eff, '^--', label=r"$ e^{-x} $", ms=10, markerfacecolor="none",color="green")
ex_plt.plot(tau, first_hit_third_scatter_eff, 's--', label=r"$ e^{-x} $", ms=10, markerfacecolor="none",color="blue")

ex_plt.set_xlim(left=0)
ex_plt.set_ylim(bottom=0)
ex_plt.set_ylim(top=1.0)

ex_plt.set_xlabel("Tau(mils)")
ex_plt.set_ylabel(r"Hits per Scatter($\frac{Hits}{Scatters}$)")
ex_plt.set_title("Log Plot Hits per Scatter vs Tau")

add_minor_ticks(ex_plt, 5)
ex_plt.set_yscale("log")
ex_plt.set_ylim(1e-2, 1)

plt.savefig("plots/log_hits_per_scatter.png")


plt.show()

