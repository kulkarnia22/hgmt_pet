import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter
from matplotlib.ticker import LogLocator, NullFormatter

tau  = np.array([5, 4, 3, 2, 1])
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
plt.savefig("plots/lors_per_annihil_total.png", dpi=300)
plt.show()
