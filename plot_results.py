import numpy as np
import matplotlib.pyplot as plt
from os.path import exists


fig, ax = plt.subplots()
for file in ["cpp", "fortran", "c", "python"]:
    if exists(file + ".out"):
        x, y, sigma, xi = np.loadtxt(file + ".out", delimiter=",").T
        ax.errorbar(x, y, sigma, fmt=".", label=file)
        ax.fill_between(x, y - sigma, y + sigma, alpha=0.2)
ax.legend()
fig.savefig("results.png", dpi=300)
