from numpy import loadtxt
from matplotlib.pyplot import subplots
import os.path
import sys


fig, ax = subplots()
for file in sys.argv[1:]:

    x, y, sigma, xi = loadtxt(file, delimiter=",").T
    ax.errorbar(x, y, sigma, fmt=".", label=os.path.basename(file)[:-4])
    ax.fill_between(x, y - sigma, y + sigma, alpha=0.2)
ax.legend()
fig.savefig("results.png", dpi=300)
