import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
import sys

df = read_csv(sys.argv[1], sep="\s+")


fig, ax = plt.subplots()
ax.plot(df["N"], df["fortran"], ".-", label="Fortran")
ax.plot(df["N"], df["cpp"], ".-", label="C++")
ax.plot(df["N"], df["c"], ".-", label="C")
ax.plot(df["N"], df["python"], ".-", label="Python")
ax.set(xscale="log", yscale="log", xlabel="Number of particles", ylabel="time (s)")
ax.legend()
fig.savefig("times.png", dpi=300)
