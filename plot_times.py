import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
import sys

df = read_csv(sys.argv[1], sep="\s+")

fig, ax = plt.subplots()

for key in df.keys()[2:]:
    ax.plot(df["N"], df[key], ".-", label=key)

ax.set(
    xscale="log", yscale="log", xlabel="Number of particles", ylabel="time (s)"
)
ax.legend()
fig.savefig("times.png", dpi=300)
