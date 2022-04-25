import numpy as np
import matplotlib.pyplot as plt

Ns = np.loadtxt("Ns.out")
Fortran_times = np.loadtxt("times_fortran.out")
Cpp_times = np.loadtxt("times_cpp.out")
C_times = np.loadtxt("times_c.out")
Python_times = np.loadtxt("times_python.out")


fig, ax = plt.subplots()
ax.plot(Ns, Fortran_times, ".-", label="Fortran")
ax.plot(Ns, Cpp_times, ".-", label="C++")
ax.plot(Ns, C_times, ".-", label="C")
ax.plot(Ns, Python_times, ".-", label="Python")
ax.set(xscale="log", yscale="log", xlabel="Number of particles", ylabel="time (s)")
ax.legend()
fig.savefig("times.png", dpi=300)
