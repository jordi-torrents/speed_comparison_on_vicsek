# import numba as nb
import numpy as np
import sys
import numba_funcs


with open(sys.argv[1], "r") as file:
    L = int(file.readline().split()[0])
    v0 = float(file.readline().split()[0])
    rho = float(file.readline().split()[0])
    N_reset = int(file.readline().split()[0])
    N_steps = int(file.readline().split()[0])
    seed = int(file.readline().split()[0])


N_cells = L * L
N = int((L * L * rho + 2.0) / 4) * 4
with np.errstate(divide="ignore"):
    inv_2L = np.divide(2.0, L)
v0 = 0.03
th = 1.0

pos = np.empty((N, 2))
vel = np.empty((N, 2))

np.random.seed(seed)

sum_phis2 = 0.0
sum_phis = 0.0
N_observations = 0

# RANDOMIZE THE SYSTEM
pos = np.random.rand(N, 2) * L
thetas = np.random.rand(N) * 2 * np.pi
vel[:, 0] = np.cos(thetas)
vel[:, 1] = np.sin(thetas)

# SET GEOMETRY
nbr_X = np.array([-1, 0, 1, -1])
nbr_Y = np.array([1, 1, 1, 0])
nbr_indx = np.zeros((N_cells, 4), np.int32)
for cell_Y in range(L):
    for cell_X in range(L):
        cell = cell_X + L * cell_Y
        nbr_x = (cell_X + nbr_X) % L
        nbr_y = (cell_Y + nbr_Y) % L
        nbr_indx[cell] = nbr_x + L * nbr_y


eta = 0.0
pos, vel, phi, sigma_phi, xi_phi = numba_funcs.integrate(
    N_reset, False, pos, vel, eta, L, nbr_indx, v0
)


for eta in np.linspace(0, 1, 21):
    pos, vel, phi, sigma_phi, xi_phi = numba_funcs.integrate(
        N_reset, False, pos, vel, eta, L, nbr_indx, v0
    )
    pos, vel, phi, sigma_phi, xi_phi = numba_funcs.integrate(
        N_steps, True, pos, vel, eta, L, nbr_indx, v0
    )
    print(f"{eta:.2f}, {phi:.3f}, {sigma_phi:.3f}, {xi_phi:.3f}")
