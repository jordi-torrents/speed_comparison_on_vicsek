import numba as nb
import numpy as np
import sys
from math import sqrt

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
nbr_indx = np.zeros((N_cells, 4), int)
for cell_Y in range(L):
    for cell_X in range(L):
        cell = cell_X + L * cell_Y
        nbr_x = (cell_X + nbr_X) % L
        nbr_y = (cell_Y + nbr_Y) % L
        nbr_indx[cell] = nbr_x + L * nbr_y


@nb.njit
def pbc(value):
    return value - int(inv_2L * value) * L


@nb.njit
def dist_PBC(pos1, pos2):
    dx = pbc(pos1[0] - pos2[0])
    dy = pbc(pos1[1] - pos2[1])
    return dx * dx + dy * dy


@nb.njit
def dist_simple(pos1, pos2):
    dx = pos1[0] - pos2[0]
    dy = pos1[1] - pos2[1]
    return dx * dx + dy * dy


@nb.njit
def integrate(steps, update_obs, pos, vel, eta):

    header = np.empty((N_cells), nb.int32)
    cell_list = np.empty(N, nb.int32)
    integr = np.empty_like(vel)
    sum_phis = 0.0
    sum_phis2 = 0.0

    factor1 = eta * 6.283185307179586
    factor2 = -eta * 3.14159265359

    for _ in range(steps):

        header[:] = -1
        integr[:] = vel

        for i in range(N):
            cell = int(pos[i, 0]) + L * int(pos[i, 1])
            cell_list[i] = header[cell]
            header[cell] = i

        for cell in range(N_cells):

            i = header[cell]
            while i > -1:

                j = cell_list[i]
                while j > -1:
                    if dist_simple(pos[i], pos[j]) < th:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbr_indx[cell, 0]]
                while j > -1:
                    if dist_PBC(pos[i], pos[j]) < th:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbr_indx[cell, 1]]
                while j > -1:
                    if dist_PBC(pos[i], pos[j]) < th:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbr_indx[cell, 2]]
                while j > -1:
                    if dist_PBC(pos[i], pos[j]) < th:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbr_indx[cell, 3]]
                while j > -1:
                    if dist_PBC(pos[i], pos[j]) < th:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                i = cell_list[i]

        particle_direction = (
            np.arctan2(integr[:, 1], integr[:, 0])
            + factor1 * np.random.rand(N)
            + factor2
        )

        vel[:, 0] = np.cos(particle_direction)
        vel[:, 1] = np.sin(particle_direction)

        pos = (pos + v0 * vel) % L

        if update_obs:
            sum_x, sum_y = np.sum(vel, axis=0)
            polar_i_sq = sum_x * sum_x + sum_y * sum_y

            sum_phis += sqrt(polar_i_sq)
            sum_phis2 += polar_i_sq

    if update_obs:
        phi = sum_phis / (N * N_steps)
        sigma_phi = sqrt(sum_phis2 / N_steps - (sum_phis / N_steps) ** 2) / N
        xi_phi = (sum_phis2 - sum_phis**2 / N_steps) / sum_phis

    return pos, vel, phi, sigma_phi, xi_phi


eta = 0.0
pos, vel, phi, sigma_phi, xi_phi = integrate(N_reset, False, pos, vel, eta)

for eta in np.linspace(0, 1, 20):
    pos, vel, phi, sigma_phi, xi_phi = integrate(N_reset, False, pos, vel, eta)
    pos, vel, phi, sigma_phi, xi_phi = integrate(N_steps, True, pos, vel, eta)
    # print(eta)
    # phi, sigma_phi, xi_phi = get_and_reset_observables()
    print(f"{eta:.2f}, {phi:.3f}, {sigma_phi:.3f}, {xi_phi:.3f}")
#
