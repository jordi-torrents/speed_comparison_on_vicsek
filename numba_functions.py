from numba.pycc import CC
import numpy as np
from math import sqrt, atan2, sin, cos

import numba as nb

cc = CC("numba_funcs")


# @cc.export('multf', 'f8(f8, f8)')
# @cc.export('multi', 'i4(i4, i4)')
# def mult(a, b):
#     return a * b

# @cc.export('square', 'f8(f8)')
# def square(a):
#     return a ** 2


@nb.jit
def pbc(value, inv_2L, L):
    return value - int(inv_2L * value) * L


@nb.jit
def dist_PBC(pos1, pos2, inv_2L, L):
    dx = pbc(pos1[0] - pos2[0], inv_2L, L)
    dy = pbc(pos1[1] - pos2[1], inv_2L, L)
    return dx * dx + dy * dy


@nb.jit
def dist_simple(pos1, pos2):
    dx = pos1[0] - pos2[0]
    dy = pos1[1] - pos2[1]
    return dx * dx + dy * dy


@cc.export(
    "integrate",
    "Tuple((f8[:,:], f8[:,:], f8, f8, f8))(i4, b1, f8[:,:], f8[:,:], f8, i4, i4[:,:], f8)",
)
@nb.njit()
def integrate(steps, update_obs, pos, vel, eta, L, nbr_indx, v0):
    N_cells = L * L
    inv_2L = 2.0 / L
    N = pos.shape[0]

    header = np.empty(N_cells, nb.int32)
    cell_list = np.empty(N, nb.int32)
    integr = np.empty_like(vel)
    sum_phis = 0.0
    sum_phis2 = 0.0

    fct1 = eta * 6.283185307179586
    fct2 = -eta * 3.14159265359

    for _ in range(steps):

        integr[:] = vel[:]
        header[:] = -1

        for i in range(N):
            cell = int(pos[i, 0]) + L * int(pos[i, 1])
            cell_list[i] = header[cell]
            header[cell] = i

        for cell in range(N_cells):

            i = header[cell]
            while i > -1:

                j = cell_list[i]
                while j > -1:
                    if dist_simple(pos[i], pos[j]) < 1.0:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbr_indx[cell, 0]]
                while j > -1:
                    if dist_PBC(pos[i], pos[j], inv_2L, L) < 1.0:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbr_indx[cell, 1]]
                while j > -1:
                    if dist_PBC(pos[i], pos[j], inv_2L, L) < 1.0:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbr_indx[cell, 2]]
                while j > -1:
                    if dist_PBC(pos[i], pos[j], inv_2L, L) < 1.0:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbr_indx[cell, 3]]
                while j > -1:
                    if dist_PBC(pos[i], pos[j], inv_2L, L) < 1.0:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                i = cell_list[i]

        directions = (
            np.arctan2(integr[:, 1], integr[:, 0]) + fct1 * np.random.rand(N) + fct2
        )

        vel[:, 0] = np.sin(directions)
        vel[:, 1] = np.cos(directions)

        pos = (pos + v0 * vel) % L

        if update_obs:
            sum_x, sum_y = np.sum(vel, axis=0)
            polar_i_sq = sum_x * sum_x + sum_y * sum_y

            sum_phis += sqrt(polar_i_sq)
            sum_phis2 += polar_i_sq

    if update_obs:
        phi = sum_phis / (N * steps)
        sigma_phi = sqrt(sum_phis2 / steps - (sum_phis / steps) ** 2) / N
        xi_phi = (sum_phis2 - sum_phis**2 / steps) / sum_phis
    else:
        phi = 0.0
        sigma_phi = 0.0
        xi_phi = 0.0
    # print(pos[:5])
    return pos, vel, phi, sigma_phi, xi_phi


if __name__ == "__main__":
    cc.compile()
