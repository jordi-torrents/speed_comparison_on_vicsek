from numba.pycc import CC
import numpy as np
from math import sqrt, atan2, sin, cos


# @cc.export('multf', 'f8(f8, f8)')
# @cc.export('multi', 'i4(i4, i4)')
# def mult(a, b):
#     return a * b

# @cc.export('square', 'f8(f8)')
# def square(a):
#     return a ** 2


def integrate(steps, update_obs, pos, vel, eta, L, nbr_indx, v0):
    N_cells = L * L
    inv_2L = 2.0 / L
    N = pos.shape[0]

    def pbc(value):
        return value - int(inv_2L * value) * L

    def dist_PBC(pos1, pos2):
        dx = pbc(pos1[0] - pos2[0])
        dy = pbc(pos1[1] - pos2[1])
        return dx * dx + dy * dy

    def dist_simple(pos1, pos2):
        dx = pos1[0] - pos2[0]
        dy = pos1[1] - pos2[1]
        return dx * dx + dy * dy

    # header = np.empty(N_cells, int)
    # cell_list = np.empty(N, int)
    cell_list = [0] * N
    integr = np.empty_like(vel)
    sum_phis = 0.0
    sum_phis2 = 0.0

    fct1 = eta * 6.283185307179586
    fct2 = -eta * 3.14159265359

    for _ in range(steps):

        np.copyto(integr, vel)
        header = [None] * N_cells

        cell_indxs = pos[:, 0].astype(int) + L * pos[:, 1].astype(int)
        for i in range(N):
            cell_list[i] = header[cell_indxs[i]]
            header[cell_indxs[i]] = i

        for i, nbrs in zip(header, nbr_indx):

            while i is not None:

                j = cell_list[i]
                while j is not None:
                    if dist_simple(pos[i], pos[j]) < 1.0:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbrs[0]]
                while j is not None:
                    if dist_PBC(pos[i], pos[j]) < 1.0:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbrs[1]]
                while j is not None:
                    if dist_PBC(pos[i], pos[j]) < 1.0:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbrs[2]]
                while j is not None:
                    if dist_PBC(pos[i], pos[j]) < 1.0:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                j = header[nbrs[3]]
                while j is not None:
                    if dist_PBC(pos[i], pos[j]) < 1.0:
                        integr[i] += vel[j]
                        integr[j] += vel[i]
                    j = cell_list[j]

                i = cell_list[i]

        directions = (
            np.arctan2(integr[:, 1], integr[:, 0])
            + fct1 * np.random.rand(N)
            + fct2
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

    return pos, vel, phi, sigma_phi, xi_phi
