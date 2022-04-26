import numba as nb
import numpy as np
import sys

spec = [
    ("v0", nb.float32),
    ("L", nb.int32),
    ("N_cells", nb.int32),
    ("N", nb.int32),
    ("inv_2L", nb.float32),
    ("v0", nb.float32),
    ("th", nb.float32),
    ("pos", nb.float32[:]),
    ("vel", nb.float32[:]),
    ("sum_phis2", nb.float32),
    ("sum_phis", nb.float32),
    ("N_observations", nb.int32),
    ("array", nb.float32[:]),
]


@nb.experimental.jitclass(spec)
class vicsek_system:
    def __init__(self, rho_in, v0_in, L_in, seed):
        self.v0 = v0_in
        self.L = L_in
        self.N_cells = L_in * L_in
        self.N = int((L_in * L_in * rho_in + 2.0) / 4) * 4
        self.inv_2L = 2.0 / float(L_in)
        self.v0 = 0.03
        self.th = 1.0

        self.pos = np.empty((self.N, 2))
        self.vel = np.empty((self.N, 2))

        np.random.seed(seed)

        self.randomize_system()
        self.set_geometry()

        self.sum_phis2 = 0.0
        self.sum_phis = 0.0
        self.N_observations = 0

    def pbc(self, value):
        return value - int(self.inv_2L * value) * L

    def dist_PBC(self, pos1, pos2):

        dx = self.pbc(pos1[0] - pos2[0])
        dy = self.pbc(pos1[1] - pos2[1])
        return dx * dx + dy * dy

    def dist_simple(self, pos1, pos2):

        dx = pos1[0] - pos2[0]
        dy = pos1[1] - pos2[1]
        return dx * dx + dy * dy

    def wrap3(self, value):
        return (
            value + self.L
            if value < 0
            else (value if value < self.L else value - self.L)
        )

    def get_and_reset_observables(self):

        phi = self.sum_phis / (self.N * self.N_observations)
        sigma_phi = (
            np.sqrt(
                self.sum_phis2 / self.N_observations
                - (self.sum_phis / self.N_observations) ** 2
            )
            / self.N
        )
        xi_phi = (
            self.sum_phis2 - self.sum_phis**2 / self.N_observations
        ) / self.sum_phis

        self.sum_phis2 = 0.0
        self.sum_phis = 0.0
        self.N_observations = 0
        return phi, sigma_phi, xi_phi

    def update_observables(self):

        sum_x, sum_y = np.sum(self.vel, axis=0)
        polar_i_sq = sum_x * sum_x + sum_y * sum_y

        self.sum_phis += np.sqrt(polar_i_sq)
        self.sum_phis2 += polar_i_sq
        self.N_observations += 1

    def set_geometry(self):
        nbr_X = np.array([-1, 0, 1, -1])
        nbr_Y = np.array([1, 1, 1, 0])

        self.nbr_indx = np.zeros((self.N_cells, 4), int)

        for cell_Y in range(L):
            for cell_X in range(L):
                cell = cell_X + L * cell_Y
                nbr_x = (cell_X + nbr_X) % L
                nbr_y = (cell_Y + nbr_Y) % L
                self.nbr_indx[cell] = nbr_x + L * nbr_y

    def randomize_system(self):
        self.pos = np.random.rand(self.N, 2) * L
        thetas = np.random.rand(self.N) * 2 * np.pi
        self.vel[:, 0] = np.cos(thetas)
        self.vel[:, 1] = np.sin(thetas)

    def integrate(self, steps=1, update_obs=0):

        integr = np.empty((self.N, 2))
        header = np.empty((self.N_cells), int)
        cell_list = np.empty(self.N, int)

        particle_direction = np.empty(self.N)

        factor1 = self.eta * 6.283185307179586
        factor2 = -self.eta * 3.14159265359

        for step in range(steps):
            # print(step)

            header[:] = -1
            integr = np.copy(self.vel)

            for i in range(self.N):
                cell = int(self.pos[i, 0]) + L * int(self.pos[i, 1])
                cell_list[i] = header[cell]
                header[cell] = i

            for cell in range(self.N_cells):

                i = header[cell]
                while i > -1:

                    j = cell_list[i]
                    while j > -1:
                        if self.dist_simple(self.pos[i], self.pos[j]) < self.th:
                            integr[i] += self.vel[j]
                            integr[j] += self.vel[i]
                        j = cell_list[j]

                    j = header[self.nbr_indx[cell, 0]]
                    while j > -1:
                        if self.dist_PBC(self.pos[i], self.pos[j]) < self.th:
                            integr[i] += self.vel[j]
                            integr[j] += self.vel[i]
                        j = cell_list[j]

                    j = header[self.nbr_indx[cell, 1]]
                    while j > -1:
                        if self.dist_PBC(self.pos[i], self.pos[j]) < self.th:
                            integr[i] += self.vel[j]
                            integr[j] += self.vel[i]
                        j = cell_list[j]

                    j = header[self.nbr_indx[cell, 2]]
                    while j > -1:
                        if self.dist_PBC(self.pos[i], self.pos[j]) < self.th:
                            integr[i] += self.vel[j]
                            integr[j] += self.vel[i]
                        j = cell_list[j]

                    j = header[self.nbr_indx[cell, 3]]
                    while j > -1:
                        if self.dist_PBC(self.pos[i], self.pos[j]) < self.th:
                            integr[i] += self.vel[j]
                            integr[j] += self.vel[i]
                        j = cell_list[j]

                    i = cell_list[i]

            particle_direction = (
                np.arctan2(integr[:, 1], integr[:, 0])
                + factor1 * np.random.rand(self.N)
                + factor2
            )

            self.vel[:, 0] = np.cos(particle_direction)
            self.vel[:, 1] = np.sin(particle_direction)

            for i in range(self.N):
                self.pos[i, 0] = self.wrap3(self.pos[i, 0] + v0 * self.vel[i, 0])
                self.pos[i, 1] = self.wrap3(self.pos[i, 1] + v0 * self.vel[i, 1])

            if update_obs:
                self.update_observables()


with open(sys.argv[1], "r") as file:
    L = int(file.readline().split()[0])
    v0 = float(file.readline().split()[0])
    rho = float(file.readline().split()[0])
    N_reset = int(file.readline().split()[0])
    N_steps = int(file.readline().split()[0])
    seed = int(file.readline().split()[0])


system = vicsek_system(rho, v0, L, seed)


system.eta = 0.0
system.integrate(N_reset, False)

for system.eta in np.linspace(0, 1, 20):
    system.integrate(N_reset, False)
    system.integrate(N_steps, True)
    phi, sigma_phi, xi_phi = system.get_and_reset_observables()
    print(f"{system.eta:.2f}, {phi:.3f}, {sigma_phi:.3f}, {xi_phi:.3f}")
