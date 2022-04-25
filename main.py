import numba as nb
import numpy as np
import sys


with open(sys.argv[1], "r") as file:
    L = int(file.readline().split()[0])
    v0 = float(file.readline().split()[0])
    rho = float(file.readline().split()[0])
    N_reset = int(file.readline().split()[0])
    N_steps = int(file.readline().split()[0])
    seed = int(file.readline().split()[0])

np.random.seed(seed)

N = int(L * L * rho)
N_cells = L * L


thetas = np.random.rand(N) * 2 * np.pi
x = np.random.rand(N) * L
y = np.random.rand(N) * L
vx = np.cos(thetas)
vy = np.sin(thetas)


nbr_X = np.array([-1, 0, 1, -1, 0, 1, -1, 0, 1])
nbr_Y = np.array([1, 1, 1, 0, 0, 0, -1, -1, -1])

near_nbr_cells = np.zeros((N_cells, 9), int)

for cell_Y in range(L):
    for cell_X in range(L):
        cell = cell_X + L * cell_Y
        nbr_x = (cell_X + nbr_X) % L
        nbr_y = (cell_Y + nbr_Y) % L
        near_nbr_cells[cell] = nbr_x + L * nbr_y


@nb.njit
def neighbours_direction(cell_direction, where_is):

    cell_vel_x = np.zeros(N_cells)
    cell_vel_y = np.zeros(N_cells)
    occupied = np.zeros(N_cells)

    # occupied[where_is] = True
    for i in range(N):
        occupied[where_is[i]] = True
        cell_vel_x[where_is[i]] = cell_vel_x[where_is[i]] + vx[i]
        cell_vel_y[where_is[i]] = cell_vel_y[where_is[i]] + vy[i]

    for cell in range(N_cells):

        if occupied[cell]:

            cell_mean_vel_x = np.sum(cell_vel_x[near_nbr_cells[cell]])
            cell_mean_vel_y = np.sum(cell_vel_y[near_nbr_cells[cell]])

            cell_direction[cell] = np.arctan2(cell_mean_vel_y, cell_mean_vel_x)

    # print(cell_direction)


@nb.njit
def integrate_and_mesure(eta, steps, x, y, vx, vy):
    cell_direction = np.zeros(N_cells)
    count_polar = 0.0
    count_polar2 = 0.0
    noise_constant = 2 * eta * np.pi
    theta = np.zeros(N)

    for step in range(steps):

        where_is = x.astype(np.int32) + L * y.astype(np.int32)
        neighbours_direction(cell_direction, where_is)
        # print(cell_direction)

        randoms = np.random.rand(N) - 0.5
        for i in range(N):
            theta[i] = cell_direction[where_is[i]] + noise_constant * randoms[i]

        vx = np.cos(theta)
        vy = np.sin(theta)
        x = (x + v0 * vx) % L
        y = (y + v0 * vy) % L

        polar_i = np.sqrt(np.sum(vx) ** 2 + np.sum(vy) ** 2)
        count_polar = count_polar + polar_i
        count_polar2 = count_polar2 + polar_i * polar_i
        # if step % 5 == 0:
        #     x = x % L
        #     y = y % L

    polar = count_polar / (float(N * steps))
    polar2 = count_polar2 / (float(N * N * steps))

    return polar, polar2, x, y, vx, vy


eta = 0.0
polar, polar2, x, y, vx, vy = integrate_and_mesure(eta, N_reset, x, y, vx, vy)

for eta in np.linspace(0, 1, 20):
    polar, polar2, x, y, vx, vy = integrate_and_mesure(eta, N_reset, x, y, vx, vy)
    polar, polar2, x, y, vx, vy = integrate_and_mesure(eta, N_steps, x, y, vx, vy)
    print(f"{eta:.2f}, {polar:.3f}, {np.sqrt(polar2 - (polar * polar)):.3f}")


# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation

# fig, ax = plt.subplots()
# line = ax.plot(x, y, "k.")[0]
# ax.set(xlim=(0, L), ylim=(0, L), aspect=1, xticks=np.arange(L), yticks=np.arange(L))
# ax.grid(1)


# def animate(i):
#     global x, y, vx, vy
#     polar, polar2, x, y, vx, vy = integrate_and_mesure(0.65, 1, x, y, vx, vy)
#     line.set_data(x, y)
#     return (line,)


# ani = FuncAnimation(fig, animate, frames=10, blit=True, interval=1)
# plt.show()
