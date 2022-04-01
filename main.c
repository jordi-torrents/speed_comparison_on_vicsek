#include <math.h>
// #include <cstring>
#include <string.h>
// #include "mt19937-64.c"
// #define RAN (genrand64_real2())

// #include "mt19937ar.c"
// #define RAN (genrand_real2())

#include "SFMT-src-1.5.1/SFMT.h"
sfmt_t sfmt;
#define RAN (sfmt_genrand_real2(&sfmt))
#define x 0
#define y 1
#define vx 2
#define vy 3

float rho, L, eta, v0;
int N, N_cells, int_L;

void set_geometry(int (*near_nbr_cells)[9])
{
    int cell, cell_X, cell_Y, nbr_cell_X, nbr_cell_Y, nbr_cell, k;
    int nbr_X[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
    int nbr_Y[9] = {1, 1, 1, 0, 0, 0, -1, -1, -1};

    for (cell_Y = 0; cell_Y < int_L; cell_Y++)
    {
        for (cell_X = 0; cell_X < int_L; cell_X++)
        {
            cell = cell_X + int_L * cell_Y;
            for (k = 0; k < 9; k++)
            {
                nbr_cell_X = cell_X + nbr_X[k];
                nbr_cell_Y = cell_Y + nbr_Y[k];

                if (nbr_cell_X < 0)
                    nbr_cell_X += int_L;

                if (nbr_cell_X >= int_L)
                    nbr_cell_X -= int_L;

                if (nbr_cell_Y < 0)
                    nbr_cell_Y += int_L;

                if (nbr_cell_Y >= int_L)
                    nbr_cell_Y -= int_L;

                nbr_cell = nbr_cell_X + int_L * nbr_cell_Y;

                near_nbr_cells[cell][k] = nbr_cell;
            }
        }
    }
}

void neighbours_direction(float (*data)[4], int (*near_nbr_cells)[9], float **particle_direction, float *cell_direction)
{
    float cell_vel_x[N_cells]; // = {0.0};
    float cell_vel_y[N_cells]; // = {0.0};
    float cell_mean_vel_x, cell_mean_vel_y;
    int cell;
    int occupied[N_cells]; // = {0};
    // memset(occupied, 0, N_cells * sizeof(int));
    // memset(cell_vel_x, 0, N_cells * sizeof(float));
    // memset(cell_vel_y, 0, N_cells * sizeof(float));

    for (cell = 0; cell < N_cells; cell++)
    {
        cell_vel_x[cell] = 0.0;
        cell_vel_y[cell] = 0.0;
        occupied[cell] = 0;
    }

    for (int i = 0; i < N; i++)
    {
        cell = (int_L + (int)data[i][x]) % int_L + int_L * ((int_L + (int)data[i][y]) % int_L);

        occupied[cell] = 1;
        cell_vel_x[cell] += data[i][vx];
        cell_vel_y[cell] += data[i][vy];
        particle_direction[i] = &cell_direction[cell];
    }

    for (cell = 0; cell < N_cells; cell++)
    {
        if (occupied[cell] == 1)
        {
            cell_mean_vel_x = 0.0;
            cell_mean_vel_y = 0.0;

            for (int i = 0; i < 9; i++)
            {
                int nbr_cell = near_nbr_cells[cell][i];
                cell_mean_vel_x += cell_vel_x[nbr_cell];
                cell_mean_vel_y += cell_vel_y[nbr_cell];
            }

            cell_direction[cell] = atan2(cell_mean_vel_y, cell_mean_vel_x);
        }
    }
}

void integrate_and_mesure(
    float (*data)[4],
    int (*near_nbr_cells)[9],
    int steps,
    float *polar,
    float *polar2)
{
    float polar_i, sum_x = 0.0, sum_y = 0.0, cell_direction[N_cells];
    float *particle_direction[N], theta, count_polar = 0.0, count_polar2 = 0.0;

    for (int step = 0; step < steps; step++)
    {
        sum_x = 0.0;
        sum_y = 0.0;
        neighbours_direction(data, near_nbr_cells, particle_direction, cell_direction);
        for (int i = 0; i < N; i++)
        {
            theta = *particle_direction[i] + eta * 6.283185307179586 * (RAN - 0.5);
            data[i][vx] = cos(theta);
            data[i][vy] = sin(theta);
            data[i][x] += v0 * data[i][vx];
            data[i][y] += v0 * data[i][vy];
            sum_x += data[i][vx];
            sum_y += data[i][vy];
        }

        polar_i = sqrt(sum_x * sum_x + sum_y * sum_y);
        count_polar += polar_i;
        count_polar2 += (polar_i * polar_i);
        if (step % 5000 == 0)
            for (int i = 0; i < N; i++)
            {
                data[i][x] = fmod(data[i][x], L);
                data[i][y] = fmod(data[i][y], L);
            }
    }
    *polar = count_polar / ((double)(N * steps));
    *polar2 = count_polar2 / ((double)(N * N * steps));
}

void integrate(float (*data)[4], int (*near_nbr_cells)[9], int steps)
{
    float cell_direction[N_cells];
    float *particle_direction[N];
    for (int step = 0; step < steps; step++)
    {
        neighbours_direction(data, near_nbr_cells, particle_direction, cell_direction);
        for (int i = 0; i < N; i++)
        {
            float theta = *particle_direction[i] + eta * 6.283185307179586 * (RAN - 0.5);
            data[i][vx] = cos(theta);
            data[i][vy] = sin(theta);
            data[i][x] += v0 * data[i][vx];
            data[i][y] += v0 * data[i][vy];
        }

        if (step % 5000 == 0)
            for (int i = 0; i < N; i++)
            {
                data[i][x] = fmod(data[i][x], L);
                data[i][y] = fmod(data[i][y], L);
            }
    }
}

int main(int argc, char *argv[])
{
    float polar = 0.0, polar2 = 0.0;
    int N_reset, N_steps, seed;

    // init_genrand(123);
    sfmt_init_gen_rand(&sfmt, 12345);
    // printf("%i %i\n", sfmt_get_min_array_size32(&sfmt), sfmt_get_min_array_size64(&sfmt));

    FILE *in_file = fopen("input.dat", "r"); // read only
    if (fscanf(in_file, "%d %*s\n%f %*s\n%f %*s\n%i %*s\n%i %*s\n%i %*s\n", &int_L, &v0, &rho, &N_reset, &N_steps, &seed) != 6)
        printf("ERROR on input!");

    fclose(in_file);

    L = (double)int_L;
    N = (int)L * L * rho;
    N_cells = int_L * int_L;

    float data[N][4];
    int near_nbr_cells[N_cells][9];

    for (int i = 0; i < N; i++)
    {
        float theta = RAN * 6.283185307;
        data[i][x] = RAN * L;
        data[i][y] = RAN * L;
        data[i][vx] = cos(theta);
        data[i][vy] = sin(theta);
    }

    set_geometry(near_nbr_cells);

    eta = 0.0;

    integrate(data, near_nbr_cells, N_reset);

    for (eta = 0.0; eta < 1.01; eta += 0.05)
    {
        integrate(data, near_nbr_cells, N_reset);
        integrate_and_mesure(data, near_nbr_cells, N_steps, &polar, &polar2);

        printf("%f, %f, %f\n", eta, polar, sqrt(polar2 - (polar * polar)));
    }
}
