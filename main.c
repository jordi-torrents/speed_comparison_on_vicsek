#include <math.h>
// #include <cstring>
#include <string.h>

#include "SFMT-src-1.5.1/SFMT.h"
sfmt_t sfmt;

float rho, eta, v0;
int N, N_cells, L, N_observations, N4;
float v0, inv_2L, eta, sum_phis, sum_phis2;
const float th = 1.0, rand_const = (1.0 / 4294967296.0);
float *pos, *vel;
int *near_nbr_cells;
sfmt_t sfmt;

inline float pbc(float value)
{
    return value - (int)(inv_2L * value) * L;
}

inline float dist_PBC(float *pos1, float *pos2)
{
    float dx = pbc(*(pos1++) - *(pos2++));
    float dy = pbc(*pos1 - *pos2);
    return dx * dx + dy * dy;
}

inline float dist_simple(float *pos1, float *pos2)
{
    float dx = *(pos1++) - *(pos2++);
    float dy = *pos1 - *pos2;
    return dx * dx + dy * dy;
}

inline float wrap3(float value)
{
    return (value > 0 ? (value > L ? value - L : value) : value + L);
}

void get_and_reset_observables(float *phi, float *sigma_phi, float *xi_phi)
{
    *phi = sum_phis / ((float)(N) * (float)(N_observations));
    *sigma_phi = sqrt(sum_phis2 / (float)(N_observations) - (sum_phis * sum_phis) / (float)(N_observations * N_observations)) / (float)(N);
    *xi_phi = (sum_phis2 - (sum_phis * sum_phis) / (float)(N_observations)) / sum_phis;

    sum_phis2 = 0.0;
    sum_phis = 0.0;
    N_observations = 0;
}

void update_observables()
{

    float sum_x = 0, sum_y = 0;

    for (int i = 0; i < 2 * N;)
    {
        sum_x += vel[i++];
        sum_y += vel[i++];
    }

    float polar_i_sq = sum_x * sum_x + sum_y * sum_y;

    sum_phis += sqrt(polar_i_sq);
    sum_phis2 += polar_i_sq;
    N_observations++;
}

void randomize_system()
{
    int random_lenght = fmax(3 * N4, sfmt_get_min_array_size32(&sfmt));
    uint32_t random[random_lenght];
    sfmt_fill_array32(&sfmt, random, random_lenght);

    int rand_indx = 0;
    float theta;
    for (int i = 0; i < N; i++)
    {
        pos[2 * i] = rand_const * random[rand_indx++] * L;
        pos[2 * i + 1] = rand_const * random[rand_indx++] * L;
        theta = rand_const * random[rand_indx++] * 6.283185307;
        vel[2 * i] = cos(theta);
        vel[2 * i + 1] = sin(theta);
    }
}

void integrate(int steps, int update_obs)
{
    float *integr = (float *)malloc(2 * N * sizeof(float));
    int *header = (int *)malloc(N_cells * sizeof(int));
    int *cell_list = (int *)malloc(2 * N * sizeof(int));

    // float integr[2 * N];
    // int header[N_cells];
    // int cell_list[2 * N];

    float particle_direction[N];
    int random_lenght = fmax(N4, sfmt_get_min_array_size32(&sfmt));
    uint32_t *random = (uint32_t *)malloc(random_lenght * sizeof(uint32_t));
    const float factor1 = eta * 2.0 * 3.14159265359 * rand_const, float_L = (float)L;
    const float factor2 = -eta * 3.14159265359;

    // TODO: Create an array of pointers do every Y element (another for X) of integr (or vel os pos...)
    // Objective. Clean up the code without (2*i + 1) while having contiguous x/y components.

    for (int step = 0; step < steps; step++)
    {
        for (int i = 0; i < N_cells; i++)
            header[i] = -1;

        for (int i = 0; i < 2 * N; i++)
            integr[i] = vel[i];

        for (int i = 0; i < 2 * N; i += 2)
        {
            int cell = (int)(pos[i]) % L + L * ((int)(pos[i + 1]) % L);
            // printf("%f, %i\n", pos[i], (int)(pos[i]) % L);
            cell_list[i] = header[cell];
            header[cell] = i;
        }

        for (int cell = 0; cell < N_cells; cell++)
        {
            int *nbr_indx = &near_nbr_cells[cell * 4];
            for (int i = header[cell]; i > -1; i = cell_list[i])
            {

                for (int j = cell_list[i]; j > -1; j = cell_list[j])
                {
                    if (dist_simple(&pos[i], &pos[j]) < th)
                    {
                        integr[i] += vel[j];
                        integr[i + 1] += vel[j + 1];
                        integr[j] += vel[i];
                        integr[j + 1] += vel[i + 1];
                    }
                }

                for (int j = header[*nbr_indx++]; j > -1; j = cell_list[j])
                {
                    if (dist_PBC(&pos[i], &pos[j]) < th)
                    {
                        integr[i] += vel[j];
                        integr[i + 1] += vel[j + 1];
                        integr[j] += vel[i];
                        integr[j + 1] += vel[i + 1];
                    }
                }

                for (int j = header[*nbr_indx++]; j > -1; j = cell_list[j])
                {
                    if (dist_PBC(&pos[i], &pos[j]) < th)
                    {
                        integr[i] += vel[j];
                        integr[i + 1] += vel[j + 1];
                        integr[j] += vel[i];
                        integr[j + 1] += vel[i + 1];
                    }
                }

                for (int j = header[*nbr_indx++]; j > -1; j = cell_list[j])
                {
                    if (dist_PBC(&pos[i], &pos[j]) < th)
                    {
                        integr[i] += vel[j];
                        integr[i + 1] += vel[j + 1];
                        integr[j] += vel[i];
                        integr[j + 1] += vel[i + 1];
                    }
                }

                for (int j = header[*nbr_indx]; j > -1; j = cell_list[j])
                {
                    if (dist_PBC(&pos[i], &pos[j]) < th)
                    {
                        integr[i] += vel[j];
                        integr[i + 1] += vel[j + 1];
                        integr[j] += vel[i];
                        integr[j + 1] += vel[i + 1];
                    }
                }
                nbr_indx -= 3;
            }
        }

        sfmt_fill_array32(&sfmt, random, random_lenght);

        for (int i = 0; i < N; i++)

            particle_direction[i] = atan2(integr[2 * i + 1], integr[2 * i]) + factor1 * random[i] + factor2;

        for (int i = 0; i < N; i++)
        {
            vel[2 * i] = cos(particle_direction[i]);
            vel[2 * i + 1] = sin(particle_direction[i]);
        }

        for (int i = 0; i < 2 * N; i++)
        {
            pos[i] = wrap3(pos[i] + v0 * vel[i]);
        }

        if (update_obs)
            update_observables();
    }
    free(header);
    free(cell_list);
    free(integr);
    free(random);
}

void set_geometry()
{
    int cell, cell_X, cell_Y, nbr_cell_X, nbr_cell_Y, nbr_cell, k;
    int nbr_X[] = {-1, 0, 1, -1};
    int nbr_Y[] = {1, 1, 1, 0};

    for (cell_Y = 0; cell_Y < L; cell_Y++)
    {
        for (cell_X = 0; cell_X < L; cell_X++)
        {
            cell = cell_X + L * cell_Y;
            for (k = 0; k < 4; k++)
            {
                nbr_cell_X = cell_X + nbr_X[k];
                nbr_cell_Y = cell_Y + nbr_Y[k];

                // nbr_cell_X = (nbr_cell_X > 0 ? (nbr_cell_X >= L ? nbr_cell_X - L : nbr_cell_X) : nbr_cell_X + L);
                // nbr_cell_Y = (nbr_cell_Y > 0 ? (nbr_cell_Y >= L ? nbr_cell_Y - L : nbr_cell_Y) : nbr_cell_Y + L);

                if (nbr_cell_X < 0)
                    nbr_cell_X += L;

                if (nbr_cell_X >= L)
                    nbr_cell_X -= L;

                if (nbr_cell_Y < 0)
                    nbr_cell_Y += L;

                if (nbr_cell_Y >= L)
                    nbr_cell_Y -= L;

                nbr_cell = nbr_cell_X + L * nbr_cell_Y;

                near_nbr_cells[4 * cell + k] = nbr_cell;
            }
        }
    }
}

int main(int argc, char *argv[])
{
    float phi, sigma_phi, xi_phi;
    int N_reset, N_steps, seed;

    FILE *in_file = fopen(argv[1], "r"); // read only
    if (fscanf(in_file, "%d %*s\n%f %*s\n%f %*s\n%i %*s\n%i %*s\n%i %*s\n", &L, &v0, &rho, &N_reset, &N_steps, &seed) != 6)
        printf("ERROR on input!");
    fclose(in_file);

    // printf("%d \n%f \n%f \n%i \n%i \n%i \n", L, v0, rho, N_reset, N_steps, seed);

    sfmt_init_gen_rand(&sfmt, seed);

    N = (int)((float)L * (float)L * rho);
    N4 = 4 * (N / 4) + 4;
    N_cells = L * L;
    inv_2L = (2.0 / (float)(L));

    pos = (float *)malloc(2 * N * sizeof(float));
    vel = (float *)malloc(2 * N * sizeof(float));
    near_nbr_cells = (int *)malloc(N_cells * 4 * sizeof(int));

    randomize_system();
    set_geometry();
    get_and_reset_observables(&phi, &sigma_phi, &xi_phi);

    eta = 0.0;

    integrate(N_reset, 0);

    for (eta = 0.0; eta < 1.01; eta += 0.05)
    {
        integrate(N_reset, 0);
        integrate(N_steps, 1);
        get_and_reset_observables(&phi, &sigma_phi, &xi_phi);

        printf("%f, %f, %f, %f\n", eta, phi, sigma_phi, xi_phi);
    }
}
