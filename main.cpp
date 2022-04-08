#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <numeric>
#include <vector>
#include <algorithm>

// #include "mt19937-64.c"
// #define RAN (genrand64_real2())

// #include "mt19937ar.c"
// #define RAN (genrand_real2())

#include "SFMT-src-1.5.1/SFMT.h"

using namespace std;

class vicsek_system
{
private:
    float v0, inv_2L, th = 1.0, rand_const = (1.0 / 4294967296.0);
    const int L, N_cells, N;
    // float *x, *y, *vx, *vy;
    // vector<float> x, y;
    // vector<vector<float>> vel;, pos;
    float *pos, *vel;
    int *near_nbr_cells;
    sfmt_t sfmt;

    inline float pbc(float value)
    {
        return value - int(inv_2L * value) * L;
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

public:
    float eta;
    float sum_phis, sum_phis2;
    int N_observations;

    void get_and_reset_observables(float *phi, float *sigma_phi, float *xi_phi)
    {
        *phi = sum_phis / (float(N) * float(N_observations));
        *sigma_phi = sqrt(sum_phis2 / float(N_observations) - (sum_phis * sum_phis) / float(N_observations * N_observations)) / float(N);
        *xi_phi = (sum_phis2 - (sum_phis * sum_phis) / float(N_observations)) / sum_phis;

        sum_phis2 = 0.0;
        sum_phis = 0.0;
        N_observations = 0;
    }

    vicsek_system(float rho_in,
                  float v0_in,
                  int L_in,
                  int seed) : v0(v0_in),
                              L(L_in),
                              N_cells(L_in * L_in),
                              N(int((L_in * L_in * rho_in + 2.0) / 4) * 4),
                              inv_2L(2.0 / float(L_in))
    {
        pos = (float *)malloc(2 * N * sizeof(float));
        vel = (float *)malloc(2 * N * sizeof(float));

        // vel.resize(N);
        // for (int i = 0; i < N; i++)
        //     vel[i].resize(2);

        // pos.resize(N);
        // for (int i = 0; i < N; i++)
        //     pos[i].resize(2);

        // vy.resize(N);
        // cout << sfmt_get_min_array_size32(&sfmt) << endl;
        near_nbr_cells = (int *)malloc(N_cells * 4 * sizeof(int));
        sfmt_init_gen_rand(&sfmt, seed);
        randomize_system();
        set_geometry();
        float junk;
        get_and_reset_observables(&junk, &junk, &junk);
    }

    ~vicsek_system()
    {
        // free(x);
        // free(y);
        // free(vx);
        // free(vy);
        free(pos);
        free(vel);
        free(near_nbr_cells);
    }

    void update_observables()
    {
        // float sum_x = accumulate(vx, vx + N, 0.0);
        // float sum_y = accumulate(vy, vy + N, 0.0);

        float sum_x = 0, sum_y = 0;

        for (int i = 0; i < 2 * N;)
        {
            sum_x += vel[i++];
            sum_y += vel[i++];
        }
        // float sum_x = accumulate(begin(vx), end(vx), 0.0);
        // float sum_y = accumulate(begin(vy), end(vy), 0.0);
        float polar_i_sq = sum_x * sum_x + sum_y * sum_y;

        sum_phis += sqrt(polar_i_sq);
        sum_phis2 += polar_i_sq;
        N_observations++;
    }

    void set_geometry()
    {
        int cell, cell_X, cell_Y, nbr_cell_X, nbr_cell_Y, nbr_cell, k;
        int nbr_X[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
        int nbr_Y[9] = {1, 1, 1, 0, 0, 0, -1, -1, -1};

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

                    near_nbr_cells[cell * 4 + k] = nbr_cell;
                }
            }
        }
    }

    void randomize_system()
    {

        int random_lenght = max(3 * N, sfmt_get_min_array_size32(&sfmt));

        uint32_t random[random_lenght];

        sfmt_fill_array32(&sfmt, random, random_lenght);
        int rand_indx = 0;
        for (int i = 0; i < N; i++)
        {
            pos[2 * i] = rand_const * random[rand_indx++] * L;
            pos[2 * i + 1] = rand_const * random[rand_indx++] * L;
            float theta = rand_const * random[rand_indx++] * 6.283185307;
            vel[2 * i] = cos(theta);
            vel[2 * i + 1] = sin(theta);
        }
    }

    void integrate(int steps = 1, int update_obs = 0)
    {

        // vector<vector<float, 2>, N> integr;
        float *integr = (float *)malloc(2 * N * sizeof(float));
        int *header = (int *)malloc(N_cells * sizeof(int));
        int *cell_list = (int *)malloc(2 * N * sizeof(int));

        // float integr[2 * N];
        // int header[N_cells];
        // int cell_list[2 * N];
        // vector<float> integr(2 * N);
        // vector<int> header(N_cells);
        // vector<int> cell_list(2 * N);

        float particle_direction[N];
        int random_lenght = max(N, sfmt_get_min_array_size32(&sfmt));
        uint32_t *random = (uint32_t *)malloc(random_lenght * sizeof(uint32_t));
        float factor1 = eta * 6.283185307179586 * rand_const;
        float factor2 = -eta * 3.14159265359;

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
                int cell = int(pos[i]) % L + L * (int(pos[i + 1]) % L);
                cell_list[i] = header[cell];
                header[cell] = i;
            }

            for (int cell = 0; cell < N_cells; cell++)
            {
                int *nbr_indx = near_nbr_cells + cell * 4;
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
};

int main()
{
    float rho, v0;
    int L, N_reset, N_steps, seed;
    float phi, sigma_phi, xi_phi;

    fstream myfile("input.dat");
    string line;
    myfile >> L >> line;
    myfile >> v0 >> line;
    myfile >> rho >> line;
    myfile >> N_reset >> line;
    myfile >> N_steps >> line;
    myfile >> seed >> line;
    myfile.close();

    vicsek_system system(rho, v0, L, seed);

    system.eta = 0.0;

    system.integrate(N_reset);

    for (system.eta = 0.0; system.eta < 1.01; system.eta += 0.05)
    {
        system.integrate(N_reset);

        system.integrate(N_steps, 1);

        system.get_and_reset_observables(&phi, &sigma_phi, &xi_phi);

        cout << system.eta << ',' << phi << ',' << sigma_phi << ',' << xi_phi << endl;
    }
}
