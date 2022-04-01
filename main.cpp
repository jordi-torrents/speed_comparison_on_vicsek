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
#define RAN (sfmt_genrand_real2(&sfmt))
private:
    const float v0, inv_2L, th = 1.0, rand_const = (1.0 / 4294967296.0);
    // float *x, *y, *vx, *vy;
    vector<float> x, y;
    vector<float> vx, vy;
    int *near_nbr_cells;
    sfmt_t sfmt;

    const float pbc(float value)
    {
        return value - int(inv_2L * value) * L;
    }

    const float dist_PBC(float x1, float y1, float x2, float y2)
    {
        float dx = pbc(x1 - x2);
        float dy = pbc(y1 - y2);
        return dx * dx + dy * dy;
    }

    const float dist_simple(float x1, float y1, float x2, float y2)
    {
        float dx = x1 - x2;
        float dy = y1 - y2;
        return dx * dx + dy * dy;
    }

    const float wrap3(float value)
    {
        return (value > 0 ? (value > L ? value - L : value) : value + L);
    }

public:
    const int L, N_cells, N;
    float eta;
    float sum_phis, sum_phis2;
    int N_observations;

    auto get_and_reset_observables()
    {
        struct retVals
        {
            float f1, f2, f3;
        };
        float phi = sum_phis / (float(N) * float(N_observations));
        float sigma_phi = sqrt(sum_phis2 / float(N_observations) - (sum_phis * sum_phis) / float(N_observations * N_observations)) / float(N);
        float xi_phi = (sum_phis2 - (sum_phis * sum_phis) / float(N_observations)) / sum_phis;

        sum_phis2 = 0.0;
        sum_phis = 0.0;
        N_observations = 0;

        return retVals{phi, sigma_phi, xi_phi};
    }

    vicsek_system(float rho_in, float v0_in, int L_in, int seed) : v0(v0_in),
                                                                   L(L_in),
                                                                   N_cells(L_in * L_in),
                                                                   N(int(L_in * L_in * rho_in)),
                                                                   inv_2L(2.0 / float(L))
    {
        // x = (float *)malloc(N * sizeof(float));
        // y = (float *)malloc(N * sizeof(float));
        // vx = (float *)malloc(N * sizeof(float));
        // vy = (float *)malloc(N * sizeof(float));

        x.resize(N);
        y.resize(N);
        vx.resize(N);
        vy.resize(N);

        near_nbr_cells = (int *)malloc(N_cells * 4 * sizeof(int));
        sfmt_init_gen_rand(&sfmt, seed);
        randomize_system();
        set_geometry();
        get_and_reset_observables();
    }

    ~vicsek_system()
    {
        // free(x);
        // free(y);
        // free(vx);
        // free(vy);
        free(near_nbr_cells);
    }

    void update_observables()
    {
        // float sum_x = accumulate(vx, vx + N, 0.0);
        // float sum_y = accumulate(vy, vy + N, 0.0);
        float sum_x = accumulate(begin(vx), end(vx), 0.0);
        float sum_y = accumulate(begin(vy), end(vy), 0.0);
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
        // uint32_t random1[N];
        // uint32_t random2[N];
        // uint32_t random3[N];
        // sfmt_fill_array32(&sfmt, random1, N);
        // sfmt_fill_array32(&sfmt, random2, N);
        // sfmt_fill_array32(&sfmt, random3, N);

        for (int i = 0; i < N; i++)
        {
            x[i] = RAN * L;
            y[i] = RAN * L;
            float theta = RAN * 6.283185307;
            vx[i] = cos(theta);
            vy[i] = sin(theta);
        }
    }

    void neighbours_direction(float *nbr_direction)
    {

        vector<float> integr_vx(vx), integr_vy(vy);
        vector<int> header(N_cells, -1), cell_list(N);

        // copy(vx.begin(), vx.end(), integr_vx);
        // copy(vx.begin(), vx.end(), integr_vx);
        fill(header.begin(), header.end(), -1);

        for (int i = 0; i < N; i++)
        {
            int cell = int(x[i]) % L + L * (int(y[i]) % L);
            cell_list[i] = header[cell];
            header[cell] = i;
            integr_vx[i] = vx[i];
            integr_vy[i] = vy[i];
        }

        for (int cell = 0; cell < N_cells; cell++)
        {
            int *nbr_indx = near_nbr_cells + cell * 4;
            for (int i = header[cell]; i > -1; i = cell_list[i])
            {

                for (int j = cell_list[i]; j > -1; j = cell_list[j])
                {
                    if (dist_simple(x[i], y[i], x[j], y[j]) < th)
                    {
                        integr_vx[i] += vx[j];
                        integr_vy[i] += vy[j];
                        integr_vx[j] += vx[i];
                        integr_vy[j] += vy[i];
                    }
                }

                for (int j = header[*nbr_indx++]; j > -1; j = cell_list[j])
                {
                    if (dist_PBC(x[i], y[i], x[j], y[j]) < th)
                    {
                        integr_vx[i] += vx[j];
                        integr_vy[i] += vy[j];
                        integr_vx[j] += vx[i];
                        integr_vy[j] += vy[i];
                    }
                }

                for (int j = header[*nbr_indx++]; j > -1; j = cell_list[j])
                {
                    if (dist_PBC(x[i], y[i], x[j], y[j]) < th)
                    {
                        integr_vx[i] += vx[j];
                        integr_vy[i] += vy[j];
                        integr_vx[j] += vx[i];
                        integr_vy[j] += vy[i];
                    }
                }

                for (int j = header[*nbr_indx++]; j > -1; j = cell_list[j])
                {
                    if (dist_PBC(x[i], y[i], x[j], y[j]) < th)
                    {
                        integr_vx[i] += vx[j];
                        integr_vy[i] += vy[j];
                        integr_vx[j] += vx[i];
                        integr_vy[j] += vy[i];
                    }
                }

                for (int j = header[*nbr_indx]; j > -1; j = cell_list[j])
                {
                    if (dist_PBC(x[i], y[i], x[j], y[j]) < th)
                    {
                        integr_vx[i] += vx[j];
                        integr_vy[i] += vy[j];
                        integr_vx[j] += vx[i];
                        integr_vy[j] += vy[i];
                    }
                }
                nbr_indx -= 3;
            }
        }

        // float random[N];
        // for (int i = 0; i < N; i++)
        // {
        //     random[i] = RAN;
        // }
        float factor1 = eta * 6.283185307179586;
        float factor2 = -eta * 3.14159265359;

        // uint32_t random[N];
        // sfmt_fill_array32(&sfmt, random, N);
        // float factor1 = eta * 6.283185307179586 / 4294967296.0;
        // float factor2 = -eta * 3.14159265359;

        for (int i = 0; i < N; i++)
        {
            nbr_direction[i] = atan2(integr_vy[i], integr_vx[i]) + factor1 * RAN + factor2;
        }
    }

    void integrate(int steps = 1, int update_obs = 0)
    {
        float particle_direction[N];

        for (int step = 0; step < steps; step++)
        {
            neighbours_direction(particle_direction);

            for (int i = 0; i < N; i++)
            {
                vx[i] = cos(particle_direction[i]);
                vy[i] = sin(particle_direction[i]);
                x[i] = wrap3(x[i] + v0 * vx[i]);
                y[i] = wrap3(y[i] + v0 * vy[i]);
            }

            if (update_obs)
                update_observables();
        }
    }
};

int main()
{
    float rho, v0;
    int L, N_reset, N_steps, seed;

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

        auto [phi, sigma_phi, xi_phi] = system.get_and_reset_observables();

        cout << system.eta << ',' << phi << ',' << sigma_phi << ',' << xi_phi << endl;
    }
}
