#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

// #include "mt19937-64.c"
// #define RAN (genrand64_real2())

// #include "mt19937ar.c"
// #define RAN (genrand_real2())

#include "SFMT-src-1.5.1/SFMT.h"
sfmt_t sfmt;
#define RAN (sfmt_genrand_real2(&sfmt))
#define sq(value) (value) * (value)

using namespace std;

// const float rho{1.0}, v0{0.03};
// const int L{15}, N_reset{1000}, N_steps{10000}, seed{1234};
// const int N_cells = L * L, N{int(float(L) * float(L) * rho)};
// const float inv_L = 1.0 / float(L), inv_2L = 2.0 / float(L);
// float eta;

// const float rho{1.0}, v0{0.03};
// const int L{15}, N_reset{1000}, N_steps{10000}, seed{1234};
// const int N_cells = L * L, N{int(float(L) * float(L) * rho)};
// const float inv_L = 1.0 / float(L), inv_2L = 2.0 / float(L);
// float eta;

float pbc_dist(float x, int L, float inv_2L)
{
    return x - int(inv_2L * x) * L;
}

float dist_PBC(float x1, float y1, float x2, float y2, int L, float inv_2L)
{
    return sq(pbc_dist(x1 - x2, L, inv_2L)) + sq(pbc_dist(y1 - y2, L, inv_2L));
}

float dist_simple(float x1, float y1, float x2, float y2)
{
    return sq(x1 - x2) + sq(y1 - y2);
}

// float dist_PBC(float *pos1, float *pos2)
// {
//     // cout << *pos1 << endl;
//     float dx = pbc_dist((*pos1++) - (*pos2++));
//     float dy = pbc_dist(*pos1 - *pos2);

//     return dx * dx + dy * dy;
// }

// float dist_simple(float *pos1, float *pos2)
// {
//     float dx = (*pos1++) - (*pos2++);
//     float dy = *pos1 - *pos2;

//     return dx * dx + dy * dy;
// }
void set_geometry(int L, int (*near_nbr_cells)[4])
{
    int cell, cell_X, cell_Y, nbr_cell_X, nbr_cell_Y, nbr_cell, k;
    int nbr_X[4] = {-1, 0, 1, -1};
    int nbr_Y[4] = {1, 1, 1, 0};

    for (cell_Y = 0; cell_Y < L; cell_Y++)
    {
        for (cell_X = 0; cell_X < L; cell_X++)
        {
            cell = cell_X + L * cell_Y;
            for (k = 0; k < 4; k++)
            {
                nbr_cell_X = cell_X + nbr_X[k];
                nbr_cell_Y = cell_Y + nbr_Y[k];

                if (nbr_cell_X < 0)
                    nbr_cell_X += L;

                if (nbr_cell_X >= L)
                    nbr_cell_X -= L;

                if (nbr_cell_Y < 0)
                    nbr_cell_Y += L;

                if (nbr_cell_Y >= L)
                    nbr_cell_Y -= L;

                nbr_cell = nbr_cell_X + L * nbr_cell_Y;

                near_nbr_cells[cell][k] = nbr_cell;
            }
        }
    }
}

// int where_is(float x, float y)
// {
//     return;
// }

void neighbours_direction(float inv_2L, int L, int N, int N_cells, float eta, float *x, float *y, float *vx, float *vy, int (*near_nbr_cells)[4], float *nbr_direction)
{
    float integr_vx[N], integr_vy[N];
    int header[N_cells], cell_list[2 * N];

    fill_n(header, N_cells, -1);

    // for (int i = 0; i < 2 * N; i++)
    // {
    //     integr_vx[i] = unit_vel[i];
    // }

    // copy(vx, vx + N, integr_vx);
    // copy(vy, vy + N, integr_vy);

    for (int i = 0; i < N; i++)
    {
        int cell = int(x[i]) % L + L * (int(y[i]) % L);
        cell_list[i] = header[cell];
        header[cell] = i;
    }

#pragma clang loop vectorize(assume_safety)
    for (int i = 0; i < N; i++)
    {
        integr_vx[i] = vx[i];
        integr_vy[i] = vy[i];
    }

    for (int cell = 0; cell < N_cells; cell++)
    {
        for (int i = header[cell]; i > -1; i = cell_list[i])
        {
            for (int j = cell_list[i]; j > -1; j = cell_list[j])
            {
                if (dist_simple(x[i], y[i], x[j], y[j]) < 1.0)
                {
                    integr_vx[i] = integr_vx[i] + vx[j];
                    integr_vy[i] = integr_vy[i] + vy[j];
                    integr_vx[j] = integr_vx[j] + vx[i];
                    integr_vy[j] = integr_vy[j] + vy[i];
                }
            }

            for (int nbr_cell : near_nbr_cells[cell])
            {
                for (int j = header[nbr_cell]; j > -1; j = cell_list[j])
                {
                    if (dist_PBC(x[i], y[i], x[j], y[j], L, inv_2L) < 1.0)
                    {
                        integr_vx[i] = integr_vx[i] + vx[j];
                        integr_vy[i] = integr_vy[i] + vy[j];
                        integr_vx[j] = integr_vx[j] + vx[i];
                        integr_vy[j] = integr_vy[j] + vy[i];
                    }
                }
            }
        }
    }

    for (int i = 0; i < N; i++)
    {
        nbr_direction[i] = atan2(integr_vy[i], integr_vx[i]) + eta * 6.283185307179586 * (RAN - 0.5);
    }
}

// float wrap1(float x)
// {
//     return (x > 0 ? 0 : L) + fmodf(x, L);
// }

// float wrap2(float x)
// {
//     return x - floor(x * inv_L) * L;
// }

float wrap3(float x, int L)
{
    return (x > 0 ? (x > L ? x - L : x) : x + L);
}

// float wrap4(float x)
// {
//     return x - lrint(x * inv_L) * L;
// }

void integrate(float inv_2L, int N, float v0, int L, float eta, float *x, float *y, float *vx, float *vy, int (*near_nbr_cells)[4])
{
    float particle_direction[N];

    neighbours_direction(inv_2L, L, N, L * L, eta, x, y, vx, vy, near_nbr_cells, particle_direction);

    // cout << "here\n";
    for (int i = 0; i < N; i++)
    {
        vx[i] = cos(particle_direction[i]);
        vy[i] = sin(particle_direction[i]);
        x[i] = wrap3(x[i] + v0 * vx[i], L);
        y[i] = wrap3(y[i] + v0 * vy[i], L);
    }
}

int main(int argc, char *argv[])
{
    float polar, polar2, sum_x, sum_y;
    // cout << N << endl;
    // int N_reset, N_steps, seed;

    // init_genrand(123);
    sfmt_init_gen_rand(&sfmt, 12345);
    // cout << sfmt_get_min_array_size32(&sfmt) << ' ' << sfmt_get_min_array_size64(&sfmt) << endl;
    // uniform_real_distribution<double> dis(0.0, 1.0);

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

    const int N_cells = L * L, N = int(float(L) * float(L) * rho);
    const float inv_L = 1.0 / float(L), inv_2L = 2.0 / float(L);
    float eta;
    // cout << N << endl;

    // float xy[N * 2], unit_vel[N * 2];
    float x[N], y[N], vx[N], vy[N];
    int near_nbr_cells[N_cells][4];

    for (int i = 0; i < N; i++)
    {
        x[i] = RAN * L;
        y[i] = RAN * L;
        float theta = RAN * 6.283185307;
        vx[i] = cos(theta);
        vy[i] = sin(theta);
    }

    set_geometry(L, near_nbr_cells);

    eta = 0.0;

    for (int step = 0; step < N_reset; step++)
        integrate(inv_2L, N, v0, L, eta, x, y, vx, vy, near_nbr_cells);

    for (eta = 0.0; eta < 1.01; eta += 0.05)
    {

        for (int step = 0; step < N_reset; step++)
            integrate(inv_2L, N, v0, L, eta, x, y, vx, vy, near_nbr_cells);

        polar = 0.0;
        polar2 = 0.0;

        for (int step = 0; step < N_steps; step++)
        {
            sum_x = 0.0;
            sum_y = 0.0;

            integrate(inv_2L, N, v0, L, eta, x, y, vx, vy, near_nbr_cells);

            for (int i = 0; i < N; i++)
            {
                sum_x += vx[i];
                sum_y += vy[i];
            }

            float polar_i_sq = sum_x * sum_x + sum_y * sum_y;
            polar += sqrt(polar_i_sq);
            polar2 += polar_i_sq;
        }

        cout << eta << ','
             << polar / (float(N) * float(N_steps)) << ','
             << sqrt(polar2 / float(N_steps) - (polar * polar) / float(N_steps * N_steps)) / float(N) << ','
             << (polar2 - (polar * polar) / float(N_steps)) / polar << endl;
    }
}

// extension BinaryFloatingPoint
// {
//     /// Returns the value after restricting it to the periodic boundary
//     /// condition [0, 1).
//     /// See https://forums.swift.org/t/why-no-fraction-in-floatingpoint/10337
//     @_transparent
//         func
//         wrappedToUnitRange()
//             ->Self
//     {
//         let fract = self - self.rounded(.down)
//                            // Have to clamp to just below 1 because very small negative values
//                            // will otherwise return an out of range result of 1.0.
//                            // Turns out this:
//                            if fract >=
//                     1.0
//         {
//             return Self(1).nextDown
//         }
//         else { return fract }
//         // is faster than this:
//         // return min(fract, Self(1).nextDown)
//     }
//     @_transparent
//         func
//         wrapped(to range
//                 : Range<Self>)
//             ->Self
//     {
//         let measure = range.upperBound
//                           let recipMeasure = Self(1) / measure
//                                                            let scaled = (self)*recipMeasure return scaled.wrappedToUnitRange() * measure
//     }
//     @_transparent
//         func
//         wrappedIteratively(to range
//                            : Range<Self>)
//             ->Self
//     {
//         var v = self
//             let measure = range.upperBound while v >= range.upperBound { v = v - measure }
//         while
//             v < 0 { v = v + measure }
//         return v
//     }
// }
