#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// Mersenne twister

#include "mt19937-64.c"
#define RAN (genrand64_real2())

#define V0 0.02

int main(int argc, char *argv[])
{
    int iarg, N, L;
    double eta, v0, rho;
    char filename[200];
    FILE *fp;
    double phi, phi2, phi4, variance;
    double G, var;
    long int seed;
    int therm, simul_time, max_time;

    for (iarg = 1; iarg < argc; iarg++)
    {
        if (!strncmp(argv[iarg], "N", 1))
            sscanf(argv[iarg], "%*[^=]%*1c%d", &N);
        if (!strncmp(argv[iarg], "seed", 4))
            sscanf(argv[iarg], "%*[^=]%*1c%ld", &seed);
        if (!strncmp(argv[iarg], "eta", 3))
            sscanf(argv[iarg], "%*[^=]%*1c%lf", &eta);
        if (!strncmp(argv[iarg], "rho", 3))
            sscanf(argv[iarg], "%*[^=]%*1c%lf", &rho);
        if (!strncmp(argv[iarg], "v0", 2))
            sscanf(argv[iarg], "%*[^=]%*1c%lf", &v0);
        if (!strncmp(argv[iarg], "therm", 5))
            sscanf(argv[iarg], "%*[^=]%*1c%d", &therm);
        if (!strncmp(argv[iarg], "time", 4))
            sscanf(argv[iarg], "%*[^=]%*1c%d", &simul_time);
        if (!strncmp(argv[iarg], "output", 6))
            sscanf(argv[iarg], "%*[^=]%*1c%s", filename);
    }

    init_genrand64(seed);

    max_time = simul_time + therm;

    L = (int)sqrt((double)N / (double)rho);

    fprintf(stderr, "System size %d(%.4f), number of particles %d\n", L, sqrt((double)N / (double)rho), N);

    vicsek_lattice_angular_verlet_optim_averages(eta, max_time, therm,
                                                 L, N, v0, &phi, &phi2, &phi4);

    variance = phi2 - phi * phi;
    var = variance * (double)L * (double)L;
    G = 1.0 - phi4 / (3.0 * phi2 * phi2);

    if ((fp = fopen(filename, "a")) == NULL)
        exit(1);

    fprintf(fp, "%.5f  %.12f  %.12f  %.12f  %.12f  %.12f\n", v0, eta, phi, var, var / phi, G);

    fclose(fp);
}

/*
 *    Computes first, second and fourth moment of the order parameter, to enable
 *    computation of the ratio to measure a first order transition
 *
 */

void vicsek_lattice_angular_verlet_optim_averages(
    double eta,
    int max_time,
    int therm,
    int L,
    int N,
    double v_0,
    double *tot_phi,
    double *tot_phi2,
    double *tot_phi4)
{
    int i, t, j, pos;
    double *vx, *vy, *x, *y, theta, posx, posy, xn, yn, dL;
    double phi_x, phi_y;
    int M, *myBox;
    double *vbar_x, *vbar_y;
    int xp, yp;
    int *nnpos;
    double norm, tp = 0.0, tv = 0.0, t4 = 0.0;
    int num_stats = 0;

    nnpos = malloc(NNBOXES * sizeof(int));

    dL = (double)L;

    // velocity
    vx = malloc(N * sizeof(double));
    vy = malloc(N * sizeof(double));

    // position
    x = malloc(N * sizeof(double));
    y = malloc(N * sizeof(double));

    // box position
    myBox = malloc(N * sizeof(double));

    // box size with R=1
    M = L * L;

    // vbar velocities of each box
    vbar_x = malloc(M * sizeof(double));
    vbar_y = malloc(M * sizeof(double));

    fprintf(stderr, "Set up finished\n");

    for (i = 0; i < M; i++)
        vbar_x[i] = vbar_y[i] = 0.0;

    // initial conditions
    for (i = 0; i < N; i++)
    {
        // velocities: initial conditions in -pi : +pi
        theta = M_PI * (2.0 * RAN - 1.0);
        vx[i] = sin(theta);
        vy[i] = cos(theta);

        // positions: random in [0,L] x [0, L]
        x[i] = RAN * L;
        y[i] = RAN * L;

        // allocate particle in one box
        xp = (int)x[i];
        yp = (int)y[i];

        pos = xp + yp * L;

        myBox[i] = pos;

        neighboring_boxes(pos, nnpos, L);

        for (j = 0; j < NNBOXES; j++)
        {
            pos = nnpos[j];
            vbar_x[pos] += vx[i];
            vbar_y[pos] += vy[i];
        }
    }

    for (t = 1; t < max_time; t++)
    {

        // update first velocities
        phi_x = phi_y = 0.0;
        for (i = 0; i < N; i++)
        {
            // retrieve the box of the particle
            pos = myBox[i];

            theta = atan2(vbar_x[pos], vbar_y[pos]) + eta * M_PI * (2.0 * RAN - 1.0);

            posx = sin(theta);
            vx[i] = posx;
            phi_x += posx;

            posy = cos(theta);
            vy[i] = posy;
            phi_y += posy;
        }

        if (t > therm)
        {
            norm = sqrt(phi_x * phi_x + phi_y * phi_y) / (double)N;
            tp += norm;
            tv += (norm * norm);
            t4 += (norm * norm * norm * norm);
            num_stats++;
        }

        // update positions using the new velocity

        for (i = 0; i < M; i++)
            vbar_x[i] = vbar_y[i] = 0.0;

        for (i = 0; i < N; i++)
        {
            xn = x[i] + v_0 * vx[i];
            yn = y[i] + v_0 * vy[i];

            // periodic boundary conditions
            if (xn < 0.0)
                xn += dL;
            else if (xn > dL)
                xn -= dL;

            if (yn < 0.0)
                yn += dL;
            else if (yn > dL)
                yn -= dL;

            // Check for position inside the [0,L] x [0,L] box

            if (xn < 0.0 || xn > dL || yn < 0.0 || yn > dL)
            {
                fprintf(stderr, "Particles out of the box\n");
                exit(0);
            }

            // update box velocities
            xp = (int)xn;
            yp = (int)yn;

            pos = xp + yp * L;

            myBox[i] = pos;

            neighboring_boxes(pos, nnpos, L);
            for (j = 0; j < NNBOXES; j++)
            {
                pos = nnpos[j];
                vbar_x[pos] += vx[i];
                vbar_y[pos] += vy[i];
            }

            x[i] = xn;
            y[i] = yn;
        }
    }

    tp /= (double)num_stats;
    tv /= (double)num_stats;
    t4 /= (double)num_stats;

    *tot_phi = tp;
    *tot_phi2 = tv;
    *tot_phi4 = t4;

    free(vx);
    free(vy);
    free(x);
    free(y);
    free(vbar_x);
    free(vbar_y);
    free(myBox);
}
