/*
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 *
 * Modified and  restructured by Scott B. Baden, UCSD
 *
 */
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "time.h"
#include "apf.h"
#include "types.h"

int solve(DOUBLE ***_E, DOUBLE ***_E_prev, DOUBLE **R, int m, int n, DOUBLE T, int iterations, DOUBLE alpha, DOUBLE dt, int do_stats, int plot_freq, int bx, int by) 
{
    // Simulated time is different from the integer timestep number
    DOUBLE t = 0.0;
    // Integer timestep number
    int niter=0;

    DOUBLE **E = *_E, **E_prev = *_E_prev;

    // We continue to sweep over the mesh until the simulation has reached
    // the desired simulation Time
    // This is different from the number of iterations
    while (niter < iterations || (iterations == -1 && t < T)) {
        #ifdef DEBUG
        printMat(E_prev, m, n);
        repNorms(E_prev, t, dt, m, n, niter);
        if (plot_freq) {
            splot(E_prev, t, niter, m + 1, n + 1, WAIT);
        }
        #endif

        t += dt;
        niter++;

        /*
        * Copy data from boundary of the computational box to the
        * padding region, set up for differencing computational box's boundary
        *
        */
        int i, j;
        int ii, jj;// for blocking.

        #pragma ivdep
        for (j = 1; j <= m + 1; j++) {
            E_prev[j][0] = E_prev[j][2];
            E_prev[j][n + 2] = E_prev[j][n];
        }

        #pragma ivdep
        for (i = 1; i <= n + 1; i++) {
            E_prev[0][i] = E_prev[2][i];
            E_prev[m + 2][i] = E_prev[m][i];
        }


        // Solve for the excitation, a PDE
        // Also make sure that each CPU works on a block on a time instead of a larger unit of data.
        #pragma omp parallel for private(i,j, ii, jj)
        for (jj = 1; jj <= m + 1; jj += by) {
            for (ii = 1; ii <= n + 1; ii += bx) {

                for (j = jj; j < jj+by && j <= m+1; j++) {
                    #pragma ivdep
                    for (i = ii; i < ii+bx && i <= n+1; i++) {
                        E[j][i] = E_prev[j][i] + alpha * (E_prev[j][i + 1]+
                                                  E_prev[j][i - 1]-
                                                  4 * E_prev[j][i]+
                                                  E_prev[j + 1][i]+
                                                  E_prev[j - 1][i]);
                    }
                }
            }
        }

        /*
        * Solve the ODE, advancing excitation and recovery variables
        * to the next timtestep
        */
        #pragma omp parallel for private(i, j)
        for (j = 1; j <= m + 1; j++) 
        {
            for (i = 1; i <= n + 1; i++) 
            {
                E[j][i] *=
                  1 - dt * (kk * (E[j][i] - a) * (E[j][i] - 1) + R[j][i]);

                R[j][i] +=
                    dt * (
                        epsilon + M1 * R[j][i] /
                        (E[j][i] + M2)) *
                        (-R[j][i] - kk * E[j][i] * (E[j][i] - b - 1)
                    );
            }
        }

        if (do_stats) {
            repNorms(E, t, dt, m, n, niter);
        }

        if (plot_freq) {
            int k = (int)(t / plot_freq);
            if ((t - k * plot_freq) < dt) 
            {
                splot(E, t, niter, m + 1, n + 1, WAIT);
            }
        }

        // Swap current and previous
        DOUBLE **tmp = E;
        E = E_prev;
        E_prev = tmp;
    }

    // Store them into the pointers passed in
    *_E = E;
    *_E_prev = E_prev;

    return niter;
}
