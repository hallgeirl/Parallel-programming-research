/*
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 *
 * Modified and  restructured by Scott B. Baden, UCSD
 *
 */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "time.h"
#include "apf.h"
#include "types.h"

int solve(DOUBLE ***_E, DOUBLE ***_E_prev, DOUBLE **R, int m, int n, DOUBLE T, int iters, DOUBLE alpha, DOUBLE dt, int tx, int ty, int do_stats, int plot_freq) 
{
    // Simulated time is different from the integer timestep number
    DOUBLE t = 0.0;
    // Integer timestep number
    int niter=0;

    DOUBLE **E = *_E, **E_prev = *_E_prev;
    int *chunksizes = malloc(sizeof(int)*ty);
    int i, j, ti;

    for (i = 0; i < ty; i++)
    {
        chunksizes[i] = (m+1)/ty;
        if (i == ty-1) chunksizes[i] += (m+1)%ty;
#ifdef NO_GHOST_CELLS
        chunksizes[i]--;
#endif
        printf("Thread %i workload: %i\n", i, chunksizes[i]);
    }


    // We continue to sweep over the mesh until the simulation has reached
    // the desired simulation Time
    // This is different from the number of iterations
    while ((iters < 0 && t < T) || niter < iters) {
#ifdef DEBUG
        //printMat(E_prev, m+3, n+3);
        repNorms(E_prev, t, dt, m, n, niter);
        if (plot_freq) {
            splot(E_prev, t, niter, m + 1, n + 1, WAIT);
        }
#endif

        t += dt;
        niter++;

        // Solve for the excitation, a PDE
        // Also make sure that each CPU works on a block on a time instead of a larger unit of data.

#pragma omp parallel for private(i,j) schedule(static, 1)
        for (ti = 0; ti < ty; ti++)
        {
            //Let top and bottom thread copy top and bottom ghost cells
            if (ti == 0)
                memcpy(&E_prev[0][0], &E_prev[2][0], sizeof(DOUBLE)*(n+3));
            if (ti == ty-1)
                memcpy(&E_prev[m+2][0], &E_prev[m][0], sizeof(DOUBLE)*(n+3));

            int ii = ((m+1)/ty)*ti + 1;
#ifdef NO_GHOST_CELLS
            ii++;
#endif
            for (i = ii; i < ii+chunksizes[ti]; i++) 
            {
                //Copy left and right ghost cells
                E_prev[i][0] = E_prev[i][2];
                E_prev[i][n + 2] = E_prev[i][n];

#pragma ivdep
                for (j = 1; j <= n+1; j++) {
                    E[i][j] = E_prev[i][j] + alpha * (E_prev[i][j + 1]+
                            E_prev[i][j - 1]-
                            4 * E_prev[i][j]+
                            E_prev[i + 1][j]+
                            E_prev[i - 1][j]);
                }
            }
        }

#ifdef DEBUG
        printMat(E_prev, m+3, n+3);
        printf("\n");
#endif

        /*
         * Solve the ODE, advancing excitation and recovery variables
         * to the next timtestep
         */
#pragma omp parallel for private(i,j) schedule(static, 1)
        for (ti = 0; ti < ty; ti++)
        {
            int ii = ((m+1)/ty)*ti + 1;
            for (i = ii; i < ii+chunksizes[ti]; i++)
            {
#pragma ivdep
                for (j = 1; j <= n + 1; j++) 
                {
                    E[i][j] *=
                        1 - dt * (kk * (E[i][j] - a) * (E[i][j] - 1) + R[i][j]);

                    R[i][j] +=
                        dt * (
                                epsilon + M1 * R[i][j] /
                                (E[i][j] + M2)) *
                        (-R[i][j] - kk * E[i][j] * (E[i][j] - b - 1)
                        );
                }
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
#ifndef NO_GHOST_CELLS
        DOUBLE **tmp = E;
        E = E_prev;
        E_prev = tmp;
#endif
    }

    // Store them into the pointers passed in
    *_E = E;
    *_E_prev = E_prev;

    return niter;
}
