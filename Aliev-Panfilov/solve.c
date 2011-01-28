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
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "time.h"
#include "apf.h"
#include "types.h"

int solve(block_t ** blocks, int m, int n, DOUBLE T, DOUBLE alpha, DOUBLE dt, int do_stats, int plot_freq, int tx, int ty, bool enableGhostCells, int iterations) 
{
    // Simulated time is different from the integer timestep number
    DOUBLE t = 0.0;
    // Integer timestep number
    int niter=0;
    int i, j, k;//, ti, tj;

    // We continue to sweep over the mesh until the simulation has reached
    // the desired simulation Time
    // This is different from the number of iterations
    while ((iterations == -1 && t < T) || niter < iterations) {
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

        //Transfer ghost cells
        if (enableGhostCells)
        {
            #pragma omp parallel for schedule(static, 1) private(i,j,k,tx,ty) shared(blocks) 
            for (k = 0; k < ty*tx; k++)
            {
                block_t* b = blocks[k];
                int bx = b->bx, by = b->by;
                
                //Left ghost cells, thread at left edge
                if (b->tj == 0)
                {
                    for (i = 1; i < by+1; i++)
                    {
                        b->E_prev[i][0] = b->E_prev[i][2];
                    }
                }
                //Left ghost cells, block exists to the left
                else 
                {
                    block_t* left = blocks[k-1];
                    for (i = 1; i < by+1; i++)
                    {
                        b->E_prev[i][0] = left->E_prev[i][left->bx];
                    }
                }
                
                //Right ghost cells, thread at right edge
                if (b->tj == tx-1)
                {
                    for (i = 1; i < by+1; i++)
                    {
                        b->E_prev[i][bx+1] = b->E_prev[i][bx-1];
                    }
                }
                //Right ghost cells, block exists to the right
                else 
                {
                    block_t* right = blocks[k+1];
                    for (i = 1; i < by+1; i++)
                    {
                        b->E_prev[i][b->bx+1] = right->E_prev[i][1];
                    }
                }
                
                //Top ghost cells, thread at top edge
                if (b->ti == 0)
                {
                    memcpy(&b->E_prev[0][1], &b->E_prev[2][1], bx*sizeof(DOUBLE));
                }
                //Top ghost cells, block exists above
                else 
                {
                    block_t *top = blocks[k-tx];
                    memcpy(&b->E_prev[0][1], &top->E_prev[top->by][1], bx*sizeof(DOUBLE));
                }
                
                //Bottom ghost cells, thread at bottom edge
                if (b->ti == ty-1)
                {
                    memcpy(&b->E_prev[b->by+1][1], &b->E_prev[b->by-1][1], bx*sizeof(DOUBLE));
                }
                //Bottom ghost cells, block exists below
                else 
                {
                    block_t* bottom = blocks[k+tx];
                    memcpy(&b->E_prev[b->by+1][1], &bottom->E_prev[1][1], bx*sizeof(DOUBLE));
                }
            }
        }

        // Solve for the excitation, a PDE
        #pragma omp parallel for schedule(static, 1) private(i,j,k,tx,ty) shared(blocks) 
        for (k = 0; k < tx*ty; k++) {
            DOUBLE **E_local = blocks[k]->E,
                   **E_prev_local = blocks[k]->E_prev;
            int bx = blocks[k]->bx, by = blocks[k]->by;

            for (i = 1; i < by+1; i++) {
                #pragma ivdep
                for (j = 1; j < bx+1; j++) {
                    E_local[i][j] = E_prev_local[i][j] + alpha * (E_prev_local[i][j + 1] +
                              E_prev_local[i][j - 1] -
                              4 * E_prev_local[i][j] +
                              E_prev_local[i + 1][j] +
                              E_prev_local[i - 1][j]);
                }
            }
        }
        
        /*
        * Solve the ODE, advancing excitation and recovery variables
        * to the next timtestep
        */
        #pragma omp parallel for schedule(static, 1) private(i,j,k) shared(blocks) 
        for (k = 0; k < tx*ty; k++) 
        {
            DOUBLE **E_local = blocks[k]->E,
                   **R_local = blocks[k]->R;
            int bx = blocks[k]->bx, by = blocks[k]->by;
            
            for (i = 1; i < by+1; i++) 
            {
                #pragma ivdep
                for (j = 1; j < bx+1; j++) 
                {
                    E_local[i][j]  +=
                        -dt * (kk * E_local[i][j]*(E_local[i][j]-a)*(E_local[i][j] - 1) + E_local[i][j]*R_local[i][j]);

                    R_local[i][j] +=
                        dt * (
                            epsilon + M1 * R_local[i][j] /
                            (E_local[i][j] + M2)) *
                            (-R_local[i][j] - kk * E_local[i][j] * (E_local[i][j] - b - 1)
                        );
                }
            }
            
            //Swap arrays
            DOUBLE** tmp = blocks[k]->E;
            blocks[k]->E = blocks[k]->E_prev;
            blocks[k]->E_prev = tmp;
        }
        
      
        /*
        if (do_stats) {
            repNorms(E, t, dt, m, n, niter);
        }

        if (plot_freq) {
            int k = (int)(t / plot_freq);
            if ((t - k * plot_freq) < dt) 
            {
                splot(E, t, niter, m + 1, n + 1, WAIT);
            }
        }*/
    }
    

    return niter;
}
