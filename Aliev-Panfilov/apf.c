/*
 * Driver for a cardiac elecrophysioly simulatin that uses the
 * Aliev-Panfilov model
 * We use an explicit method
 *
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 *
 * Modified and  restructured by Scott B. Baden, UCSD
 */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "time.h"
#include "apf.h"
#include "types.h"

// Utilities
// Allocate a 2D array
DOUBLE **alloc2D(int m, int n) {
    DOUBLE **E;
    int nx = n, ny = m;
    E = (DOUBLE **)malloc(sizeof(DOUBLE*) * ny + sizeof(DOUBLE) * nx * ny);
    assert(E);
    int j;
    for (j = 0; j < ny; j++)
        E[j] = (DOUBLE *)(E + ny) + j * nx;
    return E;
}
    
// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem
DOUBLE stats(DOUBLE **E, int m, int n, DOUBLE *_mx) {
  DOUBLE mx = -1;
  DOUBLE l2norm = 0;
  int i, j;
  for (j = 1; j <= m + 1; j++)
    for (i = 1; i <= n + 1; i++) {
      l2norm += E[j][i] * E[j][i];
      DOUBLE fe = fabs(E[j][i]);
      if (fe > mx)
        mx = fe;
    }

    l2norm /= (DOUBLE) ((m + 1) * (n + 1));
    l2norm = sqrt(l2norm);

    *_mx = mx;
    return l2norm;
}

void init (DOUBLE **E, DOUBLE **E_prev, DOUBLE **R, int m, int n) {
  int i, j;
  // Initialization
  for (j = 1; j <= m + 1; j++) {
    for (i = 1; i <= n + 1; i++) {
      E_prev[j][i] = R[j][i] = 0;
    }
  }

  for (j = 1; j <= m + 1; j++) {
    for (i = n / 2 + 2; i <= n + 1 ; i++) {
      E_prev[j][i] = 1.0;
    }
  }

  for (j = m / 2 +2; j <= m + 1; j++) {
    for (i = 1; i <= n + 1; i++) {
        R[j][i] = 1.0;
    }
  }
}

// Report statistics periodically
// Reduce the value of FREQ to increase the frequency,
// increase the value to raise the frequency
void repNorms(DOUBLE **E,DOUBLE t, DOUBLE dt, int m,int n, int niter) {
    const int FREQ = 100;
    int k = (int)(t / FREQ);
    
    if ((t - k * 100) < dt) 
    {
        DOUBLE mx;
        DOUBLE l2norm = stats(E, m, n, &mx);
        printf("iteration %d, t=%g\n", niter, t);
        printf("Max norm: %13.6e, L2norm: %13.6e\n", mx, l2norm);
    }
}

void printTOD(const char* mesg)
{
    time_t tim = time(NULL);
    char* s = ctime(&tim);
    if (strlen(mesg) ==  0) 
        printf("Time of day: %s\n", s);
    else {
        printf("[%s] ", mesg);
        printf("Time of day: %s\n", s);
    }
    printf("\n");
}


// External functions
void cmdLine(int argc, char *argv[], double* T, int* m, int* n, int* tx, int* ty, int *iterations, bool *enableGhostCells, int* do_stats, int* plot_freq);

int solve(block_t ** blocks, int m, int n, DOUBLE T, DOUBLE alpha, DOUBLE dt, int do_stats, int plot_freq, int tx, int ty, bool enableGhostCells, int iterations);

// Main program
int main(int argc, char** argv) {
    /*
    *  Solution arrays
    *   E is the "Excitation" variable, a voltage
    *   R is the "Recovery" variable
    *   E_prev is the Excitation variable for the previous timestep,
    *      and is used in time integration
    */
    DOUBLE **E, **R, **E_prev;
    DOUBLE mx, l2norm;

    // Command line arguments
    // Default value of the domain sizes
    double T = 1500.0;
    int m = 100, n = -1;
    int do_stats = 0;
    int plot_freq = 0;
    int tx = 1, ty = 1;
    int bx = m / 4, by = n / 4;
    bool enableGhostCells = true;
    int iterations = -1;

    cmdLine(argc, argv, &T, &m, &n, &tx, &ty, &iterations, &enableGhostCells, &do_stats, &plot_freq);
    if (n == -1) n = m;
    if (m == -1) m = n;

    omp_set_num_threads(tx*ty);

    printTOD("Run begins");
    printf("Max number of threads: %d\n", omp_get_max_threads());
    // Allocate contiguous memory for solution arrays
    // The computational box is defined on [1:m+1,1:n+1]
    // We pad the arrays in order to facilitate differencing on the 
    // boundaries of the computation box
    E = alloc2D(m + 3, n + 3);
    E_prev = alloc2D(m + 3, n + 3);
    R = alloc2D(m + 3, n + 3);

    init(E, E_prev, R, m, n);
    
    DOUBLE dx = 1.0 / n;

    // For time integration, these values shouldn't change 
    double rp= kk * (b + 1) * (b + 1) / 4;
    double dte= (dx * dx) / (d * 4 + ((dx * dx)) * (rp + kk));
    double dtr= 1 / (epsilon + ((M1 / M2) * rp));
    double ddt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
    DOUBLE dt = (DOUBLE) ddt;
    DOUBLE alpha = d * dt / (dx * dx);

    #ifdef DEBUG
    printMat(E_prev, m, n);
    repNorms(E_prev, -1, dt, m, n, -1);
    if (plot_freq) 
    {
        splot(E_prev, -1, -1, m + 1, n + 1, WAIT);
    }
    #endif


    #ifdef FLOAT
    printf("Using SINGLE precision arithmetic (element size: %d)\n", sizeof(DOUBLE));
    #else
    printf("Using DOUBLE precision arithmetic (element size: %d)\n", sizeof(DOUBLE));
    #endif

    printf("dt=%f, T = %lf\n", dt, T);
    printf("m x n = %d x %d\n", m, n);
    printf("Thread geometry: %d x %d\n", tx, ty);
    printf("Cache blocking factors: %d x %d\n", bx, by);
    printf("\n");

    //Initialize blocks
    block_t ** blocks = (block_t**)malloc(sizeof(block_t*)*tx*ty);

    int i,j,k;
    
    #pragma omp parallel for schedule(static, 1) private(k,i,j) // We want each thread to do one block only.
    for (k = 0; k < ty*tx; k++)
    {
        int ti = k/tx, tj = k%tx;

        //Allocate a block
        block_t* block = malloc(sizeof(block_t));
        block->ti = ti; block->tj = tj;
        block->by = ((m+1)/ty) + (ti == ty-1 ? (m+1) % ty : 0);
        block->bx = ((n+1)/tx) + (tj == tx-1 ? (n+1) % tx : 0);
        block->E = alloc2D(block->by+2, block->bx+2);
        block->E_prev = alloc2D(block->by+2, block->bx+2);
        block->R = alloc2D(block->by+2, block->bx+2);
        int bx = (n+1)/tx;
        int by = (m+1)/ty;

        //Copy data from the large array
        for (i = 1; i < block->by+1; i++)
        {
            for (j = 1; j < block->bx+1; j++)
            {
                block->E_prev[i][j] = E_prev[i+ti*by][j+tj*bx];
                block->R[i][j] = R[i+ti*by][j+tj*bx];
            }
        }
        
        blocks[k] = block;
        
        printf("Thread %d ti=%d tj=%d bx=%d by=%d\n", k, block->ti, block->tj, block->bx, block->by);
    }

    // Start the timer
    double t0 = -getTime();
    int niter = solve(blocks, m, n, T, alpha, dt, do_stats, plot_freq, tx, ty, enableGhostCells, iterations);
    t0 += getTime();


    //Copy the data back to the original arrays
    #pragma omp parallel for schedule(static, 1) private(k,i,j) // We want each thread to do one block only.
    for (k = 0; k < ty*tx; k++)
    {
        int ti = k/tx, tj = k%tx;
        block_t *block = blocks[k];
        
        for (i = 1; i < block->by+1; i++)
        {
            for (j = 1; j < block->bx+1; j++)
            {
                E_prev[i+ti*((m+1)/ty)][j+tj*((n+1)/tx)] = block->E_prev[i][j];
            }
        }
    }

    printTOD("Run completes");
    printf("End at time: %g, iteration %d\n", T, niter);
    l2norm = stats(E_prev, m, n, &mx);
    printf("Max: %13.5e, L2norm: %13.5e\n", mx, l2norm);
    printf("Running Time: %f sec.\n", t0);
    if (plot_freq) {
    printf("\n\nEnter any input to close the program and the plot...");
    int resp;
    scanf("%d", &resp);
    }

    free (E);
    free (E_prev);
    free (R);

    return 0;
}
