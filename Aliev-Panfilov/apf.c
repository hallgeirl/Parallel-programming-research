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
#include <string.h>
#include <math.h>
#include <omp.h>
#include "time.h"
#include "apf.h"
#include "types.h"

// Utilities
// Allocate a 2D array
DOUBLE **alloc2D(int m, int n, int *pitch) {
    DOUBLE **E;
    int leftover = n % (4096/sizeof(DOUBLE));
    *pitch = n + (((4096/sizeof(DOUBLE))-leftover) % (4096/sizeof(DOUBLE)));
    E = (DOUBLE **)malloc(sizeof(DOUBLE*) * m + sizeof(DOUBLE) * (*pitch) * m + 4096);
    assert(E);
   
    //Align to page boundary
    unsigned long lowerBits = ((unsigned long)(E+m) & 4095);
    int shifted = (4096-lowerBits)+4096; //Shift the whole array by this much to align each row to the 4096 border
    
    int i;
    for (i = 0; i < m; i++)
    {
        E[i] = (DOUBLE *)(((DOUBLE**)(((char*)E)+shifted)) + m) + i * (*pitch);
        #ifdef DEBUG
        printf("E[%i] mod 4096: %d\n", i, ((int)E[i] & 4095));
        #endif
    }
        
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

void init (DOUBLE **E_prev, DOUBLE **R, int m, int n) {
  int i, j;
  // Initialization
  for (i = 1; i <= m + 1; i++) {
    for (j = 1; j <= n + 1; j++) {
      E_prev[i][j] = R[i][j] = 0;
    }
  }

  for (i = 1; i <= m + 1; i++) {
    for (j = n / 2 + 2; j <= n + 1 ; j++) {
      E_prev[i][j] = 1.0;
    }
  }

  for (i = m / 2 +2; i <= m + 1; i++) {
    for (j = 1; j <= n + 1; j++) {
        R[i][j] = 1.0;
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

void cmdLine(int argc, char *argv[], double* T, int *m, int* n, int* tx, int* ty, int *iters, int* do_stats, int* plot_freq);
int solve(DOUBLE ***_E, DOUBLE ***_E_prev, DOUBLE **R, int m, int n, DOUBLE T, int iters, DOUBLE alpha, DOUBLE dt, int tx, int ty, int do_stats, int plot_freq);

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
    int iters = -1;
    int do_stats = 0;
    int plot_freq = 0;
    int tx = 1, ty = 1;
    int pitch;

    cmdLine(argc, argv, &T, &m, &n, &tx, &ty, &iters, &do_stats, &plot_freq);
    if (m == -1) m = n;
    if (n == -1) n = m;

    printTOD("Run begins");
    printf("Max number of threads: %d\n", omp_get_max_threads());
    // Allocate contiguous memory for solution arrays
    // The computational box is defined on [1:m+1,1:n+1]
    // We pad the arrays in order to facilitate differencing on the 
    // boundaries of the computation box
    E = alloc2D(m + 3, n + 3, &pitch);
    E_prev = alloc2D(m + 3, n + 3, &pitch);
    R = alloc2D(m + 3, n + 3, &pitch);

    init(E_prev, R, m, n);
    
    DOUBLE dx = 1.0 / n;

    // For time integration, these values shouldn't change 
    double rp= kk * (b + 1) * (b + 1) / 4;
    double dte= (dx * dx) / (d * 4 + ((dx * dx)) * (rp + kk));
    double dtr= 1 / (epsilon + ((M1 / M2) * rp));
    double dt = (dte < dtr) ? 0.95 * dte : 0.95 * dtr;
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
    printf("Using SINGLE precision arithmetic (element size: %u)\n", sizeof(DOUBLE));
    #else
    printf("Using DOUBLE precision arithmetic (element size: %u)\n", sizeof(DOUBLE));
    #endif

    printf("dt=%f, T = %lf\n", dt, T);
    printf("m x n = %d x %d\n", m, n);
    printf("Thread geometry: %d x %d\n", tx, ty);
    printf("\n");

    // Start the timer
    double t0 = -getTime();
    int niter = solve(&E, &E_prev, R, m, n, T, iters, alpha, dt, tx, ty, do_stats, plot_freq);
    t0 += getTime();

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
