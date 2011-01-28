#ifndef apf_h
#define apf_h

#include "types.h"
// Various constants - these definitions shouldn't change
static const DOUBLE a=0.1, b=0.1, kk=8.0, M1= 0.07, M2=0.3, epsilon=0.01, d=5e-5;
// Wait for response after plotting
#ifdef DEBUG
static const int WAIT = 1;
#else
static const int WAIT = 0;
#endif

typedef struct
{
    DOUBLE** E, ** E_prev, ** R;
    int by, bx;
    int ti, tj;
} block_t;

// External functions
void splot(DOUBLE **E, DOUBLE T, int niter, int m, int n, int WAIT);
// Uses gettimeofday() to collect timings on the host side
double getTime();

void printMat(DOUBLE **U, int m, int n);
void printMatLocal(DOUBLE **U, int m, int n);
void repNorms(DOUBLE **E, DOUBLE t, DOUBLE dt, int m,int n, int niter);

DOUBLE** alloc2D(int m, int n);



#endif
