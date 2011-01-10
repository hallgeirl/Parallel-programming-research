/*
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 *
 * Modified and  restructured by Scott B. Baden, UCSD
 *
 */
#define _XOPEN_SOURCE 600
#include <pthread.h>
#include <semaphore.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
//#include <omp.h>
#include "time.h"
#include "apf.h"
#include "types.h"

typedef struct thread_args_s
{
    int thread_id;
    DOUBLE *** E; DOUBLE *** E_prev; DOUBLE ** R;
    DOUBLE alpha;
    int offset_x; int offset_y;
    int bx; int by; 
    DOUBLE dt;
} thread_args_t;

pthread_barrier_t barr;
int threadcount;
int state = 0;
//int done = 0;
int done_count;

pthread_cond_t ready = PTHREAD_COND_INITIALIZER;
pthread_cond_t done = PTHREAD_COND_INITIALIZER;
pthread_cond_t swapped = PTHREAD_COND_INITIALIZER;

pthread_mutex_t ready_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t done_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t swapped_mutex = PTHREAD_MUTEX_INITIALIZER;

//Solves the PDE part for one block of the array
//void solve_pde(DOUBLE ** E, DOUBLE ** E_prev, int offset_x, int offset_y, int bx, int by, int m, int n, DOUBLE alpha)
void solve_block(void* _args)
{
    thread_args_t* args = (thread_args_t*) _args;

    int thread_id = args->thread_id;
    DOUBLE ** R = args->R;
    DOUBLE alpha = args->alpha;
    int bx = args->bx,  by = args->by;
    int offset_x = args->offset_x,  offset_y = args->offset_y;
    DOUBLE dt = args->dt;
    int i,j;
    //printf("Thread %d says hello.\n", thread_id);

    while (true)
    {
        //Wait for border padding
	pthread_barrier_wait(&barr);
        DOUBLE ** E = *args->E;
        DOUBLE ** E_prev = *args->E_prev;
       
        for (i = offset_y; i < offset_y + by; i++)
        {
            #pragma ivdep
            for (j = offset_x; j < offset_x + bx; j++)
            {
                E[i][j] = E_prev[i][j] + alpha * (E_prev[i][j + 1]+
                                          E_prev[i][j - 1]-
                                          4 * E_prev[i][j]+
                                          E_prev[i + 1][j]+
                                          E_prev[i - 1][j]);
            }
        }
        
        for (i = offset_y; i < offset_y + by; i++)
        {
            #pragma ivdep
            for (j = offset_x; j < offset_x + bx; j++)
            {
                E[i][j] += -dt * (kk * E[i][j]*(E[i][j]-a)*(E[i][j] - 1) + E[i][j]*R[i][j]);
                R[i][j] +=
                    dt * (
                        epsilon + M1 * R[i][j] /
                        (E[i][j] + M2)) *
                        (-R[i][j] - kk * E[i][j] * (E[i][j] - b - 1)
                    );
            }
        }
        
        //Barrier to make sure no thread waits for ready signal
        #ifdef DEBUG
        printf("Thread %d waits at barrier.\n", thread_id);
        #endif
        pthread_barrier_wait(&barr);
        #ifdef DEBUG
        printf("Thread %d past barrier.\n", thread_id);
        #endif
    }
    printf("Thread %d exits.\n", thread_id);
    pthread_exit(NULL);
}

int solve(DOUBLE ***_E, DOUBLE ***_E_prev, DOUBLE **R, int m, int n, DOUBLE T, DOUBLE alpha, DOUBLE dt, int do_stats, int plot_freq, int tx, int ty) 
{
    // Simulated time is different from the integer timestep number
    DOUBLE t = 0.0;
    // Integer timestep number
    int niter=0;
    
    //These values are used to determine the block size for the last threads, so cache the results for performance.
    int m_mod_ty = (m+1) % ty;
    int n_mod_tx = (n+1) % tx;
    int ti, tj;

    DOUBLE **E = *_E, **E_prev = *_E_prev;
    threadcount = tx*ty;
    pthread_barrier_init(&barr, NULL, tx*ty+1);
    
    //Allocate threads
    pthread_t *threads = (pthread_t*)malloc(sizeof(pthread_t)*tx*ty);
 
    //The arguments are fixed for each iteration, so initialize them only once
    thread_args_t *thread_args = (thread_args_t*)malloc(sizeof(thread_args_t)*tx*ty);
    for (ti = 0; ti < ty; ti++)
    {
        int offset_y = ti*(m+1)/ty + 1, by = (m+1)/ty;
        
        //Let the last thread take care of any leftovers
        if (ti == ty-1) by += m_mod_ty;
          
        for (tj = 0; tj < tx; tj++)
        {
            //Determine boundaries and block size
            int offset_x = tj*(n+1)/tx + 1, bx = (n+1)/tx;
            //Let the last thread take care of any leftovers
            if (tj == tx-1) bx += n_mod_tx;

            thread_args[ti*tx+tj].thread_id = ti*tx+tj;
            thread_args[ti*tx+tj].E = &E; 
            thread_args[ti*tx+tj].E_prev = &E_prev;
            thread_args[ti*tx+tj].R = R;
            thread_args[ti*tx+tj].alpha = alpha;
            thread_args[ti*tx+tj].offset_x = offset_x;
            thread_args[ti*tx+tj].offset_y = offset_y;
            thread_args[ti*tx+tj].bx = bx;
            thread_args[ti*tx+tj].by = by;
            thread_args[ti*tx+tj].dt = dt;
            
            
            //Start the thread
            pthread_create(&threads[ti*tx+tj], NULL, &solve_block, &thread_args[ti*tx+tj]);
        }
    }

    // We continue to sweep over the mesh until the simulation has reached
    // the desired simulation Time
    // This is different from the number of iterations
    while (t < T) 
    {
        #ifdef DEBUG
        printf("Main thread starts new loop.\n");
        #endif
        
        //printf("Iteration %d\n", niter);
        #ifdef DEBUG
        printMat(E_prev, m, n);
        repNorms(E_prev, t, dt, m, n, niter);
        if (plot_freq) {
            splot(E_prev, t, niter, m + 1, n + 1, WAIT);
        }
        printf("\n");
        #endif

        t += dt;
        niter++;

        /*
        * Copy data from boundary of the computational box to the
        * padding region, set up for differencing computational box's boundary
        *
        */
        int i, j;
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
	pthread_barrier_wait(&barr);


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

	pthread_barrier_wait(&barr);

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
