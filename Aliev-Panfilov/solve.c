/*
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 *
 * Modified and  restructured by Scott B. Baden, UCSD
 * 
 * Modified further by Hallgeir Lien
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
int done_count;

pthread_cond_t ready = PTHREAD_COND_INITIALIZER;
pthread_cond_t done = PTHREAD_COND_INITIALIZER;

pthread_mutex_t ready_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t done_mutex = PTHREAD_MUTEX_INITIALIZER;

//Solves the PDE and ODE part for one block of the array.
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
    
    #ifdef DEBUG
    printf("Thread %d params: bx=%d by=%d offx=%d offy=%d\n", thread_id, bx, by, offset_x, offset_y);
    fflush(stdout);
    #endif
    
    while (true)
    {
        //Wait for border padding
        pthread_mutex_lock(&ready_mutex);
        #ifdef DEBUG
        printf("Thread %d waits for ready.\n", thread_id);
        fflush(stdout);
        #endif
        while (state == 0)
            pthread_cond_wait(&ready, &ready_mutex);
            
        DOUBLE ** E = *args->E;
        DOUBLE ** E_prev = *args->E_prev;

        #ifdef DEBUG
        printf("Thread %d got ready.\n", thread_id);
        fflush(stdout);
        #endif
        pthread_mutex_unlock(&ready_mutex);
        
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
        fflush(stdout);
        #endif
        pthread_barrier_wait(&barr);
        #ifdef DEBUG
        printf("Thread %d past barrier.\n", thread_id);
        fflush(stdout);
        #endif
        
        pthread_mutex_lock(&done_mutex);
        
        //Let the first thread entering this lock reset the state variable
        if (done_count == 0)
        {
            #ifdef DEBUG
            printf("Thread %d sets state to 0.\n", thread_id);
            #endif
            state = 0;
        }
        
        //Increment count    
        done_count++;
        
        //Are we done? If so, signal the main thread to continue.
        if (done_count == threadcount)
        {
            #ifdef DEBUG
            printf("Thread %d broadcasts done.\n", thread_id);
            #endif
            pthread_cond_broadcast(&done);
        }
        pthread_mutex_unlock(&done_mutex);
    }
    printf("Thread %d exits.\n", thread_id);
    pthread_exit(NULL);
}

int solve(DOUBLE ***_E, DOUBLE ***_E_prev, DOUBLE **R, int m, int n, DOUBLE T, int iterations, DOUBLE alpha, DOUBLE dt, int do_stats, int plot_freq, int tx, int ty) 
{
    // Simulated time is different from the integer timestep number
    DOUBLE t = 0.0;
    // Integer timestep number
    int niter=0;
    
    //These values are used to determine the block size for the last threads, so cache the results for performance.
    int m_mod_ty = (m+1) % ty;
    int n_mod_tx = (n+1) % tx;
    int blocks_x = (n+1) / tx;
    int blocks_y = (m+1) / ty;
    int ti, tj;

    DOUBLE **E = *_E, **E_prev = *_E_prev;
    threadcount = tx*ty;
    pthread_barrier_init(&barr, NULL, tx*ty);
    
    //Allocate threads
    pthread_t *threads = (pthread_t*)malloc(sizeof(pthread_t)*tx*ty);
 
    //The arguments are fixed for each iteration, so initialize them only once
    thread_args_t *thread_args = (thread_args_t*)malloc(sizeof(thread_args_t)*tx*ty);
    for (ti = 0; ti < ty; ti++)
    {
        int offset_y = ti*(blocks_y) + 1, 
            by = (m+1)/ty;
        
        //Let the last thread take care of any leftovers
        if (ti == ty-1) by += m_mod_ty;
          
        for (tj = 0; tj < tx; tj++)
        {
            //Determine boundaries and block size
            int offset_x = tj*(blocks_x) + 1, bx = (n+1)/tx;
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

/*            //Set CPU affinity*/
/*            cpu_set_t cpuset;*/
/*            CPU_ZERO(&cpuset);*/
/*            CPU_SET((ti*tx+tj)%2, &cpuset);*/
/*            pthread_setaffinity_np(threads[ti*tx+tj], sizeof(cpu_set_t), &cpuset);*/
            
            //Start the thread
            pthread_create(&threads[ti*tx+tj], NULL, &solve_block, &thread_args[ti*tx+tj]);
        }
    }

    // We continue to sweep over the mesh until the simulation has reached
    // the desired simulation Time
    // This is different from the number of iterations
    while ((iterations < 0 && t < T) || niter < iterations) 
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

        //Wake up threads to do one iteration worth of work
        pthread_mutex_lock(&ready_mutex);
        #ifdef DEBUG
        printf("Main thread inited ghost cells, sets state to 1 and broadcasts ready.\n");
        #endif
        state = 1;
        pthread_cond_broadcast(&ready);
        pthread_mutex_unlock(&ready_mutex);

        //Wait for threads to finish
        pthread_mutex_lock(&done_mutex);
        #ifdef DEBUG
        printf("Main waits for done signal.\n");
        #endif
        while (done_count < tx*ty)
            pthread_cond_wait(&done, &done_mutex);
        #ifdef DEBUG
        printf("Main received done signal, sets done_count to 0.\n");
        #endif
        done_count = 0;
        pthread_mutex_unlock(&done_mutex);

        //Main loop
        /*for (ti = 0; ti < ty; ti++)
        {
            for (tj = 0; tj < tx; tj++)
            {
                //Create thread and execute solver for sub-problem
                pthread_create(&threads[ti*tx+tj], NULL, &solve_block, &thread_args[ti*tx+tj]);
            }
        }*/
        
        //Join the threads
        /*for (ti = 0; ti < ty; ti++)
        {
            for (tj = 0; tj < tx; tj++)
            {
                pthread_join(threads[ti*tx+tj], NULL);
            }
        }*/

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
        
        //Tell the threads that the arrays are swapped
        /*pthread_mutex_lock(&swapped_mutex);
        swap = 1;
        pthread_cond_broadcast(&swapped);
        pthread_mutex_unlock(&swapped_mutex);*/
    }

    // Store them into the pointers passed in
    *_E = E;
    *_E_prev = E_prev;

    return niter;
}
