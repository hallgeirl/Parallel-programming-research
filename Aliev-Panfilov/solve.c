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

pthread_t *threads;

typedef struct thread_args_s
{
    int thread_id, thread_x, thread_y;
    DOUBLE ** E; DOUBLE ** E_prev; DOUBLE ** R;            //Local arrays for this thread
    DOUBLE * edgeTop, * edgeBottom, * edgeLeft, * edgeRight; //Border elements
    DOUBLE alpha;
    //int offset_x; int offset_y; //Might not be needed anymore...
    int bx, by, tx, ty; 
    DOUBLE dt;
} thread_args_t;
thread_args_t * thread_args;

typedef struct thread_init_args_s
{
    thread_args_t * args;
    DOUBLE ** E_prev;
    DOUBLE ** R;
    int n, m, tx, ty;
} thread_init_args_t;

pthread_barrier_t barr;
int threadcount;
int state = 0;
int done_count;

pthread_cond_t ready = PTHREAD_COND_INITIALIZER;
pthread_cond_t done = PTHREAD_COND_INITIALIZER;

pthread_mutex_t ready_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t done_mutex = PTHREAD_MUTEX_INITIALIZER;

//Solves the PDE and ODE part for one block of the array.
void * solve_block(void* _args)
{
    thread_args_t* args = (thread_args_t*) _args;

    int thread_id = args->thread_id;
    DOUBLE ** R = args->R;
    DOUBLE alpha = args->alpha;
    int bx = args->bx,  by = args->by;
    //int offset_x = args->offset_x,  offset_y = args->offset_y;
    DOUBLE dt = args->dt;
    int i,j, tj = args->thread_x, ti = args->thread_y, tx = args->tx, ty = args->ty;
    
    #ifdef DEBUG
    printf("Thread %d params: bx=%d by=%d ti=%d tj=%d\n", thread_id, bx, by, ti, tj);
    fflush(stdout);
    #endif
    
    //References to the thread arguments on all four sides.
    thread_args_t *t_left   = (tj == 0    ? 0 : &thread_args[ti*tx+tj-1]);
    thread_args_t *t_right  = (tj == tx-1 ? 0 : &thread_args[ti*tx+tj+1]);
    thread_args_t *t_top    = (ti == 0    ? 0 : &thread_args[(ti-1)*tx+tj]);
    thread_args_t *t_bottom = (ti == ty-1 ? 0 : &thread_args[(ti+1)*tx+tj]);
    
    //Border elements for neighboring blocks.
    DOUBLE* left = 0, *right = 0, *top = 0, *bottom = 0;


    //If we're not on the edge of the array, point the border element pointers to the neighbor's border element array.
    if (tj == 0) left = (DOUBLE*)malloc(sizeof(DOUBLE)*by);
    else left = t_left->edgeRight;
    if (tj == tx-1) right = (DOUBLE*)malloc(sizeof(DOUBLE)*by);
    else right = t_right->edgeLeft;
    if (ti != 0) top = t_top->edgeBottom;
    if (ti != ty-1) bottom = t_bottom->edgeTop;
   
    //if (ti == 0) top = (DOUBLE*)malloc(sizeof(DOUBLE)*bx);
    //if (ti == by-1) bottom = (DOUBLE*)malloc(sizeof(DOUBLE)*bx);

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
            
        #ifdef DEBUG
        printf("Thread %d got ready.\n", thread_id);
        fflush(stdout);
        #endif
        pthread_mutex_unlock(&ready_mutex);
        
        DOUBLE ** E = args->E;
        DOUBLE ** E_prev = args->E_prev;
        
        if (ti == 0) top = E_prev[2];
        if (ti == ty-1) bottom = E_prev[by-3];
        
        if (tj == 0)
        {
            for (i = 0; i < by; i++)
                left[i] = E_prev[i][2];
        }
        if (tj == tx-1)
        {
            for (i = 0; i < by; i++)
                right[i] = E_prev[i][bx-3];
        }

        //Solve the PDE. First the inner area.
        for (i = 1; i < by-1; i++)
        {
            #pragma ivdep
            for (j = 1; j < bx-1; j++)
            {
                E[i][j] = E_prev[i][j] + alpha * (E_prev[i][j + 1]+
                                          E_prev[i][j - 1]-
                                          4 * E_prev[i][j]+
                                          E_prev[i + 1][j]+
                                          E_prev[i - 1][j]);
            }
        }

        
        //and the borders.
        //Left and right borders
        #pragma ivdep
        for (i = 0; i < by; i++)
        {
            E[i][0] = E_prev[i][0] + alpha * (E_prev[i][1] +
                                      left[i] -
                                      4 * E_prev[i][0] +
                                      (i > 0 ? E_prev[i - 1][0] : top[0]) +
                                      (i < by-1 ? E_prev[i + 1][0] : bottom[0]));
                                      
            E[i][bx-1] = E_prev[i][bx-1] + alpha * (right[i] +
                                      E_prev[i][bx-2] -
                                      4 * E_prev[i][bx-1]+
                                      (i > 0 ? E_prev[i - 1][bx-1] : top[bx-1]) +
                                      (i < by-1 ? E_prev[i + 1][bx-1] : bottom[bx-1]));
        }
        
        #pragma ivdep
        for (j = 0; j < bx; j++)
        {
            E[0][j] = E_prev[0][j] + alpha * ((j < bx-1 ? E_prev[0][j+1] : right[0]) +
                                      (j > 0 ? E_prev[0][j-1] : left[0]) -
                                      4 * E_prev[0][j] +
                                      E_prev[1][j] +
                                      top[j]);
                                      
            E[by-1][j] = E_prev[by-1][bx-1] + alpha * ((j < bx-1 ? E_prev[by-1][j+1] : right[by-1]) +
                                      (j > 0 ? E_prev[by-1][j-1] : left[by-1]) -
                                      4 * E_prev[by-1][j]+
                                      bottom[j] +
                                      E_prev[bx-1][j-1]);
        }
        
        //Solve the ODE for one time step
        for (i = 1; i < by-1; i++)
        {
            #pragma ivdep
            for (j = 0; j < bx; j++)
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
        
        //Update the edge arrays
        for (j = 0; j < bx; j++)
        {
            args->edgeTop[j] = E_prev[0][j];
            args->edgeBottom[j] = E_prev[by-1][j];
        }
        
        for (i = 0; i < by; i++)
        {
            args->edgeLeft[i] = E_prev[i][0];
            args->edgeRight[i] = E_prev[i][bx-1];
        }
        
        #ifdef DEBUG
        printf("Thread %d past barrier.\n", thread_id);
        fflush(stdout);
        #endif
        
        //Swap arrays
        DOUBLE **tmp = args->E;
        args->E = args->E_prev;
        args->E_prev = tmp;

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
    if (tj == 0) free(left);
    if (tj == tx-1) free(right);
    
    printf("Thread %d exits.\n", thread_id);
    pthread_exit(NULL);
}

void * init_thread(void * _args)
{
    thread_init_args_t * args = (thread_init_args_t*)_args;
    int i, j;
    int ti = args->args->thread_y, tj = args->args->thread_x;
    int n = args->n, m = args->m;
    int tx = args->tx, ty = args->ty;
    
    //These values are used to determine the block size for the last threads, so cache the results for performance.
    int blocks_x = (n+1) / tx;
    int blocks_y = (m+1) / ty;

    //Offsets into the global array and block size
    int offset_x = tj*(blocks_x) + 1, bx = (n+1)/tx;
    int offset_y = ti*(blocks_y) + 1, by = (m+1)/ty;

    args->args->E = alloc2D(by-1, bx-1);
    args->args->E_prev = alloc2D(by-1, bx-1);
    args->args->R = alloc2D(by-1, bx-1);

    args->args->edgeTop = (DOUBLE*)malloc(sizeof(DOUBLE)*bx);
    args->args->edgeBottom = (DOUBLE*)malloc(sizeof(DOUBLE)*bx);
    args->args->edgeLeft = (DOUBLE*)malloc(sizeof(DOUBLE)*by);
    args->args->edgeRight = (DOUBLE*)malloc(sizeof(DOUBLE)*by);

    //Let the last thread take care of any leftovers
    if (tj == tx-1) bx += (n+1) % tx;
    if (ti == ty-1) by += (m+1) % ty;

    args->args->bx = bx;
    args->args->by = by;
    
    //Copy data from the global array to the local one
    for (i = 0; i < bx; i++)
    {
        #pragma ivdep
        for (j = 0; j < by; j++)
        {
            args->args->E_prev[i][j] = args->E_prev[offset_y+i][offset_x+j];
            args->args->R[i][j] = args->R[offset_y+i][offset_x+j];
        }
    }
    
    return 0;
}

int solve(DOUBLE ***_E, DOUBLE ***_E_prev, DOUBLE **R, int m, int n, DOUBLE T, DOUBLE alpha, DOUBLE dt, int do_stats, int plot_freq, int tx, int ty) 
{
    // Simulated time is different from the integer timestep number
    DOUBLE t = 0.0;
    // Integer timestep number
    int niter=0;
    
    int ti, tj;

    DOUBLE **E = *_E, **E_prev = *_E_prev;
    threadcount = tx*ty;
    pthread_barrier_init(&barr, NULL, tx*ty);
    
    //Allocate threads
    threads = (pthread_t*)malloc(sizeof(pthread_t)*tx*ty);
    pthread_t* init_threads = (pthread_t*)malloc(sizeof(pthread_t)*tx*ty);
 
    //The arguments are fixed for each iteration, so initialize them only once
    thread_args = (thread_args_t*)malloc(sizeof(thread_args_t)*tx*ty);
    for (ti = 0; ti < ty; ti++)
    {
        for (tj = 0; tj < tx; tj++)
        {
            thread_args_t* ta = &thread_args[ti*tx+tj];
            ta->thread_id = ti*tx+tj;
            ta->thread_x = tj; ta->thread_y = ti;
            ta->alpha = alpha; ta->dt = dt;
            ta->tx = tx; ta->ty = ty;

            thread_init_args_t init_args;
            init_args.R =  R;  init_args.E_prev = E_prev;
            init_args.m =  m;  init_args.n      = n;
            init_args.tx = tx; init_args.ty     = ty;
            init_args.args = ta;

            //Start the initialization thread
            pthread_create(&threads[ti*tx+tj], NULL, &init_thread, &init_args);
        }
    }
    
    for (ti = 0; ti < ty; ti++)
    {
        for (tj = 0; tj < tx; tj++)
        {
            //Wait for initialization to finish
            pthread_join(threads[ti*tx+tj], NULL);
            //then create the solver thread
            pthread_create(&threads[ti*tx+tj], NULL, &solve_block, &thread_args[ti*tx+tj]);
        }
    }
    
    free(init_threads);

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
        //DOUBLE **tmp = E;
        //E = E_prev;
        //E_prev = tmp;
    }

    // Store them into the pointers passed in
    *_E = E;
    *_E_prev = E_prev;

    return niter;
}
