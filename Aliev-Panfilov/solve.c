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

#define STATE_NOTREADY 0
#define STATE_WORKING 1
#define STATE_DONE 2
#define STATE_COPIED = 3

typedef struct thread_args_s
{
    int thread_id, thread_x, thread_y;
    DOUBLE **E_prev_global, **R_global;  //Global arrays, for initialization
    DOUBLE ***E_prev;                    //For ghost cells
    int iterations;
    //DOUBLE * edgeTop, * edgeBottom, * edgeLeft, * edgeRight; //Border elements
    DOUBLE alpha;
    int tx, ty, m, n, bx, by; 
    DOUBLE dt;
    DOUBLE T;
} thread_args_t;

//int m, n, tx, ty, iterations;
//DOUBLE **E_prev_global, **R_global;


pthread_t *threads;
thread_args_t ** thread_args; // Need this to be global in order for threads to access neighboring blocks.

pthread_barrier_t barr;
int threadcount;
int state = STATE_NOTREADY;
int done_count;

pthread_cond_t ready = PTHREAD_COND_INITIALIZER;
pthread_cond_t done = PTHREAD_COND_INITIALIZER;

pthread_mutex_t ready_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t done_mutex = PTHREAD_MUTEX_INITIALIZER;

#define MT_PRINT(...) printf("Thread %d\t", thread_id); printf(__VA_ARGS__); printf("\n"); fflush(stdout); 

//Solves the PDE and ODE part for one block of the array.
void * solve_block(void* _args)
{
    thread_args_t* args = (thread_args_t*) _args;

    int thread_id = args->thread_id;
    DOUBLE ** E, **E_prev, **R;
    DOUBLE alpha = args->alpha;
    DOUBLE dt = args->dt;
    DOUBLE T = args->T;
    int iterations = args->iterations;
    int i,j, tj = args->thread_x, ti = args->thread_y, tx = args->tx, ty = args->ty;
    int n = args->n, m = args->m;
    double t = 0;
    int niter = 0;

    //Initialize the blocks
    //Block sizes
    int bx = (n+1) / tx;
    int by = (m+1) / ty;
    args->bx = bx; args->by = by;

    //Offsets into the global array (for initialization and copying back)
    int offset_x = tj*bx + 1;
    int offset_y = ti*by + 1;

    //Let the last thread take care of any leftovers
    if (tj == tx-1) bx += (n+1) % tx;
    if (ti == ty-1) by += (m+1) % ty;

    //Allocate local arrays
    E = alloc2D(by+2, bx+2);
    E_prev = alloc2D(by+2, bx+2); 
    R = alloc2D(by+2, bx+2);
    args->E_prev = &E_prev;

    #ifdef DEBUG
    MT_PRINT("Initializing...");
    #endif
    //Copy data from the global array to the local one
    for (i = 0; i < by; i++)
    {
        #pragma ivdep
        for (j = 0; j < bx; j++)
        {
            E_prev[i+1][j+1] = args->E_prev_global[offset_y+i][offset_x+j];
            R[i+1][j+1] = args->R_global[offset_y+i][offset_x+j];
        }
    }

    #ifdef DEBUG
    MT_PRINT("Initialization done.");
    printf("Thread %d params: bx=%d by=%d ti=%d tj=%d\n", thread_id, bx, by, ti, tj);
    fflush(stdout);
    #endif
    
    //References to the thread arguments on all four sides.
    thread_args_t *t_left   = (tj == 0    ? 0 : thread_args[ti*tx+tj-1]);
    thread_args_t *t_right  = (tj == tx-1 ? 0 : thread_args[ti*tx+tj+1]);
    thread_args_t *t_top    = (ti == 0    ? 0 : thread_args[(ti-1)*tx+tj]);
    thread_args_t *t_bottom = (ti == ty-1 ? 0 : thread_args[(ti+1)*tx+tj]);
    
    //Border elements for neighboring blocks.
    //DOUBLE* left = 0, *right = 0, *top = 0, *bottom = 0;

    //If we're not on the edge of the array, point the border element pointers to the neighbor's border element array.
    /*if (tj == 0) left = (DOUBLE*)malloc(sizeof(DOUBLE)*by);
    else left = t_left->edgeRight;
    if (tj == tx-1) right = (DOUBLE*)malloc(sizeof(DOUBLE)*by);
    else right = t_right->edgeLeft;
    if (ti != 0) top = t_top->edgeBottom;
    if (ti != ty-1) bottom = t_bottom->edgeTop;*/
   
   
    while ((iterations < 0 && t < T) || niter < iterations) 
    {
        t += dt;
        niter++;
        
        //Wait until all are initialized and we are sure that the arrays are swapped
        //pthread_barrier_wait(&barr);
        //Copy ghost cells
        if (ti == 0)
            memcpy(E_prev[0], E_prev[2], sizeof(DOUBLE)*(bx+2));
        else
            memcpy(E_prev[0], (*(t_top->E_prev))[t_top->by], sizeof(DOUBLE)*(bx+2));
        if (ti == ty-1)
            memcpy(&E_prev[by+1][0], &E_prev[by-1][0], sizeof(DOUBLE)*(bx+2));
        else
            memcpy(&E_prev[by+1][0], (*(t_bottom->E_prev))[1], sizeof(DOUBLE)*(bx+2));
            
        //MT_PRINT("PING")
        if (tj == 0)
        {
            for (i = 1; i < by+1; i++)
                E_prev[i][0] = E_prev[i][2];
        }
        else
        {
            for (i = 1; i < by+1; i++)
                E_prev[i][0] = (*(t_left->E_prev))[i][t_left->bx];
        }

        if (tj == tx-1)
        {
            for (i = 1; i < by+1; i++)
                E_prev[i][bx+1] = E_prev[i][bx-1];
        }
        else
        {
            for (i = 1; i < by+1; i++)
                E_prev[i][bx+1] = (*(t_right->E_prev))[i][1];
        }
        
        #ifdef DEBUG
        if (DEBUG >=3)
        {
            printMat(E_prev, by+2, bx+2);
            printf("\n");
        }
        #endif
        
        
        
        //Solve the PDE.
        for (i = 1; i < by+1; i++)
        {
            #pragma ivdep
            for (j = 1; j < bx+1; j++)
            {
                E[i][j] = E_prev[i][j] + alpha * (E_prev[i][j + 1]+
                                          E_prev[i][j - 1]-
                                          4 * E_prev[i][j]+
                                          E_prev[i + 1][j]+
                                          E_prev[i - 1][j]);
            }
        }
        
        //MT_PRINT("PONG3")
        //Solve the ODE for one time step
        for (i = 1; i < by+1; i++)
        {
            #pragma ivdep
            for (j = 1; j < bx+1; j++)
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
        /*for (j = 0; j < bx; j++)
        {
            args->edgeTop[j] = E_prev[0][j];
            args->edgeBottom[j] = E_prev[by-1][j];
        }
        
        for (i = 0; i < by; i++)
        {
            args->edgeLeft[i] = E_prev[i][0];
            args->edgeRight[i] = E_prev[i][bx-1];
        }*/
        
        #ifdef DEBUG
        printf("Thread %d past barrier.\n", thread_id);
        fflush(stdout);
        #endif
        
        
        //Swap arrays
        DOUBLE **tmp = E;
        E = E_prev;
        E_prev = tmp;
    }

    #ifdef DEBUG
    MT_PRINT("Started copying back...");
    #endif

    args->iterations = niter;

    //Copy data from the global array to the local one
    for (i = 0; i < by; i++)
    {
        #pragma ivdep
        for (j = 0; j < bx; j++)
        {
            args->E_prev_global[offset_y+i][offset_x+j] = E_prev[i+1][j+1];
            args->R_global[offset_y+i][offset_x+j] = R[i+1][j+1];
        }
    }  
    
    printf("Thread %d exits.\n", thread_id);

    return NULL;
}

int solve(DOUBLE ***_E, DOUBLE ***_E_prev, DOUBLE **R, int m, int n, DOUBLE T, int iterations, DOUBLE alpha, DOUBLE dt, int do_stats, int plot_freq, int tx, int ty) 
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
 
    //The arguments are fixed for each iteration, so initialize them only once
    thread_args = (thread_args_t**)malloc(sizeof(thread_args_t*)*tx*ty);
    
    for (ti = 0; ti < ty; ti++)
    {
        for (tj = 0; tj < tx; tj++)
        {
            thread_args[ti*tx+tj] = (thread_args_t*)malloc(sizeof(thread_args_t));
            thread_args_t* ta = thread_args[ti*tx+tj];
            ta->thread_id = ti*tx+tj;
            ta->thread_x = tj; ta->thread_y = ti;
            ta->alpha = alpha; ta->dt = dt;
            ta->tx = tx; ta->ty = ty;
            ta->R_global = R; ta->E_prev_global = E_prev;
            ta->m = m; ta->n = n;
            ta->T = T; ta->iterations = iterations;
            

            //Start the threads
            pthread_create(&threads[ti*tx+tj], NULL, &solve_block, thread_args[ti*tx+tj]);
        }
    }
    
    //Join the threads when we're done
    for (ti = 0; ti < ty; ti++)
    {
        for (tj = 0; tj < tx; tj++)
        {
            pthread_join(threads[ti*tx+tj], NULL);
        }
    }
    niter = thread_args[0]->iterations;
    

    // We continue to sweep over the mesh until the simulation has reached
    // the desired simulation Time
    // This is different from the number of iterations
    /*while ((iterations < 0 && t < T) || niter < iterations) 
    {
        #if DEBUG
        printf("Main thread starts new loop.\n");
        #endif
        
        //printf("Iteration %d\n", niter);
        #ifdef DEBUG
        if (DEBUG >=3)
        {
            //printMat(E_prev, m, n);
            repNorms(E_prev, t, dt, m, n, niter);
            if (plot_freq) {
                splot(E_prev, t, niter, m + 1, n + 1, WAIT);
            }
            printf("\n");
        }
        #endif

        t += dt;
        niter++;

        //Wake up threads to do one iteration worth of work
        pthread_mutex_lock(&ready_mutex);
        #ifdef DEBUG
        printf("Main thread inited ghost cells, sets state to 1 and broadcasts ready.\n");
        #endif
        state = STATE_WORKING;
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
    }*/


    //Tell the threads that we're done and copy back the data
    /*pthread_mutex_lock(&ready_mutex);
    done_count = 0;
    state = STATE_DONE;
    pthread_cond_broadcast(&ready);
    pthread_mutex_unlock(&ready_mutex);

    pthread_mutex_lock(&done_mutex);
    #ifdef DEBUG
    printf("Main thread waits for copy back.\n");
    fflush(stdout);
    #endif

    while (done_count < tx*ty)
    {
        pthread_cond_wait(&done, &done_mutex);
    }

    #ifdef DEBUG
    printf("Main thread got signal indicating copy is done. Exiting.\n");
    fflush(stdout);
    #endif
    pthread_mutex_unlock(&done_mutex);*/

    // Store them into the pointers passed in
    *_E = E;
    *_E_prev = E_prev;

    return niter;
}
