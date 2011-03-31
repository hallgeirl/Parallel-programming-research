/*
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 *
 * Modified and  restructured by Scott B. Baden, UCSD
 * 
 * Modified further by Hallgeir Lien
 */
#include "includes.h"
#include <pthread.h>
#include "apf.h"
#include "types.h"
#include "barrier.h"

typedef struct thread_args_s
{
    int thread_id, thread_x, thread_y;
    DOUBLE **E_prev_global, **R_global;  //Global arrays, for initialization
    DOUBLE ***E_prev;                    //For ghost cells
    int iterations;
    DOUBLE alpha;
    int tx, ty, m, n, bx, by; 
    DOUBLE dt;
    DOUBLE T;
} thread_args_t;

pthread_t *threads;
thread_args_t ** thread_args; // Need this to be global in order for threads to access neighboring blocks.

//pthread_barrier_t barr;
tour_barrier_t tour_barr;
int threadcount;
int done_count;

pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;

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
    //Synchronize before starting the timer
    //pthread_barrier_wait(&barr);
    tour_barrier(&tour_barr, thread_id);
    
    double t0 = -getTime();

    //References to the thread arguments on all four sides.
    thread_args_t *t_left   = (tj == 0    ? 0 : thread_args[ti*tx+tj-1]);
    thread_args_t *t_right  = (tj == tx-1 ? 0 : thread_args[ti*tx+tj+1]);
    thread_args_t *t_top    = (ti == 0    ? 0 : thread_args[(ti-1)*tx+tj]);
    thread_args_t *t_bottom = (ti == ty-1 ? 0 : thread_args[(ti+1)*tx+tj]);
    
    while ((iterations < 0 && t < T) || niter < iterations) 
    {
        t += dt;
        niter++;
        
        //Wait until all are initialized and we are sure that the arrays are swapped
#ifndef DISABLE_SYNC
        tour_barrier(&tour_barr, thread_id);
        //pthread_barrier_wait(&barr);
#endif

        //Copy ghost cells
#ifndef DISABLE_GHOST
        if (ti == 0)
            memcpy(E_prev[0], E_prev[2], sizeof(DOUBLE)*(bx+2));
        else
            memcpy(E_prev[0], (*(t_top->E_prev))[t_top->by], sizeof(DOUBLE)*(bx+2));
        if (ti == ty-1)
            memcpy(&E_prev[by+1][0], &E_prev[by-1][0], sizeof(DOUBLE)*(bx+2));
        else
            memcpy(&E_prev[by+1][0], (*(t_bottom->E_prev))[1], sizeof(DOUBLE)*(bx+2));
            
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
#endif
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
#ifndef DISABLE_SYNC
        tour_barrier(&tour_barr, thread_id);
        //pthread_barrier_wait(&barr);
#endif
        #ifdef DEBUG
        printf("Thread %d past barrier.\n", thread_id);
        fflush(stdout);
        #endif
        
        
        //Swap arrays
#ifndef DISABLE_GHOST
        DOUBLE **tmp = E;
        E = E_prev;
        E_prev = tmp;
#endif
    }
    fprintf(stderr, "Thread %d completes at t=%f.\n", thread_id, t0+getTime());
//#ifndef DISABLE_SYNC
    tour_barrier(&tour_barr, thread_id);
    //pthread_barrier_wait(&barr);
//#endif
    t0 += getTime();
    

    pthread_mutex_lock(&print_mutex);
    if (done_count == 0) 
    {
        printf("Running Time: %f sec.\n", t0);
        fflush(stdout);
    }
    done_count++;
    pthread_mutex_unlock(&print_mutex);

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
    
    fprintf(stderr, "Thread %d exits.\n", thread_id);

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
    //pthread_barrier_init(&barr, NULL, tx*ty);
    tour_barrier_init(&tour_barr, tx*ty);

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

    *_E = E;
    *_E_prev = E_prev;

    return niter;
}
