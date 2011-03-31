#include "barrier.h"
#include "includes.h"

//#define BARRIER_DEBUG

#ifdef BARRIER_DEBUG
#define BARRIER_TRACE(s) printf s; printf("\n"); fflush(stdout); 
#define THREAD_PRINT(s, tid) printf("%s [tid=%d]\n", s, tid); fflush(stdout);
#else
#define BARRIER_TRACE(s)
#define THREAD_PRINT(s,tid)
#endif

//Returns the id of the parent node in the tree given a node id
inline int node_parent(int id)
{
    return (id-1)/2;
}

void tour_barrier_init(tour_barrier_t* barrier, int threads)
{
    BARRIER_TRACE(("tour_barrier_init(barrier=%p, threads=%d)", barrier, threads));

    int i, j;
    barrier->n_threads = threads;

    if (threads > 1)
    {
        //We have n-1 inner nodes ("rounds") for n threads
        barrier->tree = (tour_barrier_node_t**)malloc(sizeof(tour_barrier_node_t*)*(threads-1));
        for (i = 0; i < threads-1; i++)
        {
            barrier->tree[i] = (tour_barrier_node_t*)malloc(sizeof(tour_barrier_node_t));
            memset(barrier->tree[i], 0, sizeof(tour_barrier_node_t));
        }

        //Each thread has its own sense variable, initially true (1)
        barrier->sense = (bool**)malloc(sizeof(bool*)*threads); 
        for (i = 0; i < threads; i++)
        {
            barrier->sense[i] = (bool*)malloc(sizeof(bool));
            *(barrier->sense[i]) = true;
        }

        int logP = (int)ceil(log2((double)threads));

        //Set up the tree 
        for (i = 0; i < logP; i++)
        {
            for (j = 0; j < (1<<i); j++)
            {
                int nid = (1<<i)-1+j; // node ID
                barrier->tree[nid]->winner = j*(1<<(logP-i));
                if (nid != 0) barrier->tree[nid]->parent = barrier->tree[node_parent(nid)];
                else barrier->tree[nid]->parent = NULL;
            }
        }
    }
    BARRIER_TRACE(("tour_barrier_init exits."));
}

//void barrier_helper(barrier_t* barrier, int tid, int depth)
void tour_barrier_helper(tour_barrier_t* barrier, tour_barrier_node_t* my_node, int tid)
{
    BARRIER_TRACE(("tour_barrier_helper(barrier=%p, my_node=(addr: %p, winner: %d, arrive=%d), tid=%d)", barrier, my_node, my_node->winner, my_node->arrive, tid));

    //One thread will spin on the flags, the other one will set the arrival flag and traverse up the tree.
    if (tid == my_node->winner)
    {
        THREAD_PRINT("winner", tid);
        //Winner must wait for the second thread to arrive
        while (my_node->arrive != *(barrier->sense[tid])) { /*printf("arrive=%d, sense=%d\n", my_node->arrive, barrier->sense[tid]); fflush(stdout);*/}
        THREAD_PRINT("parent check", tid);

        //then either set the release flag (if on top => all done) or traverse up the tree
        if (my_node->parent == NULL)
        {
            THREAD_PRINT("release", tid);
            barrier->release = *(barrier->sense[tid]);
        }
        else
        {
            THREAD_PRINT("recurse", tid);
            tour_barrier_helper(barrier, my_node->parent, tid);
        }
    }
    else
    {
        THREAD_PRINT("loser", tid);
        //Loser sets the arrival flag, then waits for the release flag
        my_node->arrive = *(barrier->sense[tid]);

        THREAD_PRINT("wait for release", tid);
       
        while (barrier->release != *(barrier->sense[tid])) {/*printf("release=%d, sense=%d\n", my_node->arrive, barrier->sense[tid]); fflush(stdout);*/}
    }
    BARRIER_TRACE(("tour_barrier_helper exits (tid=%d)", tid));
}

void tour_barrier(tour_barrier_t* barrier, int tid)
{
    BARRIER_TRACE(("tour_barrier(barrier=(addr: %p, release: %d), tid=%d, sense=%d)", barrier, barrier->release, tid, barrier->sense[tid]));
    if (barrier->n_threads > 2)
    {
        tour_barrier_node_t* my_node = barrier->tree[node_parent(barrier->n_threads - 1 + tid)];
        tour_barrier_helper(barrier, my_node, tid);
        *(barrier->sense[tid]) = !(*(barrier->sense[tid]));
    }
    BARRIER_TRACE(("tour_barrier exists (tid=%d)", tid));
}

