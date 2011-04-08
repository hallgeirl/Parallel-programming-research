#ifndef __BARRIER_H__
#define __BARRIER_H__

#include <stdbool.h>

#define BARRIER_TOURNAMENT 0
#define BARRIER_PTHREAD 1

#ifndef BARRIER
#define BARRIER 0
#endif

//Tournament barrier
typedef struct tour_barrier_node_s
{
    bool volatile arrive;
    int winner;
    struct tour_barrier_node_s* parent;
} tour_barrier_node_t;

typedef struct tour_barrier_s
{
    int n_threads;
    bool** sense;
    bool volatile release;
    tour_barrier_node_t **tree;
} tour_barrier_t;

void tour_barrier_init(tour_barrier_t* barrier, int threads);
void tour_barrier();

#endif
