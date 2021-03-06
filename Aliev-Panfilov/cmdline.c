//#include <assert.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "types.h"

void cmdLine(int argc, char *argv[], double* T, int* n, int* tx, int* ty, int* bx, int* by, int* do_stats, int* plot_freq){
/// Command line arguments
 // Default value of the domain sizes
 static struct option long_options[] = {
        {"n", required_argument, 0, 'n'},
        {"tx", required_argument, 0, 'x'},
        {"ty", required_argument, 0, 'y'},
        {"bx", required_argument, 0, 'i'},
        {"by", required_argument, 0, 'j'},
        {"tfinal", required_argument, 0, 't'},
        {"stats", no_argument, 0, 's'},
        {"plot", required_argument, 0, 'p'}
 };
    // Process command line arguments
 int ac;
 for(ac=1;ac<argc;ac++) {
    int c;
    while ((c=getopt_long(argc,argv,"n:x:y:i:j:t:sp:",long_options,NULL)) != -1){
        switch (c) {

	    // Size of the computational box
            case 'n':
                *n = atoi(optarg);
                break;

	    // X processor geometry
            case 'x':
                *tx = atoi(optarg);
                break;

	    // X processor geometry
            case 'y':
                *ty = atoi(optarg);
                break;

	    // X blocking factor
            case 'i':
                *bx = atoi(optarg);
                break;

	    // Y processor geometry
            case 'j':
                *by = atoi(optarg);
                break;

	    // Length of simulation, in simulated time units
            case 't':
                *T = atof(optarg);
                break;

	    // Print various statistics
            case 's':
                *do_stats = 1;
                break;

	    // Plot the excitation variable
            case 'p':
                *plot_freq = atoi(optarg);
                break;

	    // Error
            default:
                printf("Usage: a.out [-n <domain size>] [-t <final time >]\n\t [-s print statistics ] [-p <plot frequency>]\n\t[-tx <x thread geometry> [-ty <y thread geometry] [-bx <x blocking factor>] [-by <y blocking factor]\n");
                exit(-1);
            }
    }
 }
}
