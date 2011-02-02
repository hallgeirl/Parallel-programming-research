#include <assert.h>
#include <stdio.h>
#include "types.h"
#include "apf.h"

void printMat(DOUBLE **U, int m, int n){
    int i,j;
    for (j=0; j<m; j++){
        for (i=0; i<n; i++) {
            printf("%6.3f ", U[j][i]);
        }
        printf("\n");
    }
}
