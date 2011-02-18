#include <assert.h>
#include <stdio.h>
#include "types.h"
//#include "apf.h"

void printMat(DOUBLE **U, int m, int n){
    int i,j;
    for (i=0; i<m; i++){
        for (j=0; j<n; j++) {
            printf("%6.3f ", U[i][j]);
        }
        printf("\n");
    }
}
