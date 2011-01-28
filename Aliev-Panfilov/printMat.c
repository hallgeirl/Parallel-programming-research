#include <assert.h>
#include <stdio.h>
#include "types.h"

void printMat(DOUBLE **U, int m, int n){
    int i,j;
    for (j=0; j<=m+1; j++){
        for (i=0; i<=n+1; i++) {
            printf("%6.50f ", U[j][i]);
        }
        printf("\n");
    }
}
