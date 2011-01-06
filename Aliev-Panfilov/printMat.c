#include <assert.h>
#include <stdio.h>
#include "types.h"

void printMat(DOUBLE **U, int m, int n){
    int i,j;
    for (j=1; j<=m; j++){
        for (i=1; i<=n; i++) {
            printf("%6.3f ", U[j][i]);
        }
        printf("\n");
    }
}
