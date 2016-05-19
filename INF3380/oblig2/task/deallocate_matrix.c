#include <stdio.h>
#include <stdlib.h>
#include "matmult.h"


void deallocate_matrix (double **mat, int m, int n){
    int i;
    for (i = 0; i < n; i++){
        free(mat[i]);
        }
    free(mat);

}
