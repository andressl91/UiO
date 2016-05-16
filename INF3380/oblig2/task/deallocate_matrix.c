#include <stdio.h>
#include <stdlib.h>
#include "matmult.h"


void deallocate_matrix (matrix * matr)
{
    int i, j;
        for (i = 0; i < matr->rows; i++){
            free(matr->mat[i]);
        }
        free(matr->mat);
/*
        for (j = 0; j < mat->cols; j++){
            free(mat->matB[j]);
        }
        free(mat->matB);
        */
}
