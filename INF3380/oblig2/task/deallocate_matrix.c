#include <stdio.h>
#include <stdlib.h>
#include "matmult.h"


void deallocate_matrix (matrix * matrix)
{
    int i;
        for (i = 0; i < matrix->rows; i++){
            free(matrix->mat[i]);
        }
        free(matrix->mat);
}
