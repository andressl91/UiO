#include <stdio.h>
#include <stdlib.h>
#include "matmult.h"



double ** allocate_matrix(int rows, int cols)
{
    int i;
    double *data = malloc(rows*cols*sizeof(double));
    double **array = malloc(rows*sizeof(double *));
    for (i=0; i<rows; i++) {
        array[i] = &(data[i*cols]);
    }
    return array;
}
