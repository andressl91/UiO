#include <stdio.h>
#include <stdlib.h>
#include "matmult.h"

double ** get_my_share(matrix * A, double ** my_matrix, int start, int row, int cols){
    int i, j;

    for (i = start; i < start + row; ++i) {
        for (j = 0; j < cols; ++j) {
            my_matrix[i-start][j] = A->mat[i][j];
            }
        }
    return my_matrix;
}
