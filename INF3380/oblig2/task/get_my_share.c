#include <stdio.h>
#include <stdlib.h>
#include "matmult.h"

double ** get_my_share(matrix * A, double ** my_matrix, int start, int row, int cols){
    int i, j;

    for (i = start; i < row; ++i) {
        for (j = 0; j < cols; ++j) {
            my_matrix[i-start][j] = A->mat[i][j];
            }
        }
    my_matrix[40][49] = 1;
    return my_matrix;
}
