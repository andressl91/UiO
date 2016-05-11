#include <stdio.h>
#include <stdlib.h>
#include "matmult.h"

void get_my_share(matrix * A, matrix * my_matrix, int start){
    int i, j, k;
    my_matrix->mat = malloc(my_matrix->rows*sizeof(float*)); /*Y-LENGTH ROWS*/
    for (k = 0; k < my_matrix->rows; k++) {
        my_matrix->mat[k] = malloc(my_matrix->cols*sizeof(float)); /*X-LENGTH COLUMNS*/
    }

    for (i = start; i < start + my_matrix->rows; i++) {
        for (j = 0; j < start + my_matrix->rows; j++) {
            my_matrix->mat[i][j] = A->mat[i][j];
            }
        }

}
