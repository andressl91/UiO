#include <stdio.h>
#include <stdlib.h>
#include "matmult.h"


void find_sum(double ** C, double **A, double **B, int my_rows, int num_rows, int num_cols){
    int i, j, k;
    int count = 0;
    int sum;
    for(k=0; k < my_rows; k++){
        for(i=0; i < num_rows; i++){
            for(j = 0; j < num_cols; j++){
                sum = sum + A[k][j]*B[j][i];
                }
            C[k][count] = sum;
            count += 1;
            }
            count = 0;

    }
}
