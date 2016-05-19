#include <stdio.h>
#include <stdlib.h>
#include "matmult.h"
#include <omp.h>


void find_sum(double ** C, double **A, double **B, int my_rows, int num_rows, int num_cols){
    int i, j, k;
    int count = 0;
    double sum = 0;

    #pragma omp parallel default (private) shared (A, B, C, num_rows, num_cols) num_threads(4)
    {
    #pragma omp for schedule (static)
    for(k=0; k < my_rows; k++){
        for(i=0; i < num_rows; i++){
            for(j = 0; j < num_cols; j++){
                C[k][j] = C[k][j] + A[k][i]*B[i][j];

                }
            }
        }
    }
}
