#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "matrix.h"

void matrix_multiply(int m, int n, double *A, double *B, double *C) {
    int i, j, k;    
    
    #pragma omp parallel default (private) shared (A, B, C, m, n) num_threads(4)
    {
        #pragma omp for schedule (static)
        for ( i = 0; i < m; ++i ) {
            for ( j = 0; j < n; ++j ) {
                for ( k = 0; k < n; ++k ) {
                    C[i*m + j] += A[i*m + k] * B[k*m + j];
                }
            }
        }
    }
}
