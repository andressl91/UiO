#include <stdio.h>

void block_mat_mult(int m, int n, double **A, double **B, double **C) {
    //A is m x n
    //B is n x m
    //C is m x m
    int i = 0, j = 0, k = 0;

    /* c_{i,j} = \sum_{k=0}^{l-1}a_{i,k}b_{k,j} */

    for ( i = 0; i < m; ++i ) {
        for ( j = 0; j < n; ++j ) {
            C[i][j] = 0;
            for ( k = 0; k < n; ++k ) {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}
