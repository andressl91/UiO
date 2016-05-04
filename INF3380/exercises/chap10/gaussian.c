#include <stdio.h>
#include <stdlib.h>

void Gaussian(double **A, double *b, double *y){
    int m = 3;
    int k, j, i;
    for(k = 0; k < m; k++){
        for(j = k+1; j <m; j++){
            A[k][j] = A[k][j]/A[k][k];
        y[k] = b[k]/A[k][k];
        A[k][k] = 1;
        }
        for(i = k + 1; i < m; i++ ){
            for(j = k +1; j < m; j++){
                A[i][j] = A[i][j] - A[i][k]*A[k][j];
            }
            b[i] = b[i] -A[i][k]*y[k];
            A[i][k] = 0;
        }
    }
}

void main(){
    int i, j;
    int m = 3;
    int n = 3;
    double **A = malloc(m*sizeof(double*));
        for (i = 0; i < m; i++){
            A[i] = malloc(n*sizeof(double));
        }
    double *b = malloc(m*sizeof(double));
    double *y = malloc(m*sizeof(double));
    b[0] = 8; b[1] = 11; b[2] = -3;

    A[0][0] = 2; A[0][1] = 1; A[0][2] = -1;
    A[1][0] = -3; A[1][1] = -1; A[1][2] = 2;
    A[2][0] = -2; A[2][1] = 1; A[2][2] = 2;

    Gaussian(A, b, y);

    for (i = 0; i < m; i++){
        printf("b vecotr %f \n", y[i]);
        for(j = 0; j < n; j++){
            printf("A[%d][%d] = %.f ", i, j, A[i][j]);
        }
        printf("\n");
    }
}
