#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void array(int ***matrix ,int *n){
    int i, k;
    int j = 3;
    int ** A;
    int *x;
    int *my_x;
    A = malloc(j*sizeof(float*)); /*Y-LENGTH ROWS*/
    for (i = 0; i < j; i++) {
        A[i] = malloc(j*sizeof(float)); /*X-LENGTH COLUMNS*/
    }
    x = malloc(j*sizeof(float*));
    for (i = 0; i < j; i++) {
        for (k = 0; k < j; k++) {
        A[i][k] = i;
        }
    }
    A[0][0] = 1; A[0][1] = 2; A[0][2] = 3;
    A[1][0] = 4; A[1][1] = 5; A[1][2] = 6;
    n = &j;
    /*
    matrix = &A;
    printf("%d\n", *matrix[1][1]); */
}

int main(int nargs, char **args){
    int size, my_rank, num_procs, n;
    int num_rows, num_cols;
    int *my_vec, *my_row;
    int **matrix;


    MPI_Status status;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    /*
    if (my_rank == 0){
        array(&matrix, &n);
    }*/
    array(&matrix, &n);
    printf("%d\n", n);
    /*printf("%d\n", **matrix);*/
    MPI_Finalize();
    return 0;
}
