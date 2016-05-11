#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "matmult.h"

void main(int nargs, char **args) {
    int size, my_rank, num_procs, n;
    int start;
    int j, row_tmp;
    int num_rows, num_cols;
    int my_rows, my_cols;
    matrix A, B, my_matrix;
    /*double ** matA;
    double **matB;*/



    MPI_Status status;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if(my_rank == 0){
        //Reads matrix binary file and applies row, column size in
        //typedef struct matrix
        read_matrix_binaryformat ("small_matrix_a.bin", &A.mat,
                                &A.rows, &A.cols);
        read_matrix_binaryformat ("small_matrix_b.bin", &B.mat,
                                &B.rows, &B.cols);
        num_rows = A.rows;
        num_cols = A.cols;
        }
    MPI_Bcast(&num_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /*printf("%g\n", full_mat.A[1][1]);*/


    my_rows = num_rows / num_procs;
    if ( my_rank < num_rows % num_procs) {
        ++my_rows;
    }
    my_cols = num_cols;
    my_matrix.rows = my_rows;
    my_matrix.cols = my_cols;

/*Idea allocate_and write own matrix */
    if(my_rank == 0){
        start = my_rows;
        for (j = 1; j < num_procs; j++){
            MPI_Recv(&row_tmp, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
            MPI_Send(&start, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
            start = start + row_tmp;

        }
        //Need to reset Thread 0 start value before
        //assigning matrix to local matrix
        start = 0;
    }

    else {
        MPI_Send(&my_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(&start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    printf("Thread %d, with rows%d, start at %d \n", my_rank, my_rows, start);
    /*
    if(my_rank == 0){
        int a, b;
        get_my_share(&A, &my_matrix, start);
        for (a = 0; a < my_matrix.rows; a++) {
            for (b = 0; b < my_matrix.cols; b++) {
                printf("%f\n", my_matrix.mat[a][b]);
                }
            }
    }
    */
    //deallocate_matrix(&my_matrix);

    MPI_Finalize();
}
