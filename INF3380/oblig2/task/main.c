
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
    matrix A, B;
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

/*Idea allocate_and write own matrix */
    matrix *my_matrix;
    my_matrix = malloc(sizeof(my_matrix));
    my_matrix->mat = allocate_matrix(my_rows, my_cols);

    if(my_rank == 0){
    /*    matrix *send;
        send = malloc(sizeof(matrix));
        send->mat = allocate_matrix(my_rows, my_cols);*/
        start = 0;
        for (j = 1; j < num_procs; j++){
            MPI_Recv(&row_tmp, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
            start = start + row_tmp;
            my_matrix->mat = get_my_share(&A, my_matrix->mat, start, row_tmp, my_cols);
            MPI_Send(&(my_matrix->mat[0][0]), row_tmp*my_cols, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            //MPI_Send(&start, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
            /*
            MPI_Recv(&row_tmp, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
            start = start + row_tmp;
            send->mat = get_my_share(&A, send->mat, start, row_tmp, my_cols);
            MPI_Send(&(send->mat[0][0]), row_tmp*my_cols, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            //MPI_Send(&start, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
            free(send);*/

        }
        //Need to reset Thread 0 start value before
        //assigning matrix to local matrix
        start = 0;
    }

    else {
        MPI_Send(&my_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(&(my_matrix->mat[0][0]), my_rows*my_cols, MPI_DOUBLE, 0, 0,
                MPI_COMM_WORLD, &status);
        int l,m;
        //for (l=0; l < my_rows; l++)
            //for (m=0; m<my_cols; m++)
                //printf("FOR A[%d][%d]%f\n" , l, m, my_matrix->mat[l][m]);
    }
    printf("Thread %d, with rows%d, start at %d \n", my_rank, my_rows, start);
    //deallocate_matrix(&my_matrix);
    
    MPI_Finalize();
}
