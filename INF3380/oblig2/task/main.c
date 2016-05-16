
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "matmult.h"

double **alloc2d(int m, int n, int a) {
    int i,j;
    double *data = malloc(n*m*sizeof(double));
    double **array = malloc(m*sizeof(double *));
    for (i=0; i<m; i++) {
        array[i] = &(data[i*n]);
    }
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            array[i][j] = i*a;
        }
    }
    return array;
}


void main(int nargs, char **args) {
    int size, my_rank, num_procs, n;
    int start_row, start_rowB;
    int i, j, row_tmp, rowB_tmp;
    //ROWS AND COLUMNS FOR THE A MATRIX
    int num_rows, num_cols;
    int my_rows, my_rowsB, my_cols;
    matrix A, B;
    matrix *my_matrix;
    /*double ** matA;
    double **matB;*/

    int a = 3;
    int b = 5;
    int c = 2;
    int f, g;

    MPI_Status status;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if(my_rank == 0){
        //Reads matrix binary file and applies row, column size in
        //typedef struct matrix
        /*read_matrix_binaryformat ("small_matrix_a.bin", &A.mat,
                                &A.rows, &A.cols);
        read_matrix_binaryformat ("small_matrix_b.bin", &B.mat,
                                &B.rows, &B.cols);
        num_rows = A.rows;
        num_cols = A.cols;
        */
        num_rows = a;
        num_cols = b;
        A.mat = alloc2d(a, b, c);
        B.mat = alloc2d(b, a, c+1);

        printf("INITIAL B MATRIX\n");
        for (f=0; f<b; f++) {
         printf("%f ",B.mat[f][0]);
            for (g=1; g<a; g++) {
                printf("%f ", B.mat[f][g]);
            }
             printf("\n");
        }

        }

    MPI_Bcast(&num_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /*printf("%g\n", full_mat.A[1][1]);*/


    my_rows = num_rows / num_procs;
    if ( my_rank < num_rows % num_procs) {
        ++my_rows;
    }
    //Notive inverted choice of cols and rows for the B matrix
    my_rowsB = num_cols / num_procs;
    if ( my_rank < num_cols % num_procs) {
        ++my_rowsB;
    }

    my_matrix = malloc(sizeof(my_matrix));
    my_matrix->mat = allocate_matrix(my_rows, num_cols);
    //Notive more memory allocated due to every process need a full copy of B
    //Each process needs full array of B, gets after alltoall Broadcast
    my_matrix->matB = allocate_matrix(my_rowsB, num_rows);
    my_matrix->wholeB = allocate_matrix(num_cols, num_rows);



    //Distribute chunks of matrix A and B
    if(my_rank == 0){
        start_row = my_rows;
        start_rowB = my_rowsB;
        printf("thread %d with row %d is given %d as start\n", my_rank, my_rowsB, start_rowB);
        for (j = 1; j < num_procs; j++){
            MPI_Recv(&row_tmp, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&rowB_tmp, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);

            //SOMETHING WRONG WITH PARTITIONS OF B


            printf("thread %d with row %d is given %d as start\n", j, rowB_tmp, start_rowB);

            my_matrix->mat = get_my_share(&A, my_matrix->mat, start_row, row_tmp, num_cols);
            MPI_Send(&(my_matrix->mat[0][0]), row_tmp*num_cols, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            //The rows of matrix A is the cols of matri B
            my_matrix->matB = get_my_share(&B, my_matrix->matB, start_rowB, rowB_tmp, num_rows);
            MPI_Send(&(my_matrix->matB[0][0]), rowB_tmp*num_rows, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            start_row = start_row + row_tmp;
            start_rowB = start_rowB + rowB_tmp;
        }
        //Need to reset Thread 0 start value before
        //assigning matrix to local matrix
        start_row = 0;
        start_rowB = 0;
        my_matrix->mat = get_my_share(&A, my_matrix->mat, start_row, my_rows, num_cols);
        my_matrix->matB = get_my_share(&B, my_matrix->matB, start_rowB, my_rowsB, num_rows);

        printf("THREAD %d B matrix\n", my_rank);
        for (f=0; f<my_rowsB; f++) {
         printf("%f ",my_matrix->matB[f][0]);
            for (g=1; g<num_rows; g++) {
                printf("%f ", my_matrix->matB[f][g]);
            }
             printf("\n");
        }
    }

    else {
        MPI_Send(&my_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&my_rowsB, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

        MPI_Recv(&(my_matrix->mat[0][0]), my_rows*num_cols, MPI_DOUBLE, 0, 0,
                MPI_COMM_WORLD, &status);
        MPI_Recv(&(my_matrix->matB[0][0]), my_rowsB*num_rows, MPI_DOUBLE, 0, 0,
                MPI_COMM_WORLD, &status);

        printf("THREAD %d B matrix\n", my_rank);
        for (f=0; f<my_rowsB; f++) {
         printf("%f ",my_matrix->matB[f][0]);
            for (g=1; g<num_rows; g++) {
                printf("%f ", my_matrix->matB[f][g]);
            }
             printf("\n");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("Thread %d, with %d rows of matrix B \n", my_rank, my_rowsB);

/*
    MPI_Alltoall(&(my_matrix->matB[0][1]), 3, MPI_DOUBLE,
                  &(my_matrix->wholeB[0][0]), 3, MPI_DOUBLE,
                  MPI_COMM_WORLD);
                  */
    if(my_rank == 0){
        int e,r;
      printf("THREAD %d WHOLE B matrix AFTER ALLTOALL\n", my_rank);
      for (e=0; e< num_cols; e++) {
       printf("%f ",my_matrix->wholeB[e][0]);
          for (r=1; r < num_rows; r++) {
              printf("%f ", my_matrix->wholeB[e][r]);
          }
           printf("\n");
      }
  }



    //deallocate_matrix(my_matrix);

    MPI_Finalize();
}
