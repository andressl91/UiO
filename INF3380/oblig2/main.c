
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "matmult.h"
//http://stackoverflow.com/questions/31890523/how-to-use-mpi-gatherv-for-collecting-strings-of-diiferent-length-from-different
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
    int sum;
    //ROWS AND COLUMNS FOR THE A MATRIX
    int num_rows, num_cols;
    int my_rows, my_rowsB, my_cols;
    matrix A, B;
    matrix *my_matrix;

    /*double ** matA;
    double **matB;*/

    int a = 5;
    int b = 3;
    int c = 1;
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


        printf("INITIAL A MATRIX\n");
        for (f=0; f<a; f++) {
         printf("%f ",A.mat[f][0]);
            for (g=1; g<b; g++) {
                printf("%f ", A.mat[f][g]);
            }
             printf("\n");
        }

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
    /*
     allocate_array(&(my_matrix->mat), my_rows, num_cols);
     allocate_array(&(my_matrix->matB) , my_rows, num_cols);
     allocate_array(&(my_matrix->wholeB), num_cols, num_rows);
*/

//SJEKK ALLOKERING AV C
    //Notive more memory allocated due to every process need a full copy of B
    //Each process needs full array of B, gets after alltoall Broadcast

    my_matrix = malloc(sizeof(my_matrix));
        my_matrix->mat = allocate_matrix(my_rows, num_cols);
    my_matrix->matB = allocate_matrix(my_rowsB, num_rows);
    my_matrix->wholeB = allocate_matrix(num_cols, num_rows);



    //Distribute chunks of matrix A and B
    if(my_rank == 0){
        start_row = my_rows;
        start_rowB = my_rowsB;
        for (j = 1; j < num_procs; j++){
            MPI_Recv(&row_tmp, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&rowB_tmp, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
            my_matrix->mat = get_my_share(&A, my_matrix->mat, start_row, row_tmp, num_cols);
            MPI_Send(&(my_matrix->mat[0][0]), row_tmp*num_cols, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            //The rows of matrix A is the cols of matri B
            //CHEATING, IMPLEMENT A ALLTOALL LAITER
            my_matrix->matB = get_my_share(&B, my_matrix->wholeB, 0, num_cols, num_rows);
            MPI_Send(&(my_matrix->wholeB[0][0]), num_cols*num_rows, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            //my_matrix->matB = get_my_share(&B, my_matrix->matB, start_rowB, rowB_tmp, num_rows);
            //MPI_Send(&(my_matrix->matB[0][0]), rowB_tmp*num_rows, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            start_row = start_row + row_tmp;
            start_rowB = start_rowB + rowB_tmp;
        }
        //Need to reset Thread 0 start value before
        //assigning matrix to local matrix
        start_row = 0;
        start_rowB = 0;
        my_matrix->mat = get_my_share(&A, my_matrix->mat, start_row, my_rows, num_cols);
        my_matrix->wholeB = get_my_share(&B, my_matrix->wholeB, 0, num_cols, num_rows);
        //my_matrix->matB = get_my_share(&B, my_matrix->matB, start_rowB, my_rowsB, num_rows);

        printf("THREAD %d A matrix\n", my_rank);
        for (f=0; f<my_rows; f++) {
         printf("%f ",my_matrix->mat[f][0]);
            for (g=1; g<num_cols; g++) {
                printf("%f ", my_matrix->mat[f][g]);
            }
             printf("\n");
        }
    }

    else {
        MPI_Send(&my_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&my_rowsB, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

        MPI_Recv(&(my_matrix->mat[0][0]), my_rows*num_cols, MPI_DOUBLE, 0, 0,
                MPI_COMM_WORLD, &status);
        MPI_Recv(&(my_matrix->wholeB[0][0]), num_rows*num_cols, MPI_DOUBLE, 0, 0,
                MPI_COMM_WORLD, &status);
        //MPI_Recv(&(my_matrix->matB[0][0]), my_rowsB*num_rows, MPI_DOUBLE, 0, 0,
        //        MPI_COMM_WORLD, &status);
        printf("THREAD %d A matrix\n", my_rank);
        for (f=0; f<my_rows; f++) {
         printf("%f ",my_matrix->mat[f][0]);
            for (g=1; g<num_cols; g++) {
                printf("%f ", my_matrix->mat[f][g]);
            }
             printf("\n");
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double **y_vec = NULL;
    y_vec = allocate_matrix(my_rows, num_rows);
    find_sum(y_vec, my_matrix, my_rows, num_rows, num_cols);


  printf("THREAD %d result matrix\n", my_rank);
  for (f=0; f<my_rows; f++) {
   printf("%f ",y_vec[f][0]);
      for (g=1; g<num_rows; g++) {
          printf("%f ", y_vec[f][g]);
      }
       printf("\n");
  }



    int my_len = my_rows*num_rows;
    int *recvcounts = NULL;
    if (my_rank == 0)
       recvcounts = malloc( num_procs * sizeof(int)) ;
     MPI_Gather(&my_len, 1, MPI_INT,
          recvcounts, 1, MPI_INT,
          0, MPI_COMM_WORLD);



  /*
    * Figure out the total length of string,
    * and displacements for each rank
    */
  int totlen = 0;
   int *displs = NULL;
   double ** y_vector = NULL;


   if (my_rank == 0) {
       displs = malloc( num_procs * sizeof(int) );
       y_vector = allocate_matrix(num_rows, num_rows);
       displs[0] = 0;
       totlen += recvcounts[0];//+1;

       for (i=1; i<num_procs; i++) {
          totlen += recvcounts[i];//+1;   /* plus one for space or \0 after words */

          displs[i] = displs[i-1] + recvcounts[i-1];//+ 1;
       }

       printf("TOTLEN %d and %d \n", totlen, num_rows*num_rows);
    }


    // MPI_Gatherv(&(y_vec[0][0]), my_len, MPI_DOUBLE,
                // &(y_vector[0][0]), recvcounts, displs, MPI_DOUBLE,
                // 0, MPI_COMM_WORLD);

/*
if(my_rank == 0){
        printf("RESULT OF MATRIX \n");
        for (f=0; f<num_rows-1; f++) {
         printf("%f ", y_vector[f][0]);
            for (g=1; g<num_rows; g++) {
                printf("%f ", y_vector[f][g]);
            }
             printf("\n");
         }
     }

*/
    //printf("Thread %d, with %d rows of matrix B \n", my_rank, my_rowsB);

/*
    for (i=0; i < num_procs; i++) {
        if(i == my_rank) {

        printf("THREAD %d B matrix\n", my_rank);
        for (f=0; f<my_rows; f++) {
         printf("%f ",my_matrix->matC[f][0]);
            for (g=1; g<num_cols; g++) {
                printf("%f ", my_matrix->matC[f][g]);
            }
             printf("\n");
         }
        }


        }
        */

/*
    MPI_Alltoall(&(my_matrix->matB[0][1]), 3, MPI_DOUBLE,
                  &(my_matrix->wholeB[0][0]), 3, MPI_DOUBLE,
                  MPI_COMM_WORLD);
                  */



    //deallocate_matrix(my_matrix);

    MPI_Finalize();
}
