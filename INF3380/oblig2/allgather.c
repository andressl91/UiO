#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stddef.h>

typedef struct
{
  double** matrix; /* a 2D array of floats */
  int m;               /* # m-rows */
  int n;               /* # n-columns */

} mat;
/*
void getarray(mat ** Ax){
    int i, j = 3;
    Ax->matrix = (int**)malloc(j*sizeof(int*)); /*Y-LENGTH ROWS*/
    /*for (i = 0; i < j; i++) {
        Ax->matrix[i] = (int*)malloc(j*sizeof(int)); /*X-LENGTH COLUMNS*/
    /*}
}*/

double **alloc2d(int m, int n) {
    int i;
    double *data = malloc(n*m*sizeof(double));
    double **array = malloc(m*sizeof(double *));
    for (i=0; i<m; i++) {
        array[i] = &(data[i*n]);
    }
    return array;
}

void main(int nargs, char **args){
    int size, my_rank, num_procs;
    int start;
    int num_rows, num_cols;
    int my_rows, my_cols;
    mat A, B, my_matrix;
/*
    typedef struct car_s {
        int shifts;
        int topSpeed;
    } car;*/

    int m = 6;
    int n = 4;

    MPI_Status status;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Datatype matr;
    MPI_Type_contiguous(m*n, MPI_INT, &matr);
    MPI_Type_commit(&matr);

    mat *sendA, *sendB;
    mat *recvA, *recvB;
    recvA = malloc(sizeof(recvA));
    recvA->matrix = alloc2d(m, n);

    sendA = malloc(sizeof(sendA));
    sendA->matrix = alloc2d(2, n);

    sendB = malloc(sizeof(sendB));
    sendB->matrix = alloc2d(2, n);

    recvB = malloc(sizeof(recvB));
    recvB->matrix = alloc2d(m,n);

    if(my_rank == 0){
        int i, j;
        int a = 3;

        for (i=0; i<2; i++)
            for (j=0; j<n; j++)
                sendA->matrix[i][j] = i*a;

        printf("MATRIX PROCESS 0\n");
        for (i=0; i<2; i++) {
        printf("%f  ",sendA->matrix[i][0]);
            for (j=1; j<n; j++){
                printf("%f " , sendA->matrix[i][j]);
            }
            printf("\n");
        }
    }

      if(my_rank == 1){
          int i, j;
          int a = 2;

          for (i=0; i<2; i++)
              for (j=0; j<n; j++)
                  sendA->matrix[i][j] = i*a;
                  printf("MATRIX PROCESS 1\n");
          for (i=0; i<2; i++) {
          printf("%f  ",sendA->matrix[i][0]);
              for (j=1; j<n; j++){
                  printf("%f " , sendA->matrix[i][j]);
              }
              printf("\n");
          }

              }





    MPI_Barrier(MPI_COMM_WORLD);
    printf("BARRIER\n" );
    MPI_Alltoall(&(sendA->matrix[0][0]), 2*n, MPI_DOUBLE,
                  &(recvB->matrix[0][0]), 2*n, MPI_DOUBLE,
                  MPI_COMM_WORLD);


                  /*
                  MPI_Alltoall(&(sendA->matrix[0][0]), 2*n, MPI_DOUBLE,
                                &(recvB->matrix[0][0]), 2*n, MPI_DOUBLE,
                                MPI_COMM_WORLD);
                  */
    /*
    MPI_Allgather(&(sendA->matrix[1][0]), 1, MPI_DOUBLE,
                  &(recvB->matrix[0][0]), 2, MPI_DOUBLE,
                  MPI_COMM_WORLD);
                  */

    if(my_rank == 0){
        int i, j;

        printf("MATRIX PROCESS B AFTER ALLTOALL\n");
        for (i=0; i<m; i++) {
        printf("%f  ",recvB->matrix[i][0]);
            for (j=1; j<n; j++){
                printf("%f " , recvB->matrix[i][j]);
            }
            printf("\n");
        }
    }



    MPI_Finalize();
}
