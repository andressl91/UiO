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
    printf("ALLOCATED \n");
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

    int m = 4;
    int n = 2;

    MPI_Status status;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Datatype matr;
    MPI_Type_contiguous(m*n, MPI_INT, &matr);
    MPI_Type_commit(&matr);

    /*mat = Ax->A;*/
    if(my_rank == 0){
        int i, j;
        int a = 3;
        mat *sendA, *sendB;
        sendA = malloc(sizeof(sendA));
        sendA->matrix = alloc2d(m, n);

        sendB = malloc(sizeof(sendB));
        sendB->matrix = alloc2d(n, m);
        for (i=0; i<m; i++)
            for (j=0; j<n; j++)
                sendA->matrix[i][j] = i*a;

        for (i=0; i<n; i++)
            for (j=0; j<m; j++)
                sendB->matrix[i][j] = j*a;

        MPI_Send(&(sendA->matrix[0][0]), m*n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&(sendB->matrix[0][0]), m*n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);

        printf("SEND\n");
    }

    else{
        int i, j;

        mat *recvA, *recvB;
        recvA = malloc(sizeof(recvA));
        recvA->matrix = alloc2d(m, n);

        recvB = malloc(sizeof(recvB));
        recvB->matrix = alloc2d(n, m);


        //MPI_Recv(&(recv->offset), 1, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&(recvA->matrix[0][0]), m*n, MPI_DOUBLE, 0, 0,
                MPI_COMM_WORLD, &status);
                //MPI_Recv(&(send->matrix[0][0]), m*n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&(recvB->matrix[0][0]), m*n, MPI_DOUBLE, 0, 0,
                MPI_COMM_WORLD, &status);
            printf("REC\n");
            
        for (i=0; i<m; i++) {
        printf("%f  ",recvA->matrix[i][0]);
            for (j=1; j<n; j++){
                printf("%f " , recvA->matrix[i][j]);
            }
            printf("\n");
        }

        for (i=0; i<n; i++) {
        printf("%f  ",recvB->matrix[i][0]);
            for (j=1; j<m; j++){
                printf("%f " , recvB->matrix[i][j]);
            }
            printf("\n");
        }
    }
    MPI_Finalize();
}
