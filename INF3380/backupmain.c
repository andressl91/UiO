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

double **alloc2d(int n, int m) {
    int i;
    double *data = malloc(n*m*sizeof(double));
    double **array = malloc(n*sizeof(double *));
    for (i=0; i<n; i++) {
        array[i] = &(data[i*m]);
    }
    return array;
}

void main(int nargs, char **args){
    int size, my_rank, num_procs;
    int start;
    int num_rows, num_cols;
    int my_rows, my_cols;
    mat A, B, my_matrix;

    int n = 3;
    int m = 3;
/*
    typedef struct car_s {
        int shifts;
        int topSpeed;
    } car;*/


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
        printf("SEND\n");
        mat *send;
        send = malloc(sizeof(mat));
        send->n = 1;
        send->matrix = alloc2d(3, m);
        for (i=0; i<n; i++)
            for (j=0; j<m; j++)
                send->matrix[i][j] = i*j;
        send->matrix[4][2] = 1;
        MPI_Send(&(send->matrix[0][0]), m*n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);

    }

    else{
        int i, j;
        printf("REC\n");
        mat *recv;
        recv = malloc(sizeof(mat));
        recv->matrix = alloc2d(n, m);

        //MPI_Recv(&(recv->offset), 1, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&(recv->matrix[0][0]), n*m, MPI_DOUBLE, 0, 0,
                MPI_COMM_WORLD, &status);
                //MPI_Recv(&(send->matrix[0][0]), m*n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        for (i=0; i<n; i++)
            for (j=0; j<m; j++)
                printf("FOR A[%d][%d]%f\n" , i, j, recv->matrix[i][j]);
        //printf("MATRIX SEND value %d\n", recv.matrix[1][1]);
    }

    MPI_Finalize();
}
