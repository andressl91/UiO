#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef struct
{
  int** matrix; /* a 2D array of floats */
  int** x;
  int** C;
  int m;               /* # m-rows */
  int n;               /* # n-columns */

} mat;

void getarray(mat * Ax){
    int i, j = 3;
    Ax->matrix = malloc(j*sizeof(float*)); /*Y-LENGTH ROWS*/
    for (i = 0; i < j; i++) {
        Ax->matrix[i] = malloc(j*sizeof(float)); /*X-LENGTH COLUMNS*/
    }

}

void main(int nargs, char **args){
    int size, my_rank, num_procs;
    int start;
    int num_rows, num_cols;
    int my_rows, my_cols;
    mat A, B, my_matrix;

    int n = 3;
    int m = 3;

    //int A[m][n];
    //int B[m][n];
    //A[1][1] = 1;

    MPI_Status status;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Datatype mat;
    MPI_Type_contiguous(m*n, MPI_INT, &mat);
    MPI_Type_commit(&mat);

    /*mat = Ax->A;*/
    if(my_rank == 0){
        getarray(&A);
        getarray(&B);
        A.matrix[0][0] = 1;
        A.matrix[1][0] = 2;

        B.matrix = A.matrix;
        printf("%d\n", B.matrix[0][0]);
        printf("%d\n", B.matrix[1][0]);
        //A.matrix[0][1] = 1; A.matrix[1][1] = 3; A.matrix[2][1] = 2;
        //MPI_Bcast(&arr[0][0], m*n, MPI_INT, 0, MPI_COMM_WORLD);
        //printf("%d\n", A.A[2][1]);
        //MPI_Send(&A.matrix, 1, mat, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&A.matrix[0][0], 1, mat, 1, 0, MPI_COMM_WORLD);

    }

    else{
        getarray(&B);
        //MPI_Recv(&B.matrix, 1, mat, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&B.matrix[0][0], 1, mat, 0, 0, MPI_COMM_WORLD, &status);
        printf("REC\n");
        printf("%d \n", B.matrix[1][0]);

        //printf("%d\n", A.matrix[1][1]);
    }

    MPI_Finalize();
}
