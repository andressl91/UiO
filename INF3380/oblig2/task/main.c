#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "matmult.h"

void main(int nargs, char **args) {
    int size, my_rank, num_procs, n;
    matsys Ax;
    double ** mat;
    int num_rows, num_cols;


    MPI_Status status;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    read_matrix_binaryformat ("small_matrix_a.bin", &Ax.A,
                            &Ax.rows, &Ax.cols);
    /*mat = Ax->A;*/
    printf("%g\n", Ax.A[1][1]);
    printf("%d\n", Ax.cols);
    MPI_Finalize();
}
