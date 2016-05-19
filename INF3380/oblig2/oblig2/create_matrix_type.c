#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void create_matrix_type(MPI_Datatype *resizedmatrix, int m, int n, int sub_m, int sub_n) {
    MPI_Datatype matrix;
    int sizes[2] = {m, n};
    int subsizes[2] = {sub_m, sub_n};
    int starts[2] = {0, 0};

    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix);
    MPI_Type_create_resized(matrix, 0, sub_n*sizeof(double), resizedmatrix);
    MPI_Type_commit(resizedmatrix);

}
