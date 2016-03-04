#include <stdio.h>
#include <mpi.h>

int main (int nargs, char **args) {
    int size, my_rank;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    printf("Hello World! I am rank %d out of %d processes", my_rank, size);
    MPI_Finalize();
    return 0;
}
