#include <stdio.h>
#include <mpi.h>

int main(int nargs, char **args) {
    int size, my_rank;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int i;
    for (i = 0; i < size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (i==my_rank){
        printf("Hello World! I am rank %d out of %d processes \n", my_rank, size);
        fflush(stdout); /*Only affects outgoing output stream */
        }
    }
    MPI_Finalize();
    return 0;
}
