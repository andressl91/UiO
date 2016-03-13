#include <stdio.h>
#include <mpi.h>
#include <string.h> /*Imports strlen which counts string */

int main(int nargs, char **args) {

    char my_line[40];
    char flag[20];
    int i, size, my_rank;
    MPI_Status status;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    sprintf(my_line, "Hello world form pricess %d \n", my_rank);
    int length = strlen(my_line) + 1;
    if (my_rank == 0) {
        printf("%s\n", my_line);
        for(i = 1; i < size; i++){
            MPI_Recv(my_line, length, MPI_CHAR, i,
            10, MPI_COMM_WORLD, &status);
            printf("%s\n", my_line);
        }
    }

    if (my_rank > 0) {
    MPI_Send(my_line, length, MPI_CHAR, 0,
    10, MPI_COMM_WORLD);
        }

    MPI_Finalize();
    return 0;

}
