#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int nargs, char **args) {

    int my_rank, size, i;
    int loc_sum, glob_sum;
    time_t t;
    int n = 2;
    int list[n];
    /* Intializes random number generator */
    srand((unsigned) time(&t));
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    for(i == 0; i < n; i++){
        int numb = rand() % 20;
        printf("numb %d \n", numb);
        loc_sum +=  numb;
    }

    MPI_Reduce(&loc_sum, &glob_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    printf("my rank %d \n", my_rank);

    if(my_rank == 0){
        printf("Global sum %d \n", glob_sum);
    }

    MPI_Finalize();
    return 0;
}
