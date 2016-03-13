#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int nargs, char **args) {

    int size, my_rank, count;
    double t1 = MPI_Wtime();
    int M = 3;
    int i, m;
    double *x, *y, *local_x, *local_y;
    double c;
    double my_c = 0.0;
    MPI_Status status;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0){
        x = (double *) calloc(M, sizeof(double));
        y = (double *) calloc(M, sizeof(double));
        for(i = 0; i < M; i++){x[i] = 1; y[i] = 1;}
    }
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /*Broad cast from ROOT to all processes in communicator */
    m = M/size; /*Size of segment length of array to sum */
    int my_start = M*my_rank/size;
    int my_stop = M*(my_rank + 1)/size;

    local_x = (double *)calloc(m, sizeof(double));
    local_y = (double *)calloc(m, sizeof(double));

    MPI_Scatter(x, m, MPI_DOUBLE, local_x, m, MPI_DOUBLE, 0,
    MPI_COMM_WORLD);

    MPI_Scatter(y, m, MPI_DOUBLE, local_y, m, MPI_DOUBLE, 0,
    MPI_COMM_WORLD);

    for ( i = my_start; i < my_stop; i++ ) {
        my_c = my_c + (local_x[i] *local_y[i]);
        }
    free(local_x); free(local_y);

    MPI_Allreduce(&my_c, &c, 1, MPI_DOUBLE,
    MPI_SUM, MPI_COMM_WORLD);

    if(my_rank == 0){
        printf("Sum of array x, y is %f \n", c);
    }

    MPI_Finalize();
    return 0;
}
