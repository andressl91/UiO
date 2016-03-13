#include <stdio.h>
#include <mpi.h>

int main(int nargs, char **args) {


   int M = 10;
   int u[M], v[M]; /* n is an array of 10 integers */
   int i,j;

   /* initialize elements of array n to 0 */
   for ( i = 0; i < 10; i++ ) {
      u[i] = i + 10; /* set element at location i to i + 100 */
      v[i] = i;
   }
    int size, my_rank, count;
    int P = 2; /*Number of segments */
    int m = M/P;
    double t1 = MPI_Wtime();
    MPI_Status status;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int my_start = M*my_rank/size;
    int my_stop = M*(my_rank + 1)/size;
    int my_c = 0.;

    if (my_rank > 0) {
        MPI_Recv(&count, 1, MPI_INT, my_rank-1,
        10, MPI_COMM_WORLD, &status);
        my_c = count;
    }

    if (my_rank < size-1){
        for (i=my_start; i < my_stop; i++) {
            my_c = my_c + (u[i] + v[i]);
        MPI_Send(&my_c, 1, MPI_INT, my_rank+1,
        10, MPI_COMM_WORLD);
        }
    }

    if (my_rank == size-1){
        printf("The sum is %d \n", my_c);
        double t2 = MPI_Wtime();
        printf( "Elapsed time is %f\n", t2 - t1 );
    }

    MPI_Finalize();

    return 0;
}
