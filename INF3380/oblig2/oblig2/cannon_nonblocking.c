#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "matrix.h"

#define TRUE 1
#define FALSE 0

void cannon_nonblocking(int m, int n, double *A, double *B, double *C, MPI_Comm comm) {
    int i, j, mlocal, nlocal;
    double *A_buffers[2], *B_buffers[2];
    int size, dims[2], periods[2];
    int rank, rank_2d, mycoords[2];
    int uprank, downrank, leftrank, rightrank, coords[2];
    int shiftsource, shiftdest;
    MPI_Status status;
    MPI_Comm comm_2d;
    MPI_Request reqs[4];

    /* Get communicator related information */
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    /* Set up the Cartesian topology */
    dims[0] = dims[1] = sqrt(size);

    /* Set hte periods for wraparound connections */
    periods[0] = periods[1] = 1;
    
    /* Create the Cartesian topology with rank reordering */
    MPI_Cart_create(comm, 2, dims, periods, TRUE, &comm_2d);

    /* Get the rank and coordinates in the new topology */
    MPI_Comm_rank(comm_2d, &rank_2d);
    MPI_Cart_coords(comm_2d, rank_2d, 2, mycoords);

    /* Compute rank of up and left shifts */
    MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
    MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

    /* Determine dimensions of local matrix block */
    mlocal = m/dims[0];
    nlocal = n/dims[1];

    /* Set up the A_buffers and B_buffers arrays */
    A_buffers[0] = A;
    A_buffers[1] = (double *)malloc(mlocal*nlocal*sizeof(double));
    B_buffers[0] = B;
    B_buffers[1] = (double *)malloc(mlocal*nlocal*sizeof(double));

    /* Perform the initial matrix alignment. First for A, then for B */
    MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(A_buffers[0], mlocal*mlocal, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, \
            comm_2d, &status);

    MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(B_buffers[0], nlocal*mlocal, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, \
            comm_2d, &status);


    /* Get into the main computation loop */
    for ( i = 0; i < dims[0]; ++i ) {
        MPI_Isend(A_buffers[i % 2], mlocal*nlocal, MPI_DOUBLE, leftrank, 1, comm_2d, &reqs[0]);
        MPI_Isend(B_buffers[i % 2], mlocal*nlocal, MPI_DOUBLE, uprank, 1, comm_2d, &reqs[1]);
        MPI_Irecv(A_buffers[(i + 1) % 2], mlocal*nlocal, MPI_DOUBLE, rightrank, \
                1, comm_2d, &reqs[2]);
        MPI_Irecv(B_buffers[(i + 1) % 2], mlocal*nlocal, MPI_DOUBLE, downrank, \
                1, comm_2d, &reqs[3]);

        matrix_multiply(mlocal, nlocal, A_buffers[i % 2], B_buffers[i % 2], C);

        for ( j = 0; j < 4; ++j ) {
            MPI_Wait(&reqs[j], &status);
        }
    }

    /* Restore the original distribution of A and B */
    MPI_Cart_shift(comm_2d, 1, +mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(A_buffers[i % 2], mlocal*nlocal, MPI_DOUBLE, shiftdest, 1, shiftsource, \
            1, comm_2d, &status);

    MPI_Cart_shift(comm_2d, 0, +mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(B_buffers[i % 2], mlocal*nlocal, MPI_DOUBLE, shiftdest, 1, shiftsource, \
            1, comm_2d, &status);

    MPI_Comm_free(&comm_2d); /* Free up the ommunicator */
    
    free(A_buffers[1]);
    free(B_buffers[1]);
}
