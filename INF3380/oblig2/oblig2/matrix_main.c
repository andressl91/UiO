#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "matrix.h"

int main(int argc, char *argv[]) {
    double **matrix_a, **matrix_b, **my_matrix_a, **my_matrix_b, **matrix_out, **my_matrix_out; 
    int A_m, A_n, B_m, B_n;
    int rank, size;
    int my_Am, my_An, my_Bm, my_Bn;
    int i, j, k, disp_Am, disp_Bm, disp;
    int numbers[4];  //holds dimensions of matrices A and B, i.e. 4 ints
    MPI_Datatype MPImatrix_a, MPImatrix_b, MPImatrix_out;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int sqrtsize = (int)sqrt(size);

    /* Read matrices, only Master does this */
    double *sendMatrix_a = NULL, *sendMatrix_b = NULL, *recvMatrix = NULL;
    if (rank == 0) {
        read_matrix_binaryformat ("small_matrix_a.bin", &matrix_a, &A_m, &A_n);
        read_matrix_binaryformat ("small_matrix_b.bin", &matrix_b, &B_m, &B_n);
        /* NB! */
        allocate_array(&matrix_out, A_m, B_n); //Fix Dimensions 

        sendMatrix_a = &(matrix_a[0][0]);
        sendMatrix_b = &(matrix_b[0][0]);
        recvMatrix = &(matrix_out[0][0]);

        numbers[0] = A_m;
        numbers[1] = A_n;
        numbers[2] = B_m;
        numbers[3] = B_n;
    }

    // Broadcast matrix dimensions
    MPI_Bcast(&numbers, 4, MPI_INT, 0, MPI_COMM_WORLD); // 0 = master is root

    A_m = numbers[0];
    A_n = numbers[1];
    B_m = numbers[2];
    B_n = numbers[3];

    my_Am = A_m / sqrtsize;
    my_An = A_n / sqrtsize;
    my_Bm = B_m / sqrtsize;
    my_Bn = B_n / sqrtsize;

    if ( rank < A_m % sqrtsize ) {
        ++my_Am;
    }
    if ( rank < A_n % sqrtsize ) {
        ++my_An;
    }
    if ( rank < B_n % sqrtsize ) {
        ++my_Bn;
    }
    if ( rank < B_n % sqrtsize ) {
        ++my_Bm;
    }

    allocate_array(&my_matrix_a, my_Am, my_An);
    allocate_array(&my_matrix_b, my_Bm, my_Bn);
    allocate_array(&my_matrix_out, my_Am, my_Bn);

    /* MPI_Datatype subarray used for purposes of scatterv */
    create_matrix_type(&MPImatrix_a, A_m, A_n, my_Am, my_An);
    create_matrix_type(&MPImatrix_b, B_m, B_n, my_Bn, my_Bm);
    create_matrix_type(&MPImatrix_out, A_m, B_n, my_Am, my_Bn);

    int displs_A[size], displs_B[size], counts[size];
    if ( rank == 0 ) {
        /* How many pieces of data everyone has, in units of blocks */
        for ( i = 0; i < size; ++i ) {
            counts[i] = 1;
        }

        /* Starting point of everyone's data in matrix_a in block size */
        disp = 0;
        for ( i = 0; i < sqrtsize; ++i) {
            disp_Am = A_m / sqrtsize; 
            if ( i < A_m % sqrtsize ) {
                disp_Am += 1;
            }
            for ( j = 0; j < sqrtsize; ++j) {
                displs_A[i*sqrtsize + j] = disp + j;
            }
            disp += disp_Am*sqrtsize;
        }

        /* Then the same for B */
        disp = 0;
        for ( i = 0; i < sqrtsize; ++i ) {
            disp_Bm = B_m % sqrtsize;            
            if ( i < B_m % sqrtsize ) {
                disp_Bm += 1;
            }
            for ( j = 0; j < sqrtsize; ++j ) {
                displs_B[i*sqrtsize + j] = disp + j;
            }
            disp += disp_Bm*sqrtsize;
        }
    }

    MPI_Scatterv(sendMatrix_a, counts, displs_A, MPImatrix_a, &(my_matrix_a[0][0]), my_Am*my_An,\
             MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(sendMatrix_b, counts, displs_B, MPImatrix_b, &(my_matrix_b[0][0]), my_Bm*my_Bn,\
             MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    cannon_nonblocking(my_Am, my_An, &(my_matrix_a[0][0]), &(my_matrix_b[0][0]),\
            &(my_matrix_out[0][0]), MPI_COMM_WORLD); 

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gatherv(&(my_matrix_a[0][0]), my_Am*my_Bn, MPI_DOUBLE, recvMatrix, counts, displs_A, \
            MPImatrix_out, 0, MPI_COMM_WORLD);
    
    /* Free MPI_Datatype */
    MPI_Type_free(&MPImatrix_a);
    MPI_Type_free(&MPImatrix_b);

    if ( rank == 0 ) {
        write_matrix_binaryformat ("small_matrix_out.bin", matrix_out, A_m, B_n);
    }

    MPI_Finalize();
    return 0;
}
