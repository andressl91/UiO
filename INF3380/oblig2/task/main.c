#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include "matmult.h"

void main(int argc, char *argv[]){
    double **mat_A, **mat_B, **my_mat_A, **my_mat_B, **mat_C, **my_mat_C;
    int A_m, A_n, B_m, B_n;
    int my_Am, my_An, my_Bm, my_Bn;
    int dimentions[4];

    int i, j, l;

    int my_rank, num_procs;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


    double *send_mat_A = NULL, *send_mat_B = NULL, *recv_mat = NULL;
    if(my_rank == 0){
        read_matrix_binaryformat("small_matrix_a.bin", &mat_A, &A_m, &A_n);
        read_matrix_binaryformat("small_matrix_b.bin", &mat_B, &B_m, &B_n);

        allocate_array(&mat_C, A_m, B_n);
        send_mat_A = &(mat_A[0][0]);
        send_mat_B = &(mat_B[0][0]);
        recv_mat = &(mat_C[0][0]);

        dimentions[0] = A_m;
        dimentions[1] = A_n;
        dimentions[2] = B_m;
        dimentions[3] = B_n;
    }

    MPI_Bcast(&dimentions, 4, MPI_INT, 0, MPI_COMM_WORLD);
    A_m = dimentions[0];
    A_n = dimentions[1];
    B_m = dimentions[2];
    B_n = dimentions[3];

    my_Am = A_m / num_procs;
    if ( my_rank < A_m % num_procs) {
        ++my_Am;
    }

    my_Bm = B_m / num_procs;
    if ( my_rank <B_m % num_procs) {
        ++my_Bm;
    }
    my_An = A_n;
    my_Bn = B_n;

    /* Allocates memory for each processes array
    * and each calculated part of the result matrix
    */
    allocate_array(&my_mat_A, my_Am, my_An);
    allocate_array(&my_mat_B, my_Bm, my_Bn);
    allocate_array(&my_mat_C, my_Am, my_Bn);

    int my_len_A = my_Am*my_An;
    int *sendcount_A = NULL;
    int my_len_B = my_Bm*my_Bn;
    int *sendcount_B = NULL;
    int my_len_C = my_Am*B_n;
    int *sendcount_C = NULL;
    if (my_rank == 0){
       sendcount_A = malloc( num_procs * sizeof(int)) ;
       sendcount_B = malloc( num_procs * sizeof(int)) ;
       sendcount_C = malloc( num_procs * sizeof(int)) ;
   }
     MPI_Gather(&my_len_A, 1, MPI_INT, sendcount_A, 1, MPI_INT,
          0, MPI_COMM_WORLD);
     MPI_Gather(&my_len_B, 1, MPI_INT, sendcount_B, 1, MPI_INT,
               0, MPI_COMM_WORLD);
     MPI_Gather(&my_len_C, 1, MPI_INT, sendcount_C, 1, MPI_INT,
                 0, MPI_COMM_WORLD);
        /*
        * Calculate total length of matrix A,B,C
        * and displacements for each ranks matrix
        */

        int totlen_A = 0;
        int *displs_A = NULL;

        int totlen_B = 0;
        int *displs_B = NULL;

        int totlen_C = 0;
        int *displs_C = NULL;

       if (my_rank == 0) {

           displs_A = malloc( num_procs * sizeof(int) );
           displs_A[0] = 0;
           totlen_A += sendcount_A[0];//+1;

           displs_B = malloc( num_procs * sizeof(int) );
           displs_B[0] = 0;
           totlen_B += sendcount_B[0];//+1;

           displs_C = malloc( num_procs * sizeof(int) );
           displs_C[0] = 0;
           totlen_C += sendcount_C[0];//+1;

           for (i=1; i<num_procs; i++) {
              totlen_A += sendcount_A[i];//+1;
              displs_A[i] = displs_A[i-1] + sendcount_A[i-1]+ 1;

              totlen_B += sendcount_B[i];//+1;
              displs_B[i] = displs_B[i-1] + sendcount_B[i-1]+ 1;

              totlen_C += sendcount_C[i];//+1;
              displs_C[i] = displs_C[i-1] + sendcount_C[i-1]+ 1;
           }
        }

        MPI_Scatterv(send_mat_A, sendcount_A, displs_A, MPI_DOUBLE,
            &(my_mat_A[0][0]), my_Am*my_An,
                 MPI_DOUBLE, 0, MPI_COMM_WORLD);

        ///MPI_Scatterv(send_mat_B, sendcount_B, displs_B, MPI_DOUBLE,
        ///     &(my_mat_B[0][0]), my_Bm*my_Bn,
        //          MPI_DOUBLE, 0, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);

        if(my_rank != 0){
            allocate_array(&mat_B, B_m, B_n);
        }

        MPI_Bcast(&mat_B[0][0], B_m*B_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        /*
        if(my_rank == 0){
                printf("MY RANK %d BEFIRE SUM\n", my_rank);
            for(i = 0; i < my_Am; i++){
                printf("%f ", my_mat_A[i][0]);
                for(j = 1; j < my_An; j++){
                    printf("%f ",my_mat_A[i][j]);
                }
                printf("\n");
            }printf("END LOOp\n" );
        } */

        //FIX POINTER ARGUMENT TO
        MPI_Barrier(MPI_COMM_WORLD);
        find_sum(my_mat_C, my_mat_A, mat_B, my_Am, B_m, B_n);


            MPI_Barrier(MPI_COMM_WORLD);
            /*
        if(my_rank == 0){
                printf("MY RANK %d AFTER SUM\n", my_rank);
            for(i = 0; i < my_Am; i++){
                //printf("%f ",my_mat_C[i][0]);
                for(j = 1; j < B_n; j++){
                    //printf("%f ",my_mat_C[i][j]);
                }
                printf("\n");
            }printf("END LOOp\n" );
        } */

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gatherv(&(my_mat_C[0][0]), my_Am*B_n, MPI_DOUBLE,
        recv_mat, sendcount_C, displs_C, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if(my_rank == 0){
        double * exact = NULL;
        double **C_ex;
        double sum;
        write_binaryformat ("small_matrix_result.bin", &recv_mat, A_m, B_n);
        //read_matrix_binaryformat("small_matrix_result.bin", &C_ex, &A_m, &B_n);
        read_matrix_binaryformat("small_matrix_c.bin", &C_ex, &A_m, &B_n);

        for(i = 0; i < my_Am; i++){
            for(j = 0; j < B_n; j++){
                //printf("Exact %f\n", C_ex[i][j]);//-C_ex[i][j]);
                if(my_mat_C[i][j]- C_ex[i][j] !=0)
                printf("FAIL AT LIFE\n");
                //printf("My calc %f\n", my_mat_C[i][j]- C_ex[i][j]);
            }
        }
                    printf("THREAD %d has my_Am as %d\n", my_rank, my_Am);

    }
        //double **mat_A, **mat_B, **my_mat_A, **my_mat_B, **mat_C, **my_mat_C;

        //if(my_rank == 0)
        //deallocate_matrix(mat_A, A_m, A_n);

        if(my_rank == 0){
        free(mat_B); free(mat_A);
            //free(mat_C);
        }
        //f
        free(my_mat_A); free(my_mat_B); free(my_mat_C);
    MPI_Finalize();
}
