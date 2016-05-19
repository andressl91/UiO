#ifndef GRANDFATHER_H
#define GRANDFATHER_H
#include <mpi.h>
typedef struct
{
  double** mat; /* a 2D array of floats */
  double** matB; /* a 2D array of floats */
  double** wholeB;
  double** matC;
} matrix;

void read_matrix_binaryformat (char* filename, double *** matrix,
                                int* num_rows, int* num_cols);

double ** allocate_matrix(int rows, int cols);
double ** get_my_share(matrix * A, double ** my_matrix, int start, int row, int cols);
void deallocate_matrix (double **mat, int m, int n);
void find_sum(double ** C, double **A, double **B, int my_rows, int num_rows, int num_cols);
void allocate_array(double ***array, int m, int n);
void create_matrix_type(MPI_Datatype *resizedmatrix, int m, int n, int sub_m, int sub_n);
void write_binaryformat (char* filename, double** matrix,
int num_rows, int num_cols);
#endif


/*
void read_matrix_binaryformat (char* filename, double*** matrix,
                                int* num_rows, int* num_cols);
*/
