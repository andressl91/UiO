#ifndef GRANDFATHER_H
#define GRANDFATHER_H

typedef struct
{
  double** mat; /* a 2D array of floats */
  int rows;               /* # m-rows */
  int cols;               /* # n-columns */

} matrix;

void read_matrix_binaryformat (char* filename, double *** matrix,
                                int* num_rows, int* num_cols);

void get_my_share(matrix * A, matrix * my_mat, int start);
void deallocate_matrix (matrix * mat);

#endif


/*
void read_matrix_binaryformat (char* filename, double*** matrix,
                                int* num_rows, int* num_cols);
*/
