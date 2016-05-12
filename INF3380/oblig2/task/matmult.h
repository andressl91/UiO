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

double ** allocate_matrix(int rows, int cols);
double ** get_my_share(matrix * A, double ** my_matrix, int start, int row, int cols);
void deallocate_matrix (matrix * mat);


#endif


/*
void read_matrix_binaryformat (char* filename, double*** matrix,
                                int* num_rows, int* num_cols);
*/
