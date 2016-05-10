#ifndef GRANDFATHER_H
#define GRANDFATHER_H

typedef struct
{
  double** A; /* a 2D array of floats */
  double** x;
  double** C;
  int rows;               /* # m-rows */
  int cols;               /* # n-columns */

} matsys;

void read_matrix_binaryformat (char* filename, double*** matrix,
                                int* num_rows, int* num_cols);

#endif
