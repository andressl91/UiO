#include <stdio.h>
#include <stdlib.h>  /*For the malloc function */
/*
int array(double ** ptr, int m) {
  int i;
  double *A = (double*)malloc(m*sizeof(double*));
  ptr = &A;
  return 1;
}
*/
int main() {
    //double * ar;
    int m = 6;
    int i;
    //array(&ar, m);

    double *A = (double*)malloc(m*sizeof(double*));
      for (i = 0; i < m; i++)
        A[i] = (double*)malloc(n*sizeof(double*));

    double **A = (double**)malloc(m
    *sizeof(double
    *));
    for (i=0; i<m; i++)
    A[i] = (double
    *)malloc(n
    *sizeof(double));

    for (i = 0; i < m; i++) {
      A[i] = i;
      printf("%f \n", A[i]);
      }
    return 0;
}
