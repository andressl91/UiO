#include <stdio.h>
#include <stdlib.h>

typedef struct
{
  int** A; /* a 2D array of floats */
  int** x;
  int** C;
  int m;               /* # m-rows */
  int n;               /* # n-columns */

} linsys;

void read_matrix_binaryformat (char* filename, double*** matrix,
                                int* num_rows, int* num_cols)
{

    int i;
    FILE* fp = fopen (filename,"rb");
    fread (num_rows, sizeof(int), 1, fp);
    fread (num_cols, sizeof(int), 1, fp);

    /* storage allocation of the matrix */
    *matrix = (double**)malloc((*num_rows)*sizeof(double*));
    (*matrix)[0] = (double*)malloc((*num_rows)*(*num_cols)*sizeof(double));
    for (i=1; i<(*num_rows); i++)
    (*matrix)[i] = (*matrix)[i-1]+(*num_cols);

    /* read in the entire matrix */
    fread ((*matrix)[0], sizeof(double), (*num_rows)*(*num_cols), fp);
    fclose (fp);
}

void getarray(linsys * Ax){
    int i, j = 3;
    Ax->A = malloc(j*sizeof(float*)); /*Y-LENGTH ROWS*/
    for (i = 0; i < j; i++) {
        Ax->A[i] = malloc(j*sizeof(float)); /*X-LENGTH COLUMNS*/
    }
    Ax->A[1][1] = 3;
}

void main(){
    linsys Ax;
    int ** mat;
    getarray(&Ax);
    /*mat = Ax->A;*/
    printf("%d\n", Ax.A[1][1]);
}
