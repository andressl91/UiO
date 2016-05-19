#include <stdio.h>

void write_binaryformat (char* filename, double** matrix,
int num_rows, int num_cols)
{
    FILE *fp = fopen (filename,"wb");
    fwrite (&num_rows, sizeof(int), 1, fp);
    fwrite (&num_cols, sizeof(int), 1, fp);
    
    fwrite (matrix[0], sizeof(double), num_rows*num_cols, fp);
    fclose (fp);
}
