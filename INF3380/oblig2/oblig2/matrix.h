void read_matrix_binaryformat (char* filename, double*** matrix, int* num_rows, int* num_cols);
void write_matrix_binaryformat (char* filename, double** matrix, int num_rows, int num_cols);
void block_mat_mult(int m, int n, double **A, double **B, double **C);
void allocate_array(double ***A, int m, int n);
void create_matrix_type(MPI_Datatype *resizedmatrix, int m, int n, int sub_m, int sub_n);
void cannon_nonblocking(int m, int n, double *A, double *B, double *C, MPI_Comm comm);
void matrix_multiply(int mlocal, int nlocal, double *A, double *B, double *C);
