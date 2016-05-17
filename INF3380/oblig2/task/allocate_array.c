#include <stdio.h>
#include <stdlib.h>

void allocate_array(double ***array, int m, int n) {
    int i;

    /* allocate the n*m contiguous items */
    double *p = (double *)malloc(n*m*sizeof(double));

    /* allocate the row pointers into the memory */
    (*array) = (double **)malloc(m*sizeof(double*));
    if (!(*array)) {
       free(p);
       printf("Something went wrong mallocing\n");
       fflush(stdout);
    }
    fflush(stdout);

    /* set up the pointers into the contiguous memory */
    for ( i=0; i<m; i++) {
       (*array)[i] = &(p[i*n]);
    }
}
