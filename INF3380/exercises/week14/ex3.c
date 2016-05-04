#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "time.h"

void array(int i, int *a, int *b, int *c){
    int j;
    for(j = 0; j < i; j++){
        a[j] = b[j] = c[j] = 1;
    }
}

void main(int argc, char *argv[]) {

    clock_t begin, end;
    double time_spent;
    begin = clock();
    /* here, do your time-consuming job */

    int j, i = 1000;
    int chunk = 100;
    int *a = malloc(i*sizeof(int));
    int *b = malloc(i*sizeof(int));
    int *c = malloc(i*sizeof(int));
    array(i, a, b, c);

    #pragma omp parallel shared(a, b, c, chunk) private(i) num_threads(5)
    {
        #pragma omp for schedule(dynamic, chunk) nowait
            for(j = 0; j < i; j++) {
                c[j] = a[j]*b[j];
            }
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%f \n", time_spent);

}
