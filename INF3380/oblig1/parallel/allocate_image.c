#include <stdio.h>
#include <stdlib.h>
#include "parallel_main.h"
#include <mpi.h>

void allocate_image (image *u ,int j, int k) {
    int i;
    u->image_data = (float**)malloc(j*sizeof(float*)); /*Y-LENGTH ROWS*/
    for (i = 0; i < j; i++) {
        u->image_data[i] = (float*)malloc(k*sizeof(float)); /*X-LENGTH COLUMNS*/
    }
    u->m = j; u->n = k;
    return;

}
