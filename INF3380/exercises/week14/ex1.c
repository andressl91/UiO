#include <omp.h>
#include <stdio.h>
#include <math.h>

void main(int argc, char *argv[]) {
    int i;
    int x = sqrt(16);
    int a[5], b[5];
    for (i=0; i < x; i++) {
        a[i] = 2.3 * x;
        if (i < 10) b[i] = a[i];
        }
}
