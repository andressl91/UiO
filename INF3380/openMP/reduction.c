#include <omp.h>
#include <stdio.h>

void main(int argc, char *argv[]) {
    int th_id, nthreads;
    int sum = 0;

    #pragma omp parallel reduction(+: sum) num_threads(4)  private(th_id)
    {
        th_id = omp_get_thread_num();
        sum = 2;
        printf("The sum is %d \n", sum);

    #pragma omp barrier


    if (th_id == 0) {
        nthreads = omp_get_num_threads();
        printf("There are %d threads \n", nthreads);

        }
    }
    printf("The sum is %d \n", sum);
}
