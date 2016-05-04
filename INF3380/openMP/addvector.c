#include <omp.h>
#include <stdio.h>
#define CHUNKSIZE 4
#define N 10

void main () {
    int i, chunk, th_id;
    float a[N], b[N], c[N];

    /*Build arrays */
    for(i = 0; i < N; i++){
        a[i] = b[i] = i * 1.0;
    }
    chunk = CHUNKSIZE;

    #pragma omp parallel shared(a,b,c, chunk) private(i)
    {
        th_id = omp_get_thread_num();
        printf("Hello from thread %d\n", th_id);
        #pragma omp for schedule(dynamic, chunk) nowait
            for(i=0; i < N; i++)
                c[i] = a[i] + b[i];
            }//Ending parallel section of code
    int sum = 0;
    for(i = 0; i < N; i++){
        sum += c[i];
    }
    printf("Sum of array is %d \n", sum);
}
