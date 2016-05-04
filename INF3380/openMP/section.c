#include <omp.h>
#include <stdio.h>

void main(int argc, char *argv[]) {

    int plus, minus;
    int a = 2;
    int b = 4;
    #pragma omp parallel num_threads(3)
    {
        #pragma omp sections
        {
            #pragma omp section
                {
                    plus = a + b;
                    printf("a + b = %d \n", plus);
                }
            #pragma omp section
                {
                    minus = b - a;
                    printf("b - a = %d \n", minus);
                }
        }

    }

}
