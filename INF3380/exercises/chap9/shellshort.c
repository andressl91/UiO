#include <stdio.h>
#include <stdlib.h>
#define MAX 10
int* makearray(){
    int i, n = 10;
    int *A = malloc(n*sizeof(int));
    for(i = 0; i < n; i++){
        A[i] = rand() % 20;
    }
    printf("ARRAY TO BE SORTED \n");
    for(i = 0; i < n; i++){
        printf("%d, ", A[i]);
    }
    printf("\n");
    return A;
}

void shellsort(){
    int inner, outer;
    int valueToInsert;
    int interval = 1;
    int count = 0;
    int elements = MAX;

    while(interval <= elements/3){
        interval = interval*3 + 1;
        printf("%s\n", );
    }

    /*while(interval > 0){
        printf("This is iteration number %d \n", count);
*/




    }


}

void main() {
    int *B = makearray();
    shellsort();
}
