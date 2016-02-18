#include "mkarray.h"
#include <stdio.h>
#include <stdlib.h> /*For malloc function */

int *** array3d(int i, int j, int k) {
    /*Unitcube(0,0,0) as bottom corner */
    int m = 10;
    int del_h = 1/(float)m;
    int start = 0, stop = 1;

    int ***u = malloc(m*sizeof(int**));
    int ***u_prev = malloc(m*sizeof(int**));

   for (i = 0; i < m; i++) {
     u[i] = malloc(m*sizeof(int*));
     for (j = 0; j < m; j++) {
      u[i][j] = malloc(m*sizeof(int));
     }
   }
   for (i = 0; i < m; i++) {
     u_prev[i] = malloc(m*sizeof(int*));
     for (j = 0; j < m; j++) {
      u_prev[i][j] = malloc(m*sizeof(int));
     }
   }
/*
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < m; k++)
        printf("%d \n", i*m);
        */
    return u;
}
