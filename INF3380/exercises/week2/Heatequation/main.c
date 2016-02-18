#include "mkarray.h"
#include <stdio.h>
#include <stdlib.h> /*For malloc function */


int getarray() {
  int i = 10;
  int *** ptr = array3d(i, i ,i);
  ptr[1][1][1] = 1;
  printf("%d \n", ptr[1][1][1]);
  printf("This is working \n");
  return 0;
}

int main(){
  getarray();
  return 0;

}
