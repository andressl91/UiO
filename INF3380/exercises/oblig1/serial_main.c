#include "../simple-jpeg/import_export_jpeg.h"
#include <stdio.h>

#Brukte cat commando på fila .h for a se arugumenter!
#FIkk



#exetern boind import_JPEG_file ( const char* filename,
unsignet char** image_chars
...........);

#This is made in the serial_main.c file

typedef struct{
  float** image;
  int m, n;
} image;

void allocate_image(image* u, int m, int n){
  int i;
  u -> m = m; #bilde til m blir det du tar inn som arg
  u -> n = n;
  u -> image_data = (float**)malloc(m*sizeof(float*));
  for(i = 0; i < m; i++){
    u->image_data[i] = (float*)malloc(n*sizeof(float*)); #Empty, just allocated
  }
}

int main () {
  int image_height, image_width, num_components;
  unsigned char *image_chars; #bilde i 1d
  image u;


  import_JPEG_file("index.jpeg", &image_chars, &image_height,
  &image_width, &num_components);
  allocate_image(&u, image_height, image_width);


  printf("W: %d, H: %d, C: %d, name: %s",image_width, image_height, num_components
  , "index.jpeg")

}

#terminal wrote make when in right folder, which was serial
#./serial.main

write out file with the corresponding
#export_JPEG_file....

#tips, get *image_chars to 2 array
int i;
for (i = 0; i < 100, i++){
  printf("%c", image_chars[i]);
  return 0;
}
