/*Author: Andreas Slyngstad
Project : Serial program for smoothing noise on picture
read from command line: kappa, iters, input_jpeg_filename, output_jpeg_filename
*/

#include <stdio.h>
#include <stdlib.h>
/*#include "../simple-jpeg/import_export_jpeg.h"*/
// make use of two functions from the simplejpeg library

typedef struct {
float** image_data;  /* a 2D array of floats */
int m;             /* # pixels in x-direction */
int n;             /* # pixels in y-direction */
}
image;

void import_JPEG_file(const char *filename, unsigned char **image_chars,
int *image_height, int *image_width, int *num_components);

void export_JPEG_file(const char *filename, unsigned char *image_chars,
int image_height, int image_width, int num_components, int quality);

/*Function allocate image is supposed to allocate the 2D array image data inside u
when m and n are given as input. */
void allocate_image (image *u ,int j, int k) {
    int i;
    float **A = malloc(j*sizeof(float*));
    for (i = 0; i < k; i++) {
        A[i] = malloc(k*sizeof(float));
    }
    &u.j = 2;
    return;
}

int main(int argc, char *argv[])  {
    int m, n, c;
    float kappa = atoi(argv[1]);
    int iters = atoi(argv[2]);

    unsigned char *image_chars;
    char *input_jpeg_filename = argv[3];
    char *output_jpeg_filename = argv[4];

    image u, u_bar;


    import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
    allocate_image (&u, m, n);

    /*allocate_image (&u, m, n);*/
    /*
    allocate_image (&u_bar, m, n);
    convert_jpeg_to_image (image_chars, &u);
    iso_diffusion_denoising (&u, &u_bar, kappa, iters);
    convert_image_to_jpeg (&u_bar, image_chars);
    export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
    deallocate_image (&u);
    deallocate_image (&u_bar); */
    return 0;
}
