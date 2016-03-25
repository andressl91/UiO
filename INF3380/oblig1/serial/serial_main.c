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

/*Function allocate image allocates memory the 2D array image data inside u
when m and n are given as input. */
void allocate_image (image *u ,int j, int k) {
    int i;

    float **A = malloc(j*sizeof(float*)); /*X-LENGTH ROWS*/
    for (i = 0; i < k; i++) {
        A[i] = malloc(k*sizeof(float)); /*Y-LENGTH COLUMNS */
    }
    u->image_data = A;
    u->m = k; u->n = j;
    return;
    }

/*Function convert jpeg to image is supposed to convert a 1D array of
unsigned char values into an image struct.*/
void convert_jpeg_to_image (unsigned char * chars, image *u) {
    int i, j;
    for (i = 0; i < u->n; i++){ /*ROWS */
        for(j = 0; j < u->m; j++){ /*COLUMS */
            u->image_data[i][j] = chars[i+j];
        }
    }
    return;
    }

void iso_diffusion_denoising (image *u, image *u_bar, float kappa, int iters) {
    int i, j;
    /*u_bar->image_data = u->image_data;*/
    for (i = 1; i < u->n - 1; i++){
        for(j = 1; j < u->m - 1; j++){
            u_bar->image_data[i][j] = kappa * (u->image_data[i-1][j] + u->image_data[i][j-1]
                - 4*u->image_data[i][j] + u->image_data[i+1][j] + u->image_data[i][j+1] );
        }
    }
    return;
    }

void convert_image_to_jpeg (unsigned char ** chars, image *u_bar){
    int i, j;
    for (i = 0; i < u_bar->n; i++){
        for(j = 0; j < u_bar->m; j++){
            (*chars)[i+j] = u_bar->image_data[i][j];
        }
    }
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
    import_JPEG_file(input_jpeg_filename, &image_chars, &n, &m, &c);
    allocate_image (&u, m, n);
    allocate_image (&u_bar, m, n);

    /*Value 1 should also be given to num components before invoking export
    JPEG file to export a grey-scale JPEG image*/
    /*
    convert_jpeg_to_image (image_chars, &u);

    printf("%f \n", u_bar.image_data[20][10]);
    iso_diffusion_denoising (&u, &u_bar, kappa, iters);
    printf("%f \n", u_bar.image_data[20][10]);

    convert_image_to_jpeg (&image_chars, &u_bar);*/
    /*export_JPEG_file(output_jpeg_filename, image_chars, n, m, c, 75);*/
    /*deallocate_image (&u);
    deallocate_image (&u_bar); */
    return 0;
}
