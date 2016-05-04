/*Author: Andreas Slyngstad
Project : Serial program for smoothing noise on picture
read from command line: kappa, iters, input_jpeg_filename, output_jpeg_filename
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*#include "../simple-jpeg/import_export_jpeg.h"*/
// make use of two functions from the simplejpeg library



typedef struct
{
  float** image_data;  /* a 2D array of floats */
  int m;               /* # pixels in height */
  int n;               /* # pixels in width */

} image;


void import_JPEG_file(const char *filename, unsigned char **image_chars,
int *image_height, int *image_width, int *num_components);

void export_JPEG_file(const char *filename, unsigned char *image_chars,
int image_height, int image_width, int num_components, int quality);

/*Function allocate image allocates memory the 2D array image data inside u
when m and n are given as input. */
void allocate_image (image *u ,int j, int k) {
    int i;
    u->image_data = (float**)malloc(j*sizeof(float*)); /*Y-LENGTH ROWS*/
    for (i = 0; i < j; i++) {
        u->image_data[i] = (float*)malloc(k*sizeof(float)); /*X-LENGTH COLUMNS*/
    }

    u->m = j; u->n = k;
    return;
    }

/*Function convert jpeg to image is supposed to convert a 1D array of
unsigned char values into an image struct.*/
void convert_jpeg_to_image (unsigned char *image_chars, image *u) {

    int i, j, count = 0;
    for (i = 0; i < u->m; i++){ /*ROWS */
        for(j = 0; j < u->n; j++){ /*COLUMS */
            u->image_data[i][j] = (float)image_chars[count];
            count += 1;
            }
        }
    }

void iso_diffusion_denoising (image *u, image *u_bar, float kappa, int iters) {
    int i, j, it;
    float ** u_p, **copy, u_upt ;
    for(it = 0; it < iters; it++){
        u_p = u->image_data;

        for (i = 1; i < u->m - 1; i++){
            for(j = 1; j < u->n - 1; j++){
                u_upt = u_p[i][j];
                u_upt += kappa * (u_p[i-1][j] + u_p[i][j-1] - 4* u_upt + u_p[i+1][j] + u_p[i][j+1]);
                u_bar->image_data[i][j] = u_upt;
            }
        }
        copy = u->image_data;
        u->image_data = u_bar->image_data;
        u_bar->image_data = copy;
        }
    }

void convert_image_to_jpeg (image *u, unsigned char *image_chars){
    int i, j;
    int count = 0;
    for (i = 0; i < u->m; i++){
        for(j = 0; j < u->n; j++){
            image_chars[count] = u->image_data[i][j];
            count+= 1;
            }
        }
    }

void deallocate_image (image *u) {
    int i;
    for (i = 0; i < u->m; i++){
        free(u->image_data[i]);
        }
    free(u->image_data);
    }

int main(int argc, char *argv[])  {
    double kappa = atof(argv[2]);
    int m, n, c, iters = atoi(argv[1]);

    unsigned char *image_chars;
    char *input_jpeg_filename = argv[3], *output_jpeg_filename = argv[4];
    image u, u_bar;
    import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
    allocate_image (&u, m, n);
    allocate_image (&u_bar, m, n);

    /*Value 1 should also be given to num components before invoking export
    JPEG file to export a grey-scale JPEG image*/
    convert_jpeg_to_image (image_chars, &u);

    iso_diffusion_denoising (&u, &u_bar, kappa, iters);

    convert_image_to_jpeg (&u_bar, image_chars);
    export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 100);
    deallocate_image (&u);
    deallocate_image (&u_bar);
    return 0;
}
