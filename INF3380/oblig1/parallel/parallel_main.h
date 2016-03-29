#ifndef PARALLEL_MAIN
#define PARALLEL_MAIN

typedef struct
{
    float** image_data;
    int m;
    int n;
}
image;

// make use of two functions from the simplejpeg library
void import_JPEG_file(const char *filename, unsigned char **image_chars,
    int *image_height, int *image_width, int *num_components);
void export_JPEG_file(const char *filename, unsigned char *image_chars,
    int image_height, int image_width, int num_components, int quality);


void allocate_image (image *u ,int j, int k);
void convert_jpeg_to_image (unsigned char *image_chars, image *u);
void iso_diffusion_denoising (image *u, image *u_bar, float kappa, int iters);
void convert_image_to_jpeg (image *u, unsigned char *image_chars);
void deallocate_image (image *u);

#endif
