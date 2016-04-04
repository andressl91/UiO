#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef struct
{
    float** image_data;
    int m;
    int n;
}
image;

void import_JPEG_file(const char *filename, unsigned char **image_chars,
                    int *image_height, int *image_width,
                     int *num_components);

void export_JPEG_file(const char *filename, unsigned char *image_chars,
                    int image_height, int image_width,
                    int num_components, int quality);

void allocate_image (image *u ,int j, int k) {
    int i;
    u->image_data = (float**)malloc(j*sizeof(float*)); /*Y-LENGTH ROWS*/
    for (i = 0; i < j; i++) {
        u->image_data[i] = (float*)malloc(k*sizeof(float)); /*X-LENGTH COLUMNS*/
    }
    u->m = j; u->n = k;
}

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
    int m = u->m, n = u->n;
    float ** u_p, **copy, u_upt ;
    int my_rank, num_procs;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


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


        /*Boundary operations for the different proesses */
        if ( my_rank == 0 ) {
            MPI_Request reqs[2];
            MPI_Isend(u->image_data[m - 2], n, MPI_FLOAT, my_rank + 1, 0, MPI_COMM_WORLD, reqs);
            MPI_Irecv(u->image_data[m-1], n, MPI_FLOAT, my_rank + 1, 0, MPI_COMM_WORLD, reqs + 1);
            MPI_Waitall(2, reqs, MPI_STATUS_IGNORE);
        }
        else if  ( my_rank == num_procs - 1 ) {
            MPI_Request reqs[2];
            MPI_Isend(u->image_data[1], n, MPI_FLOAT, my_rank - 1, 0, MPI_COMM_WORLD, reqs);
            MPI_Irecv(u->image_data[0], n, MPI_FLOAT, my_rank - 1, 0, MPI_COMM_WORLD, reqs + 1);
            MPI_Waitall(2, reqs, MPI_STATUS_IGNORE);
        }
        else {
            MPI_Request reqs[4];
            //Bottom
            MPI_Isend(u->image_data[1], n, MPI_FLOAT, my_rank - 1, 0, MPI_COMM_WORLD, reqs);
            MPI_Irecv(u->image_data[0], n, MPI_FLOAT, my_rank - 1, 0, MPI_COMM_WORLD, reqs + 1);

            //Top
            MPI_Isend(u->image_data[m - 2], n, MPI_FLOAT, my_rank + 1, 0, MPI_COMM_WORLD, reqs + 2);
            MPI_Irecv(u->image_data[m - 1], n, MPI_FLOAT, my_rank + 1, 0, MPI_COMM_WORLD, reqs + 3);
            MPI_Waitall(4, reqs, MPI_STATUS_IGNORE);
            }
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


int main(int argc, char *argv[]) {
    int i, j;
    int start, len, m_tmp;
    int m, n, c, iters = atoi(argv[1]);
    int my_m, my_n, my_rank, num_procs;
    double kappa = atof(argv[2]);
    char *input_jpeg_filename = argv[3], *output_jpeg_filename = argv[4];
    image u, u_bar, whole_image;
    unsigned char *image_chars, *my_image_chars;

    MPI_Status status;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

    if (my_rank==0) {
        import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
        allocate_image (&whole_image, m, n);
}
    /*Broadcast from ROOT to all processes in communicator */
    MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* divide the m x n pixels evenly among the MPI processes */
    my_m = m / num_procs;
    if ( my_rank < m % num_procs) {
        ++my_m;
    }
    my_n = n;



    /*Allocate memory to neighbour process for access of values */
    if(my_rank == 0 || my_rank == num_procs - 1) {
        allocate_image (&u, my_m + 1, my_n);
        allocate_image (&u_bar, my_m + 1, my_n);
    }
    else {
        allocate_image (&u, my_m + 2, my_n);
        allocate_image (&u_bar, my_m + 2, my_n);
    }

    /* each process asks process 0 for a partitioned region */
    /* of image_chars and copy the values into u */

    if (my_rank == 0) {
        /*Process 0 starts at the top of the picture */
        start = my_m*n;
        my_image_chars = image_chars;

        for(j = 1; j < num_procs - 1; j++){
            MPI_Recv(&m_tmp, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
            len = n*(2 + m_tmp);

            MPI_Send(image_chars + start - n, len, MPI_UNSIGNED_CHAR, j, 0, MPI_COMM_WORLD);
            start = start + m_tmp*n; /*Update new partition */
        }
        //Last process has one neighbour image
        MPI_Recv(&m_tmp, 1, MPI_INT, num_procs - 1, 0 , MPI_COMM_WORLD, &status);
        len = n*(1 + m_tmp);
        MPI_Send(image_chars + start - n, len, MPI_UNSIGNED_CHAR, num_procs -1, 0, MPI_COMM_WORLD);
        start = start + m_tmp*n;
    }

    else if (my_rank != num_procs -1){
        my_image_chars = (unsigned char*)malloc((my_m + 3)*n*sizeof(unsigned char));
        MPI_Send(&my_m, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        len = (my_m + 2)*n;
        MPI_Recv(my_image_chars, len, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &status);
    }

    else {
        my_image_chars = (unsigned char*)malloc((my_m + 2)*n*sizeof(unsigned char));
        MPI_Send(&my_m, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        len = (my_m + 1)*n;
        MPI_Recv(my_image_chars, len, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        }

    convert_jpeg_to_image(my_image_chars, &u);
    iso_diffusion_denoising (&u, &u_bar, kappa, iters);
    convert_image_to_jpeg(&u_bar, my_image_chars);

    /* each process sends its resulting content of u_bar to process 0 */
    /* process 0 receives from each process incoming values and */
    /* copy them into the designated region of struct whole_image */

    if (my_rank == 0) {
        image_chars = my_image_chars;
        start = my_m*n;
        for(j = 1; j < num_procs; j++) {

            MPI_Recv(&m_tmp, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
            len = m_tmp*n;
            MPI_Recv(image_chars + start, len, MPI_UNSIGNED_CHAR, j, 0, MPI_COMM_WORLD, &status);
            start = start + m_tmp*n;
        }
    }

    else {
        MPI_Send(&my_m, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        len = my_m*n;
        MPI_Send(my_image_chars + n, len, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
        free(my_image_chars);
    }

    if (my_rank==0) {
        export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
        deallocate_image (&whole_image);
    }

    deallocate_image (&u);
    deallocate_image (&u_bar);
    MPI_Finalize ();
    return 0;
}
