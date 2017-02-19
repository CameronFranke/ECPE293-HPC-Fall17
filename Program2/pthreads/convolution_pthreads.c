#include "image_template.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

// compilation:  g++ -pthread convolution_pthreads.c -o ExecConvolve
// execution: 	 ./ExecConvolve Lenna_org_256.pgm 1.2 8

void create_gaussian_kernels(float sigma, float **gaussian_filter,float **gaussian_d, float *w);
void convolve(float *source_image, int source_x, int source_y, float *filter, int filter_width, int filter_height, float **output_image, int start_y, int my_num_rows);
void *create_gradients(void *thread_args);

struct thread_data
{
    int group_id, threads_in_group, imageX, imageY, w;
    char type;
    float *image, *temp_image, *final_image, *gaussian, *gaussian_d;
    pthread_barrier_t *barrier_a, *barrier_b;
};

int main(int argc, char *argv[])
{
    if(argc != 4)
    {
        printf("Arguments should be: exec <file> <sigma> <processors>\n");
        exit(0); 
    }   
    
    //	// start timer //	// 
    struct timeval start, end;  
    gettimeofday(&start, NULL);
    
    char *filename = argv[1]; // full path to file 
    float sigma = atof(argv[2]);
    int processors = atoi(argv[3]);
    int imageX, imageY, ret; 
    float *gaussian, *gaussian_d, width, *temp_horizontal, *temp_vertical, *image, *horizontal_gradient, *vertical_gradient;
    pthread_attr_t attr; 
    pthread_barrier_t barrier_va, barrier_ha, barrier_vb, barrier_hb;
    
    pthread_barrier_init(&barrier_va, NULL, processors/2);  // barrier a and barrier b for vertical and horizontal thread groups 
    pthread_barrier_init(&barrier_ha, NULL, processors/2);
    pthread_barrier_init(&barrier_vb, NULL, processors/2);
    pthread_barrier_init(&barrier_hb, NULL, processors/2);
    
    if (processors != 2 && processors != 4 && processors != 8 && processors != 16 )
    {
	printf("\nERROR: invalid processor count.\nMust be 2,4,8, or 16.\n");
	exit(0);
    }
    
    printf("Image file: %s\n", filename);
    printf("Sigma: %f\n", sigma);
    printf("Processors %d\n", processors);

    // read in image     
    read_image_template<float>(filename, &image, &imageX, &imageY);
   
    // create gaussian distributions 
    create_gaussian_kernels(sigma, &gaussian, &gaussian_d, &width);
    
    //allocate space for temporary images
    temp_horizontal = (float *)malloc(sizeof(float)*imageX*imageY);
    temp_vertical = (float *)malloc(sizeof(float)*imageX*imageY);
    
    //allocate space for gradient images
    horizontal_gradient = (float *)malloc(sizeof(float)*imageX*imageY);
    vertical_gradient = (float *)malloc(sizeof(float)*imageX*imageY);
    
    pthread_t threads[processors];				// create thread objects
    struct thread_data thread_args[processors];			// create aregument structures
    
    
    for (int x = 0; x < processors; x++)
    {
	if (x < processors/2)                     // init group specific arguments 
	{
	    thread_args[x].type = 'v';
	    thread_args[x].group_id = x;
	    thread_args[x].barrier_a = &barrier_va;
	    thread_args[x].barrier_b = &barrier_vb;
	    thread_args[x].temp_image = temp_vertical;
	    thread_args[x].final_image = vertical_gradient;
	}
	else
	{
	    thread_args[x].type = 'h';
	    thread_args[x].group_id = x - processors/2;
	    thread_args[x].barrier_a = &barrier_ha;
	    thread_args[x].barrier_b = &barrier_hb;
	    thread_args[x].temp_image = temp_horizontal;
	    thread_args[x].final_image = horizontal_gradient;
	}
	thread_args[x].imageX = imageX;                         // init general arguments 
	thread_args[x].imageY = imageY;
	thread_args[x].threads_in_group = processors/2;
	thread_args[x].image = image;
	thread_args[x].gaussian = gaussian;
	thread_args[x].gaussian_d = gaussian_d;
	thread_args[x].w = width;
	
	ret = pthread_create(&threads[x], NULL, create_gradients, (void *)&thread_args[x]); // id, attr, function, args
    }
    
    for (int x = 0; x < processors; x++) // wait for all threads to complete 
    {
	pthread_join(threads[x], NULL);
    }
    
    // stop timer //	//
    gettimeofday(&end, NULL);
    printf("Parallel execution time(microseconds): %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
    return 0; 
}


void *create_gradients(void *thread_args)
{
    // reclaim data from struct 
    struct thread_data *mydata;
    mydata = (struct thread_data *)thread_args;
    int group_id = mydata->group_id;
    int imageX = mydata->imageX;
    int imageY = mydata->imageY;
    int group_size = mydata->threads_in_group;
    int w = mydata->w;
    char type = mydata->type;
    float *source_image = mydata->image;
    float *temp_image = mydata->temp_image;
    float *final_image = mydata->final_image;
    float *gaussian = mydata->gaussian;
    float *gaussian_d = mydata->gaussian_d;
    pthread_barrier_t *barrier_a = mydata->barrier_a;
    pthread_barrier_t *barrier_b = mydata->barrier_b;
    
    int my_start, rows_to_compute;
    rows_to_compute = (imageY/group_size);	// number of rows of the image to compute 
    my_start = rows_to_compute*group_id; 	// this is the starting row
    
    if (type == 'v')
    {
      convolve(source_image, imageX, imageY, gaussian, w, 1, &temp_image, my_start, rows_to_compute);         // is silultanious writing to different regions of an array safe?
      pthread_barrier_wait(barrier_a);      
      convolve(temp_image, imageX, imageY, gaussian_d, 1, w, &final_image, my_start, rows_to_compute);
      pthread_barrier_wait(barrier_b);
   
      if (group_id == 0)
      {
	char nameh[30] = "Vertical_gradient.pgm";
	write_image_template<float>(nameh, final_image, imageX, imageY);
      }
    }
   
    if (type == 'h')
    {
      convolve(source_image, imageX, imageY, gaussian, 1, w, &temp_image, my_start,rows_to_compute);
      pthread_barrier_wait(barrier_a);
      convolve(temp_image, imageX, imageY, gaussian_d, w, 1, &final_image, my_start,rows_to_compute);
      pthread_barrier_wait(barrier_b);
    
      if (group_id == 0)
      {
	char nameh[30] = "Horizontal_gradient.pgm";
	write_image_template<float>(nameh, final_image, imageX, imageY);
      }
    }
    pthread_exit(NULL);
}


void convolve(float *source_image, int source_x, int source_y, float *filter, int filter_width, int filter_height, float **output_image, int start_y, int my_num_rows)
{
    int offset, coordX, coordY;
    float pixelTotal;
    
    for(int y = start_y; y<(start_y + my_num_rows); y++) // iterate over rows of image
    {
        for(int x = 0; x<source_x; x++) //iterate over columns of image 
        {
            pixelTotal = 0;
            for(int w = 0; w<filter_width; w++)
	    {
		for(int z = 0; z<filter_height; z++)
		{
		    coordX = x - floor(filter_width/2) + w;
		    coordY = y - floor(filter_height/2) + z;
		    
		    if ((coordX < 0) || (coordX > source_x) || (coordY > source_y) || (coordY < 0))
		    {
		      continue; 
		    }
		    pixelTotal = pixelTotal + (source_image[(source_x*(coordY)) + coordX] * filter[w+z]);
		}
	    }
	    (*output_image)[(source_x*y) + x] = pixelTotal;
        }
    }
}


void create_gaussian_kernels(float sigma, float **gaussian_filter,float **gaussian_d_filter, float *w)
{
    float a = round(2.5*sigma-0.5);                            // this function implements the pseudocode from lecture 
    float sum = 0;
    *w = 2*a+1;
    *gaussian_filter = (float *)malloc(sizeof(float)*(*w));
    *gaussian_d_filter = (float *)malloc(sizeof(float)*(*w));
    
    // Gaussian filter calculation
    for (int i = 0; i<(*w); i++)
    {
        (*gaussian_filter)[i] =exp((-1*(i-a)*(i-a))/(2*sigma*sigma));
        sum = sum + (*gaussian_filter)[i];
    }
    printf("Gaussian Filter: \n");
    for (int i = 0; i<(*w); i++)
    {
        (*gaussian_filter)[i] = (*gaussian_filter)[i]/sum;
        printf("\t%f\n", (*gaussian_filter)[i]);
    }
    
    // Gaussian_d_filter calculations 
    for (int i = 0; i<(*w); i++)
    {
        (*gaussian_d_filter)[i] = -1*(i-a)*exp((-1*(i-a)*(i-a))/(2*sigma*sigma));
        sum = sum - i*(*gaussian_d_filter)[i];
    }
    printf("Gaussian Derivitive Filter: \n");
    for (int i = 0; i<(*w); i++)
    {
        (*gaussian_d_filter)[i] = (*gaussian_d_filter)[i]/sum;
        printf("\t%f\n", (*gaussian_d_filter)[i]);
    }
}
























