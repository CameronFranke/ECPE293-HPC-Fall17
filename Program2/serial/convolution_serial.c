#include "image_template.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

void create_gaussian_kernels(float sigma, float **gaussian_filter,float **gaussian_d, float *w);
void convolve(float *source_image, int source_x, int source_y, float *filter, int filter_width, int filter_height, float **output_image);

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        printf("Arguments should be: exec <file> <sigma>\n");
        exit(0); 
    }   
    
    //	// start timer //	// 
    struct timeval start, end;  
    gettimeofday(&start, NULL);
    
    char *filename = argv[1]; // full path to file 
    float sigma = atof(argv[2]);	
    int imageX, imageY; 
    float *gaussian, *gaussian_d, width, *temp_horizontal, *temp_vertical, *image, *horizontal_gradient, *vertical_gradient;
    printf("Image file: %s\n", filename);
    printf("Sigma: %f\n", sigma);

    // read in image     
    read_image_template<float>(filename, &image, &imageX, &imageY);
    
    printf("Image dimensions: %d x %d\n", imageX, imageY);
    
    // create gaussian distributions 
    create_gaussian_kernels(sigma, &gaussian, &gaussian_d, &width);
    
    // allocate space for temporary images
    temp_horizontal = (float *)malloc(sizeof(float)*imageX*imageY);
    temp_vertical = (float *)malloc(sizeof(float)*imageX*imageY);
    
    // allocate space for final images 
    horizontal_gradient = (float *)malloc(sizeof(float)*imageX*imageY);
    vertical_gradient = (float *)malloc(sizeof(float)*imageX*imageY);
    
    
    
    // Horizontal Gradient
    convolve(image, imageX, imageY, gaussian, 1, width, &temp_horizontal);
    convolve(temp_horizontal, imageX, imageY, gaussian_d, width, 1, &horizontal_gradient);
    
    // Vertical Gradient 
    convolve(image, imageX, imageY, gaussian, width, 1, &temp_vertical);
    convolve(temp_vertical, imageX, imageY, gaussian_d, 1, width, &vertical_gradient);
    
    char nameh[30] = "Horizontal_gradient.pgm";
    char namev[30] = "Vertical_gradient.pgm";
    write_image_template<float>(nameh, horizontal_gradient, imageX, imageY);
    write_image_template<float>(namev, vertical_gradient, imageX, imageY);
    
    // stop timer //	//
    gettimeofday(&end, NULL);
    printf("Serial execution time(microseconds): %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
    
    return 0; 
}


void convolve(float *source_image, int source_x, int source_y, float *filter, int filter_width, int filter_height, float **output_image)
{
    int offset, coordX, coordY;
    float pixelTotal;
    for(int y = 0; y<source_y; y++) // iterate over rows of image
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
    float a = round(2.5*sigma-0.5);
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