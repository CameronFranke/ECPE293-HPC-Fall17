#include "image.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>
// make sure to get the new image.h !!!!!!!!!!!!!!!!!!!!!!!!


void create_gaussian_kernels(float sigma, float **gaussian_filter,float **gaussian_d, float *w);

int main(int argc, char *argv[])
{
    // ARGS: Path to image		sigma 		# threads
    if(argc != 3)
    {
        printf("Arguments should be: exec <file> <sigma>");
        exit(0); 
    }
    
    char *filename = argv[1]; // full path to file 
    float sigma = atof(argv[2]);	
	
    printf("Image file: %s\n", filename);
    printf("Sigma: %f\n", sigma);

    char tempname[30];
    strcpy(tempname, filename);	

    // Get image size 
    int imageX, imageY; 
  
    // read in image 
    int *image; 
    read_image_template(filename, &image, &imageX, &imageY);
    
    float *gaussian, *gaussian_d;
    float width;
    create_gaussian_kernels(sigma, &gaussian, &gaussian_d, &width);

    int pixelTotal;
    
    // convolve 
    for(int y = 0; y<imageY; y++) // iterate over rows of image
    {
        for(int x = 0; x<imageX; x++) //iterate over columns of image 
        {
            //pixelTotal = 0;
            //for(int w = 0; w<filterSize; w++)
                //check for out of bounds left, right 
                //perform multiplication 
        }
        //divide pixel totel by filter size and assign value to result image
    }
    
	return 0; 
}

void create_gaussian_kernels(float sigma, float **gaussian_filter,float **gaussian_d_filter, float *w)
{
    float a = round(2.5*sigma-0.5);
    *w = 2*a+1;
    *gaussian_filter = (float *)malloc(sizeof(float)*(*w));
    float sum = 0;
    
    // Gaussian filter calculation
    for (int i = 0; i<(*w); i++)
    {
        (*gaussian_filter)[i] =exp((-1*(i-a)*(i-a))/(2*sigma*sigma));
        sum = sum + (*gaussian_filter)[i];
    }
    printf("sum: %f\n", sum);
    for (int i = 0; i<(*w); i++)
    {
        (*gaussian_filter)[i] = (*gaussian_filter)[i]/sum;
        printf("%f\n", (*gaussian_filter)[i]);
    }
    
    // Gaussian_d_filter calculations 
    //
    //
    //
    
    
    
}

























