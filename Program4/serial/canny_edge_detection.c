#include "image_template.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265358979323846

void create_gaussian_kernels(float sigma, float **gaussian_filter,float **gaussian_d, float *w);
void convolve(float *source_image, int source_x, int source_y, float *filter, int filter_width, int filter_height, float **output_image);
void mag_and_phase(float *horizontal_gradient, float *vertical_gradient, int source_x, int source_y, float **magnitude, float **phase);
void suppression(float *magnitude, float *phase, float **suppressed, int source_x, int source_y);
void hysteresis(float *suppressed, float **hyst, int source_x, int source_y);
void finalize_edges(float **edges, float *hyst, int source_x, int source_y);
int floatcomp(const void* elem1, const void* elem2);



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
    float *gaussian, *gaussian_d, width, *temp_horizontal, *temp_vertical, *image, *horizontal_gradient, *vertical_gradient, *magnitude, *phase;
	float *suppressed, *hyst, *edges;
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

	// allocate space for phase and magnitude images
	magnitude = (float *)malloc(sizeof(float)*imageX*imageY);
	phase = (float *)malloc(sizeof(float)*imageX*imageY);

	// allocate space for suppressed image
	suppressed = (float *)malloc(sizeof(float)*imageX*imageY);

	// allocate space for hysteresis image
	hyst = (float *)malloc(sizeof(float)*imageX*imageY);
   
	// allocate space for edges image
	edges = (float *)malloc(sizeof(float)*imageX*imageY);
 
    // Horizontal Gradient
    convolve(image, imageX, imageY, gaussian, 1, width, &temp_horizontal);
    convolve(temp_horizontal, imageX, imageY, gaussian_d, width, 1, &horizontal_gradient);
    
    // Vertical Gradient 
    convolve(image, imageX, imageY, gaussian, width, 1, &temp_vertical);
    convolve(temp_vertical, imageX, imageY, gaussian_d, 1, width, &vertical_gradient);
    
	// Magnitude and phase 
	mag_and_phase(horizontal_gradient, vertical_gradient, imageX, imageY, &magnitude, &phase);

	// Suppression 
	suppression(magnitude, phase, &suppressed, imageX, imageY);
	
	// Hysteresis
	hysteresis(suppressed, &hyst, imageX, imageY);

	// Edges
	finalize_edges(&edges, hyst, imageX, imageY); 


    char name_h[30] = "Horizontal_gradient.pgm";
    char name_v[30] = "Vertical_gradient.pgm";
	char name_m[30] = "Magnitude.pgm";
	char name_ph[30] = "Phase.pgm";
	char name_s[30] = "Suppressed.pgm";
	char name_hy[30] = "Hysteresis.pgm";
	char name_e[30] = "Edges.pgm";

    write_image_template<float>(name_h, horizontal_gradient, imageX, imageY);
    write_image_template<float>(name_v, vertical_gradient, imageX, imageY);
	write_image_template<float>(name_m, magnitude, imageX, imageY);
	write_image_template<float>(name_ph, phase, imageX, imageY);
	write_image_template<float>(name_s, suppressed, imageX, imageY);	   
	write_image_template<float>(name_hy, hyst, imageX, imageY);
	write_image_template<float>(name_e, edges, imageX, imageY);
 
    // stop timer //	//
    gettimeofday(&end, NULL);
    printf("Serial execution time(microseconds): %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
   
	// free resources
	free(temp_horizontal);
	free(temp_vertical);
	free(horizontal_gradient);
	free(vertical_gradient);
	free(image);
 
    return 0; 
}

void mag_and_phase(float *horizontal_gradient, float *vertical_gradient, int source_x, int source_y, float **magnitude, float **phase)
{
	for(int i=0; i<source_y; i++)
	{
		for(int j=0; j<source_x; j++)
		{
			(*magnitude)[source_x*i+j]=sqrt(pow(vertical_gradient[source_x*i+j],2) + pow(horizontal_gradient[source_x*i+j],2));
      		(*phase)[source_x*i+j]=atan2(vertical_gradient[source_x*i+j], horizontal_gradient[source_x*i+j]);
			
		}
	} 
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
		    
		    		if (((coordX >= 0) && (coordX < source_x) && (coordY < source_y) && (coordY >= 0)))
		    		{
		    			pixelTotal = pixelTotal + (source_image[(source_x*(coordY)) + coordX] * filter[w+z*filter_width]);
		    		}
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

void suppression(float *magnitude, float *phase, float **suppressed, int source_x, int source_y)
{
	float left_pixel, right_pixel, topleft_pixel, bottomright_pixel, top_pixel, bottom_pixel, topright_pixel, bottomleft_pixel;
	float theta;
	for(int y = 0; y<source_y; y++) // iterate over rows of image
    {
        for(int x = 0; x<source_x; x++) //iterate over columns of image 
        {
			(*suppressed)[(source_x*y) + x] = magnitude[(source_x*y) + x];
			theta = phase[(source_x*y) + x];
			if (theta < 0)
			{
				theta = theta + PI; 
			}
			theta = (180/PI)*theta; // will division like this work 
			
			if (theta <= 22.5 || theta > 157.5)
			{
				if (x - 1 == 0) left_pixel = 0;
				else left_pixel =  magnitude[(source_x*y) + x - 1];
				if (x + 1 == source_x) right_pixel = 0;
				else right_pixel = magnitude[(source_x*y) + x + 1];				

				if (left_pixel > magnitude[(source_x*y) + x] || right_pixel > magnitude[(source_x*y) + x])
				{
					(*suppressed)[(source_x*y) + x] = 0;
				}
			}
			else if (theta > 22.5 && theta <= 67.5)
			{
				if (x - 1 == 0 || y - 1 == 0) topleft_pixel = 0;
                else topleft_pixel =  magnitude[(source_x*(y-1)) + x - 1]; 	
				if (x + 1 == source_x || y + 1 == source_y) bottomright_pixel = 0;
                else bottomright_pixel = magnitude[(source_x*(y + 1)) + x + 1];    

                if (topleft_pixel > magnitude[(source_x*y) + x] || bottomright_pixel > magnitude[(source_x*y) + x]) 
                {   
                    (*suppressed)[(source_x*y) + x] = 0;
                }
			}
			else if (theta > 67.5 && theta <= 112.5)
			{
				if (y - 1 == 0) top_pixel = 0;
                else top_pixel =  magnitude[(source_x*(y-1)) + x];    
                if (y + 1 == source_y) bottom_pixel = 0;
                else bottom_pixel = magnitude[(source_x*(y + 1)) + x];       

                if (top_pixel > magnitude[(source_x*y) + x] || bottom_pixel > magnitude[(source_x*y) + x]) 
                {
                    (*suppressed)[(source_x*y) + x] = 0;
                }   
			}
			else if (theta > 112.5 && theta <= 157.5)
			{
				if (x + 1 == source_x || y - 1 == 0) topright_pixel = 0;
                else topright_pixel =  magnitude[(source_x*(y-1)) + x + 1];    
                
				if (x - 1 == 0 || y + 1 == source_y) bottomleft_pixel = 0;
                else bottomleft_pixel = magnitude[(source_x*(y + 1)) + x - 1];       

                if (bottomleft_pixel > magnitude[(source_x*y) + x] || topright_pixel > magnitude[(source_x*y) + x]) 
                {
                    (*suppressed)[(source_x*y) + x] = 0;
                }   
			} 			
		}
	}
}

void hysteresis(float *suppressed, float **hyst, int source_x, int source_y)
{
	float *sorted_hyst, t_high, t_low;
	int index; 
	sorted_hyst = (float *)malloc(sizeof(float)*source_x*source_y);
	memcpy(sorted_hyst, suppressed, sizeof(float)*source_x*source_y);
	qsort(sorted_hyst, source_x*source_y, sizeof(float), floatcomp);
	memcpy((*hyst), suppressed, sizeof(float)*source_x*source_y);	
	index = source_x*source_y*.9;	
	t_high = sorted_hyst[index];
	t_low = t_high/5;
	
	for(int y = 0; y<source_y; y++) // iterate over rows of image
    {   
        for(int x = 0; x<source_x; x++) //iterate over columns of image 
        {
			if ((*hyst)[(source_x*y) + x] > t_high)
			{
				(*hyst)[(source_x*y) + x] = 255;
			}
			else if (t_low < (*hyst)[(source_x*y) + x] && t_high > (*hyst)[(source_x*y) + x])
			{
				(*hyst)[(source_x*y) + x] = 125;
			}
			else
			{
				(*hyst)[(source_x*y) + x] = 0;
			}
		}
	}		
}

int floatcomp(const void* elem1, const void* elem2)
{
    if(*(const float*)elem1 < *(const float*)elem2)
        return -1;
    return *(const float*)elem1 > *(const float*)elem2;
}

void finalize_edges(float **edges, float *hyst, int source_x, int source_y)
{
	memcpy((*edges), hyst, sizeof(float)*source_x*source_y);
	for(int y = 0; y<source_y; y++) // iterate over rows of image
    {   
        for(int x = 0; x<source_x; x++) //iterate over columns of image 
        {
			if (hyst[(source_x*y) + x] == 125)
			{
				(*edges)[(source_x*y) + x] = 0;	
				for(int w = 0; w<3; w++)
				{
					for(int z = 0; z<3; z++)
					{
						if (x - 1 + w == -1 || x - 1 + w == source_x || y - 1 + z == -1 || y -1 + z == source_y) continue;
						else
						{
							if (hyst[(source_x*(y-1+z)) + (x-1+w)] == 255)
							{
								(*edges)[(source_x*y) + x] = 255;
							}
						}
					}
				}
			}
		}
	}
} 










