#include "image.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>

void write_segment(int row_start, int row_stop, int image_size, int *host_image, int x_dimension, int y_dimension, int pid, int processors); 

int main(int argc, char *argv[])
{
	// ARGS: Path to image		sigma 		# threads
	char *filename = argv[1]; // full path to file 
	int sigma = atoi(argv[2]);	
	int p = atoi(argv[3]); // number of processors  
	
	if (p != 1 && p != 2 && p != 4 && p != 8 && p != 16)
	{
		printf("\nERROR: invalid processor count.\nMust be 1,2,4,8,16 or 32.\n");
		exit(0);
	}
	
	printf("Processors: %d\n", p);
	printf("Image file: %s\n", filename);
	printf("Sigma: %d\n", p);

	char tempname[30];
    strcpy(tempname, filename);	

    // Get image size 
    int imageX, imageY; 
    printf("Image width: ");
    scanf("%d", &imageX);
    printf("Image height: ");
    scanf("%d", &imageY);
    
	printf("Image dimensions: %d x %d\n", imageX, imageY);	
	
	// calculate the number of rows per image segment
	int rows_per_segment = imageY / p;  
	printf("Pixel rows per image segment: %d\n", rows_per_segment); 
	
	// calculate the size of the image segments including ghost pixels
	int segment_x, segment_y;
	segment_x = imageX + 2;
	segment_y = rows_per_segment + 2; 
	printf("Sub-image dimensions: %d x %d\n", segment_x, segment_y);

    
    
	// read in image 
	int *image; 
	//read_image_template(filename, &image, &image_size, &image_size); 
		

	// call write segment for each segemnt of the image 
	//int segment_start, segment_stop; 
	//for(int x = 0; x < p; x++)
	//{
	//	segment_start = (x*rows_per_segment);
	//	segment_stop = (x*rows_per_segment) + (rows_per_segment); 
	//	//write_segment(segment_start, segment_stop, image_size, image, segment_x, segment_y, x, p); 	
	//}
	return 0; 
}

