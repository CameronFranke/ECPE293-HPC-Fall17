#include "image.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>

void write_segment(int row_start, int row_stop, int image_size, int *host_image, int x_dimension, int y_dimension, int pid, int processors); 

int main(int argc, char *argv[])
{
	int p = atoi(argv[1]); 
	char *filename = argv[2];
	
	if (p != 1 && p != 2 && p != 4 && p != 8 && p != 16 && p != 32)
	{
		printf("\nERROR: invalid processor count.\nMust be 1,2,4,8,16 or 32.\n");
		exit(0);
	}
	
	//	// start timer //	// 
	struct timeval start, end;  
	gettimeofday(&start, NULL);


	printf("Processors: %d\n", p);
	printf("Filename: %s\n", filename);
	char tempname[30];
    strcpy(tempname, filename);	

	// parse out image size 
	char *parse; 
	parse = strtok(tempname, ".");
	int image_size = atoi(strtok(parse,"Lenna_org_"));
	printf("Image size: %d\n", image_size);	
	
	// calculate the number of rows per image segment
	int rows_per_segment = image_size / p;  
	printf("Pixel rows per image segment: %d\n", rows_per_segment); 
	
	// calculate the size of the image segments including ghost pixels
	int segment_x, segment_y;
	segment_x = image_size + 2;
	segment_y = rows_per_segment + 2; 
	printf("Sub-image dimensions: %d x %d\n", segment_x, segment_y);

	// read in image 
	int *image; 
	read_image_template(filename, &image, &image_size, &image_size); 
		

	// call write segment for each segemnt of the image 
	int segment_start, segment_stop; 
	for(int x = 0; x < p; x++)
	{
		segment_start = (x*rows_per_segment);
		segment_stop = (x*rows_per_segment) + (rows_per_segment); 
		write_segment(segment_start, segment_stop, image_size, image, segment_x, segment_y, x, p); 	
	}



		// stop timer //	//
	gettimeofday(&end, NULL);
	printf("Execution time(microseconds): %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

	return 0; 
}

void write_segment(int row_start, int row_stop, int image_size, int *host_image, int x_dimension, int y_dimension, int pid, int processors)
{
	printf("Writing section for rows %d - %d\n", row_start, row_stop);
	int *image_segment =malloc(sizeof(int)*(x_dimension)*(y_dimension)); // allocate space for image segment 
	
	//copy from host image
	int size_of_row = x_dimension; 
	for (int j=0; j<row_stop-row_start+1; j++)
	{
		for(int i=0; i<image_size; i++)
		{
			image_segment[((j + 1)*size_of_row ) + (i + 1)] = host_image[(row_start*image_size) + (j*image_size) + i];
		}
	}

	// fill in top/ bottom 
	// if not first or last segment 
	if (pid == 0)
	{
		// copy top row
		for(int i=0; i<x_dimension; i++)
    	{
        	image_segment[i] = image_segment[size_of_row + i];
		}
	}
	else
	{
		// get pixels for top row from host image
		for(int i=0; i<image_size; i++)
        {
			image_segment[i + 1] = host_image[((row_start - 1) *image_size) + i];
        }
 
	}

	if (pid == processors - 1) // create ghost bottom row
	{
		for(int i=0; i<x_dimension; i++)
    	{
        	image_segment[(y_dimension-1) * size_of_row + i] = image_segment[(y_dimension-2) * size_of_row + i];
    	} 

	}
	else
	{
		// get bottom row of pixels from host image 
		for(int i=0; i<x_dimension; i++)
        {
            image_segment[(y_dimension-1) * size_of_row + 1] = host_image[((row_stop + 1) *image_size) + i];
        }

	}	


	// fill in values for left and right ghost columns 
	for (int j=0; j<y_dimension; j++)
	{
		image_segment[j*size_of_row] = image_segment[j*size_of_row + 1]; // fill in left side 
		image_segment[j*size_of_row + x_dimension - 1] = image_segment[j*size_of_row + x_dimension - 2]; //fill in right side 
	}

	
	// assemble filename string 	
	char my_segment_filename[20] = "op_";
	char temp[3]; 
	sprintf(temp, "%d", pid);
	strcat(my_segment_filename, temp);
	strcat(my_segment_filename, "_");
	sprintf(temp, "%d", processors);
	strcat(my_segment_filename, temp);
	strcat(my_segment_filename, ".pgm");
       		
	write_image_template(my_segment_filename, image_segment, x_dimension, y_dimension); 
}
