/* This program was originally written by
Sumedh Naik (now at Intel) at Clemson University
as a part of his thesis titled, "Connecting Architectures,
fitness, Optimizations and Performance using an Anisotropic
Diffusion Filter. This header was also used
in Dr. Pallipuram's dissertation work. */

/*
	Additional comments made by Cameron Franke 
*/

#include "stdio.h" 
#include "math.h"
#include "stdlib.h"
#include "string.h"

#define BUFFER 512

// Function Declaration


void read_image(char *name, unsigned char **image, int *im_width, int *im_height);
void read_image_template(char *name, int **image, int *im_width, int *im_height);
void write_image(char *name, unsigned char *image, int im_width, int im_height);
void write_image_template(char *name, int *image, int im_width, int im_height);

//Function Definition

/*Call this function alone to read images*/

void read_image_template(char *name, int **image, int *im_width, int *im_height)
{
        unsigned char *temp_img;

	int i;

        read_image(name, &temp_img, im_width, im_height);
			// temp_image now poins to the image 

        *image=(int *)malloc(sizeof(int)*(*im_width)*(*im_height)); // allocate space for image

        for(i=0;i<(*im_width)*(*im_height);i++) // iterate over each all pixels in image 
        {
                (*image)[i]=(int)temp_img[i];   // convert the unsigned chars from temp+img to integers 
        }
        free(temp_img); // free up temp_img memory space 
}

void read_image(char *name, unsigned char **image, int *im_width, int *im_height)
{
	FILE *fip;
	char buf[BUFFER];
	char *parse;
	int im_size;

	fip=fopen(name,"rb"); // open name in read mode as a binary file 
	if(fip==NULL)
	{
		fprintf(stderr,"ERROR:Cannot open %s\n",name); // throw error if there is a read problem 
		exit(0);
	}
	
	fgets(buf,BUFFER,fip); // read line from file 
	do
	{
		fgets(buf,BUFFER,fip);
	}
	while(buf[0]=='#'); // read entire file into buf

	parse=strtok(buf," ");
	(*im_width)=atoi(parse); // parse out width 

	parse=strtok(NULL,"\n");
	(*im_height)=atoi(parse); // parse out height 

	fgets(buf,BUFFER,fip);
	parse=strtok(buf," ");
	
	im_size=(*im_width)*(*im_height); // calculate pixel count 
	(*image)=(unsigned char *)malloc(sizeof(unsigned char)*im_size); // allocate memory for image, insigned char so that we can use the chars as numeric values (because we are using greyscale) 
	fread(*image,1,im_size,fip); // fread file into image pointer 
	
	fclose(fip); // close file 
	// the critical part of this function is that the image * now points to the image data 
}

/* Call this function alone to write an image*/

void write_image_template(char *name, int *image, int im_width, int im_height)
{
	int i;
        unsigned char *temp_img=(unsigned char*)malloc(sizeof(unsigned char)*im_width*im_height); //allocate space for image as unsigned chars
	for(i=0;i<(im_width*im_height);i++)
        {
                temp_img[i]=(unsigned char)image[i]; // convert image to array of chars from array of ints 
        }
        write_image(name,temp_img,im_width,im_height); // call function to write the image 

        free(temp_img); // free up allocated space 
}


void write_image(char *name, unsigned char *image, int im_width, int im_height)
{
	FILE *fop; // create file pointer 
	int im_size=im_width*im_height;

	fop=fopen(name,"w+"); //open file 
	fprintf(fop,"P5\n%d %d\n255\n",im_width,im_height); // display image information 
	fwrite(image, sizeof(unsigned char),im_size,fop); // write file
	fclose(fop); // close file
}






