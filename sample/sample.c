/* sample.c
/*
/* The enclosed function reads in a one-byte/pixel unsigned char image */
/* with filename in variable fn, and of size xsize*ysize  */ 
/* the array data is declared with the statement:    */ 
/* data=(unsigned char*)calloc(xsize*ysize, sizeof(char));*/ 
/* upon exiting cdatainput, data will contain the input image */ 

#include <stdio.h>

void cdatainput(fn,xsize,ysize,data)
char *fn;
int xsize,ysize;
unsigned char *data;
{
    int i;
    FILE *fp_inp;
    unsigned char pixel;
    unsigned char max, min;

    max=0;
    min=255;

    if ((fp_inp=fopen(fn,"r")) == NULL) exit(0);
    for (i=0; i<xsize*ysize; i++) 
        {
        fread(&pixel, sizeof(char),1,fp_inp);
        data[i]=pixel;
        if (pixel>max) max=pixel;
        if (pixel<min) min=pixel;
        }

(void)fclose(fp_inp);
printf("Max and min of image are: %d %d \n", (int)max, (int)min);
}

