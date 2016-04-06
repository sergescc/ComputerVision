//////////////////////// homework3b.c /////////////////////////////////////////
/*
By:   Sergio Coronado
        16.484 Computer Vision
        Assignemnt #2
        Part 1

PURPOSE:
  Program reads in an image performs a FFT outputs the spectrum and then
  performs the inverse FFT

USAGE:

  HW2a <input-file> <ospectrum-file> <reverse-file> <xSize> <ySize>

*/

///////////////////////// Includes /////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include "ImageProcessing.h"
#include "CursorCntl.h"

///////////////////////// Constants ////////////////////////////////////////////


#define NUM_ARGS 5


///////////////////////// Main /////////////////////////////////////////////////


int main ( int argc, char ** argv)
{
	unsigned char ** img;

	unsigned char ** sobel;
					           // Matrix Holding Original Image Values
  	unsigned xSize;                 // NUmber of Horizontal pixels
  	unsigned ySize;
  	unsigned tiers;                 // Number of Vertical Pixels
 

  	if (argc <  NUM_ARGS)
  	{
    	printError("Usage: HW3b <Outputfile>  <columns> <rows> <tiers> \n");
    	exit(0);
  	}

  	xSize = atoi(argv[2]);

  	ySize = atoi(argv[3]);

  	tiers = atoi(argv[4]);

  	printOK("Creating  Image \n");

  	img = MakeStandard(xSize, ySize, tiers);

  	printOK("Applying Sobel");

  	sobel = ApplySobel(img, xSize, ySize);

  	printOK("Outputting Image \n");

  	//MakeBinary(sobel, xSize, ySize, 30);

  	OutputImage ( argv[1], img, xSize, ySize);

  	printOK("Clean up\n");

  	DestroyImage(img, xSize, ySize);

	exit(0);

}
