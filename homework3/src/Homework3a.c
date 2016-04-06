//////////////////////// homework2a.c /////////////////////////////////////////
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


#define NUM_ARGS 6
#define FFT 1
#define RFFT -1


///////////////////////// Main /////////////////////////////////////////////////


int main ( int argc, char ** argv)
{
  unsigned char ** img;           // Matrix Holding Original Image Values
  unsigned char ** filteredImg;           // Matrix holding vectorized image
  unsigned xSize;                 // NUmber of Horizontal pixels
  unsigned ySize;                 // Number of Vertical Pixels
  unsigned nRows;                 // Number of Rows Read in
  unsigned threshold;

  if (argc <  NUM_ARGS)
  {
    printError("Usage: HW3a <inputFile> <Outputfile>  <columns> <rows> <threshold>\n");
    exit(0);
  }


  xSize = atoi(argv[3]);

  ySize = atoi(argv[4]);

  threshold = atoi(argv[5]);


  printOK("Reading Image \n");

  img = ReadImage( argv[1], xSize, ySize, &nRows);

  if (img == NULL)
  {
    exit(-1);
  }

  if ( nRows != ySize)
  {
    ySize = nRows;
  }


  printOK("Applying Sobel\n");

  filteredImg = ApplySobel(img, xSize, ySize);

  printOK("Outputing Filtered Image \n");

  MakeBinary (filteredImg, xSize, ySize, threshold);

  OutputImage(argv[2], filteredImg, xSize, ySize);

  printOK ("Cleanup\n");

  DestroyImage (filteredImg , xSize, ySize);

  DestroyImage (img, xSize, ySize);

  exit(0);

}
