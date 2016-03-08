#include <stdlib.h>
#include <stdio.h>
#include "ImageProcessing.h"
#include "CursorCntl.h"

#define BUFFER_SIZE  50
#define NUM_ARGS 7
#define FFFT 1
#define RFFT -1


int main ( int argc, char ** argv)
{
  unsigned char ** img;           //Matrix Holding Image Values
  float ** vectoredImg;
  unsigned char ** result;
  unsigned xSize;                 //NUmber of horizontal pixels
  unsigned ySize;
  unsigned nRows;
  int cutOff;
  char message[BUFFER_SIZE];


  if (argc <  NUM_ARGS)
  {
    printError("Usage: HW2a <inputFile> <spectrumOut> <filteredOut> <cutoff> <rows> <columns>\n");
    exit(0);
  }

  cutOff = atoi(argv[4]);

  xSize = atoi(argv[5]);

  ySize = atoi(argv[6]);

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

  printOK("Vectorizing Image \n");

  vectoredImg = VectorizeImage( img, xSize, ySize);

  printOK("Centering Spectrum\n");

  CenterSpectrum(vectoredImg, xSize, ySize);

  printOK("Performing Fourier Transform \n");

  Fourier2D(vectoredImg , xSize, ySize, FFFT );

  printOK("Normalizing Image \n");

  result = NormalizeImage(vectoredImg, xSize, ySize);

  printOK("Outputing Spectrum Image \n");

  OutputImage(argv[2], result, xSize, ySize);

  printOK("Destroying Spectrum Image \n");

  DestroyImage (result , xSize, ySize);

  sprintf(message, "Applying Low Pass Filter Cutoff = %d\n", cutOff);
  printOK(message);

  ApplyLowPass(vectoredImg, xSize, ySize, cutOff);

  printOK("Performing Reverse Fourier \n");

  Fourier2D(vectoredImg, xSize, ySize, RFFT);

  printOK("Normalizing Image \n");

  result = NormalizeImage(vectoredImg, xSize, ySize);

  printOK("Outputing Reverse Image \n");

  OutputImage(argv[3], result, xSize, ySize);

  printOK ("cleanup\n");

  DestroyImage (result , xSize, ySize);

  DestroyFloatImage (vectoredImg , xSize, ySize);

  DestroyImage (img, xSize, ySize);

  exit(0);

}

