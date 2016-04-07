//////////////////////// ImageProcessing.h /////////////////////////////////////////
/*
By:   Sergio Coronado
        16.484 Computer Vision
        Assignemnt #2
        Part 1

PURPOSE:
  Function Protoypes for image processing

*/
#ifndef IMAGE_PROC
#define IMAGE_PROC



unsigned char ** ReadImage(
  char * filename,
  unsigned xSize,
  unsigned ySize,
  unsigned * numRowsi);

void Fourier2D( float ** img, unsigned xSize, unsigned ySize, int iSign);

float ** VectorizeImage( unsigned char ** img, unsigned xSize, unsigned ySize);

unsigned char ** NormalizeImage( float ** img, unsigned xsize, unsigned ySize);

void OutputImage ( char * filename, unsigned char ** img, unsigned xSize, unsigned ySize);

void DestroyFloatImage (float ** img, unsigned xSize, unsigned ySize);

void DestroyImage ( unsigned char ** img, unsigned xSize, unsigned ySize);

void CenterSpectrum ( float ** img, unsigned xSize, unsigned ySize);

void ApplyLowPass ( float ** img, unsigned xSize, unsigned ySize, int cutOff);

unsigned char ** ApplySobel ( unsigned char ** img, unsigned xSize, unsigned ySize);

void MakeBinary ( unsigned char ** img, unsigned xSize, unsigned ySize, unsigned threshold);

void ApplyNoise ( unsigned char ** img, unsigned xSize, unsigned ySize, unsigned char intensity);

unsigned char **  MakeStandard (unsigned xSize, unsigned ySize, unsigned tiers);

float ** ApplyConvolution ( unsigned char ** img,unsigned xSize, unsigned ySize, float **filter ,  unsigned size);

unsigned char ** ApplyGaussian (unsigned char ** img, unsigned xSize, unsigned ySize, unsigned scale);

unsigned FindOptimalThreshold(unsigned char ** standard, unsigned char ** img, unsigned xSize, unsigned ySize, int * errorValue);

unsigned CompareImage ( unsigned char ** standard, unsigned char ** img, unsigned xSize, unsigned ySize);

unsigned char ** NormalizeFloatImage ( float ** img, unsigned xSize, unsigned ySize);

#endif
