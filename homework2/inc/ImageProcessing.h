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

#endif
