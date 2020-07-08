#include <stdlib.h>
#include <tiffio.h>
#include "paramsdata.h"

int** read_data(TIFF* tif, paramsdata *params)
{
	int* raster = NULL;
	int **raster2D = NULL;

  if (tif)
  {
      int w, h;
  		size_t npixels;
  		int i=0, j=0;

  		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
  		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
      params->Mx = w;
      params->My = h;


  		npixels = w * h;
  		raster = (int*) _TIFFmalloc(npixels * sizeof (int));
  		raster2D = (int **)malloc(w * sizeof(int *));

  		if (raster != NULL)
  		{
  		    if (TIFFReadRGBAImage(tif, w, h, raster, 0))
          {

  				for (i = 0; i<w; i++)
  				{
  				      raster2D[i] = (int *)malloc(h * sizeof(int));
                for (j =0; j<h; j++)
                {
  					       int pixel = raster[i*w+j];
                   int R = (pixel & 0x000000FF);
                   int G = (pixel & 0x0000FF00) >> 8;
                   int B = (pixel & 0x00FF0000) >> 16;
                   int A = (pixel & 0xFF000000) >> 24;
                   int grey_scale = (R+G+B)/3;
                   raster2D[i][j] = grey_scale;
									 //printf("%d\t", grey_scale);
                 }
  					//printf("\n");
  				}
  		   }
  		}
      }
	return raster2D;
}


void get_PSFs(int ****allPSFs_adres, paramsdata *params)
{
  TIFF* tif = TIFFOpen("spots.tif", "r");
  int stack_size = 0;
  int dircount = 0;
  if (tif)
  {

  do
  {
      dircount++;
  } while (TIFFReadDirectory(tif));
  stack_size = dircount;


  TIFFClose(tif);
  }

  tif = TIFFOpen("spots.tif", "r");
  int*** image_vector = (int***)malloc(stack_size * sizeof(int**));

  if (tif)
  {
  int dircount = 0;
  do {
    image_vector[dircount] = read_data(tif, params);
      dircount++;
  } while (TIFFReadDirectory(tif));
  TIFFClose(tif);

  }
  params->Ncfg = dircount / params->K;
  *allPSFs_adres = image_vector;
}
