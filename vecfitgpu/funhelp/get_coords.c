#include "paramsdata.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void get_coords(paramsdata *params)
{


//defining the needed variables
double Npupil = params -> Npupil;
int Mx = params -> Mx;
int My = params -> My;
int ImageSizex = params -> xrange;
int ImageSizey = params -> yrange;

/// computing [YPupil,XPupil]

double PupilSize = 1.0;
double DxYPupil = 2*PupilSize/Npupil;

double *XYPupil;

int XYPupil_size = 0;
double stop = -PupilSize+DxYPupil/2;

while(stop <= PupilSize)
{
  XYPupil_size += 1;
  stop += DxYPupil;

}

XYPupil = (double *) malloc(XYPupil_size * sizeof(double));

stop = -PupilSize+DxYPupil/2;

for(int i = 0; i< XYPupil_size; i++)
{
  XYPupil[i] = stop;
  stop += DxYPupil;
}
double ***res, **YPupil, **XPupil;

res = meshgrid(XYPupil, XYPupil_size, XYPupil, XYPupil_size);
YPupil = res[0];
XPupil = res[1];
/////////////////////////////////
///computing [YImage,XImage]
double DxImage = 2*ImageSizex/Mx;
double DyImage = 2*ImageSizey/My;

double *ximagelin;
double *yimagelin;

int ximagelin_size = 0;
int yimagelin_size = 0;

double stop1 = -ximagelin_size+DxImage/2;
double stop2 = -yimagelin_size+DyImage/2;

while(stop1 <= ximagelin_size)
{
  ximagelin_size += 1;
  stop1 += DxImage;

}

while(stop2 <= yimagelin_size)
{
  yimagelin_size += 1;
  stop2 += DyImage;

}

ximagelin_size = (double *) malloc(ximagelin_size * sizeof(double));
yimagelin_size = (double *) malloc(yimagelin_size * sizeof(double));

stop1 = -ximagelin_size+DxImage/2;
stop2 = -yimagelin_size+DyImage/2

for(int i = 0; i< ximagelin_size; i++)
{
  ximagelin[i] = stop1;
  stop1 += DxImage;
}
for(int i = 0; i< yimagelin_size; i++)
{
  yimagelin[i] = stop2;
  stop2 += DyImage;
}
double ***res1, **YImage, **XImage;

res = meshgrid(yimagelin, yimagelin_size, ximagelin, ximagelin_size);
YImage = res1[0];
XImage = res1[1];
}
