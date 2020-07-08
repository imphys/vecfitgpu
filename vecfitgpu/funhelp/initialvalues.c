#include "paramsdata.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/// this function to compute the dot product ( you may add it to the commonfunctions.h)
void dot_product(double v[], double u[], int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
         result += v[i]*u[i];
    }
}
///This function provides initial values for the fit parameters by means of a centroid fit.

void initialvalues (TIFF* tif,paramsdata *params)
{

//defining the needed variables
int  K = params -> K;
double  m = params -> m;
int Mx = params -> Mx;
int My = params -> My;
int Ncfg = params -> Ncfg;
int numparams = params -> numparams;
char fitmodel = params -> fitmodel;
double allalpha = params -> allalpha;
int ImageSizex = params -> xrange;
int ImageSizey = params -> yrange;
char doetype = params -> doetype;
//sampling coordinates in image plane

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

thetainit = (double _Complex **)malloc(numparams * sizeof(double _Complex *));
//background estimate from the median value of the rim pixels,
for(int i = 0; i < XYPupil_size; i++)
{
  thetainit[i] = (double _Complex *)malloc(Ncfg * sizeof(double _Complex ));
}

for(int i = 0; i < numparams; i++)
{
  for(int j = 0; j < Ncfg; j++)
  {
    thetainit[i][j] = 0;
  }
}

errorfun = (double _Complex **)malloc(Ncfg * sizeof(double _Complex *));

for(int i = 0; i < XYPupil_size; i++)
{
  errorfun[i] = (double _Complex *)malloc(1 * sizeof(double _Complex ));
}

for(int i = 0; i < Ncfg; i++)
{
  for(int j = 0; j < 1; j++)
  {
    errorfun[i][j] = 0;
  }
}
alpha = (double _Complex **)malloc(Ncfg * sizeof(double _Complex *));

for(int i = 0; i < Ncfg; i++)
{
  alpha[i] = (double _Complex *)malloc(1 * sizeof(double _Complex ));
}

dummat_all = (double _Complex **)malloc(Ncfg * sizeof(double _Complex *));
for(int i = 0; i < Ncfg; i++)
{
  dummat_all[i] = (double _Complex *)malloc(1 * sizeof(double _Complex ));
}

for (int jcfg = 1; jcfg < Ncfg; i++)
{
    alpha [jcfg][1] = allalpha [jcfg][1];

}
double dummat [dummat_all_size][dummat_all_size];
for (i =0; i<dummat_all_size;i++)
{
    for (j =0; j<dummat_all_size;j++)
    {
    dummat1[i][j] = dummat_all[i][j];
    }
}
double Nph =  dummat_all[1]+dummat_all[2];
double x0 = ((dummat[1]+dummat[2])*XImage)/Nph);
double yo = ((dummat[1]+dummat[2])*YImage)/Nph);

    ///centroid estimate of axial position

if((contains (params->fitmodel, "xyz")))
{
    double mask[My][Mx];
    for (int i =0; i <My;i++){
        for (int j =0; j <Mx;j++){
            mask[My][Mx] = 1;
        }
    }

  if((My>Mx)
  {
    for (int i =0; i <My;i++){
        for (int j =1; j <((My-Mx)/2);i++){
            mask[My][Mx] = 0;
        }
    }
    for (int i =0; i <My;i++){
        for (int j =(My-(My-Mx)/2)+1; j <My;i++){
            mask[My][Mx] = 0;
        }
    }
  }
double dummasksum = 0;

for (int i =0; i <My;i++){
        for (int j =0; j <Mx;j++){

            dummasksum += mask[My][Mx]*dummat[My][Mx] ;
        }
    }
    double dummasksumX2 = 0;

for (int i =0; i <My;i++){
        for (int j =0; j <Mx;j++){

            dummasksumX2 += mask[My][Mx]*dummat[My][Mx]*XImage*XImage ;
        }
    }
    double dummasksumY2 = 0;

for (int i =0; i <My;i++){
        for (int j =0; j <Mx;j++){

            dummasksumY2 += mask[My][Mx]*dummat[My][Mx]*YImage*YImage ;
        }
    }
    double dummasksumXY = 0;

for (int i =0; i <My;i++){
        for (int j =0; j <Mx;j++){

            dummasksumXY += mask[My][Mx]*dummat[My][Mx]*YImage*XImage ;
        }
    }


double Momxx = dummasksumX2;
double   Momyy = dummasksumY2;
double  Momxy = dummasksumXY;
double Nphr = dummasksum;
double Axx = Momxx-Nphr*x0*x0;
double  Ayy = Momyy-Nphr*y0*y0;
double  Axy = Momxy-Nphr*x0*y0;
double  z0 = 1250*Axy/(Axx+Ayy);

  }
double pola0 = 90*M_PI/180;
double azim0 = 45*M_PI/180;

double dummatsum2 = 0;

for (int i =0; i <2;i++){
        for (int j =0; j <Mx;j++){

            dummatsum2 += dummat[My][Mx] ;
        }
    }
double dummatsumX = 0;

for (int i =0; i <2;i++){
        for (int j =0; j <Mx;j++){

            dummatsumX += dummat[My][Mx]*XImage ;
        }
    }
double dummatsumY = 0;

for (int i =0; i <2;i++){
        for (int j =0; j <Mx;j++){

            dummatsumY += dummat[My][Mx]*YImage ;
        }
    }
double dummatsumX2 = 0;

for (int i =0; i <2;i++){
        for (int j =0; j <Mx;j++){

            dummatsumX += dummat[My][Mx]*XImage*XImage ;
        }
    }
double dummatsumY2 = 0;

for (int i =0; i <2;i++){
        for (int j =0; j <Mx;j++){

            dummatsumY += dummat[My][Mx]*YImage*YImage ;
        }
    }

double dummatsumXY = 0;

for (int i =0; i <2;i++){
        for (int j =0; j <Mx;j++){

            dummatsumY += dummat[My][Mx]*YImage*YImage ;
        }
    }
double dummatsumX2Y = 0;

for (int i =0; i <2;i++){
        for (int j =0; j <Mx;j++){

            dummatsumX += dummat[My][Mx]*XImage*XImage*YImage ;
        }
    }
double dummatsumY2X = 0;

for (int i =0; i <2;i++){
        for (int j =0; j <Mx;j++){

            dummatsumY += dummat[My][Mx]*YImage*YImage*XImage ;
        }
    }
if (K == 1)
{

    switch (doetype)
        {

            case 'none':
                azim0 = M_PI/2;
                break;
            case 'vortex':
               double M00 = dummatsum2;
               double M10 = dummatsumX;
               double M01 = dummatsumY;
               double M20 = dummatsumX2;
               double M02 = dummatsumY2;
               double M11 = dummatsumXY;
               double M21 = dummatsumX2Y;
               double M12 = dummatsumY2X;
               ///centroids
               x0 = M10/M00;
               y0 = M01/M00;
               Nph = M00*0.80;
               ///central moments
               double mu11 = (M11-x0*M01)/M00;
               double mu20 = (M20-x0*M10)/M00;
               double mu02 = (M02-y0*M01)/M00;
               double mu12 = (M12-2*y0*M11-x0*M02+2*y0*y0*M10)/M00;
               double mu21 = (M21-2*x0*M11-y0*M20+2*x0*x0*M01)/M00;

                if (abs(mu21-mu12)>abs(2*mu11)){
                    azim0 = -atan2(mu12,mu21);
                }else{
                    azim0 = 1/2*atan(mu11/(mu20-mu02));
                    if (mu11<0 && mu02>mu20){
                    }else if (mu11<0 && mu02<mu20){
                        azim0 = azim0+pi/2;
                    }else if (mu11>0 && mu02<mu20){
                        azim0 = azim0+pi/2;
                    }else if (mu11>0 && mu02>mu20){

                        }
                }
                azim0 = azim0%(2*M_PI);
                break;
        }

    }else if (K>1)
        {
            double csig = 0;
            double ssig = 0;
            for (int i = 0; i<Nph_all_size;i++){
                csig += cos(2*alpha[i])*Nph_all[i];
            }
            for (int i = 0; i<Nph_all_size;i++){
                ssig += sin(2*alpha[i])*Nph_all[i];
            }
            azim0 = (((M_PI/2)-atan2(csig,ssig))/2)%M_PI;
        }

    if (K>1){
            double signal [Nph_size];
            double model [Nph_size];
            for (int i =0; i<Nph_size;i++){

             signal[i] = Nph_all[i]/Nph[i];
            }
            for (int i =0; i<Nph_size;i++){
            double model[i] = 2*(1-m*sin(azim0-alpha[i])*sin(azim0-alpha[i]))/K/(2-m);
            }
        double* normalize;
        normalize = get_normalization(signal,model);
        Cs = dot_product(signal,model,Nph_size)/normalize;
        double distance = 2*acos(Cs)/pi;
        errorfun = distance/sqrt(normalize);
    }
    switch (fitmodel)
        {

        case 'xy':
            thetainit(1,jcfg) = x0;
            thetainit(2,jcfg) = y0;
            thetainit(3,jcfg) = Nph;
            thetainit(4,jcfg) = bg;
            break;
        case 'xyz':
            thetainit(1,jcfg) = x0;
            thetainit(2,jcfg) = y0;
            thetainit(3,jcfg) = z0;
            thetainit(4,jcfg) = Nph;
            thetainit(5,jcfg) = bg;
            break;
        case 'xyz-azim':
            thetainit(1,jcfg) = x0;
            thetainit(2,jcfg) = y0;
            thetainit(3,jcfg) = z0;
            thetainit(4,jcfg) = Nph;
            thetainit(5,jcfg) = bg;
            thetainit(6,jcfg) = azim0;
            break;
        case 'xy-azim':
            thetainit(1,jcfg) = x0;
            thetainit(2,jcfg) = y0;
            thetainit(3,jcfg) = Nph;
            thetainit(4,jcfg) = bg;
            thetainit(5,jcfg) = azim0;
            break;
        case 'xy-azim-diffusion':
            thetainit(1,jcfg) = x0;
            thetainit(2,jcfg) = y0;
            thetainit(3,jcfg) = Nph;
            thetainit(4,jcfg) = bg;
            thetainit(5,jcfg) = azim0;
            thetainit(6,jcfg) = .5;
            break;
        case 'xyz-azim-pola':
            thetainit(1,jcfg) = x0;
            thetainit(2,jcfg) = y0;
            thetainit(3,jcfg) = z0;
            thetainit(4,jcfg) = Nph;
            thetainit(5,jcfg) = bg;
            thetainit(6,jcfg) = azim0;
            thetainit(7,jcfg) = pola0;
            break;
        case 'xy-azim-pola':
            thetainit(1,jcfg) = x0;
            thetainit(2,jcfg) = y0;
            thetainit(3,jcfg) = Nph;
            thetainit(4,jcfg) = bg;
            thetainit(5,jcfg) = azim0;
            thetainit(6,jcfg) = pola0;
            break;

        }
}
