#include "../paramsdata.h"
#include "get_pupil_matrix.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double *** meshgrid(double *x, int size_x, double *y, int size_y)
{

  double ***res, **X, **Y;

  res = (double ***)malloc(2 * sizeof(double **));


  X = (double **)malloc(size_y * sizeof(double *));
  Y = (double **)malloc(size_y * sizeof(double *));
  for(int i = 0; i < size_y; i++)
  {
    X[i] = (double *)malloc(size_x * sizeof(double ));
    Y[i] = (double *)malloc(size_x * sizeof(double ));
  }

  for(int i = 0; i < size_y; i++)
  {
    for(int j = 0; j < size_x; j++)
    {
      X[i][j] = x[j];
      Y[i][j] = y[i];
    }

  }

  res[0] = X;
  res[1] = Y;

  return res;

}



void get_pupil_matrix(paramsdata *params)
{
//This function calculates the puXYpupil_sizepil matrix Q_{jk}, which gives the j-th
//electric field component proportional to the k-th dipole vector
//component.

//parameters: NA, refractive indices of medium, cover slip, immersion fluid,
//wavelength (in nm), sampling in pupil

double NA = params -> NA;
double refmed = params -> refmed;
double refcov = params -> refcov;
double refimm = params -> refimm;
double refimmnom = params -> refimmnom;
double lambda = params -> lambda;
double Npupil = params -> Npupil;

//pupil radius (in diffraction units) and pupil coordinate sampling
double PupilSize = 1.0;
double DxyPupil = 2*PupilSize/Npupil;

double *XYPupil;


int XYpupil_size = 0;
double stop = -PupilSize+DxyPupil/2;

while(stop <= PupilSize)
{
  XYpupil_size += 1;
  stop += DxyPupil;

}

XYPupil = (double *) malloc(XYpupil_size * sizeof(double));

stop = -PupilSize+DxyPupil/2;

for(int i = 0; i< XYpupil_size; i++)
{
  XYPupil[i] = stop;
  stop += DxyPupil;
}
double ***res, **Ypupil, **Xpupil;

res = meshgrid(XYPupil, XYpupil_size, XYPupil, XYpupil_size);


Ypupil = res[0];
Xpupil = res[1];


double _Complex CosThetaMed[XYpupil_size][XYpupil_size], CosThetaCov[XYpupil_size][XYpupil_size], CosThetaImm[XYpupil_size][XYpupil_size], CosThetaImmnom[XYpupil_size][XYpupil_size];
double _Complex FresnelPmedcov[XYpupil_size][XYpupil_size], FresnelSmedcov[XYpupil_size][XYpupil_size], FresnelPcovimm[XYpupil_size][XYpupil_size], FresnelScovimm[XYpupil_size][XYpupil_size];
double _Complex FresnelP[XYpupil_size][XYpupil_size], FresnelS[XYpupil_size][XYpupil_size];

double Phi[XYpupil_size][XYpupil_size], CosPhi[XYpupil_size][XYpupil_size], SinPhi[XYpupil_size][XYpupil_size];
double _Complex CosTheta[XYpupil_size][XYpupil_size];
double SinTheta[XYpupil_size][XYpupil_size];

double _Complex pvec[XYpupil_size][XYpupil_size][3], svec[XYpupil_size][XYpupil_size][3];

for(int i = 0; i < XYpupil_size; i++)
{
  for(int j = 0; j < XYpupil_size; j++)
  {

    CosThetaMed[i][j] = csqrt(1 - (Xpupil[i][j]*Xpupil[i][j] + Ypupil[i][j]*Ypupil[i][j])*(NA*NA)/(refmed*refmed));
    CosThetaCov[i][j] = csqrt(1 - (Xpupil[i][j]*Xpupil[i][j] + Ypupil[i][j]*Ypupil[i][j])*(NA*NA)/(refcov*refcov));
    CosThetaImm[i][j] = csqrt(1 - (Xpupil[i][j]*Xpupil[i][j] + Ypupil[i][j]*Ypupil[i][j])*(NA*NA)/(refimm*refimm));
    CosThetaImmnom[i][j] = csqrt(1 - (Xpupil[i][j]*Xpupil[i][j] + Ypupil[i][j]*Ypupil[i][j])*(NA*NA)/(refimmnom*refimmnom));
    FresnelPmedcov[i][j] = 2.0/(refmed*CosThetaCov[i][j] + refcov*CosThetaMed[i][j]);
    FresnelSmedcov[i][j] = 2.0/(refmed*CosThetaMed[i][j] + refcov*CosThetaCov[i][j]);
    FresnelPcovimm[i][j] = 2*refcov*CosThetaCov[i][j]/(refcov*CosThetaImm[i][j]+refimm*CosThetaCov[i][j]);
    FresnelScovimm[i][j] = 2*refcov*CosThetaCov[i][j]/(refcov*CosThetaCov[i][j]+refimm*CosThetaImm[i][j]);
    FresnelP[i][j] = FresnelPmedcov[i][j]*FresnelPcovimm[i][j];
    FresnelS[i][j] = FresnelSmedcov[i][j]*FresnelScovimm[i][j];

    //setting of vectorial functions
    Phi[i][j] = atan2(Ypupil[i][j], Xpupil[i][j]);
    CosPhi[i][j] = cos(Phi[i][j]);
    SinPhi[i][j] = sin(Phi[i][j]);
    CosTheta[i][j] = csqrt(1 - (Xpupil[i][j]*Xpupil[i][j] + Ypupil[i][j]*Ypupil[i][j])*(NA*NA)/(refmed*refmed));
    SinTheta[i][j] = csqrt(1 - (CosTheta[i][j]*CosTheta[i][j]));

    pvec[i][j][0] = FresnelP[i][j] * CosTheta[i][j] * CosPhi[i][j];
    pvec[i][j][1] = FresnelP[i][j] * CosTheta[i][j] * SinPhi[i][j];
    pvec[i][j][2] = -1*FresnelP[i][j] * SinTheta[i][j];
    svec[i][j][0] = -FresnelS[i][j] * SinPhi[i][j];
    svec[i][j][1] = FresnelS[i][j]*CosPhi[i][j];
    svec[i][j][2] = 0;
  }

}

double _Complex PolarizationVector[XYpupil_size][XYpupil_size][2][3];

for(int jtel = 0; jtel<3 ; jtel++)
{
  for(int i = 0; i < XYpupil_size; i++)
  {
    for(int j = 0; j < XYpupil_size; j++)
    {
      PolarizationVector[i][j][1][jtel] = CosPhi[i][j] * pvec[i][j][jtel] - SinPhi[i][j] * svec[i][j][jtel];
      PolarizationVector[i][j][2][jtel] = SinPhi[i][j] * pvec[i][j][jtel] + CosPhi[i][j] * svec[i][j][jtel];


    }
  }
}

// definition aperture
double ApertureMask[XYpupil_size][XYpupil_size], Amplitude[XYpupil_size][XYpupil_size], Waberration[XYpupil_size][XYpupil_size];

for(int i = 0; i < XYpupil_size; i++)
{
  for(int j = 0; j < XYpupil_size; j++)
  {
    ApertureMask[i][j] = ((Xpupil[i][j]*Xpupil[i][j] + Ypupil[i][j] * Ypupil[i][j]) < 1.0);
    Amplitude[i][j] = ApertureMask[i][j] * csqrt(CosThetaImm[i][j]);


    Waberration[i][j] = 0;
  }
}

//calculation aberration function
double **orders, zernikecoefs[9], normfac[9];

orders =(double **) malloc(9 * sizeof(double*));

for(int i = 0; i < 9;i++ )
{
  orders[i] = (double *)malloc(2 * sizeof(double));
}


for(int i = 0; i < 9; i++)
{
  orders[i][0] = params -> aberrations[i][0];
  orders[i][1] = params -> aberrations[i][1];



  zernikecoefs[i] = params -> aberrations[i][2];

  normfac[i] = sqrt(2*(orders[i][0] + 1)/(1 + (orders[i][1] == 0)));

  zernikecoefs[i] = normfac[i] * zernikecoefs[i];

  printf("%lf\n", zernikecoefs[i]);
}


get_zernikefunctions(orders, 9, Xpupil,Ypupil, XYpupil_size);

//NO SQUEEZE functin used


// allzernikes = get_zernikefunctions(orders,XPupil,YPupil);
// for j = 1:numel(zernikecoefs)
//     Waberration = Waberration+zernikecoefs(j)*allzernikes(:,:,j);
// end


}
