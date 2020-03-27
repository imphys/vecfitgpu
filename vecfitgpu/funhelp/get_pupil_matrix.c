#include "../paramsdata.h"
#include "get_pupil_matrix.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
//This function calculates the puXYPupil_sizepil matrix Q_{jk}, which gives the j-th
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


double _Complex CosThetaMed[XYPupil_size][XYPupil_size], CosThetaCov[XYPupil_size][XYPupil_size], CosThetaImm[XYPupil_size][XYPupil_size], CosThetaImmnom[XYPupil_size][XYPupil_size];
double _Complex FresnelPmedcov[XYPupil_size][XYPupil_size], FresnelSmedcov[XYPupil_size][XYPupil_size], FresnelPcovimm[XYPupil_size][XYPupil_size], FresnelScovimm[XYPupil_size][XYPupil_size];
double _Complex FresnelP[XYPupil_size][XYPupil_size], FresnelS[XYPupil_size][XYPupil_size];

double Phi[XYPupil_size][XYPupil_size], CosPhi[XYPupil_size][XYPupil_size], SinPhi[XYPupil_size][XYPupil_size];
double _Complex CosTheta[XYPupil_size][XYPupil_size];
double SinTheta[XYPupil_size][XYPupil_size];

double _Complex pvec[XYPupil_size][XYPupil_size][3], svec[XYPupil_size][XYPupil_size][3];

for(int i = 0; i < XYPupil_size; i++)
{
  for(int j = 0; j < XYPupil_size; j++)
  {

    CosThetaMed[i][j] = csqrt(1 - (XPupil[i][j]*XPupil[i][j] + YPupil[i][j]*YPupil[i][j])*(NA*NA)/(refmed*refmed));
    CosThetaCov[i][j] = csqrt(1 - (XPupil[i][j]*XPupil[i][j] + YPupil[i][j]*YPupil[i][j])*(NA*NA)/(refcov*refcov));
    CosThetaImm[i][j] = csqrt(1 - (XPupil[i][j]*XPupil[i][j] + YPupil[i][j]*YPupil[i][j])*(NA*NA)/(refimm*refimm));
    CosThetaImmnom[i][j] = csqrt(1 - (XPupil[i][j]*XPupil[i][j] + YPupil[i][j]*YPupil[i][j])*(NA*NA)/(refimmnom*refimmnom));
    FresnelPmedcov[i][j] = 2.0/(refmed*CosThetaCov[i][j] + refcov*CosThetaMed[i][j]);
    FresnelSmedcov[i][j] = 2.0/(refmed*CosThetaMed[i][j] + refcov*CosThetaCov[i][j]);
    FresnelPcovimm[i][j] = 2*refcov*CosThetaCov[i][j]/(refcov*CosThetaImm[i][j]+refimm*CosThetaCov[i][j]);
    FresnelScovimm[i][j] = 2*refcov*CosThetaCov[i][j]/(refcov*CosThetaCov[i][j]+refimm*CosThetaImm[i][j]);
    FresnelP[i][j] = FresnelPmedcov[i][j]*FresnelPcovimm[i][j];
    FresnelS[i][j] = FresnelSmedcov[i][j]*FresnelScovimm[i][j];

    //setting of vectorial functions
    Phi[i][j] = atan2(YPupil[i][j], XPupil[i][j]);
    CosPhi[i][j] = cos(Phi[i][j]);
    SinPhi[i][j] = sin(Phi[i][j]);
    CosTheta[i][j] = csqrt(1 - (XPupil[i][j]*XPupil[i][j] + YPupil[i][j]*YPupil[i][j])*(NA*NA)/(refmed*refmed));
    SinTheta[i][j] = csqrt(1 - (CosTheta[i][j]*CosTheta[i][j]));

    pvec[i][j][0] = FresnelP[i][j] * CosTheta[i][j] * CosPhi[i][j];
    pvec[i][j][1] = FresnelP[i][j] * CosTheta[i][j] * SinPhi[i][j];
    pvec[i][j][2] = -1*FresnelP[i][j] * SinTheta[i][j];
    svec[i][j][0] = -FresnelS[i][j] * SinPhi[i][j];
    svec[i][j][1] = FresnelS[i][j]*CosPhi[i][j];
    svec[i][j][2] = 0;
  }

}

double _Complex PolarizationVector[XYPupil_size][XYPupil_size][2][3];

for(int jtel = 0; jtel<3 ; jtel++)
{
  for(int i = 0; i < XYPupil_size; i++)
  {
    for(int j = 0; j < XYPupil_size; j++)
    {
      PolarizationVector[i][j][1][jtel] = CosPhi[i][j] * pvec[i][j][jtel] - SinPhi[i][j] * svec[i][j][jtel];
      PolarizationVector[i][j][2][jtel] = SinPhi[i][j] * pvec[i][j][jtel] + CosPhi[i][j] * svec[i][j][jtel];


    }
  }
}

// definition aperture
double ApertureMask[XYPupil_size][XYPupil_size], Amplitude[XYPupil_size][XYPupil_size], Waberration[XYPupil_size][XYPupil_size];

for(int i = 0; i < XYPupil_size; i++)
{
  for(int j = 0; j < XYPupil_size; j++)
  {
    ApertureMask[i][j] = ((XPupil[i][j]*XPupil[i][j] + YPupil[i][j] * YPupil[i][j]) < 1.0);
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

  //printf("%lf\n", zernikecoefs[i]);
}

double ***allzernikes;

allzernikes = get_zernikefunctions(orders, 9, XPupil,YPupil, XYPupil_size);

for(int i = 0; i < XYPupil_size; i++)
{
  for(int j = 0; j < XYPupil_size; j++)
  {
    //numel not used
    for (int k = 0; k < 9;k++)
    {
      Waberration[i][j] = Waberration[i][j] + zernikecoefs[k]*allzernikes[i][j][k];
    }
  }
}

double rho[XYPupil_size][XYPupil_size], zoneindex[XYPupil_size][XYPupil_size], ZoneFunction[XYPupil_size][XYPupil_size];

if((strcmp (params->doetype, "none")))
{
  if(!(strcmp (params->doetype, "vortex")))
  {
    //int numzones = (sizeof(params->ringradius)/sizeof(params->ringradius[0])) - 1; if params ringradius is a vector use This
    int numzones = 0;
    for(int i = 0; i < XYPupil_size; i++)
    {
      for(int j = 0; j < XYPupil_size; j++)
      {
        rho[i][j] = sqrt((XPupil[i][j] * XPupil[i][j]) + (YPupil[i][j] * YPupil[i][j]));
        zoneindex[i][j] = 0;
      }
    }
    //Use if numzones => 1
    // for(int jj = 1; jj <=numzones)
    // {
    //   for(int i = 0; i < XYPupil_size; i++)
    //   {
    //     for(int j = 0; j < XYPupil_size; j++)
    //     {
    //       zoneindex[i][j] = zoneindex[i][j] + jj * (rho[i][j]>params->ringradius[jj-1]) * (rho[i][j] <= params -> ringradius[jj]);
    //     }
    //   }
    // }
    if(params->doelevels == 0)
    {
      for(int i = 0; i < XYPupil_size; i++)
      {
        for(int j = 0; j < XYPupil_size; j++)
        {
          ZoneFunction[i][j] = (2*zoneindex[i][j]-1) * Phi[i][j]/(2*M_PI)+0.5;
        }
      }
    }
    else
    {
      for(int i = 0; i < XYPupil_size; i++)
      {
        for(int j = 0; j < XYPupil_size; j++)
        {
          ZoneFunction[i][j] = round((params->doelevels) * ((2*zoneindex[i][j]-1) * Phi[i][j]/(2*M_PI)+0.5))/(params->doelevels);
        }
      }
    }
  }
  double DOEaberration[XYPupil_size][XYPupil_size];
  if(!(strcmp (params->doetype, "vortex")))
  {
    for(int i = 0; i < XYPupil_size; i++)
    {
      for(int j = 0; j < XYPupil_size; j++)
      {
        DOEaberration[i][j] = params->doephasedepth * (ZoneFunction[i][j]-floor(ZoneFunction[i][j]));
      }
    }
  }
}

///% compute effect of refractive index mismatch, in this function we set NA=refmed
//when actually NA>refmed, so it is not fully correct for TIRF-conditions
//zvals = [nominal stage position, free working distance, -image depth from cover slip]
  get_rimismatchpars(params);

}
