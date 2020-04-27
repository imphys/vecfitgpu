#include "../paramsdata.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

double * get_normalization(double _Complex ****PupilMatrix, paramsdata *parameters)
{


int Npupil, lambda, pixelsize;
double NA, pola, azim;

Npupil = parameters -> Npupil;
NA = parameters -> NA;
lambda = parameters -> lambda;
pixelsize = parameters -> pixelsize;
pola = parameters -> pola;
azim = parameters -> azim;

double dipor[3] = {(sin(pola)*cos(azim)), (sin(pola)*sin(azim)), cos(pola)};


double _Complex pupmat1[Npupil][Npupil], pupmat2[Npupil][Npupil];
double IntensityMatrix[3][3];
double sum;

for(int i = 0; i<3; i++)
{
  for(int j=0; j<3; j++)
  {
    IntensityMatrix[i][j] = 0;
  }
}

for(int itel = 0; itel < 3; itel++)
{
  for(int jtel = 0; jtel < 3; jtel++)
  {

    for(int ztel = 0; ztel < 2; ztel++)
    {
      sum = 0;
      for(int i = 0; i<Npupil; i++)
      {
        for(int j=0; j<Npupil; j++)
        {
          sum += creal(PupilMatrix[i][j][ztel][itel]*conj(PupilMatrix[i][j][ztel][jtel]));
        }
      }
      IntensityMatrix[itel][jtel] = IntensityMatrix[itel][jtel] + sum;
      sum = 0;
    }
  }
}

double DxyPupil = 2.0/Npupil;
double normfac = (DxyPupil*DxyPupil)/((pixelsize*NA/lambda)* (pixelsize*NA/lambda));

for(int i = 0; i<3; i++)
{
  for(int j=0; j<3; j++)
  {
    IntensityMatrix[i][j] = normfac * IntensityMatrix[i][j];
  }
}

double normint_free, normint_fixed;
normint_free = 0;
normint_fixed = 0;

for(int i = 0; i<3; i++)
{
    normint_free += IntensityMatrix[i][i];
}

normint_free = normint_free/3.0;




//printf("%s\n", );
double Imm[3];
for(int i = 0; i<3; i++)
{
  Imm[i] = 0;
  for(int j=0; j<3; j++)
  {
    Imm[i] += IntensityMatrix[i][j] * dipor[j];
  }
}

for (int i = 0; i < 3; i++) {
  normint_fixed += dipor[i] * Imm[i];
}

double *result = (double *)malloc(2 * sizeof(double));

result[0] = normint_free;
result[1] = normint_fixed;

return result;


}
