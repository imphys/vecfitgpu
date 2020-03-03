

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
//This function calculates the pupil matrix Q_{jk}, which gives the j-th
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
double sum = -PupilSize+DxyPupil/2;
int XYPupil_size = 0;
while (sum < 2.0)
{
	XYPupil_size ++;
	sum += DxyPupil;
	
}

XYPupil = (double*) malloc(XYPupil_size * sizeof(double));

XYPupil[0] = -PupilSize+DxyPupil/2;

for(int i = 1; i < XYPupil_size; i++)
{
	XYPupil[i] = XYPupil[i-1] + DxyPupil
}


double ***res, **YPupil, **XPupil;

res = meshgrid(XYPupil, XYPupil_size, XYPupil, XYPupil_size);

YPupil = res[1];
XPupil = res[0];

//CosThetaMed = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2)

double **XPupil2, **YPupil2;

}


int main(int argc, char const *argv[]) {
  /* code */


  return 0;
}
