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
double XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;

//[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);

}
