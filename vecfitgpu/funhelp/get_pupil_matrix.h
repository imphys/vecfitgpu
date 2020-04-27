//#include "../paramsdata.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

double *** meshgrid(double *x, int size_x, double *y, int size_y);

void get_pupil_matrix(paramsdata *params);

double *** get_zernikefunctions(double **orders, int size_orders, double **x, double **y, int XYpupil_size);

double * get_rimismatchpars(paramsdata *parameters);

double * get_normalization(double _Complex ****PupilMatrix, paramsdata *parameters);
