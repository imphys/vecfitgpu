#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <time.h>
#include "tiffio.h"

#include "set_parameters_simulation.h"
#include "funhelp/get_pupil_matrix.h"

//gcc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c FFT/fft.c -lm -g
//gcc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c get_allPSFs.c FFT/fft.c -lm -g -ltiff
//cc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c get_allPSFs.c funhelp/get_pupil_matrix.c FFT/fft.c -lm -g -ltiff

//gcc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c get_allPSFs.c funhelp/get_pupil_matrix.c funhelp/get_zernikefunctions.c FFT/fft.c -lm -g -ltiff

//gcc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c get_allPSFs.c funhelp/get_pupil_matrix.c funhelp/get_zernikefunctions.c funhelp/get_rimismatchpars.c FFT/fft.c -lm -g -ltiff


int main()
{
    paramsdata *params = (paramsdata*)malloc(sizeof(paramsdata));
    set_parameters_simulation(params);

    params -> flg_parallel = 1;
    params -> flg_writemat = 1;
    params -> flg_showplot = 1;
    params -> flg_showconv = 1;
    params -> flg_gpu = 0;

    get_pupil_matrix(params);


    return 0;
}
