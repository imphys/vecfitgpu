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

//gcc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c FFT/fft.c -lm -g




int main()
{
    paramsdata *params = (paramsdata*)malloc(sizeof(paramsdata));
    set_parameters_simulation(params);

    params -> flg_parallel = 1;
    params -> flg_writemat = 1;
    params -> flg_showplot = 1;
    params -> flg_showconv = 1;
    params -> flg_gpu = 0;




    return 0;
}
