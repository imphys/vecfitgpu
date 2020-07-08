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
#include "funhelp/get_stack.h"

//gcc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c FFT/fft.c -lm -g
//gcc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c get_allPSFs.c FFT/fft.c -lm -g -ltiff
//cc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c get_allPSFs.c funhelp/get_pupil_matrix.c FFT/fft.c -lm -g -ltiff

//gcc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c get_allPSFs.c funhelp/get_pupil_matrix.c funhelp/get_zernikefunctions.c FFT/fft.c -lm -g -ltiff

//gcc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c get_allPSFs.c funhelp/get_pupil_matrix.c funhelp/get_zernikefunctions.c funhelp/get_rimismatchpars.c FFT/fft.c -lm -g -ltiff

//gcc -o test vectorPSFfitter_simulation.c set_parameters_simulation.c get_allPSFs.c funhelp/get_pupil_matrix.c funhelp/get_zernikefunctions.c funhelp/get_rimismatchpars.c funhelp/get_normalization.c funhelp/get_stack.c FFT/fft.c -lm -g -ltiff


//long double
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

    params -> phtrue = 1e4;
    params -> bgtrue = 10;



    int Mx = params->Mx, My = params->My, Ncfg = params->Ncfg;


    int *** spots = (int ***) malloc(Mx * sizeof(int **));

    for(int i = 0; i < Mx; i++)
    {
      spots[i] = (int **) malloc(My * sizeof(int *));
      for(int j = 0; j < Ncfg; j++)
      {
        spots[i][j] = (int *) malloc(Ncfg * sizeof(int));
      }
    }


    for(int i = 0; i < Ncfg; i++)
    {
      for(int j = 0; j < My; j++)
      {
        for(int k = 0; k < Mx; k++)
        {
          spots[k][j][i] = params -> phtrue * params->allPSFs[k][j][i] + params->bgtrue/params->K;
        }
      }
    }


    for (int i = 0; i < Mx; i++)
    {
      for (int j = 0; j < Mx; j++)
      {
        printf("%d\t", params->allPSFs[0][i][j]);
      }
      printf("\n" );
    }


    int **** PSFs;
    //spots = 1e12*imnoise(spots*1e-12,'poisson');
    PSFs = get_stack(params, spots);

    //% do MLE fit
    //%%% initial values
    thetainit = initialvalues(PSFs,params);

    return 0;
}
