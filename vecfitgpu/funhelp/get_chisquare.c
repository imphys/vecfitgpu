#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <tiffio.h>
#include "paramsdata.h"
#include <time.h>

double* get_chisquare(int (*allspots)[Mx][My][Ncfg], paramsdata *params, double (*mustore)[Mx][My][Ncfg]) /// is mustore is predefined ? if yes replace it

{
clock_t t;
t = clock();
printf("\nchi-square: ");
fflush(stdout);
    //defining the used parameters.
int  K = params -> K;
int Mx = params -> Mx;
int My = params -> My;
int Ncfg = params -> Ncfg;
int keps = 0;
keps = 1000*pow(2.22,-16);
double norm = 0;
norm = K*Mx*My;

int *** spots = (int ***) malloc(Mx * sizeof(int **));

for(int i = 0; i < Mx; i++)
  {
    spots[i] = (int **) malloc(My * sizeof(int *));
    for(int j = 0; j < Ncfg; j++)
      {
        spots[i][j] = (int *) malloc(Ncfg * sizeof(int));
      }
  }
int i, j, k;
    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < Ncfg; k++) {
                    spots [i][j][k] = allspots [i][j][k];
            }
        }
    }
double mupos1 [Mx][My[Ncfg];
    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < Ncfg; k++) {
                    if (mustore [i][j][k]>0){
                    mupos1 [i][j][k] = mustore [i][j][k];
                    }else{
                        mupos1 [i][j][k]=0;
                    }
            }
        }
    }
double mupos2 [Mx][My[Ncfg];
    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < Ncfg; k++) {
                    if (mustore [i][j][k]<0){
                    mupos2 [i][j][k] = mustore [i][j][k]*keps;
                    }else{
                        mupos1 [i][j][k]=0;
                    }
            }
        }
    }

    double mupos [Mx][My[Ncfg];
    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < Ncfg; k++) {

                    mupos [i][j][k] = mupos1 [i][j][k]+mupos2 [i][j][k];
            }
        }
    }

    double chisquare [Mx][My[Ncfg];

    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < Ncfg; k++) {

                    chisquare [i][j][k] = (((mupos [i][j][k]- spots[i][j][k])*(mupos [i][j][k]- spots[i][j][k]))/mupos [i][j][k])/norm;
            }
        }
    }
double chisquare_data[];
int chisquare_size[2];
chisquare_size[0] = 1;
chisquare_size[1] = Mx*My*Ncfg;
memcpy(&chisquare_data[0], &chisquare[0], (Mx*My*Ncfg) * sizeof(double));

double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
printf("s \n", time_taken);
return chisquare_data;
}
