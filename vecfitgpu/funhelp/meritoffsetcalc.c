#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <tiffio.h>
#include "paramsdata.h"
#include <time.h>

double* meritoffsetcalc(double & allspots[Mx][My][K][Ncfg], paramsdata *params, double &dummat [Mx][My][K])
{
int  K = params -> K;
double  m = params -> m;
int Mx = params -> Mx;
int My = params -> My;
int varfit = params -> varfit;
int Ncfg = params -> Ncfg;

Mx = sizeof(allspots) / sizeof(int);
My = sizeof(allspots[0]) / sizeof(int);
K = sizeof(allspots[0][0]) / sizeof(int);
Ncfg = sizeof(allspots[0][0][0]) / sizeof(int);
double meritoffset[Ncfg];
int i, j, k,q;
for (i = 0; i<Ncfg;i++){
    meritoffset[i] = 0;
}

    for(q=0; q<Ncfg; q++) {
        meritoffset[q] = 0;
        for(j=0;x<Mx;j++) {
            for(i=0; i<My; i++) {
                for (k = 0; k < K; k++) {

                dummat[j][i][k] = (allspots [j][i][k][q]);
                dummat[j][i][k] = fmax(dummat[j][i][k], 0.1);
                meritoffset[q] = meritoffset[q] -gammln(dummat[j][i][k]+1+varfit);

                }
            }
        }

    }
    return meritoffset;

}
