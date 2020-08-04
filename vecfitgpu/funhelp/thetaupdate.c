#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <tiffio.h>
#include "paramsdata.h"
#include <time.h>

double* thetaupdate( paramsdata *params, double &dummat [Mx][My][K],double &gradlogL [numparams],double &HessianlogL[numparams][numparams],
                        double &thetaold [numparams],double &thetamax [numparams],double &thetamin [numparams],double &thetaretry [numparams])
{
int numparams = params -> numparams;
double alamda = 1.0;
double Hessian_temp[numparams][numparams];
double Bmat[numparams][numparams];
double dtheta[numparams];
double thetanew[numparams];
int i, j;
    for(i=0; i<numparams; i++) {
        for(j=0;j<numparams;j++) {
                Hessian_temp[i][j] = 0;
        }
    }
    for(i=0; i<numparams; i++) {

        Hessian_temp[i][i] = HessianlogL[i][i];
    }

    for(i=0; i<numparams; i++) {
        for(j=0;j<numparams;j++) {
           Bmat =  HessianlogL[i][j]+Hessian_temp[i][j]*alamda;
           dtheta [i]+= -Bmat[i][j]/gradlogL[i];
        }
    }
    for(i=0; i<numparams; i++) {

        thetanew[i] = thetaold[i]+dtheta[i];
    }
    if (thetanew[5]>thetamax[5]){ %%In matlab this is a for loop
        thetanew[5] = thetanew[5]-thetamax[5];
    }else if( thetanew[5]<thetamin[5]){
        thetanew[5] = thetanew[5]+thetamax[5];
    }
    if (thetanew[6]>thetamax[6]){
        thetanew[6] = thetanew[6]-thetamax[6];
    }else if( thetanew[6]<thetamin[6]){
        thetanew[6] = thetanew[5]+thetamax[6];
    }

    for(i=0; i<numparams; i++) {

        if ((thetanew[i]>thetamax[i])||(thetanew[i]<thetamin[i])){
         thetanew[i] = thetaretry[i];
        }
    }

return thetanew;
}
