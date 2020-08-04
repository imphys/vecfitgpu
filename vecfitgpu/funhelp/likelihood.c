#include "paramsdata.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
double likelihood (paramsdata *params,double &dmudtheta [Mx][My][K][numparams],double &mu [Mx][My][K],
                    double &gradlogL [numparams],double &HessianlogL[numparams][numparams])
{
int  K = params -> K;
double  m = params -> m;
int Mx = params -> Mx;
int My = params -> My;
int numparams = params -> numparams;
int varfit = params -> varfit;
int keps = 0;
keps = 1000*pow(2.22,-16);
int i, j, k,q,x;

double mupos1 [Mx][My][K];
    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < K; k++) {
                    if (mu [i][j][k]>0){
                    mupos1 [i][j][k] = mu [i][j][k];
                    }else{
                        mupos1 [i][j][k]=0;
                    }
            }
        }
    }
double mupos2 [Mx][My][K];
    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < K; k++) {
                    if (mu [i][j][k]<0){
                    mupos2 [i][j][k] = mu [i][j][k]*keps;
                    }else{
                        mupos1 [i][j][k]=0;
                    }
            }
        }
    }

    double mupos [Mx][My][K];
    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < K; k++) {

                    mupos [i][j][k] = mupos1 [i][j][k]+mupos2 [i][j][k];
            }
        }
    }
    double weight [Mx][My][K];
    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < K; k++) {

                    weight [i][j][k] = (m-mupos [i][j][k])/(varfit+ mupos [i][j][k]);
            }
        }
    }
    double dweight [Mx][My][K];
    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < K; k++) {

                    dweight [i][j][k] = (m+varfit)/((varfit+ mupos [i][j][k])*(varfit+ mupos [i][j][k]));
            }
        }
    }
    double logl = 0;
    for(i=0; i<Mx; i++) {
        for(j=0;j<My;j++) {
            for (k = 0; k < K; k++) {

                    logl  += ((m+varfit)*log(varfit+ mupos [i][j][k]))-(varfit+ mupos [i][j][k]);
            }
        }
    }
    double gradlogL[numparams];
    for (q=0; q<numparams;q++){
        for(i=0; i<Mx; i++) {
            for(j=0;j<My;j++) {
                for (k = 0; k < K; k++) {

                    gradlogL[numparams-q-1]  +=  (weight [i][j][k])*(dmudtheta [i][j][k][q]);
                }
            }
        }
    }

    for(i=0; i<numparams; i++) {
        for(j=0;j<numparams;j++) {
                HessianlogL[i][j] = 0;
        }
    }

    for(q=0; q<numparams; q++) {
        for(x=0;x<numparams;x++) {
            for(i=0; i<Mx; i++) {
                for(j=0;j<My;j++) {
                    for (k = 0; k < K; k++) {

                HessianlogL[q][x] += (-dweight [i][j][k])*(dmudtheta [i][j][k][q])*(dmudtheta [i][j][k][x]);
                    }
                }
            }
        }
    }
    return logl;

}
