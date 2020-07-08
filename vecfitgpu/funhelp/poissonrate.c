#include "paramsdata.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/// the output matrices will be passed by reference and filled inside this function, why will be defined in the main file.
// define then like this in the main
// double dmudtheta [Mx][My][K][numparams];
//double mu [Mx][My][K];
/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void poissonrate (paramsdata *params, double *theta[K], double *alpha,double &dmudtheta [Mx][My][K][numparams],double &mu [Mx][My][K],double &PSFder [Mx][My][K])
{
int  K = params -> K;
double  m = params -> m;
int Mx = params -> Mx;
int My = params -> My;
int numparams = params -> numparams;
char fitmodel = params -> fitmodel;
char excitation = params -> excitation;
double  xemit = params -> xemit;
double  yemit = params -> yemit;
double  zemit = params -> zemit;
double  azim = params -> azim;
double  pola = params -> pola;
double  alpha = params -> alpha;
double _Complex PupilMatrix = params -> PupilMatrix;
double  g2 = params -> g2;
double npg = 0;
double nbh = 0;
double phi = 0;
double p = 0;
double dp = 0;

switch (fitmodel){
    case 'xy':
        xemit = theta[0];
        yemit = theta[1];
        Nph = theta[2];
        Nbg = theta[3];
        break;
    case 'xyz':
        xemit = theta[0];
        yemit = theta[1];
        zemit = theta[2];
        Nph = theta[3];
        Nbg = theta[4];
        break;
    case 'xy-azim':
        xemit = theta[0];
        yemit = theta[1];
        Nph = theta[2];
        Nbg = theta[3];
        phi = theta[4];
        azim = phi;
        break;
    case 'xy-azim-pola':
        xemit = theta[0];
        yemit = theta[1];
        Nph = theta[2];
        Nbg = theta[3];
        phi = theta[4];
        pola = theta[5];
        azim = phi;
        break;
    case 'xyz-azim-pola':
        xemit = theta[0];
        yemit = theta[1];
        zemit = theta[2];
        Nph = theta[3];
        Nbg = theta[4];
        phi = theta[5];
        pola = theta[6];
        azim = phi;
        break;
    case 'xy-azim-diffusion':
        xemit = theta[0];
        yemit = theta[1];
        Nph = theta[2];
        Nbg = theta[3];
        phi = theta[4];
        azim = phi;
        g2 = theta[6];
        break;
}

switch (excitation)
{

    case 'constant':
        P = 1;
        dP = 0;
        break;
    case 'cos2'
        P = 2*(1-m*sin(phi-alpha)*sin(phi-alpha))/K/(2-m);
        dP = -2*m*sin(2*phi-2*alpha)/K/(2-m);
        break;
    case 'cos4':
        P = 8*(1-m*(cos(phi-alpha)*cos(phi-alpha)+1).*sin(phi-alpha).^2)/K/(8-5*m);
        dP = -4*m*(2*sin(2*(phi-alpha))+sin(4*(phi-alpha)))/K/(8-5*m);
        break;
}
get_normalization(PupilMatrix,params);
get_field_matrix_derivatives(params);
int PSF = 0;
//PSFder will be passed by reference.
PSF = get_psfs_derivatives(FieldMatrix,FieldMatrixDerivatives,params);

int i,j,k,q;

for (i = 0; i< Mx;i++){
    for (j = 0; j< My;j++){
        for (k = 0; k< K;k++){
            mu[i][j][k] = 0;
        }
    }
}
for (i = 0; i< Mx;i++){
    for (j = 0; j< My;j++){
        for (k = 0; k< K;k++){
            for (q =0; q<numparams;q++){
               dmudtheta [i][j][k][q] =0;
            }
        }
    }
}

for (i = 0; i< Mx;i++){
    for (j = 0; j< My;j++){
        for (k = 0; k< K;k++){

            mu[i][j][k] = Nph*P[k]*PSF+Nbg/K;

            if (fitmodel=='xy'){
                dmudtheta [i][j][k][0] =Nph*P[k]*PSFder[i][j][0];
                dmudtheta [i][j][k][1] =Nph*P[k]*PSFder[i][j][1];
                dmudtheta [i][j][k][2] =P[k]*PSF;
                dmudtheta [i][j][k][3] =1/K;
            }
            if (fitmodel=='xyz'){
                dmudtheta [i][j][k][0] =Nph*P[k]*PSFder[i][j][0];
                dmudtheta [i][j][k][1] =Nph*P[k]*PSFder[i][j][1];
                dmudtheta [i][j][k][2] =Nph*P[k]*PSFder[i][j][2];
                dmudtheta [i][j][k][3] =P[k]*PSF;
                dmudtheta [i][j][k][4] =1/K;
            }
            if (fitmodel=='xy-azim'){
                dmudtheta [i][j][k][0] =Nph*P[k]*PSFder[i][j][0];
                dmudtheta [i][j][k][1] =Nph*P[k]*PSFder[i][j][1];
                dmudtheta [i][j][k][2] =P[k]*PSF;
                dmudtheta [i][j][k][3] =1/K;
                dmudtheta [i][j][k][4] =Nph*dP[k]*PSF+Nph*P[k]*PSFder[i][j][2];
            }
            if (fitmodel=='xy-azim-pola'){
                dmudtheta [i][j][k][0] =Nph*P[k]*PSFder[i][j][0];
                dmudtheta [i][j][k][1] =Nph*P[k]*PSFder[i][j][1];
                dmudtheta [i][j][k][2] =P[k]*PSF;
                dmudtheta [i][j][k][3] =1/K;
                dmudtheta [i][j][k][4] =Nph*dP[k]*PSF+Nph*P[k]*PSFder[i][j][2];
                dmudtheta [i][j][k][5] =Nph*P[k]*PSFder[i][j][3];
            }
            if (fitmodel=='xyz-azim-pola'){
                dmudtheta [i][j][k][0] =Nph*P[k]*PSFder[i][j][0];
                dmudtheta [i][j][k][1] =Nph*P[k]*PSFder[i][j][1];
                dmudtheta [i][j][k][2] =Nph*P[k]*PSFder[i][j][2];
                dmudtheta [i][j][k][3] =P[k]*PSF;
                dmudtheta [i][j][k][4] =1/K;
                dmudtheta [i][j][k][5] =Nph*dP[k]*PSF+Nph*P[k]*PSFder[i][j][3];
                dmudtheta [i][j][k][6] =Nph*P[k]*PSFder[i][j][4];
            }
            if (fitmodel=='xy-azim-pola-diffusion'){
                dmudtheta [i][j][k][0] =Nph*P[k]*PSFder[i][j][0];
                dmudtheta [i][j][k][1] =Nph*P[k]*PSFder[i][j][1];
                dmudtheta [i][j][k][2] =P[k]*PSF;
                dmudtheta [i][j][k][3] =1/K;
                dmudtheta [i][j][k][4] =Nph*dP[k]*PSF+Nph*P[k]*PSFder[i][j][2];
                dmudtheta [i][j][k][5] =Nph*P[k]*PSFder[i][j][3];
                dmudtheta [i][j][k][6] =Nph*P[k]*PSFder[i][j][4];
            }
        }
    }
}

}
