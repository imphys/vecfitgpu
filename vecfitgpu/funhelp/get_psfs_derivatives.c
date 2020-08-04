#include "paramsdata.h"
#include <complex.h>
#include <tgmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void get_psfs_derivatives (paramsdata *params,double _Complex &FieldMatrix [M][M][2][3],
                                   double _Complex &FieldMatrixDerivatives [M][M][3][2][3],double &PSF [Mx][My],double &PSFder [Mx][My][100])
{                                                                                             /// the 100 just to be safe from access errors
int Npupil = params -> Npupil;
int Mx = params -> Mx;
int M = params -> cztM;
int My = params -> My;
double azim = params -> azim;
double g2 = params -> g2;
double normint_free = params -> normint_free;
double pola = params -> pola;
int templen = M*M*3*3*2;
double _Complex tempFieldMatrixDerivatives1[templen];
double _Complex tempFieldMatrixDerivatives[M][M][2][3][3];
double _Complex tmpPSFderivatives[M][M][3];




double dipor[3];
double diporDerivatives[3][3];

dipor[0]= sin(pola)*cos(azim);
dipor[1]= sin(pola)*sin(azim);
dipor[2]= cos(pola);

char *azim1 = "azim";
if (strstr(params -> fitmodel ,azim1)!= NULL){
    diporDerivatives[0][0] = -sin(pola)*sin(azim);
    diporDerivatives[0][1] = sin(pola)*cos(azim);
    diporDerivatives[0][2] = 0;
}

char *pola1 = "pola";
if (strstr(params -> fitmodel ,pola1)!= NULL){
    diporDerivatives[1][0] = cos(pola)*cos(azim);
    diporDerivatives[1][1] = cos(pola)*sin(azim);
    diporDerivatives[1][2] = -sin(pola);
}
int i, j, k,q,x,y;

if(!(strcmp (params -> dipoletype, "free"))){

    for (i = 0;i<Mx;i++){
        for (j = 0;j<My;j++){
            for (k = 0;q<2;k++){
                for (q = 0;k<3;q++){

                        PSF[i][j] += (1/3)*pow(abs(FieldMatrix[i][j][k][q]),2);
                }
            }
        }
    }
    x = 0;
    for (i = 0; i<Mx;i++){
        for (j = 0; j<My;j++){
            for (y = 0; y<3;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){


                        tempFieldMatrixDerivatives1[x] = FieldMatrixDerivatives[i][j][y][k][q] ;
			x++;

                    }
                }
            }
        }
    }
    x = 0;
    for (i = 0; i<Mx;i++){
        for (j = 0; j<My;j++){
            for (y = 0; y<2;y++){
                for (k=0;k<3;k++){
                    for (q=0;q<3;q++){


                        tempFieldMatrixDerivatives[i][j][y][k][q] = tempFieldMatrixDerivatives1[x];
			x++;	


                    }
                }
            }
        }
    }
    for (i = 0; i<Mx;i++){
        for (j = 0; j<My;j++){
            for (y = 0; y<2;y++){
                for (k=0;k<3;k++){
                    for (q=0;q<3;q++){

                            tmpPSFderivatives[i][j][q] += (2/3)*creal((creal(FieldMatrix[i][j][y][k])-cimag(FieldMatrix[i][j][y][k])*I)*tempFieldMatrixDerivatives[i][j][y][k][q]) ;
                            PSFder [i][j][q] = tmpPSFderivatives[i][j][q];
    /// this step is redundant PSFder could have been assigned directly

                    }
                }
            }
        }
    }


}
double Ex[M][M];
double Ey[Mx][My];
double Exder[Mx][My][3];
double Eyder[Mx][My][3];

if(!(strcmp (params -> dipoletype, "fixed"))){
    for (i = 0;i<M;i++){
        for (j = 0;j<M;j++){
            Ex[i][j] = dipor[0]*FieldMatrix[i][j][0][0]+dipor[1]*FieldMatrix[i][j][0][1]+dipor[2]*FieldMatrix[i][j][0][2];
            Ey[i][j] = dipor[0]*FieldMatrix[i][j][1][0]+dipor[1]*FieldMatrix[i][j][1][1]+dipor[2]*FieldMatrix[i][j][1][2];
            PSF [i][j] = pow(abs(Ex[i][j]),2)+pow(abs(Ey[i][j]),2);
        }
    }

    for (i = 0;i<Mx;i++){
        for (j = 0;j<My;j++){
            for(k= 0; k<3;k++){ // the loop here is up to 3 only to not make an access error with field matrix drevative
            Exder[i][j][k] = dipor[0]*FieldMatrixDerivatives[i][j][k][0][0]+dipor[1]*FieldMatrixDerivatives[i][j][k][0][1]+dipor[2]*FieldMatrixDerivatives[i][j][k][0][2];
            Eyder[i][j][k] = dipor[0]*FieldMatrixDerivatives[i][j][k][1][0]+dipor[1]*FieldMatrixDerivatives[i][j][k][1][1]+dipor[2]*FieldMatrixDerivatives[i][j][k][1][2];
            }
        }
    }
    if (strstr(params -> fitmodel ,azim1)!= NULL){
        for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){
                for(k= 0; k<3;k++){
			//sizeof(Exder[2]) = 3
                    Exder[i][j][(sizeof(Exder[2]) / sizeof(double))+1] = diporDerivatives[0][0]*FieldMatrix[i][j][0][0]+diporDerivatives[0][0]*FieldMatrix[i][j][0][1]+diporDerivatives[0][0]*FieldMatrix[i][j][0][2];
                    Exder[i][j][(sizeof(Exder[2]) / sizeof(double))+1] = diporDerivatives[0][0]*FieldMatrix[i][j][1][0]+diporDerivatives[0][0]*FieldMatrix[i][j][1][1]+diporDerivatives[0][0]*FieldMatrix[i][j][1][2];
                }
            }
        }
    }

    if (strstr(params -> fitmodel ,pola1)!= NULL){
        for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){
                for(k= 0; k>3;k++){
                    Exder[i][j][(sizeof(Exder[2]) / sizeof(double))+1] = diporDerivatives[1][0]*FieldMatrix[i][j][0][0]+diporDerivatives[1][0]*FieldMatrix[i][j][0][1]+diporDerivatives[1][0]*FieldMatrix[i][j][0][2];
                    Eyder[i][j][(sizeof(Eyder[2]) / sizeof(double))+1] = diporDerivatives[1][0]*FieldMatrix[i][j][1][0]+diporDerivatives[1][0]*FieldMatrix[i][j][1][1]+diporDerivatives[1][0]*FieldMatrix[i][j][1][2];
                }
            }
        }
    }
        for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){
                for(k= 0; k>3;k++){
                    PSFder[i][j][k] = 2*creal((creal(Ex[i][j])-cimag(Ex[i][j])*I)*Exder[i][j][k])+2*creal((creal(Ey[i][j])-cimag(Ey[i][j])*I)*Eyder[i][j][k]);

                }
            }
        }


}
double FreePSF [Mx][My];
double FixedPSF [Mx][My];
double FreePSFder [Mx][My][100];

if(!(strcmp (params -> dipoletype, "diffusion"))){
    for (i = 0;i<Mx;i++){
        for (j = 0;j<My;j++){
            for (k = 0;k<2;k++){
                for (q = 0;q<3;q++){

                        FreePSF[i][j] += (1/3)*pow(abs(FieldMatrix[i][j][k][q]),2);
                        FreePSF[i][j] = FreePSF[i][j]/normint_free; //
                }
            }
        }
    }
    x = 0;
    for (i = 0; i<Mx;i++){
        for (j = 0; j<My;j++){
            for (y = 0; y<3;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){


                            tempFieldMatrixDerivatives1[x] = FieldMatrixDerivatives[i][j][y][k][q] ;
				x++;
                    }
                }
            }
        }
    }
    x = 0;
    for (i = 0; i<Mx;i++){
        for (j = 0; j<My;j++){
            for (y = 0; y<2;y++){
                for (k=0;k<3;k++){
                    for (q=0;q<3;q++){


                            tempFieldMatrixDerivatives[i][j][y][k][q] = tempFieldMatrixDerivatives1[x] ;
				x++;
                    }
                }
            }
        }
    }
    for (i = 0; i<Mx;i++){
        for (j = 0; j<My;j++){
            for (y = 0; y<2;y++){
                for (k=0;k<3;k++){
                    for (q=0;q<3;q++){

                            tmpPSFderivatives[i][j][q] += (2/3)*creal((creal(FieldMatrix[i][j][y][k])-cimag(FieldMatrix[i][j][y][k])*I)*tempFieldMatrixDerivatives[i][j][y][k][q]) ;
                            FreePSFder [i][j][q] = tmpPSFderivatives[i][j][q];
    /// this step is redundant PSFder could have been assigned directly

                    }
                }
            }
        }
    }

    if((!(strcmp (params -> fitmodel, "xyz")))||(!(strcmp (params -> fitmodel, "xy-azim")))) {

        for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){

                FreePSFder[i][j][2] = 0;
            }
        }
    }else if ((!(strcmp (params -> fitmodel, "xy-azim-pola")))){

        for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){

                FreePSFder[i][j][2] = 0;
                FreePSFder[i][j][3] = 0;
                FreePSFder[i][j][4] = 0;
            }
        }

    }
        for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){
                for (k = 0;k<3;k++){

                    FreePSFder[i][j][k] = FreePSFder[i][j][k]/normint_free;
                }
            }
        }

        for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){
                Ex[i][j] = dipor[0]*FieldMatrix[i][j][0][0]+dipor[1]*FieldMatrix[i][j][0][1]+dipor[2]*FieldMatrix[i][j][0][2];
                Ey[i][j] = dipor[0]*FieldMatrix[i][j][1][0]+dipor[1]*FieldMatrix[i][j][1][1]+dipor[2]*FieldMatrix[i][j][1][2];
                FixedPSF [i][j] = pow(abs(Ex[i][j]),2)+pow(abs(Ey[i][j]),2);
                FixedPSF [i][j] = FixedPSF [i][j]/normint_free;
            }
        }

        for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){
                for(k= 0; k<3;k++){
                    Exder[i][j][k] = dipor[0]*FieldMatrixDerivatives[i][j][k][0][0]+dipor[1]*FieldMatrixDerivatives[i][j][k][0][1]+dipor[2]*FieldMatrixDerivatives[i][j][k][0][2];
                    Eyder[i][j][k] = dipor[0]*FieldMatrixDerivatives[i][j][k][1][0]+dipor[1]*FieldMatrixDerivatives[i][j][k][1][1]+dipor[2]*FieldMatrixDerivatives[i][j][k][1][2];
                }
            }
        }
    char *xyz = "xyz";
    if (strstr(params -> fitmodel ,azim1)!= NULL){
        if (strstr(params -> fitmodel ,xyz)!= NULL){
            for (i = 0;i<Mx;i++){
                for (j = 0;j<My;j++){
                    
                        Exder[i][j][3] = diporDerivatives[0][0]*FieldMatrix[i][j][0][0]+diporDerivatives[0][1]*FieldMatrix[i][j][0][1]+diporDerivatives[0][2]*FieldMatrix[i][j][0][2];
                        Eyder[i][j][3] = diporDerivatives[0][0]*FieldMatrix[i][j][1][0]+diporDerivatives[0][1]*FieldMatrix[i][j][1][1]+diporDerivatives[0][2]*FieldMatrix[i][j][1][2];
                   
                }
            }

        }else {
            for (i = 0;i<Mx;i++){
                for (j = 0;j<My;j++){
                        Exder[i][j][2] = diporDerivatives[0][0]*FieldMatrix[i][j][0][0]+diporDerivatives[0][1]*FieldMatrix[i][j][0][1]+diporDerivatives[0][2]*FieldMatrix[i][j][0][2];
                        Eyder[i][j][2] = diporDerivatives[0][0]*FieldMatrix[i][j][1][0]+diporDerivatives[0][1]*FieldMatrix[i][j][1][1]+diporDerivatives[0][2]*FieldMatrix[i][j][1][2];
                    
                }
            }

        }
    }

    if (strstr(params -> fitmodel ,pola1)!= NULL){
        if (strstr(params -> fitmodel ,xyz)!= NULL){
            for (i = 0;i<Mx;i++){
                for (j = 0;j<My;j++){
                        Exder[i][j][4] = diporDerivatives[1][0]*FieldMatrix[i][j][0][0]+diporDerivatives[1][1]*FieldMatrix[i][j][0][1]+diporDerivatives[1][2]*FieldMatrix[i][j][0][2];
                        Eyder[i][j][4] = diporDerivatives[1][0]*FieldMatrix[i][j][1][0]+diporDerivatives[1][1]*FieldMatrix[i][j][1][1]+diporDerivatives[1][2]*FieldMatrix[i][j][1][2];

                }
            }

        }else {
            for (i = 0;i<Mx;i++){
                for (j = 0;j<My;j++){
                        Exder[i][j][3] = diporDerivatives[1][0]*FieldMatrix[i][j][0][0]+diporDerivatives[1][1]*FieldMatrix[i][j][0][1]+diporDerivatives[1][2]*FieldMatrix[i][j][0][2];
                        Eyder[i][j][3] = diporDerivatives[1][0]*FieldMatrix[i][j][1][0]+diporDerivatives[1][1]*FieldMatrix[i][j][1][1]+diporDerivatives[1][2]*FieldMatrix[i][j][1][2];
                }
            }

        }
    }

    for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){
                for(k= 0; k<5;k++){
                    FixedPSFder[i][j][K] = 2*creal((creal(Ex[i][j])-cimag(Ex[i][j])*I)*Exder[i][j][k])+2*creal((creal(Ey[i][j])-cimag(Ey[i][j])*I)*Eyder[i][j][k]);
                    FixedPSFder[i][j][K] = FixedPSFder[i][j][K]/normint_fixed;
                        PSF[i][j] = (1-g2)*FreePSF[i][j]+g2*FixedPSF[i][j];
                        PSFder[i][j][k] = (1-g2)*FreePSFder[i][j][k]+g2*FixedPSFder[i][j][k];
                }
            }
        }

    if(!(strcmp (params -> fitmodel, "xy-azim-pola-diffusion"))){
        for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){

              PSFder[i][j][4] = -FreePSF[i][j]+FixedPSF[i][j];
            }
        }
    }else if(!(strcmp (params -> fitmodel, "xy-azim-diffusion"))){
        for (i = 0;i<Mx;i++){
            for (j = 0;j<My;j++){

              PSFder[i][j][3] = -FreePSF[i][j]+FixedPSF[i][j];
            }
        }
    }


}

}
