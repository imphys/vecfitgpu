#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <tiffio.h>
#include "paramsdata.h"
#include <time.h>
#include "fft.h"


double *cztfunc2D ( paramsdata *params, double &datain [][][][],double _Complex *Amt[][],double _Complex *Bmt[][],double _Complex *Dmt[][],int sizex,int sizey)
{
    int  K = sizex;
    int  N = params -> cztN;
    int  M = params -> cztM;
    int  L = params -> cztL;
    int permuteL = K*L*2*3;
    int i, j, k,q,x;
    double _Complex datain [K][L][2][3];
    double _Complex cztin  [K][L][2][3];
    double _Complex temp   [K][L][2][3];
    double _Complex temp2   [permuteL];
    double _Complex cztout [K][L][2][3];
    double _Complex dataout[K][M][2][3];
    double _Complex dataout2[M][K][2][3];
    double *refoutreal = malloc(permuteL * sizeof(double));
	double *refoutimag = malloc(permuteL * sizeof(double));

    for (i = 0; i<K;i++){
        for (j = 0; j<L;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){

                    cztin[i][j][k][q] = 0;
                }
            }
        }
    }

    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){

                    cztin[i][j][k][q] = Amt[i][j]*datain[i][j][k][q];
                }
            }
        }
    }
    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){
                    for (x=0;x<permuteL;x++){

                        refoutreal[x] = creal(cztin[i][j][k][q]);
                        refoutimag[x] = cimag(cztin[i][j][k][q]);
                    }
                }
            }
        }
    }

    for (i = 1; i<K+1;i++){
        for (k=1;k<2+1;k++){
            for (q=1;q<3+1;q++){
                Fft_transform(refoutreal,refoutimag,N*i*k*q);
            }
        }
    }

    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){
                    for (x=0;x<permuteL;x++){

                        creal(cztin[i][j][k][q])=refoutreal[x] ;
                        cimag(cztin[i][j][k][q])=refoutimag[x] ;
                    }
                }
            }
        }
    }
    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){

                    temp[i][j][k][q] = Dmt[i][j]*(creal(cztin[i][j][k][q])+cimag(cztin[i][j][k][q])*i);
                }                }
            }
        }
    }
    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){
                    for (x=0;x<permuteL;x++){

                        refoutreal[x] = creal(temp[i][j][k][q]);
                        refoutimag[x] = cimag(temp[i][j][k][q]);
                    }
                }
            }
        }
    }
    for (i = 1; i<K+1;i++){
        for (k=1;k<2+1;k++){
            for (q=1;q<3+1;q++){
                Fft_inverseTransform(refoutreal,refoutimag,N*i*k*q);
              }
        }
    }

    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){
                    for (x=0;x<permuteL;x++){

                        creal(temp[i][j][k][q])=refoutreal[x] ;
                        cimag(temp[i][j][k][q])=refoutimag[x] ;
                    }
                }
            }
        }
    }
    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){

                    cztout[i][j][k][q] = (creal(temp[i][j][k][q])+cimag(temp[i][j][k][q])*i);
                }
            }
        }
    }

    for (i = 0; i<K;i++){
        for (j = 0; j<M;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){

                    dataout[i][j][k][q] = Bmt[i][j]*cztout[i][j][k][q] ;
                }
            }
        }
    }

/// the following code is to implement the permute function

    for (i = 0; i<K;i++){
        for (j = 0; j<M;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){
                    for (x=0;x<permuteL;x++){

                    temp2[x] = dataout[i][j][k][q] ;
                    }
                }
            }
        }
    }

    for (i = 0; i<M;i++){
        for (j = 0; j<K;j++){
            for (k=0;k<2;k++){
                for (q=0;q<3;q++){
                    for (x=0;x<permuteL;x++){

                    dataout2[i][j][k][q]= temp2[x];
                    }
                }
            }
        }
    }
    return dataout2;

}
