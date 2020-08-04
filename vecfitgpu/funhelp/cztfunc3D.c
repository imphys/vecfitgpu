#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <tiffio.h>
#include "paramsdata.h"
#include <time.h>
#include "fft.h"


double *cztfunc3D ( paramsdata *params, double &datain [][][][][],double _Complex *Amt[][],double _Complex *Bmt[][],double _Complex *Dmt[][],int sizex,int sizey)
{
    int  K = sizex;
    int  N = params -> cztN;
    int  M = params -> cztM;
    int  L = params -> cztL;
    char fitmodel = params -> fitmodel;
    int dim =0;

    if (fitmodel == 'xyz'){
        dim = 3;
    }else if (fitmodel == 'xy'){
        dim = 2;
    }
    int permuteL = K*L*dim*2*3;
    int i, j, k,q,x,y;
    double _Complex datain [K][L][dim][2][3];
    double _Complex cztin  [K][L][dim][2][3];
    double _Complex temp   [K][L][dim][2][3];
    double _Complex temp2   [permuteL];
    double _Complex cztout [[K][L][dim][2][3];
    double _Complex dataout[K][L][dim][2][3];
    double _Complex dataout2[M][K][dim][2][3];
    double *refoutreal = malloc(permuteL * sizeof(double));
	double *refoutimag = malloc(permuteL * sizeof(double));

    for (i = 0; i<K;i++){
        for (j = 0; j<L;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){

                    cztin[i][j][y][k][q] = 0;
                    }
                }
            }
        }
    }

    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){

                        cztin[i][j][y][k][q] = Amt[i][j]*datain[i][j][y][k][q];
                    }
                }
            }
        }
    }
    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){
                        for (x=0;x<permuteL;x++){

                        refoutreal[x] = creal(cztin[i][j][y][k][q]);
                        refoutimag[x] = cimag(cztin[i][j][y][k][q]);
                        }
                    }
                }
            }
        }
    }
    for (i = 1; i<K+1;i++){
        for (y = 1; y<dim+1;y++){
            for (k=1;k<2+1;k++){
                for (q=1;q<3+1;q++){
                    Fft_transform(refoutreal,refoutimag,N*i*y*k*q);

                }
            }
        }
    }
    for (i = 1; i<
    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){
                        for (x=0;x<permuteL;x++){

                        creal(cztin[i][j][y][k][q])= refoutreal[x] ;
                        cimag(cztin[i][j][y][k][q])= refoutimag[x] ;
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){

                        temp[i][j][y][k][q] = Dmt[i][j]*(creal(cztin[i][j][y][k][q])+cimag(cztin[i][j][y][k][q])*i);
                    }

                }
            }
        }
    }
    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){
                        for (x=0;x<permuteL;x++){

                        refoutreal[x] = creal(temp[i][j][y][k][q]);
                        refoutimag[x] = cimag(temp[i][j][y][k][q]);
                        }
                    }
                }
            }
        }
    }

    for (i = 1; i<K+1;i++){
        for (y = 1; y<dim+1;y++){
            for (k=1;k<2+1;k++){
                for (q=1;q<3+1;q++){
                    Fft_inverseTransform(refoutreal,refoutimag,N*i*y*k*q);

                }
            }
        }
    }

for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){
                        for (x=0;x<permuteL;x++){

                        creal(temp[i][j][y][k][q])= refoutreal[x] ;
                        cimag(temp[i][j][y][k][q])= refoutimag[x] ;
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i<K;i++){
        for (j = 0; j<N;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){

                        cztout[i][j][y][k][q] = (creal(temp[i][j][y][k][q])+cimag(temp[i][j][k][q])*i);
                    }
                }
            }
        }
    }

    for (i = 0; i<K;i++){
        for (j = 0; j<M;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){

                        dataout[i][j][y][k][q] = Bmt[i][j]*cztout[i][j][y][k][q] ;
                    }
                }
            }
        }
    }

/// the following code is to implement the permute function

    for (i = 0; i<K;i++){
        for (j = 0; j<M;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){
                        for (x=0;x<permuteL;x++){

                            temp2[x] = dataout[i][j][y][k][q] ;
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i<M;i++){
        for (j = 0; j<K;j++){
            for (y = 0; y<dim;y++){
                for (k=0;k<2;k++){
                    for (q=0;q<3;q++){
                        for (x=0;x<permuteL;x++){

                            dataout2[i][j][y][k][q]= temp2[x];
                        }
                    }
                }
            }
        }
    }
    return dataout2;

}
