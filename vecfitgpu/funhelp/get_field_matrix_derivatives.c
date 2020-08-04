#include "paramsdata.h"
#include <complex.h>
#include <tgmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void get_field_matrix_derivatives (paramsdata *params,double _Complex &FieldMatrix [][][][],double _Complex &FieldMatrixDerivatives [][][][][])
{

//char fitmodel = params -> fitmodel;
//char ztype = params -> ztype;
double _Complex ***wavevector = params -> wavevector;
double _Complex ****PupilMatrix = params -> PupilMatrix;
double _Complex **wavevectorzmed = params -> wavevectorzmed;
double  xemit = params -> xemit;
double  yemit = params -> yemit;
double  zemit = params -> zemit;
int     Npupil = params -> Npupil;
int     cztM = params -> cztM;


double _Complex Ax = params -> Axmt;
double _Complex Bx = params -> Bxmt;
double _Complex Dx = params -> Dxmt;
double _Complex Ay = params -> Aymt;
double _Complex By = params -> Bymt;
double _Complex Dy = params -> Dymt;

int i, j, k,q,x;

double _Complex Wpos[Npupil][Npupil];
double _Complex PositionPhaseMask[Npupil][Npupil];
double _Complex PupilFunction[Npupil][Npupil][2][3];
double _Complex PupilFunctionDerivatives[Npupil][Npupil][3][2][3];

    for (i = 0;i<Npupil;i++){
        for (j = 0;j<Npupil;j++){

                Wpos[i][j] = wavevector[i][j][0]*xemit+wavevector[i][j][1]*yemit+wavevector[i][j][2]*zemit;
        }
    }

    for (i = 0;i<Npupil;i++){
        for (j = 0;j<Npupil;j++){

               PositionPhaseMask[i][j] =  exp(-I*Wpos[i][j]) ;
        }
    }

    for (i = 0;i<Npupil;i++){
        for (j = 0;j<Npupil;j++){
            for (k = 0;k<2;k++){
                for (q = 0;k<3;q++){

                        PupilFunction[i][j][k][q] =  PositionPhaseMask[i][j]*PupilMatrix[i][j][k][q];
                }
            }
        }
    }

    for (i = 0;i<Npupil;i++){
        for (j = 0;j<Npupil;j++){
            for (x = 0;x<2;x++){
                for (k = 0;k<2;k++){
                    for (q = 0;k<3;q++){

                        PupilFunctionDerivatives[i][j][x][k][q] = -I*wavevector[i][j][x]*PupilFunction[i][j][k][q];
                    }
                }
            }
        }
    }
    char *xyz = "xyz";
    if (strstr(params -> fitmodel ,xyz)!= NULL){
        if(!(strcmp (params -> ztype, "stage"))){

            for (i = 0;i<Npupil;i++){
                for (j = 0;j<Npupil;j++){
                    for (k = 0;k<2;k++){
                        for (q = 0;k<3;q++){

                                PupilFunctionDerivatives[i][j][2][k][q] = -I*wavevectorzmed[i][j]*PupilFunction[i][j][k][q];
                        }
                    }
                }
            }
        }else{
            for (i = 0;i<Npupil;i++){
                for (j = 0;j<Npupil;j++){
                    for (k = 0;k<2;k++){
                        for (q = 0;k<3;q++){

                                PupilFunctionDerivatives[i][j][2][k][q] = -I*wavevector[i][j][2]*PupilFunction[i][j][k][q];
                        }
                    }
                }
            }
        }
    }

double _Complex IntermediateImage [Npupil][Npupil][3][2][3];

int sizex,sizey;
sizex = Npupil;
sizey = Npupil;
IntermediateImage = cztfunc2D(PupilFunction,Ay,By,Dy,params,sizex,sizey);
sizex = cztM;
sizey = Npupil;
FieldMatrix = cztfunc2D(IntermediateImage,Ax,Bx,Dx,params,sizex,sizey);
sizex = Npupil;
sizey = Npupil;
IntermediateImage = cztfunc3D(PupilFunctionDerivatives,Ay,By,Dy,params,sizex,sizey);
sizex = cztM;
sizey = Npupil;
FieldMatrixDerivatives = cztfunc3D(IntermediateImage,Ax,Bx,Dx,params,sizex,sizey);




}
