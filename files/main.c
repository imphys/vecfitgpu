#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
//#include "set_parameters_simulation.h"

typedef struct paramsdata{
    int flg_parallel;
    int flg_writemat;
    int flg_showplot;
    int flg_showconv;
    int flg_gpu;
    int K;
    int Ncfg;
    int xemit;
    int yemit;
    int zemit;
    int Npupil;
    int pixelsize;
    int Mx;
    int My;
    int Mz;
    int xrange;
    int yrange;
    int numparams;
    int Nitermax;
    int fwd;
    int depth;
    int lambda;
    int lambdacentral;
    int aberrationcorrected;
    int ringradius;
    int doelevels;
    int doephasedepth;
    int readnoisestd;
    int readnoisevariance;
    int cztN;
    int cztM;
    int cztL;
    int debugmode;
    int phtrue;
    int bgtrue;
    int varfit;

    float azim;
    float m;
    float alpha;
    float tollim;
    float NA;
    float refmed;
    float refcov;
    float refimm;
    float refimmnom;
    float pola;
    float welldepth;
    float g2;

    float * aberrationsoffset;
    float phi[1000];
    float theta[1000];
    float allalpha[1000];
    float Duxstr[1000];
    float Duystr[1000];


    char fitmodel[50];
    char excitation[30];
    char ztype[50];
    char dipoletype[50];
    char doetype[50];

    int zrange[2];
    int zspread[2];
    int lambdaspread[2];
    int aberrations[9][3];
    int *zonefunction;
}paramsdata;

void get_pupil_matrix(paramsdata *params)
{
//This function calculates the pupil matrix Q_{jk}, which gives the j-th
//electric field component proportional to the k-th dipole vector
//component.

//parameters: NA, refractive indices of medium, cover slip, immersion fluid,
//wavelength (in nm), sampling in pupil

double NA = params -> NA;
double refmed = params -> refmed;
double refcov = params -> refcov;
double refimm = params -> refimm;
double refimmnom = params -> refimmnom;
double lambda = params -> lambda;
double Npupil = params -> Npupil;

//pupil radius (in diffraction units) and pupil coordinate sampling
double PupilSize = 1.0;
double DxyPupil = 2*PupilSize/Npupil;
double XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;

//[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);

}



void set_parameters_simulation(paramsdata *params){
    params->K = 1;
    params->m = 0.95;
    params->alpha = 3.141592653589793;

    //params.Ncfg = size(allPSFs,3)/K;
    //params.allPSFs = allPSFs;
    //params.allalpha = allalpha;
    //params.Duxstr = Duxstr;
    //params.Duystr = Duystr;

    for(int i = 0;i<1001; i++)
    {
        params->allalpha[i] = params->alpha;
    }

//    printf("Specify fitmodel\n");
//    scanf("%s",params->fitmodel);
    strcpy(params->fitmodel, "xyz-azim-pola");

    if(!strcmp(params->fitmodel, "xy"))
    {
        params->numparams = 4;
    }
    else if(!strcmp(params->fitmodel,"xy"))
    {
        params->numparams = 4;

    }
    else if(!strcmp(params->fitmodel,"xyz"))
    {
        params->numparams = 5;
    }
    else if(!strcmp(params->fitmodel,"xy-azim"))
    {
        params->numparams = 5;
    }
    else if(!strcmp(params->fitmodel, "xy-azim-pola"))
    {
        params->numparams = 6;
    }
    else if(!strcmp(params->fitmodel,"xyz-azim-pola"))
    {
        params->numparams = 7;
    }
    else if(!strcmp(params->fitmodel,"xy-azim-diffusion"))
    {
        params->numparams = 6;
    }
    else if(!strcmp(params->fitmodel,"xy-azim-pola-diffusion"))
    {
        params->numparams = 7;
    }
    else
        {
        printf("Invalid fitmodel\n");
        abort();
    }


//MLE/fitting parameters
    params->Nitermax = 50;
    params->tollim = 1e-5;
    params->varfit = 0;
    params->readnoisevariance = 0;


//PSF/optical parameters
    params->NA = 1.49;
    params->refmed = 1.518;
    params->refcov = 1.518;
    params->refimm = 1.518;
    params->refimmnom = params->refcov;
    params->fwd = 140e3;
    params->depth = 0;
    params->zrange[0] = 0;
    params->zrange[1] = 0;
    params->zspread[0] = 0;
    params->zspread[1] = 0;
    strcpy(params->ztype, "stage"); //'medium'
    params->lambda = 552;
    params->lambdacentral = 552;
    params->lambdaspread[0] = 552;
    params->lambdaspread[1] = 552;
    params->xemit = 0.0;
    params->yemit = 0.0;
    params->zemit = 0;
    params->Npupil = 32;
    params->pixelsize = 80;
    //params.Mx = size(allPSFs,1);
    //params.My = size(allPSFs,2);
    params->Mz = 1;
    //params->xrange = params.pixelsize*params.Mx/2;
    //params.yrange = params.pixelsize*params.My/2;


    //excitation pattern
    if(params->K == 1)
    {
        strcpy(params->excitation, "constant");
    }
    else if(params->K > 1)
    {
        strcpy(params->excitation, "cos2");
        //strcpy(params->excitation, "cos4"); //two-photon excitation
    }

    int zcheck;
    int zmin;
    //sanity check on position emitter w.r.t. cover slip
    if(!strcmp(params->ztype, "stage"))
    {
        zcheck = params->depth+params->zemit;
    }
    if(!strcmp(params->ztype,'medium'))
    {
        zmin = params->zrange[0];
        zcheck = zmin+params->depth+params->zemit;
    }
    if (zcheck < 0)
    {
        printf("Warning! Emitter must be above the cover slip:\n");
        printf("Adjust parameter settings for physical results.\n");
    }

    //sanity check on refractive index values
    if (params->NA > params->refimm)
    {
        printf('Warning! Refractive index immersion medium cannot be smaller than NA.\n');
    }
    if(params->NA > params->refcov)
    {
        printf('Warning! Refractive index cover slip cannot be smaller than NA.\n');
    }

    //parameters needed for fixed dipole PSF only: emitter/absorber dipole
    //orientation (characterized by angles pola and azim)
    strcpy(params->dipoletype, 'diffusion');
    printf("%s", params->dipoletype);
    //strcpy(params->dipoletype, 'free');
    //strcpy(params->dipoletype, 'fixed');
    params->pola = 90.0*M_PI/180;
    params->azim = 0.0*M_PI/180;


    //diffusion coefficient

    float g2;
    float eps = 2.2204e-16;
    printf("%f", eps);

    float welldepth = 1/eps;
    g2 = (3 + pow(welldepth, 2) - 3 * welldepth * (cosh(welldepth)/sinh(welldepth)))/pow(welldepth, 2); //coth(welldepth) replaced by cosh/sinh
    params->welldepth = welldepth;
    params->g2 = g2;

    //aberrations (Zernike orders [n1,m1,A1,n2,m2,A2,...] with n1,n2,... the
    //radial orders, m1,m2,... the azimuthal orders, and A1,A2,... the Zernike
    //coefficients in lambda rms, so 0.072 means diffraction limit)
    params->aberrationcorrected = 0; //0 instead of false

    for(int i = 0; i<9;i++)
    {
        for(int j = 0; j<3; j++)
        {
            params->aberrations[i][j] = 0;
        }
    }

    //params->aberrationsoffset = [];
    params->aberrations[3][1] = -1;
    params->aberrations[4][1] = 1;
    params->aberrations[0][0] = 2;
    params->aberrations[1][0] = 2;
    params->aberrations[2][0] = 2;
    params->aberrations[1][1] = -2;
    params->aberrations[2][1] = 2;
    params->aberrations[8][1] = -2;
    params->aberrations[3][0] = 3;
    params->aberrations[4][0] = 3;
    params->aberrations[6][0] = 3;
    params->aberrations[7][0] = 3;
    params->aberrations[6][1] = -3;
    params->aberrations[7][1] = 3;
    params->aberrations[5][0] = 4;
    params->aberrations[8][0] = 4;

    for(int i = 0; i<9; i++)
    {
        params->aberrations[i][2] = params->aberrations[i][2]*params->lambdacentral;

    }
    for(int i = 0; i<9;i++)
    {
        for(int j = 0; j<3; j++)
        {
            printf("%d\t", params->aberrations[i][j]);
        }
        printf("\n");
    }

    //DOE/SLM
    //strcpy(params->doetype, 'none');
    strcpy(params->doetype, 'vortex');
    params->ringradius = 1;
    params->doelevels = 32;
    params->zonefunction = params->aberrations;
    params->doephasedepth = 1*params->lambdacentral;

    //Fit model parameters: signal photon count, background photons/pixel, read
    //noise variance for sCMOS camera's, fit model, output labels depending on
    //fit model
    //params.signalphotoncount = 1000;
    //params.backgroundphotoncount = 10;
    params->readnoisestd = 0;

    //calculate auxiliary vectors for chirpz
    int Kx, Ky;
    float PupilSize, ImageSizex, ImageSizey;
    float Ax, Ay, Bx, By, Dx, Dy;


    Kx = params->Mx;
    Ky = params->My;
    PupilSize = 1.0;
    ImageSizex = params->xrange*params->NA/params->lambda;
    ImageSizey = params->yrange*params->NA/params->lambda;

    //[Ax,Bx,Dx] = prechirpz(PupilSize, ImageSizex, params->Npupil, Kx);
    //[Ay,By,Dy] = prechirpz(PupilSize, ImageSizey, params->Npupil, Ky);

}

float prechirpz(xsize, qsize, N, M)
{
    int L;
    float sigma;
    float _Complex Afac, Bfac, Gfac, sqW, W;
    L = N + M - 1;
    sigma = 2 * M_PI * xsize * qsize / N / M;
    Afac = cexpf(2 * _Complex_I * sigma * (1 - M));
    Bfac = cexpf(2 * _Complex_I * sigma * (1 - N));
    sqW  = cexpf(2* _Complex_I * sigma);
    W = sqW * sqW;
    Gfac = (2*xsize/N)*cexpf(_Complex_I * sigma * (1-N) * (1-M));

    float _Complex *Utmp, *Vtmp, *A, *B, *D;

    Utmp = (float _Complex*) calloc(N, sizeof(float _Complex));
    A = (float _Complex*) calloc(N, sizeof(float _Complex));

    A[0] = 1.0;
    Utmp[0] = sqW * Afac;
    for(int i = 1; i < N; i++)
    {
        A[i] = Utmp[i-1]*A[i-1];
        Utmp[i] = Utmp[i-1]*W;
    }


    free(Utmp);


    Utmp = (float _Complex*) calloc(M, sizeof(float _Complex));
    B = (float _Complex*) calloc(M, sizeof(float _Complex));
    Utmp[0] = sqW * Bfac;
    B[0] = Gfac;

    for(int i = 1; i < M; i++)
    {
        B[i] = Utmp[i-1] * B[i-1];
        Utmp[i] = Utmp[i-1] * W;
    }



    free(Utmp);

    int size = max(N,M)+1;

    Utmp = (float _Complex*) calloc(size, sizeof(float _Complex));
    Vtmp = (float _Complex*) calloc(size, sizeof(float _Complex));
    Utmp[0] = sqW;
    Vtmp[0] = 1.0;
    for(int i = 1; i < size; i++)
    {
        Vtmp[i] = Utmp[i-1]*Vtmp[i-1];
        Utmp[i] = Utmp[i-1]*W;
    }

    D = (float _Complex*) calloc(L, sizeof(float _Complex));

    for(int i = 0; i < M; i++ )
    {
        D[i] = conj(Vtmp[i]);

    }

    for(int i = 0; i < N)
    {
        D[L+1-i] = conj(Vtmp[i+1]);
    }

}

int max(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}




int main()
{
    paramsdata *params = (paramsdata*)malloc(sizeof(paramsdata));
    set_parameters_simulation(params);
    return 0;
}
