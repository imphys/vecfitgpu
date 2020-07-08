#include "paramsdata.h"
#include "FFT/fft.h"
#include "get_allPSFs.h"

#include <string.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

double _Complex ** repmat(double _Complex *src, int size_src, int N, int r)
{

  double _Complex **dst;
  if(r==1)
  {
    dst =(double _Complex **) malloc(N * sizeof(double _Complex*));

    for(int i = 0; i < N;i++ )
    {
      dst[i] = (double _Complex *)malloc(size_src * sizeof(double _Complex));
    }

    for(int i = 0; i < N;i++ )
    {
      for(int j = 0; j < size_src; j++)
      {
        dst[i][j] = src[j];
      }

    }
  }
  else
  {
      printf("r is Larger than 1 in repmat\n");
  }
  return dst;
}

void prechirpz(int xsize, int qsize, int N, int M, double _Complex *A, double _Complex *B, double _Complex *D);

int max(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

void set_parameters_simulation(paramsdata *params){
    params->K = 1;
    params->m = 0.95;
    params->alpha = 3.141592653589793;


    get_PSFs(&params->allPSFs, params);
    params->pixelsize = 80;
    params->xrange = params->pixelsize * params->Mx/2;
    params->yrange = params->pixelsize * params->My/2;
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


    params->Mz = 1;





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
    if(!strcmp(params->ztype,"medium"))
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
        printf("Warning! Refractive index immersion medium cannot be smaller than NA.\n" );
    }
    if(params->NA > params->refcov)
    {
        printf("Warning! Refractive index cover slip cannot be smaller than NA.\n");
    }

    //parameters needed for fixed dipole PSF only: emitter/absorber dipole
    //orientation (characterized by angles pola and azim)
    strcpy(params->dipoletype, "diffusion");
    //printf("%s", params->dipoletype);
    //strcpy(params->dipoletype, 'free');
    //strcpy(params->dipoletype, 'fixed');
    params->pola = 90.0*M_PI/180;
    params->azim = 0.0*M_PI/180;


    //diffusion coefficient

    double g2;
    double eps = 2.2204e-16;

    double welldepth = 1/eps;
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
            params->zonefunction[i][j] = 0;
        }
    }

    //params->aberrationsoffset = [];
    params->aberrations[3][1] = -1;
    params->zonefunction[3][1] = -1;

    params->aberrations[4][1] = 1;
    params->zonefunction[4][1] = 1;

    params->aberrations[0][0] = 2;
    params->zonefunction[0][0] = 2;

    params->aberrations[1][0] = 2;
    params->zonefunction[1][0] = 2;

    params->aberrations[2][0] = 2;
    params->zonefunction[2][0] = 2;

    params->aberrations[1][1] = -2;
    params->zonefunction[1][1] = -2;

    params->aberrations[2][1] = 2;
    params->zonefunction[2][1] = 2;

    params->aberrations[8][1] = -2;
    params->zonefunction[8][1] = -2;

    params->aberrations[3][0] = 3;
    params->zonefunction[3][0] = 3;

    params->aberrations[4][0] = 3;
    params->zonefunction[4][0] = 3;

    params->aberrations[6][0] = 3;
    params->zonefunction[6][0] = 3;

    params->aberrations[7][0] = 3;
    params->zonefunction[7][0] = 3;

    params->aberrations[6][1] = -3;
    params->zonefunction[6][1] = -3;

    params->aberrations[7][1] = 3;
    params->zonefunction[7][1] = 3;

    params->aberrations[5][0] = 4;
    params->zonefunction[5][0] = 4;

    params->aberrations[8][0] = 4;
    params->zonefunction[8][0] = 4;

    for(int i = 0; i<9; i++)
    {
        params->aberrations[i][2] = params->aberrations[i][2]*params->lambdacentral;
        params->zonefunction[i][2] = params->aberrations[i][2]*params->lambdacentral;
    }

    // for(int i = 0; i<9;i++)
    // {
    //     for(int j = 0; j<3; j++)
    //     {
    //         printf("%d\t", params->aberrations[i][j]);
    //     }
    //     printf("\n");
    // }

    //DOE/SLM
    //strcpy(params->doetype, 'none');
    strcpy(params->doetype, "vortex");
    params->ringradius = 1;
    params->doelevels = 32;
    //params->zonefunction = params->aberrations;
    params->doephasedepth = 1*params->lambdacentral;

    //Fit model parameters: signal photon count, background photons/pixel, read
    //noise variance for sCMOS camera's, fit model, output labels depending on
    //fit model
    //params.signalphotoncount = 1000;
    //params.backgroundphotoncount = 10;
    params->readnoisestd = 0;

    //calculate auxiliary vectors for chirpz
    int Kx, Ky;
    double PupilSize, ImageSizex, ImageSizey;
    //double Ax, Ay, Bx, By, Dx, Dy;


    Kx = params->Mx;
    Ky = params->My;
    PupilSize = 1.0;
    ImageSizex = params->xrange*params->NA/params->lambda;
    ImageSizey = params->yrange*params->NA/params->lambda;


    double _Complex *Ax, *Bx, *Dx;
    double _Complex *Ay, *By, *Dy;
    int Lx, Ly;

    Lx = params->Npupil + Kx - 1;
    Ly = params->Npupil + Ky - 1;

    Ax = (double _Complex*) calloc(params->Npupil, sizeof(double _Complex));
    Bx = (double _Complex*) calloc(Kx, sizeof(double _Complex));
    Dx = (double _Complex*) calloc(Lx, sizeof(double _Complex));


    prechirpz(PupilSize, ImageSizex, params->Npupil, Kx, Ax, Bx, Dx);

    Ay = (double _Complex*) calloc(params->Npupil, sizeof(double _Complex));
    By = (double _Complex*) calloc(Ky, sizeof(double _Complex));
    Dy = (double _Complex*) calloc(Ly, sizeof(double _Complex));


    prechirpz(PupilSize, ImageSizey, params->Npupil, Ky, Ay, By, Dy);

    double _Complex **Axmt = repmat(Ax, params->Npupil, params->Mx, 1);
    double _Complex **Bxmt = repmat(Bx, Kx, params->Mx, 1);
    double _Complex **Dxmt = repmat(Dx, Lx, params->Mx, 1);

    double _Complex **Aymt = repmat(Ay, params->Npupil, params->Npupil, 1);
    double _Complex **Bymt = repmat(By, Ky, params->Npupil, 1);
    double _Complex **Dymt = repmat(Dy, Ly, params->Npupil, 1);

    params->cztN = params -> Npupil;
    params->cztM = params -> Mx;
    params->cztL = params -> Npupil + params -> Mx - 1;


    params->debugmode = 0;

}

void prechirpz(int xsize, int qsize, int N, int M, double _Complex *A, double _Complex *B, double _Complex *D)
{
    int L;
    double sigma;
    double _Complex Afac, Bfac, Gfac, sqW, W;
    L = N + M - 1;
    sigma = 2 * M_PI * xsize * qsize / N / M;
    Afac = cexpf(2 * _Complex_I * sigma * (1 - M));
    Bfac = cexpf(2 * _Complex_I * sigma * (1 - N));
    sqW  = cexpf(2* _Complex_I * sigma);
    W = sqW * sqW;
    Gfac = (2*xsize/N)*cexpf(_Complex_I * sigma * (1-N) * (1-M));

    double _Complex *Utmp, *Vtmp;

    Utmp = (double _Complex*) calloc(N, sizeof(double _Complex));
    //A = (double _Complex*) calloc(N, sizeof(double _Complex));

    A[0] = 1.0;
    Utmp[0] = sqW * Afac;
    for(int i = 1; i < N; i++)
    {https://www.google.com/search?client=ubuntu&channel=fs&q=what+colour+is+infrared+light&ie=utf-8&oe=utf-8
        A[i] = Utmp[i-1]*A[i-1];
        Utmp[i] = Utmp[i-1]*W;
    }


    free(Utmp);

    Utmp = (double _Complex*) calloc(M, sizeof(double _Complex));
    //B = (double _Complex*) calloc(M, sizeof(double _Complex));
    Utmp[0] = sqW * Bfac;
    B[0] = Gfac;

    for(int i = 1; i < M; i++)
    {
        B[i] = Utmp[i-1] * B[i-1];
        Utmp[i] = Utmp[i-1] * W;
    }



    free(Utmp);

    int size = max(N,M)+1;

    Utmp = (double _Complex*) calloc(size, sizeof(double _Complex));
    Vtmp = (double _Complex*) calloc(size, sizeof(double _Complex));
    Utmp[0] = sqW;
    Vtmp[0] = 1.0;
    for(int i = 1; i < size; i++)
    {
        Vtmp[i] = Utmp[i-1]*Vtmp[i-1];
        Utmp[i] = Utmp[i-1]*W;
    }

    //D = (double _Complex*) calloc(L, sizeof(double _Complex));

    for(int i = 0; i < M; i++ )
    {
        D[i] = conj(Vtmp[i]);
    }

    for(int i = 0; i < N; i++)
    {
        D[L-1-i] = conj(Vtmp[i+1]);
    }


    double *Dreal, *Dimag;

    Dreal = (double *)malloc(L * sizeof(double));
    Dimag = (double *)malloc(L * sizeof(double));

    for(int i = 0; i < L; i++)
    {
      Dreal[i] = creal(D[i]);
      Dimag[i] = cimag(D[i]);
    }


    Fft_transform(Dreal, Dimag, L);

    for(int i = 0; i < L; i++)
    {
      D[i] = Dreal[i] + (Dimag[i] * _Complex_I);
    }



}
