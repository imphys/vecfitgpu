#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// This function computes the Zernike basis aberration functions for the
// Zernike orders given in the input argument orders.
//
double *** get_zernikefunctions(double **orders, int size_orders, double **x, double **y, int XYpupil_size)
{
  int zersize[2] = {size_orders, 2};
  int Nzer = zersize[0], radormax, azormax;
  int rmax = -1, amax = -1;
  int Nx = XYpupil_size, Ny = XYpupil_size;

  for(int i=0; i<size_orders; i++)
  {
    if(orders[i][0] > rmax)
    {
      rmax = orders[i][0];
    }
    if(abs(orders[i][1]) > amax)
    {
      amax = abs(orders[i][1]);
    }
  }

  radormax = rmax;
  azormax = amax;
  //Evaluation of the radial Zernike polynomials using the recursion relation for
  //the Jacobi polynomials.
  int sizez3 = (2*(radormax + 2))+1;
  int sizez4 = (azormax + 2);
  double zerpol[Nx][Ny][sizez3][sizez4];
  double rhosq[XYpupil_size][XYpupil_size], rho[XYpupil_size][XYpupil_size], phi[XYpupil_size][XYpupil_size];

  double *** allzernikes;

  allzernikes = (double ***)malloc(Nx * sizeof(double **));
  for(int i = 0; i< Nx; i++)
  {
    allzernikes[i] = (double **)malloc(Ny * sizeof(double *));
    for(int j = 0; j < Ny; j++)
    {
        allzernikes[i][j] = (double *)malloc(Nzer * sizeof(double));

    }
  }



  for(int i = 0; i< Nx; i++)
  {
    for(int j = 0; j < Ny; j++)
    {
      rhosq[i][j] = x[i][j] * x[i][ j]+ y[i][j] * y[i][j];
      rho[i][j] = sqrt(rhosq[i][j]);

      for(int k = 0; k < sizez3; k++)
      {
        for(int l = 0; l < sizez4; l++)
        {
          if(k == 0 && l == 0)
          {
            zerpol[i][j][k][l] = 1;
          }
          else
          {
            zerpol[i][j][k][l] = 0;
          }

        }
      }
    }

  }



  int m, mm;

  for(int i = 0; i< Nx; i++)
  {
    for(int j = 0; j < Ny; j++)
    {
      phi[i][j] = atan2(y[i][j], x[i][j]);
      for(int jm = 0; jm < (azormax + 2); jm++)
      {
        m = jm - 1;
        //mm used because m is for the index in matrix and mm is the same as m in matlab
        mm = m+1;
        if(m > -1)
        {
          zerpol[i][j][jm][jm] = rho[i][j] * zerpol[i][j][jm-1][jm - 1];
        }
        ////
        ////NO SQUEEZE USED
        zerpol[i][j][jm+2][jm] = ((mm + 2) * rhosq[i][j] - mm - 1) * zerpol[i][j][jm][jm];
        for(int p = 1; p < radormax - mm + 2; p++)
        {
          int n = mm+2*(p+1);
          int jn = n;
          zerpol[i][j][jn][jm] = (2*(n-1)*(n*(n-2)*(2*rhosq[i][j]-1)-(mm * mm))*zerpol[i][j][jn-2][jm] - n*(n+mm-2)*(n-mm-2)*zerpol[i][j][jn-4][jm])/((n-2)*(n+mm)*(n-mm));
        }

      }
    }

  }

  for(int i = 0; i< Nx; i++)
  {
    for(int j = 0; j < Ny; j++)
    {
      for(int k = 0; k < Nzer; k++)
      {
        int n = orders[k][0];
        int m = orders[k][1];
        if(m >= 0)
        {
          allzernikes[i][j][k] = zerpol[i][j][n][m] * cos(m*phi[i][j]);
        }
        else
        {
          allzernikes[i][j][k] = zerpol[i][j][n][-m] * sin(-m*phi[i][j]);
        }
      }
    }

  }

  return allzernikes;

}
