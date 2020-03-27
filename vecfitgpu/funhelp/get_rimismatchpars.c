#include "../paramsdata.h"
#include <stdio.h>
#include <math.h>

void get_rimismatchpars(paramsdata *parameters)
{
  double refins[3] = {parameters->refimm, parameters->refimmnom, parameters->refmed};
  double zvals[3] = {0, parameters->fwd, -parameters->depth};
  double NA = parameters->NA;
  //reduce NA in case of TIRF conditions
  if (NA>parameters->refmed)
  {
    NA = parameters->refmed;
  }

  double paraxiallimit = 0.2;
  int K = sizeof(refins)/sizeof(refins[0]);

  long double zvalsratio[K];
  long double Wrmsratio[K][K];

  for(int i = 0; i < K; i++)
  {
    zvalsratio[i] = 0;
    for(int j = 0; j < K; j++)
    {
      Wrmsratio[i][j] = 0;
    }
  }


  if(NA > paraxiallimit)
  {
    long double fsqav[K], fav[K], Amat[K][K];

    for(int i = 0; i < K; i++)
    {
      fsqav[i] = (refins[i]*refins[i]) - (0.5) * (NA*NA);
      fav[i] = (2.0/3/(NA*NA)) * ((refins[i]*refins[i]*refins[i]) - pow(((refins[i]*refins[i])-(NA*NA)),(3.0/2)));
      // fav[i] = (2/3)*(3*refins[i]^4 - 3 * (refins[i]*refins[i]) * (NA * NA) + NA^4 ) / (refins[i]^3 + ((refins[i] * refins[i]) - NA*NA)^(3/2);

      Amat[i][i] = fsqav[i] - (fav[i] * fav[i]);



      for(int kk = 0; kk < i; kk++)
      {
        Amat[i][kk] = (1.0/4/(NA*NA))*(refins[i]*refins[kk]*((refins[i]*refins[i])+(refins[kk]*refins[kk]))
                -((refins[i]*refins[i])+(refins[kk]*refins[kk])-2*(NA*NA))*sqrt((refins[i]*refins[i])-(NA*NA))*sqrt((refins[kk]*refins[kk])-(NA*NA))
                +pow(((refins[i]*refins[i])-(refins[kk]*refins[kk])),2)*log((sqrt((refins[i]*refins[i])-(NA*NA))+sqrt((refins[kk]*refins[kk])-(NA*NA)))/(refins[i]+refins[kk])));
        Amat[i][kk] = Amat[i][kk] - fav[i]*fav[kk];

        Amat[kk][i] = Amat[i][kk];


      }


    }
    for(int jv = 1; jv < K; jv++)
    {
      zvalsratio[jv] = Amat[0][jv]/Amat[1][1];
      for(int kv = 1; kv < K; kv++)
      {
        Wrmsratio[jv][kv] = Amat[jv][kv] -  Amat[0][jv] * Amat[0][kv]/Amat[0][0];
        printf("%.20Lf\n",  1000000 * Amat[jv][kv] -  1000000 *Amat[0][jv] * Amat[0][kv]/Amat[0][0]);
      }
    }




  }
  else
  {
    //paraxial limit, Taylor-series in NA
    for(int jv = 1; jv < K; jv++)
    {
      zvalsratio[jv] = refins[0]/refins[jv]+(NA * NA)*((refins[0] * refins[0])-(refins[jv]*refins[jv]))/(4.0*refins[0]*(refins[jv] * refins[jv] * refins[jv]));
      for(int kv = 1; kv < K; kv++)
      {
        Wrmsratio[jv][kv] = pow(NA,8) * ((refins[0]*refins[0]) - (refins[jv] * refins[jv])) * ((refins[0]*refins[0]) - (refins[kv]*refins[kv]))/(11520.0*pow(refins[0],4)*pow(refins[jv],3)*pow(refins[kv],3));
      }
    }


  }


  for(int i = 0; i < K; i++)
  {

    for(int j = 0; j < K; j++)
    {
        printf("%.20Lf\t", Wrmsratio[i][j]);
    }
    printf("\n");
  }

  float Wrms;
  zvals[0] = zvalsratio[1]*zvals[1]+zvalsratio[2]*zvals[2];

  Wrms = Wrmsratio[1][1]*(zvals[1]*zvals[1])+Wrmsratio[2][2]*(zvals[2]*zvals[2])+2*Wrmsratio[1][2]*zvals[1]*zvals[2];
  Wrms = sqrt(Wrms);


  // if(1)
  // {
  //   printf("image plane depth from cover slip = %4.0f nm\n",-zvals[2]);
  //   printf("free working distance = %6.2f mu\n",1e-3*zvals[1]);
  //   printf("nominal z-stage position = %6.2f mu\n",1e-3*zvals[0]);
  //   printf("rms aberration due to RI mismatch = %4.1f mlambda\n",1e3*Wrms/parameters->lambda);
  // }
}
