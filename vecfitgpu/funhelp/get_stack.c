#include <stdlib.h>
#include "../paramsdata.h"
#include "get_stack.h"

int **** get_stack(paramsdata *params, int ***allspots)
{

  int K = params->K;
  int numparams = params->numparams;

  int **** spots;

  int Mx = params->Mx, My = params->My, Ncfg = params->Ncfg;

  spots = (int ****)malloc(Mx * sizeof(int ***));
  for(int i = 0; i < Mx; i++)
  {
    spots[i] = (int ***)malloc(My * sizeof(int **));
    for(int j = 0; j < My; j++)
    {
      spots[i][j] = (int **)malloc(K * sizeof(int *));
      for(int k = 0; k < K; k++)
      {
        spots[i][j][k] = (int *)malloc((Ncfg / K) * sizeof(int));
      }
    }
  }

  for(int i = 0; i < Ncfg; i++)
  {
    for(int j = 0; j < Mx; j++)
    {
      for(int k = 0; k < My; k++)
      {

          spots[j][k][(i%K)][i] = allspots[j][k][i];
      }
    }
  }



  //allspots = allspots_temp;
  //% if numparams == 4
  //%     allspots = sum(allspots,4);
  //% end



return spots;


}
