#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <time.h>
#include "tiffio.h"

#include "set_parameters_simulation.h"






int main()
{
    paramsdata *params = (paramsdata*)malloc(sizeof(paramsdata));
    set_parameters_simulation(params);
    return 0;
}
