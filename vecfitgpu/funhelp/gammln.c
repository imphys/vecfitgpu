#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

double gammln( double xx)
{
    double cofs[] = {76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.12086509738666179e-2, -0.5395239384953e-5};
    double stp = 2.5066282746310005;
    double x = 0;
    double y = 0;
    x = xx;
    y =x;
    double temp = 0;
    tmp = x+5.5;
    tmp = (x+0.5)*log(tmp)-tmp;
    double ser = 1.000000000190015;
    for (int i =0; i<6;i++){
        y = y+1.0;
        ser = ser+cofs[j]/y;
    }
    double out =0;
    out = tmp+log(stp*ser/x);
    return out;
}
