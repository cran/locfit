#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "local.h"

/* crude log-gamma function. If there isn't one in the math library, define
   LGAMMA(x) = lfgamma(x) in local.h */
#define HL2PI 0.91893853320467267
double lfgamma(x)
double x;
{ double x1;
  static double ilg[] = { 0.0, 0.0, 0.69314718055994529,
    1.791759469228055, 3.1780538303479458, 4.7874917427820458, 6.5792512120101012,
    8.5251613610654147, 10.604602902745251, 12.801827480081469 };
  static double hlg[] = { 0.57236494292470008, -0.12078223763524520,
    0.28468287047291918, 1.20097360234707430,  2.45373657084244230,
    3.95781396761871650, 5.66256205985714270,  7.53436423675873360,
    9.54926725730099870, 11.68933342079726900 };

  if (x<=0.0) return(0.0);
  if (x<10)
  { if (x==(int)x) return(ilg[(int)x-1]);
    if ((x-0.5)==(int)(x-0.5)) return(hlg[(int)(x-0.5)]);
  }
  if (x<3) return(lfgamma(x+1)-log(x));

  x1 = x-1;
  return(HL2PI+(x1+0.5)*log(x1)-x1+1/(12*x1));
}

double daws(x)
double x;
{ static double val[] = {
  0,  0.24485619356002, 0.46034428261948, 0.62399959848185, 0.72477845900708,
      0.76388186132749, 0.75213621001998, 0.70541701910853, 0.63998807456541,
      0.56917098836654, 0.50187821196415, 0.44274283060424, 0.39316687916687,
      0.35260646480842, 0.31964847250685, 0.29271122077502, 0.27039629581340,
      0.25160207761769, 0.23551176224443, 0.22153505358518, 0.20924575719548,
      0.19833146819662, 0.18855782729305, 0.17974461154688, 0.17175005072385 };
  double h, f0, f1, f2, y, z, xx;
  int j, m;
  if (x<0) return(-daws(-x));
  if (x>6)
  { /* Tail series: 1/x + 1/x^3 + 1.3/x^5 + 1.3.5/x^7 + ...  */
    y = z = 1/x;
    j = 0;
    while (((f0=(2*j+1)/(x*x))<1) && (y>1.0e-10*z))
    { y *= f0;
      z += y;
      j++;
    }
    return(z);
  }
  m = (int) (4*x);
  h = x-0.25*m;
  if (h>0.125)
  { m++;
    h = h-0.25;
  }
  xx = 0.25*m;
  f0 = val[m];
  f1 = 1-xx*f0;
  z = f0+h*f1;
  y = h;
  j = 2;
  while (fabs(y)>z*1.0e-10)
  { f2 = -(j-1)*f0-xx*f1;
    y *= h/j;
    z += y*f2;
    f0 = f1; f1 = f2;
    j++;
  }
  return(z);
}

double ptail(x) /* exp(x*x/2)*int_{-\infty}^x exp(-u^2/2)du for x < -6 */
double x;
{ double y, z, f0;
  int j;
  y = z = -1.0/x;
  j = 0;
  while ((fabs(f0= -(2*j+1)/(x*x))<1) && (fabs(y)>1.0e-10*z))
  { y *= f0;
    z += y;
    j++;
  }
  return(z);
}

double logit(x)
double x;
{ return(log(x/(1-x)));
}

double expit(x)
double x;
{ if (x<0) return(1-1/(1+exp(x)));
  return(1/(1+exp(-x)));
}
