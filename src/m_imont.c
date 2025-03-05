/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *  Multivariate integration of a vector-valued function
 *  using Monte-Carlo method.
 *
 *  uses drand48() random number generator. Does not seed.
 */

#include <stdlib.h>
#include "R.h"
#include "mutil.h"
extern void setzero();

void monte(int (*f)(), double *ll, double *ur, int d, double *res, int n);

/* static int lfindex[MXIDIM];
   static double M[(1+MXIDIM)*MXIDIM*MXIDIM]; */

void monte(int (*f)(), double *ll, double *ur, int d, double *res, int n)
/*int (*f)(), d, n;
double *ll, *ur, *res;*/
{ int i, j, nr=0;
  double z, x[MXIDIM], tres[MXRESULT];

/* srand48(234L); */
  GetRNGstate(); /* Use R's RNG */

  for (i=0; i<n; i++)
  { for (j=0; j<d; j++) x[j] = ll[j] + (ur[j]-ll[j])*unif_rand(); /* drand48();*/
    nr = f(x,d,tres,NULL);
    if (i==0) setzero(res,nr);
    for (j=0; j<nr; j++) res[j] += tres[j];
  }

  z = 1;
  for (i=0; i<d; i++) z *= (ur[i]-ll[i]);
  for (i=0; i<nr; i++) res[i] *= z/n;
 
  PutRNGstate();
}
