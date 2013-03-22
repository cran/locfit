/*
 *   Copyright (c) 1998 Lucent Technologies.
 *   See README file for details.
 */

/*
  Functions implementing the adaptive bandwidth selection.
  Usual entry point is afit(des,lf,x).
  Will make the final call to nbhd() to set smoothing weights
  for selected bandwidth, But will **not** make the
  final call to locfit().
*/

#include <math.h>
#include <stdio.h>
#include "local.h"

double acri(lk,t0,t2,pen)
double lk, t0, t2, pen;
{ /* return(-2*lk/(t0*exp(pen*log(1-t2/t0)))); */
  /* return((-2*lk+pen*t2)/t0); */
  return((MAX(-2*lk,t0-t2)+pen*t2)/t0);
}

double mmse(lf,des,x)
struct tree *lf;
struct design *des;
double *x;
{ int i, ii, j, p, p1;
  double sv, sb, *l, dp;
  l = des->wd;
  makelxd(lf,des,x,l,0,NULL,0,1);
  sv = sb = 0;
  p = lf->mi[MP];
  for (i=0; i<des->n; i++)
  { sv += l[i]*l[i];
    ii = des->ind[i];
    dp = des->di[ii];
    for (j=0; j<p; j++) dp *= des->di[ii];
    sb += fabs(l[i])*dp;
  }
  p1 = factorial(lf->mi[MDEG]+1);
  return(sv+sb*sb*lf->dp[DADP]*lf->dp[DADP]/(p1*p1));
}

static double mcp, clo, cup;

/*
  Initial bandwidth will be (by default)
  k-nearest neighbors for k small, just lage enough to
  get defined estimate (unless user provided nonzero DALP
  or DFXH components)
*/

double ainitband(des,lf,x)
struct design *des;
struct tree   *lf;
double *x;
{ INT m, p, z;
  double h, t[6];
  p = des->p;
  z = (INT)(lf->mi[MN]*lf->dp[DALP]); if (z<p+3) z = p+3;
  do
  { h = nbhd(x,lf,des,1.0*z/lf->mi[MN],lf->dp[DFXH],0);
    if (h>0) m = locfit(lf,des,x,h,1);
    z++;
  } while ((h==0)||(m!=0));

  switch(lf->mi[MACRI])
  { case ACP:
      ldf(lf,des,t,1,lf->mi,NULL);
      mcp = acri(des->llk,t[0],t[2],lf->dp[DADP]);
      return(h);
    case AKAT:
      ldf(lf,des,t,1,lf->mi,NULL);
      clo = des->cf[0]-lf->dp[DADP]*t[5];
      cup = des->cf[0]+lf->dp[DADP]*t[5];
      return(h);
    case AMDI:
      mcp = mmse(lf,des,x);
      return(h);
  }
  ERROR(("aband1: unknown criterion"));
}

/*
  aband2 increases the initial bandwidth until lack of fit results,
  or the fit is close to a global fit. Increase h by 1+0.3/d at
  each iteration.
*/

double aband2(des,lf,x,h0)
struct design *des;
struct tree   *lf;
double *x, h0;
{ double t[6], h, h1, nu1, cp, ncp, tlo, tup;
  INT d, inc, n, m, p, done;
  d = lf->mi[MDIM]; n = lf->mi[MN]; p = lf->mi[MP];
  h1 = h = h0;
  done = 0; nu1 = 0.0;
  inc = 0; ncp = 0.0;
  while ((!done) & (nu1<(n-p)*0.95))
  { h = nbhd(x,lf,des,0.0,(1+0.3/d)*h,1);
    if (locfit(lf,des,x,h,1)>0) if (m>0) WARN(("aband2: failed fit"));
    ldf(lf,des,t,1,lf->mi,NULL);
    nu1 = t[0]-t[2]; /* tr(A) */
    switch(lf->mi[MACRI])
    { case AKAT:
        tlo = des->cf[0]-lf->dp[DADP]*t[5];
        tup = des->cf[0]+lf->dp[DADP]*t[5];
printf("h %8.5f  tlo %8.5f  tup %8.5f\n",h,tlo,tup);
        done = ((tlo>cup) | (tup<clo));
        if (!done)
        { clo = MAX(clo,tlo);
          cup = MIN(cup,tup);
          h1 = h;
        }
        break;
      case ACP:
        cp = acri(des->llk,t[0],t[2],lf->dp[DADP]);
printf("h %8.5f  lk %8.5f  t0 %8.5f  t2 %8.5f  cp %8.5f\n",h,des->llk,t[0],t[2],cp);
        if (cp<mcp) { mcp = cp; h1 = h; }
        if (cp>ncp) inc++; else inc = 0;
        ncp = cp;
        done = ((inc>=3) & ((t[0]-t[2])>=10) & (cp>1.5*mcp));
        break;
      case AMDI:
        cp = mmse(lf,des,x);
        if (cp<mcp) { mcp = cp; h1 = h; }
        if (cp>ncp) inc++; else inc = 0;
        ncp = cp;
        done = (inc>=3);
        break;
    }
  }
  return(h1);
}

/*
  aband3 does a finer search around best h so far. Try
  h*(1-0.2/d), h/(1-0.1/d), h*(1+0.1/d), h*(1+0.2/d)
*/
double aband3(des,lf,x,h0)
struct design *des;
struct tree   *lf;
double *x, h0;
{ double t[6], nu1, h, h1, cp, tlo, tup;
  INT i, i0, d, n;
  d = lf->mi[MDIM]; n = lf->mi[MN];

  h1 = h0;
  i0 = (lf->mi[MACRI]==AKAT) ? 1 : -2;
  for (i=i0; i<=2; i++)
  { if (i==0) i++;
    h = h0*(1+0.1*i/lf->mi[MDIM]);
    h = nbhd(x,lf,des,0.0,h,1);
    if (locfit(lf,des,x,h,1)>0) WARN(("aband3: failed fit"));
    ldf(lf,des,t,1,lf->mi,NULL);
    switch (lf->mi[MACRI])
    { case AKAT:
        tlo = des->cf[0]-lf->dp[DADP]*t[5];
        tup = des->cf[0]+lf->dp[DADP]*t[5];
        if ((tlo>cup) | (tup<clo)) /* done */
          i = 2;
        else
        { h1 = h;
          clo = MAX(clo,tlo);
          cup = MIN(cup,tup);
        }
        break;
      case ACP:
        nu1 = t[0]-t[2]; /* tr(A) */
        cp = acri(des->llk,t[0],t[2],lf->dp[DADP]);
printf("h %8.5f  lk %8.5f  t0 %8.5f  t2 %8.5f  cp %8.5f\n",h,des->llk,t[0],t[2],cp);
        if (cp<mcp) { mcp = cp; h1 = h; }
        else
        { if (i>0) i = 2; }
      case AMDI:
        cp = mmse(lf,des,x);
        if (cp<mcp) { mcp = cp; h1 = h; }
        else
        { if (i>0) i = 2; }
    }
  }
  return(h1);
}

double afit(des,lf,x)
struct design *des;
struct tree   *lf;
{ double h;
  h = ainitband(des,lf,x);
  h = aband2(des,lf,x,h);
  h = aband3(des,lf,x,h);
  return(nbhd(x,lf,des,0.0,h,1));
}
