/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

/*
  Compute minimax weights for local regression.
*/

static double alp[10];
static int debug;
#define SINGTOL 1.0e-3
#define CONVTOL 1.0e-8

double ipower(x,n) /* use for n not too large!! */
double x;
INT n;
{ if (n==0) return(1.0);
  if (n<0) return(1/ipower(x,-n));
  return(x*ipower(x,n-1));
}

double setmmwt(des,lf,a,gam)
design *des;
lfit *lf;
double *a, gam;
{ double ip, w0, w1, sw, wt;
  INT i, p;
  sw = 0.0;
  p = lf->mi[MP];
  for (i=0; i<lf->mi[MN]; i++)
  { ip = innerprod(a,&des->X[i*p],p);
    wt = prwt(lf,i);
    w0 = wt*(ip - gam*des->wd[i]);
    w1 = wt*(ip + gam*des->wd[i]);
    des->w[i] = 0.0;
    if (w0>0) { des->w[i] = w0; sw += w0*w0; }
    if (w1<0) { des->w[i] = w1; sw += w1*w1; }
  }
if (debug) printf("try: %8.5f %8.5f %8.5f\n",a[0],a[1],sw/2-a[0]);
  return(sw/2-a[0]);
}

/* compute sum_{w!=0} AA^T; e1-sum wA  */
INT mmsums(des,lf,z,y) /* y=1: don't decompose jacobian */
double *z;
design *des;
lfit *lf;
INT y;
{ INT i, j, p;
  double *A;
  A = des->xtwx.Z;
  p = lf->mi[MP];
  for (i=0; i<p*p; i++) A[i] = 0.0;
  z[0] = 1.0;
  for (i=1; i<p; i++) z[i] = 0.0;

  for (i=0; i<lf->mi[MN]; i++)
    if (des->w[i]!=0.0)
    { addouter(A,&des->X[i*p],&des->X[i*p],p,1.0);
      for (j=0; j<p; j++) z[j] -= des->w[i]*des->X[i*p+j];
    }
  if (y) return(0);

  /* now, decompose. */
  for (i=0; i<p; i++)
  { des->xtwx.dg[i] = A[i*p+i];
    if (des->xtwx.dg[i]>0) des->xtwx.dg[i] = 1/sqrt(A[i*p+i]);
  }
  for (i=0; i<p; i++)
    for (j=0; j<p; j++)
      A[i*p+j] *= des->xtwx.dg[i]*des->xtwx.dg[j];
  eigen(A,des->xtwx.Q,p,lf->mi[MMXIT]);
  des->xtwx.sm = 1;

  /* is it nearly singular? */
  for (i=0; i<p; i++)
    if (A[i*p+i]<SINGTOL) return(1);
  return(0);
}

double updatenr(des,lf,z,p,a,a0,sw0,gam)
design *des;
lfit *lf;
INT p;
double *z, *a, *a0, sw0, gam;
{ double f, sw;
  INT i, done;
  vxtwx(&des->xtwx,z,p);
  if (lf_error) return(0.0);
  done = 0;
  f = 1.0;
  while (!done)
  { for (i=0; i<p; i++) a[i] = alp[i]+f*z[i];
    sw = setmmwt(des,lf,a,gam);
    if ((sw==0.0) | (sw>sw0+CONVTOL)) f /= 2.0;
    else done = 1;
    if (f<0.000000001)
    { printf("updatenr failure\n");
/* printf("z: %8.5f %8.5f\n",z[0],z[1]);
printf("a0 %8.5f %8.5f  gam %8.5f\n",a0[0],a0[1],gam);
for (i=0; i<lf->mi[MN]; i++)
  printf("%2d %8.5f %8.5f\n",i,des->wd[i],des->w[i]); */
      return(0.0);
    }
  }
  return(sw);
}

double updatesa(des,lf,z,p,a,a0,sw0,gam)
design *des;
lfit *lf;
int p;
double *z, *a, *a0, sw0, gam;
{ double f, sw, c0, c1, lb;
  INT i, j, done;
  c0 = c1 = 0.0;
  for (i=0; i<p; i++) c0 += z[i]*z[i];

  mmsums(des,lf,z,1);
if (debug) {
printf("A: %8.5f %8.5f  %8.5f\n",des->xtwx.Z[0],des->xtwx.Z[1],z[0]);
printf("A: %8.5f %8.5f  %8.5f\n",des->xtwx.Z[2],des->xtwx.Z[3],z[1]);
}
  for (i=0; i<p; i++)
    for (j=0; j<p; j++)
      c1 += z[i]*des->xtwx.Z[i*p+j]*z[j];
  lb = c0/c1;
if (debug) printf("sa: c0 %8.5f c1 %8.5f lb %8.5f\n",c0,c1,lb);

  done = 0;
  f = 1.0;
  while (!done)
  { for (i=0; i<p; i++) a[i] = alp[i]+f*lb*z[i];
    sw = setmmwt(des,lf,a,gam);
if (debug) printf("sa: f %8.5f  sw %8.5f\n",f,sw);
    if ((sw==0.0) | (sw>sw0+CONVTOL)) f /= 2.0;
    else done = 1;
    if (f<0.001)
    { printe("updatesa failure\n");
      return(0.0);
    }
  }
  return(sw);
}

double updatesd(des,lf,z,p,a,a0,sw0,gam)
design *des;
lfit *lf;
int p;
double *z, *a, *a0, sw0, gam;
{ double f, sw, c0, c1, tmp[10];
  INT i, j, sd;

if (debug) printf("updatesd\n");
  for (i=0; i<p; i++) if (des->xtwx.Z[i*p+i]<SINGTOL) sd = i;
  if (des->xtwx.dg[sd]>0)
    for (i=0; i<p; i++) tmp[i] = des->xtwx.Q[p*i+sd]*des->xtwx.dg[i];
  else
  { for (i=0; i<p; i++) tmp[i] = 0.0;
    tmp[sd] = 1.0;
  }

  mmsums(des,lf,z,1);
  c0 = c1 = 0.0;
  for (i=0; i<p; i++)
  { c0 += tmp[i]*z[i];
    for (j=0; j<p; j++)
      c1 += tmp[i]*des->xtwx.Z[i*p+j]*tmp[j];
  }
if (debug) printf("sdir: c0 %8.5f  c1 %8.5f  z %8.5f %8.5f  tmp %8.5f %8.5f\n",c0,c1,z[0],z[1],tmp[0],tmp[1]);
  if (c0<0) for (i=0; i<p; i++) tmp[i] = -tmp[i];

  f = 1.0;
  for (i=0; i<p; i++) a[i] = a0[i]+tmp[i];
  sw = setmmwt(des,lf,a,gam);
  
  if (sw<sw0) /* double till we drop */
  { while(1)
    { f *= 2;
      sw0 = sw;
      for (i=0; i<p; i++) a[i] = a0[i]+f*tmp[i];
      sw = setmmwt(des,lf,a,gam);
      if (sw>sw0-CONVTOL) /* go back one step */
      { f /= 2;
        for (i=0; i<p; i++) a[i] = a0[i]+f*tmp[i];
        sw0 = setmmwt(des,lf,a,gam);
        return(sw0);
      }
    }
  }

  /* halve till success */
  while (1)
  { f *= 0.5;
    for (i=0; i<p; i++) a[i] = a0[i]+f*tmp[i];
    sw = setmmwt(des,lf,a,gam);
    if (sw<sw0+CONVTOL) return(sw);
  }
}

double updatea0(des,lf,z,p,a,a0,sw0,gam)
design *des;
lfit *lf;
int p;
double *z, *a, *a0, sw0, gam;
{ double alo, ahi, sw, inc, f;
  int i;
  for (i=0; i<p; i++) a[i] = a0[i];
  sw = setmmwt(des,lf,a,gam);
if (debug) printf("1 a %8.5f  sw %8.5f\n",a[0],sw);
  if (sw!=sw0) printf("updatea0: start error\n");
  mmsums(des,lf,z,1);
  
  /* z[0] = 1-sum(wi). If +ve, increase a */
  inc = (z[0]>0) ? 1.0 : -1.0;
  if (z[0]<0) WARN(("updatea0: z0<0 not implemented"));
  f = 1.0; alo = a0[0];
  while(z[0]>0)
  { a[0] = ahi = a0[0]+f*inc;
    sw = setmmwt(des,lf,a,gam);
    mmsums(des,lf,z,1);
if (debug) printf("2 a %8.5f  sw %8.5f  z %8.5f\n",a[0],sw,z[0]);

    if (sw<sw0) sw0 = sw;
    if (fabs(z[0])<1.0e-8) return(sw); /* bingo! */
    if (z[0]>0) { f *= 2.0; alo = ahi; }
  }
if (debug) printf("alo %8.5f ahi %8.5f\n",alo,ahi);
  while (1)
  { a[0] = (ahi+alo)/2;
    sw = setmmwt(des,lf,a,gam);
    mmsums(des,lf,z,1);
if (debug) printf("3 a %8.5f  sw %8.5f  z %8.5f\n",a[0],sw,z[0]);
    if ((sw<sw0+CONVTOL) & (a[0]+sw>0)) return(sw);
    if (z[0]>0) alo = a[0]; else ahi = a[0];
  }
}

double findab(gam,lf,des,a)
lfit *lf;
design *des;
double gam, *a;
{ double *A, *z, sl, sw, sw0;
  INT i, j, n, p, sing, done;
  if (debug) printf("findab: gam = %8.5f\n",gam);
  n = lf->mi[MN]; p = lf->mi[MP];
  A = des->xtwx.Z; z = des->f1;
  
  done = 0;
  for (i=0; i<p; i++) a[i] = alp[i] = 0.0;
  sw0 = sw = updatea0(des,lf,z,p,a,alp,0.0,gam);

  for (j=0; j<lf->mi[MMXIT]; j++)
  { /* store initial coefficients */
    for (i=0; i<p; i++) alp[i] = a[i];

    /* compute z = Newton-Raphson increment */
    sing = mmsums(des,lf,z,0);
    if (debug) printf("sm %2d\n",des->xtwx.sm);
    if ((debug) & (sing)) printf("SINGULAR!!!!\n");
    sw0 = sw;

    sw = (sing) ? updatesd(des,lf,z,p,a,alp,sw0,gam)
                : updatenr(des,lf,z,p,a,alp,sw0,gam);
    if (sw==0.0) return(0.0);

    if (debug)
    { for (i=0; i<p; i++) printf("%8.5f ",a[i]);
      printf(" sw %8.5f\n",sw);
    }
    for (i=0; i<p; i++) alp[i] = a[i];
    done = ((j>0) & (fabs(sw-sw0)<CONVTOL));
    if (done) j = lf->mi[MMXIT];
  }
  if (!done) { WARN(("findab not converged")); return(0.0); }
  sl = 0.0;
  for (i=0; i<n; i++) sl += fabs(des->w[i])*des->wd[i];
  if (debug) printf("sl = %8.5f\n",sl);
  return(sl);
}

double findgam(lf,des)
lfit *lf;
design *des;
{ double g[5], z[5], *a;
  INT i, n;
  n = lf->mi[MN];
  a = des->xtwx.f2;
  g[0] = 0.0;
  a[0] = 1.0/n;
  for (i=1; i<lf->mi[MP]; i++) a[i] = 0.0;
  z[0] = findab(g[0],lf,des,a);
  g[4] = 1.0;
  z[4] = findab(g[4],lf,des,a);
  while (z[4]>g[4])
  { z[0] = z[4]; g[0] = g[4];
    g[4] *= 2.0;
    z[4] = findab(g[4],lf,des,a);
  }
  g[1] = g[0]; z[1] = z[0];
  g[2] = g[4]; z[2] = z[4];
  do
  { g[3] = g[2] + (g[2]-g[1])*(z[2]-g[2])/(z[1]-g[1]-z[2]+g[2]);
    if ((g[3]<=g[0]) | (g[3]>=g[4])) g[3] = (g[0]+g[4])/2;
    z[3] = findab(g[3],lf,des,a);
    if (z[3]>g[3]) { g[0] = g[3]; z[0] = z[3]; }
        else { g[4] = g[3]; z[4] = z[3]; }
    z[1] = z[2]; z[2] = z[3];
    g[1] = g[2]; g[2] = g[3];
  } while ((fabs(g[3]-z[3])>0.0000001) & (g[1] != g[2]));
  return(g[3]);
}

double weightmm(di,ff,mi,gam)
double di, *ff, gam;
INT *mi;
{ double y1, y2, ip;
  ip = innerprod(ff,alp,mi[MP]);
  y1 = ip-gam*di; if (y1>0) return(y1/ip);
  y2 = ip+gam*di; if (y2<0) return(y2/ip);
  return(0.0);
}

double minmax(lf,des)
lfit *lf;
design *des;
{ double h, u[MXDIM], gam;
  INT i, j, m, p1;

  p1 = factorial(lf->mi[MDEG]+1);
  for (i=0; i<lf->mi[MN]; i++)
  { for (j=0; j<lf->mi[MDIM]; j++) u[j] = datum(lf,j,i);
    des->wd[i] = lf->dp[DALP]/p1*ipower(des->di[i],1+lf->mi[MDEG]);
    fitfun(lf,u,des->xev,&des->X[i*lf->mi[MP]],NULL,(INT)0);
  }
  gam = findgam(lf,des);
  h = 0.0; m = 0;
  for (i=0; i<lf->mi[MN]; i++)
  { des->w[m] = weightmm(des->wd[i],&des->X[i*lf->mi[MP]],lf->mi,gam);
    if (des->w[m]>0)
    { if (des->di[i]>h) h = des->di[i];
      des->ind[m] = i;
      m++;
    }
  }
  des->n = m;
  return(h);
}
