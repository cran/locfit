/*
 *   Copyright (c) 1998 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"
#include <stdio.h>
#include <math.h>

/*
  Defines the weight functions and related quantities used
  in LOCFIT.
*/

static double alp[10], gam;
static int debug;

double setmmwt(des,lf,a,gam)
struct design *des;
struct tree *lf;
double *a, gam;
{ double ip, w0, w1, sw, wt;
  INT i, p;
  sw = 0.0;
  p = lf->mi[MP];
  for (i=0; i<lf->mi[MN]; i++)
  { ip = innerprod(a,&des->X[i*p],p);
    wt = prwt(lf,i);
    w0 = wt*(ip - gam*des->di[i]);
    w1 = wt*(ip + gam*des->di[i]);
    des->w[i] = 0.0;
    if (w0>0) { des->w[i] = w0; sw += w0*w0; }
    if (w1<0) { des->w[i] = w1; sw += w1*w1; }
  }
  if (sw==0.0) /* disaster! no nonzero weights */
    return(0.0);
  return(sw/2-a[0]);
}

double findab(gam,lf,des,xev,a)
struct tree *lf;
struct design *des;
double gam, *a, *xev;
{ double *A, *z, sl, sw, sw0;
  INT i, j, k, n, p, p1, done;
  if (debug) printf("findab: gam = %8.5f\n",gam);
  n = lf->mi[MN]; p = lf->mi[MP];
  A = des->Z; z = des->f1;
  
  do
  { sw = sw0 = setmmwt(des,lf,a,gam);
    if (sw==0.0) a[0] *= 2.0;
  } while (sw==0.0);

  for (j=0; j<lf->mi[MMXIT]; j++)
  { /* compute z = Newton-Raphson increment */
    z[0] = 1.0;
    for (i=1; i<p; i++) z[i] = 0.0;
    for (i=0; i<p*p; i++) A[i] = 0.0;
    for (i=0; i<p; i++) alp[i] = a[i];
    for (i=0; i<n; i++) if (des->w[i]!=0.0)
    { addouter(A,&des->X[i*p],&des->X[i*p],p,1.0);
      for (k=0; k<p; k++) z[k] -= des->w[i]*des->X[i*p+k];
    }
    for (i=0; i<p; i++)
    { des->dg[i] = A[i*p+i];
      if (des->dg[i]>0) des->dg[i] = 1/sqrt(A[i*p+i]);
    }
    for (i=0; i<p; i++)
      for (k=0; k<p; k++)
        A[i*p+k] *= des->dg[i]*des->dg[k];
    eigen(A,des->Q,p,lf->mi[MMXIT]);
    des->sm = 1;
    vxtwx(des,z,p);

    /* compute new parameter vector; 
       compute new weights;
       check for decreased sw */
    done = 0;
    while (!done)
    { for (i=0; i<p; i++) a[i] = alp[i]+z[i];
      sw = setmmwt(des,lf,a,gam);
      if (sw>sw0+1.0e-6)
        for (i=0; i<p; i++) z[i] *= 0.5;
      else done = 1;
    }
    if (done)
    { for (i=0; i<p; i++) alp[i] = a[i];
      done = ((j>0) & (fabs(sw-sw0)<0.00000001));
      sw0 = sw;
    }
    if (done) j = lf->mi[MMXIT];
  }
  if (!done) { WARN(("findab not converged")); }
  sl = 0.0;
  for (i=0; i<n; i++) sl += fabs(des->w[i])*des->di[i];
  p1 = factorial(lf->mi[MDEG]+1);
  sl *= lf->dp[DALP]*lf->dp[DALP] / (p1*p1);
  if (debug) printf("sl = %8.5f\n",sl);
  return(sl);
}

void findgam(lf,des,xev)
struct tree *lf;
struct design *des;
double *xev;
{ double g[5], z[5], *a;
  INT i, n;
  n = lf->mi[MN];
  a = des->f2;
  g[0] = 0.0;
  a[0] = 1.0/n;
  for (i=1; i<lf->mi[MP]; i++) a[i] = 0.0;
  z[0] = findab(g[0],lf,des,xev,a);
  g[4] = 1.0;
  z[4] = findab(g[4],lf,des,xev,a);
  while (z[4]>g[4])
  { z[0] = z[4]; g[0] = g[4];
    g[4] *= 2.0;
    z[4] = findab(g[4],lf,des,xev,a);
  }
  g[1] = g[0]; z[1] = z[0];
  g[2] = g[4]; z[2] = z[4];
  do
  { g[3] = g[2] + (g[2]-g[1])*(z[2]-g[2])/(z[1]-g[1]-z[2]+g[2]);
    if ((g[3]<=g[0]) | (g[3]>=g[4])) g[3] = (g[0]+g[4])/2;
    z[3] = findab(g[3],lf,des,xev,a);
    if (z[3]>g[3]) { g[0] = g[3]; z[0] = z[3]; }
        else { g[4] = g[3]; z[4] = z[3]; }
    z[1] = z[2]; z[2] = z[3];
    g[1] = g[2]; g[2] = g[3];
  } while ((fabs(g[3]-z[3])>0.0000001) & (g[0]<g[4]));
  gam = g[3];
}

/* The weight functions themselves.
   Used everywhere. */
double W(u,ker)
double u;
INT ker;
{ u = fabs(u);
  switch(ker)
  { case WRECT: return((u>1) ? 0.0 : 1.0);
    case WEPCH: return((u>1) ? 0.0 : 1-u*u);
    case WBISQ: if (u>1) return(0.0);
                u = 1-u*u; return(u*u);
    case WTCUB: if (u>1) return(0.0);
                u = 1-u*u*u; return(u*u*u);
    case WTRWT: if (u>1) return(0.0);
                u = 1-u*u; return(u*u*u);
    case WQUQU: if (u>1) return(0.0);
                u = 1-u*u; return(u*u*u*u);
    case WTRIA: if (u>1) return(0.0);
                return(1-u);
    case W6CUB: if (u>1) return(0.0);
                u = 1-u*u*u; u = u*u*u; return(u*u);
    case WGAUS: return(exp(-SQR(GFACT*u)/2.0));
    case WMINM: ERROR(("WMINM in W"));
                return(0.0);
    case WPARM: return(1.0);
  }
  return(0.0);
}

double weightmm(di,ff,mi)
double di, *ff;
INT *mi;
{ double y1, y2, ip;
  ip = innerprod(ff,alp,mi[MP]);
  y1 = ip-gam*di; if (y1>0) return(y1/ip);
  y2 = ip+gam*di; if (y2<0) return(y2/ip);
  return(0.0);
}

double minmax(lf,des,xev)
struct tree *lf;
struct design *des;
double *xev;
{ double h, u[MXDIM], w;
  INT i, j, m;
  for (i=0; i<lf->mi[MN]; i++)
  { for (j=0; j<lf->mi[MDIM]; j++) u[j] = lf->x[j][i]-xev[j];
    w = des->di[i]; des->di[i] = 1.0;
    for (j=0; j<=lf->mi[MDEG]; j++) des->di[i] *= w;
    fitfun(u,&des->X[i*lf->mi[MP]],lf->sca,lf->mi[MDIM],lf->mi[MDEG],
      lf->mi[MKT],NULL,(INT)0,lf->sty);
  }
  findgam(lf,des,xev);
  h = 0.0; m = 0;
  for (i=0; i<lf->mi[MN]; i++)
  { des->w[m] = weightmm(des->di[i],&des->X[i*lf->mi[MP]],lf->mi);
    if (des->w[m]>0)
    { if (des->di[i]>h) h = des->di[i];
      des->ind[m] = i;
      m++;
    }
  }
  des->n = m;
  h = exp(log(h)/(1+lf->mi[MDEG]));
  return(h);
}

double weight(u,sc,d,ker,kt,h,sty)
double *u, *sc, h;
INT d, ker, kt, *sty;
{ INT i;
  double w, z;
  w = 1; z = 0.0;
  for (i=0; i<d; i++)
  { switch(sty[i])
    { case KLEF:
        if (u[i]>0.0) return(0);
        if (d==1) return(W(-u[0]/(h*sc[0]),ker));
        if (kt==KPROD)
          w *= W(-u[i]/(h*sc[i]),ker);
        else z += u[i]*u[i]/(sc[i]*sc[i]);
        break;
      case KRIG:
        if (u[i]<0.0) return(0);
        if (d==1) return(W(u[0]/(h*sc[0]),ker));
        if (kt==KPROD)
          w *= W(u[i]/(h*sc[i]),ker);
        else z += u[i]*u[i]/(sc[i]*sc[i]);
        break;
      case KANG:
        if (d==1) return(W(2*fabs(sin(u[0]/(2*sc[0])))/h,ker));
        if (kt==KPROD)
          w *= W(2*fabs(sin(u[i]/(2*sc[i])))/h,ker);
        else z += 4*sin(u[i]/(2*sc[i]))*sin(u[i]/(2*sc[i]));
        break;
      default:
        if (d==1) return(W(u[0]/(h*sc[0]),ker));
        if (kt==KPROD)
          w *= W(u[i]/(h*sc[i]),ker);
        else z += u[i]*u[i]/(sc[i]*sc[i]);
    }
  }
  if (kt==KPROD) return(w);
  return(W(sqrt(z)/h,ker));
}

double sgn(x)
double x;
{ if (x>0) return(1.0);
  if (x<0) return(-1.0);
  return(0.0);
}

double WdW(u,ker) /* W'(u)/W(u) */
double u;
INT ker;
{ double eps=1.0e-10;
  if (ker==WGAUS) return(-GFACT*GFACT*u);
  if (ker==WPARM) return(0.0);
  if (fabs(u)>=1) return(0.0);
  switch(ker)
  { case WRECT: return(0.0);
    case WTRIA: return(-sgn(u)/(1-fabs(u)+eps));
    case WEPCH: return(-2*u/(1-u*u+eps));
    case WBISQ: return(-4*u/(1-u*u+eps));
    case WTRWT: return(-6*u/(1-u*u+eps));
    case WTCUB: return(-9*sgn(u)*u*u/(1-u*u*fabs(u)+eps));
  }
  ERROR(("WdW: invalid kernel"));
  return(0.0);
}

/* deriv. weights .. spherical, product etc
   u, sc, sty needed only in relevant direction
   Acutally, returns (d/dx W(||x||/h) ) / W(.)
*/
double weightd(u,sc,d,ker,kt,h,sty,di)
double u, sc, h, di;
INT d, ker, kt, sty;
{ if (sty==KANG)
  { if (kt==KPROD)
      return(-WdW(2*sin(u/(2*sc)),ker)*cos(u/(2*sc))/(h*sc));
    if (di==0.0) return(0.0);
    return(-WdW(di/h,ker)*sin(u/sc)/(h*sc*di));
  }
  if (kt==KPROD)
    return(-WdW(u/(h*sc),ker)/(h*sc));
  if (di==0.0) return(0.0);
  return(-WdW(di/h,ker)*u/(h*di*sc*sc));
}

double weightdd(u,sc,d,ker,kt,h,sty,di,i0,i1)
double *u, *sc, h, di;
INT d, ker, kt, *sty, i0, i1;
{ double w;
  w = 1;
  if (kt==KPROD)
  {
    w = WdW(u[i0]/(h*sc[i0]),ker)*WdW(u[i1]/(h*sc[i1]),ker)/(h*h*sc[i0]*sc[i1]);
  }
  return(0.0);
}

/* Derivatives W'(u)/u.
   Used in simult. conf. band computations,
   and kernel density bandwidth selectors. */
double Wd(u,ker)
double u;
INT ker;
{ double v;
  if (ker==WGAUS) return(-SQR(GFACT)*exp(-SQR(GFACT*u)/2));
  if (ker==WPARM) return(0.0);
  if (fabs(u)>1) return(0.0);
  switch(ker)
  { case WEPCH: return(-2.0);
    case WBISQ: return(-4*(1-u*u));
    case WTCUB: v = 1-u*u*u;
                return(-9*v*v*u);
    case WTRWT: v = 1-u*u;
                return(-6*v*v);
    default: ERROR(("Invalid kernel %d in Wd",ker))
  }
  return(0.0);
}

/* Second derivatives W''(u)-W'(u)/u.
   used in simult. conf. band computations in >1 dimension. */
double Wdd(u,ker)
double u;
INT ker;
{ double v;
  if (ker==WGAUS) return(SQR(u*GFACT*GFACT)*exp(-SQR(u*GFACT)/2));
  if (ker==WPARM) return(0.0);
  if (u>1) return(0.0);
  switch(ker)
  { case WBISQ: return(12*u*u);
    case WTCUB: v = 1-u*u*u;
                return(-9*u*v*v+54*u*u*u*u*v);
    case WTRWT: return(24*u*u*(1-u*u));
    default: ERROR(("Invalid kernel %d in Wdd",ker))
  }
  return(0.0);
}

/* int u1^j1..ud^jd W(u) du.
   Used for local log-linear density estimation.
   Assume all j_i are even.
   Also in some bandwidth selection.
*/
double wint(d,j,nj,ker)
INT d, *j, nj, ker;
{ double I, z;
  INT k, dj;
  dj = d;
  for (k=0; k<nj; k++) dj += j[k];
  switch(ker) /* int_0^1 u^(dj-1) W(u)du  */
  { case WRECT: I = 1.0/dj; break;
    case WEPCH: I = 2.0/(dj*(dj+2)); break;
    case WBISQ: I = 8.0/(dj*(dj+2)*(dj+4)); break;
    case WTCUB: I = 162.0/(dj*(dj+3)*(dj+6)*(dj+9)); break;
    case WTRWT: I = 48.0/(dj*(dj+2)*(dj+4)*(dj+6)); break;
    case WTRIA: I = 1/(dj*(dj+1)); break;
    case WQUQU: I = 384.0/(dj*(dj+2)*(dj+4)*(dj+6)*(dj+8)); break;
    case W6CUB: I = 524880.0/(dj*(dj+3)*(dj+6)*(dj+9)*(dj+12)*(dj+15)*(dj+18)); break;
    case WGAUS: switch(d)
                { case 1: I = S2PI/GFACT; break;
                  case 2: I = 2*PI/(GFACT*GFACT); break;
                  default: I = exp(d*log(S2PI/GFACT)); /* for nj=0 */
                }
                for (k=0; k<nj; k++) /* deliberate drop */
                  switch(j[k])
                  { case 4: I *= 3.0/(GFACT*GFACT);
                    case 2: I /= GFACT*GFACT;
                  }
                return(I);
    default: ERROR(("Unknown kernel %d in exacint",ker))
  }
  if ((d==1) && (nj==0)) return(2*I); /* common case quick */
  z = (d-nj)*LOGPI/2-LGAMMA(dj/2.0);
  for (k=0; k<nj; k++) z += LGAMMA((j[k]+1)/2.0);
  return(2*I*exp(z));
}

/* taylor series expansion of weight function around x.
   0 and 1 are common arguments, so are worth programming
   as special cases.
   Used in density estimation.
*/
INT wtaylor(f,x,ker)
double *f, x;
INT ker;
{ double v;
  switch(ker)
  { case WRECT:
      f[0] = 1.0;
      return(1);
    case WEPCH:
      f[0] = 1-x*x; f[1] = -2*x; f[2] = -1;
      return(3);
    case WBISQ:
      v = 1-x*x;
      f[0] = v*v;   f[1] = -4*x*v; f[2] = 4-6*v;
      f[3] = 4*x;   f[4] = 1;
      return(5);
    case WTCUB:
      if (x==1.0)
      { f[0] = f[1] = f[2] = 0; f[3] = -27; f[4] = -81; f[5] = -108;
        f[6] = -81; f[7] = -36; f[8] = -9; f[9] = -1; return(10); }
      if (x==0.0)
      { f[1] = f[2] = f[4] = f[5] = f[7] = f[8] = 0;
        f[0] = 1; f[3] = -3; f[6] = 3; f[9] = -1; return(10); }
      v = 1-x*x*x;
      f[0] = v*v*v; f[1] = -9*v*v*x*x; f[2] = x*v*(27-36*v);
      f[3] = -27+v*(108-84*v);         f[4] = -3*x*x*(27-42*v);
      f[5] = x*(-108+126*v);           f[6] = -81+84*v;
      f[7] = -36*x*x; f[8] = -9*x;     f[9] = -1;
      return(10);
    case WTRWT:
      v = 1-x*x;
      f[0] = v*v*v; f[1] = -6*x*v*v; f[2] = v*(12-15*v);
      f[3] = x*(20*v-8); f[4] = 15*v-12; f[5] = -6; f[6] = -1;
      return(7);
    case WTRIA:
      f[0] = 1-x; f[1] = -1;
      return(2);
    case WQUQU:
      v = 1-x*x;
      f[0] = v*v*v*v; f[1] = -8*x*v*v*v; f[2] = v*v*(24-28*v);
      f[3] = v*x*(56*v-32); f[4] = (70*v-80)*v+16; f[5] = x*(32-56*v);
      f[6] = 24-28*v; f[7] = 8*x; f[8] = 1;
      return(9);
    case W6CUB:
      v = 1-x*x*x;
      f[0] = v*v*v*v*v*v;
      f[1] = -18*x*x*v*v*v*v*v;
      f[2] = x*v*v*v*v*(135-153*v);
      f[3] = v*v*v*(-540+v*(1350-816*v));
      f[4] = x*x*v*v*(1215-v*(4050-v*3060));
      f[5] = x*v*(-1458+v*(9234+v*(-16254+v*8568)));
      f[6] = 729-v*(10206-v*(35154-v*(44226-v*18564)));
      f[7] = x*x*(4374-v*(30132-v*(56862-v*31824)));
      f[8] = x*(12393-v*(61479-v*(92664-v*43758)));
      f[9] = 21870-v*(89100-v*(115830-v*48620));
      f[10]= x*x*(26730-v*(69498-v*43758));
      f[11]= x*(23814-v*(55458-v*31824));
      f[12]= 15849-v*(34398-v*18564);
      f[13]= x*x*(7938-8568*v);
      f[14]= x*(2970-3060*v);
      f[15]= 810-816*v;
      f[16]= 153*x*x;
      f[17]= 18*x;
      f[18]= 1;
      return(19);
  }
  ERROR(("Invalid kernel %d in wtaylor",ker))
  return(0);
}

/* convolution int W(x)W(x+v)dx.
   used in kde bandwidth selection.
*/
double Wconv(v,ker)
double v;
INT ker;
{ double v2;
  switch(ker)
  { case WGAUS: return(SQRPI/GFACT*exp(-SQR(GFACT*v)/4));
    case WRECT:
      v = fabs(v);
      if (v>2) return(0.0);
      return(2-v);
    case WEPCH:
      v = fabs(v);
      if (v>2) return(0.0);
      return((2-v)*(16+v*(8-v*(16-v*(2+v))))/30);
    case WBISQ:
      v = fabs(v);
      if (v>2) return(0.0);
      v2 = 2-v;
      return(v2*v2*v2*v2*v2*(16+v*(40+v*(36+v*(10+v))))/630);
  }
  ERROR(("Wconv not implemented for kernel %d",ker))
  return(0.0);
}

/* derivative of Wconv.
   1/v d/dv int W(x)W(x+v)dx
   used in kde bandwidth selection.
*/
double Wconv1(v,ker)
double v;
INT ker;
{ double v2;
  v = fabs(v);
  switch(ker)
  { case WGAUS: return(-0.5*SQRPI*GFACT*exp(-SQR(GFACT*v)/4));
    case WRECT:
      if (v>2) return(0.0);
      return(1.0);
    case WEPCH:
      if (v>2) return(0.0);
      return((-16+v*(12-v*v))/6);
    case WBISQ:
      if (v>2) return(0.0);
      v2 = 2-v;
      return(-v2*v2*v2*v2*(32+v*(64+v*(24+v*3)))/210);
  }
  ERROR(("Wconv1 not implemented for kernel %d",ker))
  return(0.0);
}

/* 4th derivative of Wconv.
   used in kde bandwidth selection (BCV, SJPI, GKK)
*/
double Wconv4(v,ker)
double v;
INT ker;
{ double gv;
  switch(ker)
  { case WGAUS:
      gv = GFACT*v;
      return(exp(-SQR(gv)/4)*GFACT*GFACT*GFACT*(12-gv*gv*(12-gv*gv))*SQRPI/16);
  }
  ERROR(("Wconv4 not implemented for kernel %d",ker))
  return(0.0);
}

/* 5th derivative of Wconv.
   used in kde bandwidth selection (BCV method only)
*/
double Wconv5(v,ker) /* (d/dv)^5 int W(x)W(x+v)dx */
double v;
INT ker;
{ double gv;
  switch(ker)
  { case WGAUS:
      gv = GFACT*v;
      return(-exp(-SQR(gv)/4)*GFACT*GFACT*GFACT*GFACT*gv*(60-gv*gv*(20-gv*gv))*SQRPI/32);
  }
  ERROR(("Wconv5 not implemented for kernel %d",ker))
  return(0.0);
}

/* 6th derivative of Wconv.
   used in kde bandwidth selection (SJPI)
*/
double Wconv6(v,ker)
double v;
INT ker;
{ double gv, z;
  switch(ker)
  { case WGAUS:
      gv = GFACT*v;
      gv = gv*gv;
      z = exp(-gv/4)*(-120+gv*(180-gv*(30-gv)))*0.02769459142;
      gv = GFACT*GFACT;
      return(z*gv*gv*GFACT);
  }
  ERROR(("Wconv6 not implemented for kernel %d",ker))
  return(0.0);
}

/* int W(v)^2 dv / (int v^2 W(v) dv)^2
   used in some bandwidth selectors
*/
double Wikk(ker,deg)
INT ker, deg;
{ switch(deg)
  { case 0:
    case 1: /* int W(v)^2 dv / (int v^2 W(v) dv)^2 */
      switch(ker)
      { case WRECT: return(4.5);
        case WEPCH: return(15.0);
        case WBISQ: return(35.0);
        case WGAUS: return(0.2820947918*GFACT*GFACT*GFACT*GFACT*GFACT);
        case WTCUB: return(34.15211105);
        case WTRWT: return(66.08391608);
      }
    case 2:
    case 3: /* 4!^2/8*int(W1^2)/int(v^4W1)^2
               W1=W*(n4-v^2n2)/(n0n4-n2n2) */
      switch(ker)
      { case WRECT: return(11025.0);
        case WEPCH: return(39690.0);
        case WBISQ: return(110346.9231);
        case WGAUS: return(14527.43412);
        case WTCUB: return(126500.5904);
        case WTRWT: return(254371.7647);
      }
  }
  ERROR(("Wikk not implemented for kernel %d",ker))
  return(0.0);
}
