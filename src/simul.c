/*
 *   Copyright (c) 1998 Lucent Technologies.
 *   See README file for details.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "local.h"

static double k0[3], pen, sig2;
static INT corr, side;
extern arstruct aru;

void goldensec(f,des,tr,eps,xm,ym,meth)
double (*f)(), eps, *xm, *ym;
INT meth;
struct design *des;
struct tree *tr;
{ double x[4], y[4], xx[11], yy[11];
  INT i, im;
  xx[0] = tr->dp[DFXH];
  if (xx[0]<=0) { ERROR(("regband: initialize h>0")) return; }
  for (i=0; i<=10; i++)
  { if (i>0) xx[i] = (1+GOLDEN)*xx[i-1];
    yy[i] = f(xx[i],des,tr,meth);
    if ((i==0) || (yy[i]<yy[im])) im = i;
  }
  if (im==0) im = 1;
  if (im==10)im = 9;
  x[0] = xx[im-1]; y[0] = yy[im-1];
  x[1] = xx[im];   y[1] = yy[im];
  x[3] = xx[im+1]; y[3] = yy[im+1];
  x[2] = GOLDEN*x[3]+(1-GOLDEN)*x[0];
  y[2] = f(x[2],des,tr,meth);
  while (x[3]-x[0]>eps)
  { if (y[1]<y[2])
    { x[3] = x[2]; y[3] = y[2];
      x[2] = x[1]; y[2] = y[1];
      x[1] = GOLDEN*x[0]+(1-GOLDEN)*x[3];
      y[1] = f(x[1],des,tr,meth);
    }
    else
    { x[0] = x[1]; y[0] = y[1];
      x[1] = x[2]; y[1] = y[2];
      x[2] = GOLDEN*x[3]+(1-GOLDEN)*x[0];
      y[2] = f(x[2],des,tr,meth);
    }
  }
  im = 0;
  for (i=1; i<4; i++) if (y[i]<y[im]) im = i;
  *xm = x[im]; *ym = y[im];
}

double qmax(x0,x1,f0,f1,d0,d1)
double x0, x1, f0, f1, d0, d1;
{ double a, b, c, h, u, z1, z2;
  h = x1-x0;
  a = h*d0;
  b = 6*(f1-f0)-h*(4*d0+2*d1);
  c = 3*h*(d0+d1)-6*(f1-f0);
  if (fabs(c)>1.0e-30)
  { u = b*b-4*a*c;
    if (u<0) return(-1.2344);
    u = sqrt(u);
    z1 = (-b+u)/(2*c);
    z2 = (-b-u)/(2*c);
    if (c<0) return(x0 + h*((z1>z2) ? z1:z2));
        else return(x0 + h*((z1<z2) ? z1:z2));
  }
  z1 = -d0/(2*b);
  return(x0+h*z1);
}

double max2(des,tr,x0,x1,f0,f1,d0,d1)
struct design *des;
struct tree *tr;
double x0, x1, f0, f1, d0, d1;
{ double x2, f2, d2, l0, l1, ll, lli;
  l0 = x0; l1 = x1; d2 = 1;
  while (fabs(d2)>1.0e-10)
  { x2 = qmax(x0,x1,f0,f1,d0,d1);
    if ((x2<l0) | (x2>l1) | (x2==-1.2344))
      x2 = (l0+l1)/2;
    tr->xev[0] = x2;
    procv(des,tr,0);
    f2 = tr->coef[0]; d2 = tr->coef[tr->nvm];
    ll = tr->nlx[0]; lli = tr->nlx[tr->nvm];
    d2 = (ll*d2-lli*f2)*SGN(f2)/(ll*ll);
    f2 = fabs(f2)/ll;
    if (d2>0) l0 = x2; else l1 = x2;
    x0 = x1; x1 = x2;
    f0 = f1; f1 = f2;
    d0 = d1; d1 = d2;
  }
  return(f2);
}

double cbbound(des,tr)
struct design *des;
struct tree *tr;
{ INT i, k, nvm, d;
  double dlt, fm, f0, f1, llx[MXDIM], nlx, *coef, *cd, max;
  d = tr->mi[MDIM]; nvm = tr->nvm;
  coef = tr->coef; cd = &coef[nvm];
  max = 0;
  for (i=0; i<tr->nv; i++)
  { nlx = tr->nlx[i];
    for (k=0; k<d; k++) llx[k] = tr->nlx[(k+1)*tr->nvm+i];
    for (k=0; k<d; k++)
      cd[k*nvm+i] = (nlx*cd[k*nvm+i]-llx[k]*coef[i])*SGN(coef[i])/(nlx*nlx);
    coef[i] = fabs(coef[i])/nlx;
    if (coef[i]>max) max = coef[i];
  }
  if ((d==1)&&(ident==1)) for (i=1; i<tr->nv; i++)
    if ((cd[i-1]>=0) && (cd[i]<=0))
    { dlt = tr->xev[i*d] - tr->xev[(i-1)*d];
      f0 = coef[i-1]+dlt*cd[i-1];
      f1 = coef[i]-dlt*cd[i];
      if (MAX(f0,f1)>max)
      { fm = max2(des,tr,tr->xev[(i-1)*d],tr->xev[i*d],coef[i-1],coef[i],cd[i-1],cd[i]);
        if (fm>max) max = fm;
    } }
  return(max);
}

void linkt(th,tg,l2,l3)
double th, *l2, *l3;
INT tg;
{ double p;
  switch (tg&63)
  { case TGAUS: *l2 = 1.0; *l3 = 0.0; return;
    case TPOIS: *l2 = *l3 = exp(th); return;
    case TLOGT: p = expit(th);
                *l2 = p*(1-p); *l3 = *l2*(1-2*p);
                return;
    case TGAMM: *l2 = *l3 = 0.0; return;
    default: ERROR(("compr: unknown tg %d",tg)) return;
  }
}

void adj(lf,des,mu,A,B)
struct tree *lf;
struct design *des;
double *mu, *A, *B;
{ INT i, j, k, l, p;
  double l2, l3, t0, ul[5];
  p = des->p;
  for (j=0; j<p; j++) mu[j] = 0.0;
  for (i=0; i<p*p; i++) B[i] = 0.0;
  for (i=0; i<p*p*p; i++) A[i] = 0.0;
  for (l=0; l<lf->mi[MN]; l++)
  { locfit(lf,des,&lf->x[0][l],0.0,0); /* careful for lf */
    linkt(des->cf[0],lf->mi[MTG],&l2,&l3);
    makelxd(lf,des,&lf->x[0][l],ul,0,NULL,0,2);
    t0 = 0;
    for (j=0; j<p; j++) t0 += ul[j]*ul[j];
    for (j=0; j<p; j++) mu[j] += l3*t0*ul[j];
    for (i=0; i<p; i++)
      for (j=0; j<p; j++)
        for (k=0; k<p; k++)
          A[(i*p+j)*p+k] += l3*ul[i]*ul[j]*ul[k];
  }
  for (l=0; l<lf->mi[MN]; l++)
  { locfit(lf,des,&lf->x[0][l],0.0,0); /* careful for lf */
    linkt(des->cf[0],lf->mi[MTG],&l2,&l3);
    makelxd(lf,des,&lf->x[0][l],ul,0,NULL,0,2);
    for (i=0; i<p; i++)
    { t0 = 0;
      for (j=0; j<p; j++)
        for (k=0; k<p; k++)
          t0 += A[(i*p+j)*p+k]*ul[i]*ul[j];
      for (j=0; j<p; j++)
        B[i*p+j] += l3*ul[j]*t0;
    }
  }
  for (j=0; j<p; j++) mu[j] *= -0.5;
}

double compr(lf,des)
struct tree *lf;
struct design *des;
{ INT i, j, p;
  double t0, l2, l3, sm, su, st, wi[5], mu[5], trc[6];
  p = des->p;
  sm = su = 0.0;
  for (i=0; i<p; i++) mu[i] = 0.0;
  for (i=0; i<lf->mi[MN]; i++)
  { locfit(lf,des,&lf->x[0][i],0.0,0); /* careful for lf */
    ldf(lf,des,trc,0,lf->mi,&t0);
    linkt(des->cf[0],lf->mi[MTG],&l2,&l3);
    makelxd(lf,des,&lf->x[0][i],wi,0,NULL,0,2);
    st = 0.0;
    for (j=0; j<p; j++) mu[j] += l3*t0*wi[j];
    su += fabs(l3)*t0*sqrt(t0);
  }
  for (i=0; i<p; i++) sm += mu[i]*mu[i];
  sm = sqrt(sm/4);
printf("su %8.5f ",su);
  su = su/3*(1+2*k0[0]/(PI*exp(1.0))+2/(exp(1.0)*sqrt(PI)));
printf("%8.5f\n",su);
printf("bias corr: sm %8.5f  su %8.5f\n",sm,su);
  return(sm+su);
}

INT nqmax(x,f,z,xl,mbe)
double *x, *f, *xl;
INT z, mbe; /* z=0, first call; z=1, subsequent calls  mbe, may be end */
{ INT im;
  double c, h1, h2, d1, d2;
/* printf("%6.4f %8.5f  %6.4f %8.5f  %6.4f %8.5f  %6.4f %8.5f\n",x[0],f[0],x[1],f[1],x[2],f[2],x[3],f[3]); */
  if (z==1)
  { if (f[3]<f[2])
    { x[3] = (x[3]+x[0])/2; 
      return(0);
    }
    im = 2;
    if (f[3]>f[1]) { im = 1; f[2] = f[1]; x[2] = x[1]; }
    if (f[3]>f[0]) { im = 0; f[1] = f[0]; x[1] = x[0]; }
    f[im] = f[3]; x[im] = x[3];
  }
  else /* sort on f */
  { if (f[0]<f[1])
    { x[3] = x[0]; x[0] = x[1]; x[1] = x[3];
      f[3] = f[0]; f[0] = f[1]; f[1] = f[3];
    }
    if (f[1]<f[2])
    { x[3] = x[2]; x[2] = x[1]; x[1] = x[3];
      f[3] = f[2]; f[2] = f[1]; f[1] = f[3];
    }
    if (f[0]<f[1])
    { x[3] = x[0]; x[0] = x[1]; x[1] = x[3];
      f[3] = f[0]; f[0] = f[1]; f[1] = f[3];
    }
  }

  d2 = f[2]-f[0]; d1 = f[1]-f[0];
  h2 = x[2]-x[0]; h1 = x[1]-x[0];
  if (d2*h1-d1*h2==0) /* linear, return x0 */
  { x[3] = x[0]; f[3] = f[0];
    return(0);
  }

  x[3] = x[0]+(d2*h1*h1-d1*h2*h2)/(d2*h1-d1*h2);
  if (mbe)
  { if (x[3]<xl[0]) { x[3] = xl[0]; return(1); }
    if (x[3]>xl[1]) { x[3] = xl[1]; return(1); }
    c = ((x[2]-x[1])*f[0]-h2*f[1]+h1*f[2])/((x[2]-x[1])*h2*h1);
    if (c>0.0) return(1); /* local minimum */
  }

  while ((x[3]<xl[0]) | (x[3]>xl[1])) x[3] = (x[3]+x[0])/2;
  if (fabs(x[3]-x[0])*10<fabs(h1)) x[3] = (x[3]+x[1])/2;
  if (fabs(x[3]-x[1])*10<fabs(h1)) x[3] = (x[3]+x[0])/2;
  if (fabs(x[3]-x[2])*10<fabs(h2)) x[3] = (x[3]+x[0])/2;
  return(0);
}

#ifdef CVERSION

double b2(th,tg,w)
double th, w;
INT tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(w);
    case TPOIS: return(w*exp(th));
    case TLOGT:
      y = expit(th);
      return(w*y*(1-y));
  }
  ERROR(("b2: invalid family %d",tg)) return(0.0);
}

double b3(th,tg,w)
double th, w;
INT tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(0.0);
    case TPOIS: return(w*exp(th));
    case TLOGT:
      y = expit(th);
      return(w*y*(1-y)*(1-2*y));
  }
  ERROR(("b3: invalid family %d",tg)) return(0.0);
}

double b4(th,tg,w)
double th, w;
INT tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(0.0);
    case TPOIS: return(w*exp(th));
    case TLOGT:
      y = expit(th); y = y*(1-y);
      return(w*y*(1-6*y));
  }
  ERROR(("b4: invalid family %d",tg)) return(0.0);
}

double cumulant(lf,des,x,w)
struct tree *lf;
struct design *des;
double *x, w;
{ double c1, c2, c3, c4, c5, c6, c7, c8, c9;
  double b2i, b3i, b3j, b4i, p2, s[10], *ui, *uj;
  double ss, si, sj, uii, uij, ujj, k1, k2, k4;
  INT i, j, *mi;
  c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = c9 = 0;
  k1 = 0;
  mi = lf->mi;

  locfit(lf,des,x,0.0,0);
  makelxd(lf,des,x,s,0,NULL,0,2);
  ss = innerprod(s,s,mi[MP]);

  for (i=0; i<mi[MN]; i++)
  { b2i = b2(des->th[i],mi[MTG],prwt(lf,i));
    b3i = b3(des->th[i],mi[MTG],prwt(lf,i));
    b4i = b4(des->th[i],mi[MTG],prwt(lf,i));
    ui = &lf->L[i*mi[MP]];
    si = innerprod(s,ui,mi[MP]);
    uii= innerprod(ui,ui,mi[MP]);
    if (lf_error) return(0.0);

    c2 += b4i*si*si*uii;
    c6 += b4i*si*si*si*si;
    c7 += b3i*si*uii;
    c8 += b3i*si*si*si;
    c9 += b2i*b2i*si*si*si*si;
    k1 += b3i*si*(si*si/ss-uii);
/* printf("b3i %8.5f si %8.5f ss %8.5f uii %8.5f k1 %8.5f\n",b3i,si,ss,uii,k1); */

    /* i=j components */
    c1 += b3i*b3i*si*si*uii*uii;
    c3 += b3i*b3i*si*si*si*si*uii;
    c4 += b3i*b3i*si*si*uii*uii;

    for (j=i+1; j<mi[MN]; j++)
    { b3j = b3(des->th[j],mi[MTG],prwt(lf,j));
      uj = &lf->L[j*mi[MP]];
      sj = innerprod(s,uj,mi[MP]);
      uij= innerprod(ui,uj,mi[MP]);
      ujj= innerprod(uj,uj,mi[MP]);

      c1 += 2*b3i*b3j*si*sj*uij*uij;
      c3 += 2*b3i*b3j*si*si*sj*sj*uij;
      c4 += b3i*b3j*uij*(si*si*ujj+sj*sj*uii);
      if (lf_error) return(0.0);
    }
  }
  c5 = c1;
  c7 = c7*c8;
  c8 = c8*c8;

  c1 /= ss; c2 /= ss; c3 /= ss*ss; c4 /= ss;
  c5 /= ss; c6 /= ss*ss; c7 /= ss*ss; c8 /= ss*ss*ss;
  c9 /= ss*ss;
/* printf("%7.5f %7.5f  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  ",x[0],w,c1,c2,c3,c4,c5,c6,c7,c8,c9); */

  k1 = k1/(2*sqrt(ss));
  k2 = 1+c1/2-c2/2-3*c3+c4/2+c6-c7/2+1.75*c8;
  k4 = -9*c3+3*c6+6*c8+3*c9;
/* printf("kappa: %8.5f %8.5f %8.5f %8.5f\n",k1,k2,sqrt(c8),k4); */
  if (k2<0) WARN(("cumul: k2<0 in cumul"))

  w = (w-k1)/sqrt(k2);
/* printf("%7.5f %7.5f %7.5f  %7.5f  ",k1,k2,k4,w); */
  p2 = -w*(k4*(w*w-3)/24 + c8*((w*w-10)*w*w+15)/72);
/* printf("%7.5f  %7.5f\n",p2,fabs(w)+p2); */
  switch(side)
  { case -1: return(-w*sqrt(k2));
    case  0: return(fabs(w)+p2);
    case  1: return(w*sqrt(k2));
  }
}

void procvscb(des,lf,v)
struct design *des;
struct tree *lf;
INT v;
{ double mean, w;
  procv(des,lf,v);
  mean = dareval(&aru,&lf->xev[v*lf->mi[MDIM]]);
  switch(lf->mi[MLINK])
  { case LIDENT: break;
    case LLOG:   mean = log(mean); break;
    case LLOGIT: mean = logit(mean); break;
    default: ERROR(("procvscb: invalid link %d",lf->mi[MLINK])) return;
  }
  w = (lf->coef[v]-mean)/lf->nlx[v];
  if (corr)
  { lf->coef[v] = cumulant(lf,des,&lf->xev[v*lf->mi[MDIM]],w);
    return;
  }
  switch(side)
  { case -1:lf->coef[v] = -w; return;
    case 0: lf->coef[v] = fabs(w); return;
    case 1: lf->coef[v] = w; return;
  }
}

void scbmax(lf,des,co,si)
struct tree *lf;
struct design *des;
INT co, si;
{ double max, xmx, x[4], f[4], kap[3];
  INT c, i, im, nv, v, dv[MXDIM], mbe;
  if ((co) & (ident==0))
  { WARN(("Correction doesn't work with locfit; ignoring."))
    co = 0;
  }
  corr = co;
  side = si;
  if (corr)
  { v = lf->mi[MEV]; lf->mi[MEV] = EDATA;
    i = calcp(lf->mi[MDEG],lf->mi[MDIM],lf->mi[MKT]);
    checkvl(&lf->L,&lf->ll,lf->mi[MN]*i);
    evaluator(des,lf,procvhatm);
    lf->mi[MEV] = v;
  }
  evaluator(des,lf,procvscb);
  if (lf_error) return;
  max = 0.0; im = 0;
  nv = lf->nv;
  for (i=0; i<nv; i++)
    if (lf->coef[i]>max) { max = lf->coef[i]; im = i; }
  if (nv>=4)
  { mbe = 0;
    xmx = lf->xev[im];
    if (im==0) { mbe = 1; im = 1; }
    if (im==nv-1) { mbe = 1; im = nv-2; }
    if (im<=1) v = 3; else v = 0;
    x[0] = lf->xev[im];   f[0] = lf->coef[im];
    x[1] = lf->xev[im-1]; f[1] = lf->coef[im-1];
    x[2] = lf->xev[im+1]; f[2] = lf->coef[im+1];
    i = 0;
    do
    { c = nqmax(x,f,(i>0),lf->fl,mbe);
      i = 1;
      if (!c)
      { lf->xev[v] = x[3]; procvscb(des,lf,v);
        f[3] = lf->coef[v];
        if (f[3]>max) { max = f[3]; xmx = x[3]; }
      }
    } while ((!c) && (f[0]-f[2]>0.000001));
  }

  /* now, compute kappa0 and tail probability */

  constants(des,lf,kap,dv,0);
  printf("xmx: %10.6f  max: %10.6f  k0 %10.6f %10.6f  pr %10.6f\n",xmx,max,kap[0],kap[1],tailp(max,kap,1+(side==0),lf->mi[MDIM],0.0));
}

#endif

double compsda(x,h,n)
double *x, h;
INT n;
/* n/(n-1) * int( fhat''(x)^2 dx ); bandwidth h */
{ INT i, j;
  double ik, sd, z;
  ik = wint(1,NULL,0,WGAUS);
  sd = 0;

  for (i=0; i<n; i++)
    for (j=i; j<n; j++)
    { z = (x[i]-x[j])/h;
      sd += (2-(i==j))*Wconv4(z,WGAUS)/(ik*ik);
    }
  sd = sd/(n*(n-1)*h*h*h*h*h);
  return(sd);
}

void xxxx(x,h,r,d,c,k,ker,n)
double *x, h, *r, *d, c;
INT k, ker, n;
{ INT i, j;
  double pen, s, r0, r1, d0, d1, ik;
  *r = 0.0; *d = 0.0;
  ik = wint(1,NULL,0,ker);
  switch(k)
  { case 1: /* aic */
      pen = 2.0;
      for (i=0; i<n; i++)
      { r0 = d0 = 0.0;
        for (j=0; j<n; j++)
        { s = (x[i]-x[j])/h;
          r0 += W(s,ker);
          d0 += s*s*Wd(s,ker);
        }
        d0 = -(d0+r0)/(n*h*h*ik);  /* d0 = d/dh fhat(xi) */
        r0 /= n*h*ik;              /* r0 = fhat(xi) */
        *r += -2*log(r0)+pen*W(0.0,ker)/(n*h*ik*r0);
        *d += -2*d0/r0-pen*W(0.0,ker)/(n*h*ik*r0)*(d0/r0+1.0/h);
      }
      return;
    case 2: /* ocv */
      for (i=0; i<n; i++)
      { r0 = 0.0; d0 = 0.0;
        for (j=0; j<n; j++) if (i!=j)
        { s = (x[i]-x[j])/h;
          r0 += W(s,ker);
          d0 += s*s*Wd(s,ker);
        }
        d0 = -(d0+r0)/((n-1)*h*h);
        r0 = r0/((n-1)*h);
        *r -= log(r0);
        *d -= d0/r0;
      }
      return;
    case 3: /* lscv */
      r0 = r1 = d0 = d1 = 0.0;
      for (i=0; i<n; i++)
      { for (j=i+1; j<n; j++)
        { s = (x[i]-x[j])/h;
          r1 += 2*W(s,ker);
          d1 += -2*s*s*Wd(s,ker);
          r0 += 2*Wconv(s,ker);
          d0 += -2*s*s*Wconv1(s,ker);
        }
      }
      d1 -= r1;
      r0 += n*Wconv(0.0,ker);
      d0 -= r0;
      *r = r0/(n*n*h*ik)   - 2*r1/(n*(n-1)*h);
      *d = d0/(n*n*h*h*ik) - 2*d1/(n*(n-1)*h*h);
      return;
    case 4: /* bcv */
      r0 = d0 = 0.0;
      for (i=0; i<n; i++)
        for (j=i+1; j<n; j++)
        { s = (x[i]-x[j])/h;
          r0 += 2*Wconv4(s,ker);
          d0 += 2*s*Wconv5(s,ker);
        }
      d0 = (-d0-r0)/(n*n*h*h*ik*ik);
      r0 = r0/(n*n*h*ik*ik);
      j = 2.0;
      d1 = wint(1,&j,1,ker);
      r1 = Wconv(0.0,ker);
      *r = (d1*d1*r0/4+r1/(n*h))/(ik*ik);
      *d = (d1*d1*d0/4-r1/(n*h*h))/(ik*ik);
      return;
    case 5: /* sjpi */
      s = c*exp(5*log(h)/7)/SQRT2;
      d0 = compsda(x,s,n);
      *r = d0; /* this is S(alpha) in SJ */
      *d = exp(log(Wikk(WGAUS,0)/(d0*n))/5)-h;
      return;
    case 6: /* gas-k-k */
      s = exp(log(1.0*n)/10)*h;
      d0 = compsda(x,s,n);
      *r = d0;
      *d = exp(log(Wikk(WGAUS,0)/(d0*n))/5)-h;
      return;
  }
  ERROR(("xxxx: what???")) return;
}

double widthsj(x,lambda,n)
double *x, lambda;
INT n;
{ double ik, a, b, td, sa, z, c, c1, c2, c3;
  INT i, j;
  a = GFACT*0.920*lambda*exp(-log((double)n)/7)/SQRT2;
  b = GFACT*0.912*lambda*exp(-log((double)n)/9)/SQRT2;
  ik = wint(1,NULL,0,WGAUS);

  td = 0;
  for (i=0; i<n; i++)
    for (j=i; j<n; j++)
    { z = (x[i]-x[j])/b;
      td += (2-(i==j))*Wconv6(z,WGAUS)/(ik*ik);
    }

  td = -td/(n*(n-1));
  j = 2.0;
  c1 = Wconv4(0.0,WGAUS);
  c2 = wint(1,&j,1,WGAUS);
  c3 = Wconv(0.0,WGAUS);  /* (2*c1/(c2*c3))^(1/7)=1.357 */
  sa = compsda(x,a,n);
  c = b*exp(log(c1*ik/(c2*c3*GFACT*GFACT*GFACT*GFACT)*sa/td)/7)*SQRT2;
  return(c);
}

double esolve(x,j,h0,h1,k,c,ker,n)
double *x, h0, h1, c;
INT j, k, ker, n;
{ double h[7], d[7], r[7], min, minh, fact;
  INT i, nc;
  min = 1.0e30; minh = 0.0;
  fact = 1.00001;
  h[6] = h0; xxxx(x,h[6],&r[6],&d[6],c,j,ker,n);
  if (lf_error) return(0.0);
  nc = 0;
  for (i=0; i<k; i++)
  { h[5] = h[6]; r[5] = r[6]; d[5] = d[6];
    h[6] = h0*exp((i+1)*log(h1/h0)/k);
    xxxx(x,h[6],&r[6],&d[6],c,j,ker,n);
    if (lf_error) return(0.0);
    if (d[5]*d[6]<0)
    { h[2] = h[0] = h[5]; d[2] = d[0] = d[5]; r[2] = r[0] = r[5];
      h[3] = h[1] = h[6]; d[3] = d[1] = d[6]; r[3] = r[1] = r[6];
      while ((h[3]>fact*h[2])|(h[2]>fact*h[3]))
      { h[4] = h[3]-d[3]*(h[3]-h[2])/(d[3]-d[2]);
        if ((h[4]<h[0]) | (h[4]>h[1])) h[4] = (h[0]+h[1])/2;
        xxxx(x,h[4],&r[4],&d[4],c,j,ker,n);
        if (lf_error) return(0.0);
        h[2] = h[3]; h[3] = h[4];
        d[2] = d[3]; d[3] = d[4];
        r[2] = r[3]; r[3] = r[4];
        if (d[4]*d[0]>0) { h[0] = h[4]; d[0] = d[4]; r[0] = r[4]; }
                    else { h[1] = h[4]; d[1] = d[4]; r[1] = r[4]; }
      }
      if (j>=4) return(h[4]); /* first min for BCV etc */
      if (r[4]<=min) { min = r[4]; minh = h[4]; }
      nc++;
    }
  }
  if (nc==0) minh = (r[5]<r[6]) ? h0 : h1;
  return(minh);
}

double dnk(x,k)
double x;
INT k;
{ double f;
  switch(k)
  { case 0: f = 1; break;
    case 1: f = -x; break;
    case 2: f = x*x-1; break;
    case 3: f = x*(x*x-3); break;
    case 4: f = 3-x*x*(6-x*x); break;
    case 5: f = -x*(15-x*x*(10-x*x)); break;
    case 6: f = -15+x*x*(45-x*x*(15-x*x)); break;
    default: ERROR(("dnk: k=%d too large",k)) return(0.0);
  }
  return(f*exp(-x*x/2)/S2PI);
}

/*
double ise(da,h,k,mu,sg,l)
struct data *da;
double h, *mu, sg;
INT k, l;
{ double i0, i1, i2, z, h2s;
  INT fact[]={1,1,2,6,24,120,720}, m4[]={1,-4,16,-64,256,-1024},
      m2[]={1,-2,4,-8,16,-32,64,-128};
  INT i, j, r, s, n;
  n = da->n;
  i0 = 0;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
    { z = (da->x[0][i]-da->x[0][j])/(h*SQRT2);
      for (r=0; r<l; r++)
        for (s=0; s<l; s++)
          i0 += dnk(z,2*(r+s))/(fact[s]*fact[r]*m4[s+r]);
    }
  i0 /= (n*n*SQRT2*h);
  i1 = 0;
  for (i=0; i<n; i++)
    for (j=0; j<k; j++)
    { z = (da->x[0][i]-mu[j])/sqrt(sg*sg+h*h);
      for (s=0; s<l; s++)
      { h2s = (s==0) ? 1.0 : h2s*h*h/(sg*sg+h*h);
        i1 += dnk(z,2*s)*h2s/(fact[s]*m2[s]);
      }
    }
  i1 /= n*k*sqrt(sg*sg+h*h);
  i2 = 0;
  for (i=0; i<k; i++)
    for (j=0; j<k; j++)
      i2 += dnk((mu[i]-mu[j])/(SQRT2*sg),0);
  i2 /= k*k*SQRT2*sg;
  return(i0-2*i1+i2);
}
*/

double locai(h,des,tr)
double h;
struct design *des;
struct tree *tr;
{ double cp;
  tr->dp[DALP] = h;
  evaluator(des,tr,procv);
  ressumm(tr,des);
  cp = -2*tr->dp[DLK]+pen*tr->dp[DT0];
  return(cp);
}

void kdescore(x,n,h,nh,z,m)
double *x, *h, *z;
INT *n, *nh, *m;
{ double d;
  INT i;
  for (i=0; i<*nh; i++)
    switch(*m)
    { case 1: xxxx(x,h[i],&z[i],&d,0.0,3,WGAUS,*n); break;
      case 2: xxxx(x,h[i],&z[i],&d,0.0,4,WGAUS,*n); break;
      case 3: d = compsda(x,h[i],*n);
              z[i] = exp(log(Wikk(WGAUS,0)/(*n*d))/5);
              break;
    }
}

void kdeselect(band,x,ind,h0,h1,meth,nm,ker,n)
double h0, h1, *band, *x;
INT *ind, *meth, nm, ker, n;
{ double scale, c;
  INT i, k;
  k = n/4;
  for (i=0; i<n; i++) ind[i] = i;
  scale = kordstat(x,n+1-k,n,ind) - kordstat(x,k,n,ind);
  c = widthsj(x,scale,n);
  for (k=0; k<nm; k++)
    band[k] = esolve(x,meth[k],h0,h1,10,c,ker,n);
}

double loccp(h,des,tr,m) /* m=1: cp    m=2: gcv */
double h;
struct design *des;
struct tree *tr;
{ double cp;
  INT dg;
  tr->dp[DALP] = 0;
  tr->dp[DFXH] = h;
  dg = tr->mi[MDEG]; tr->mi[MDEG] = tr->mi[MDEG0];
  evaluator(des,tr,procv);
  ressumm(tr,des);
  if (m==1)
    cp = -2*tr->dp[DLK]/sig2 - tr->mi[MN] + 2*tr->dp[DT0];
  else cp = -2*tr->mi[MN]*tr->dp[DLK]/((tr->mi[MN]-tr->dp[DT0])*(tr->mi[MN]-tr->dp[DT0]));
  printf("h %8.5f  deg %2d  rss %8.5f  trl %8.5f  cp: %8.5f\n",h,tr->mi[MDEG],-2*tr->dp[DLK],tr->dp[DT0],cp);
  tr->mi[MDEG0] = tr->mi[MDEG]; tr->mi[MDEG] = dg;
  return(cp);
}

double cp(des,tr,meth)
struct design *des;
struct tree *tr;
INT meth;
{ double hm, ym;
  goldensec(loccp,des,tr,0.001,&hm,&ym,meth);
  return(hm);
}

double gkk(des,tr)
struct design *des;
struct tree *tr;
{ double h, h5, nf, th;
  INT i, j, n, dg0, dg1;
  tr->mi[MEV]=EDATA;
  tr->dp[DALP] = 0;
  n = tr->mi[MN];
  dg0 = tr->mi[MDEG0];     /* target degree */
  dg1 = dg0+1+(dg0%2==0);  /* pilot degree */
  nf = exp(log(1.0*n)/10); /* bandwidth inflation factor */
  h = tr->dp[DFXH];        /* start bandwidth */
  for (i=0; i<=10; i++)
  { tr->mi[MDEG] = dg1;
    tr->dp[DFXH] = h*nf;
    evaluator(des,tr,procv);
    th = 0;
    for (j=10; j<n-10; j++)
      th += tr->coef[dg1*n+j]*tr->coef[dg1*n+j];
th *= n/(n-20.0);
    h5 = sig2*Wikk(tr->mi[MKER],dg0)/th;
    h = exp(log(h5)/(2*dg1+1));
/* printf("pilot %8.5f  sel %8.5f\n",tr->dp[DFXH],h); */
  }
  return(h);
}

double rsw(des,tr,kk)
struct design *des;
struct tree *tr;
INT *kk;
{ INT i, j, k, nmax, nvm, n, mk, ev, dg0, dg1;
  double rss[6], cp[6], th22, dx, d2, hh;
  nmax = 5;
  ev = tr->mi[MEV];  tr->mi[MEV] = EGRID;
  mk = tr->mi[MKER]; tr->mi[MKER]= WRECT;
  dg0 = tr->mi[MDEG0];
  dg1 = 1 + dg0 + (dg0%2==0);
  tr->mi[MDEG]= 4;
  for (k=nmax; k>0; k--)
  { tr->mg[0] = k;
    tr->fl[0] = 1.0/(2*k); tr->fl[1] = 1-1.0/(2*k);
    tr->dp[DALP] = 0; tr->dp[DFXH] = 1.0/(2*k);
    evaluator(des,tr,procv);
    nvm = tr->nvm;
    rss[k] = 0;
    for (i=0; i<k; i++) rss[k] += -2*tr->lik[i];
  }
  n = tr->mi[MN]; k = 1;
  for (i=1; i<=nmax; i++)
  { /* cp[i] = (n-5*nmax)*rss[i]/rss[nmax]-(n-10*i); */
    cp[i] = rss[i]/sig2-(n-10*i);
    if (cp[i]<cp[k]) k = i;
  }
  *kk = k;
  tr->mg[0] = k;
  tr->fl[0] = 1.0/(2*k); tr->fl[1] = 1-1.0/(2*k);
  tr->dp[DALP] = 0; tr->dp[DFXH] = 1.0/(2*k);
  evaluator(des,tr,procv);
  tr->mi[MKER] = mk; tr->mi[MEV] = ev;
  nvm = tr->nvm;
  th22 = 0;
  for (i=10; i<n-10; i++)
  { j = floor(k*tr->x[0][i]);
    if (j>=k) j = k-1;
    dx = tr->x[0][i]-tr->xev[j];
    if (dg1==2)
      d2 = tr->coef[2*nvm+j]+dx*tr->coef[3*nvm+j]+dx*dx*tr->coef[4*nvm+j]/2;
    else d2 = tr->coef[4*nvm+j];
    th22 += d2*d2;
  }
  hh = Wikk(mk,dg0)*sig2/th22*(n-20.0)/n;
  return(exp(log(hh)/(2*dg1+1)));
}

void rband(des,tr,hhat,meth,nmeth,kk)
struct design *des;
struct tree *tr;
double *hhat;
INT *meth, *nmeth, *kk;
{ INT i, deg;
  double h0;

  /* first, estimate sigma^2 */
  deg = tr->mi[MDEG]; tr->mi[MDEG] = 2;
  h0 = tr->dp[DFXH];  tr->dp[DFXH] = 0.05;
printf("alp: %8.5f  h: %8.5f  deg %2d  ev %2d\n",tr->dp[DALP],tr->dp[DFXH],tr->mi[MDEG],tr->mi[MEV]);
  evaluator(des,tr,procv);
  ressumm(tr,des);
  tr->mi[MDEG] = deg; tr->dp[DFXH] = h0;
  sig2 = tr->dp[DRV]; 
  printf("sd est: %8.5f\n",sqrt(tr->dp[DRV]));

  for (i=0; i<*nmeth; i++)
  { switch(meth[i])
    { case 1: hhat[i] = cp(des,tr,1);
              break;
      case 2: hhat[i] = cp(des,tr,2);
              break;
      case 3: hhat[i] = gkk(des,tr);
              break;
      case 4: hhat[i] = rsw(des,tr,kk);
              break;
      default: hhat[i] = 0;
    }
    tr->dp[DFXH] = h0;
    tr->mi[MDEG] = deg;
  }
}
