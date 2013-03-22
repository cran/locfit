/*
 *   Copyright (c) 1998 Lucent Technologies.
 *   See README file for details.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "local.h"

static double *wd;
INT ident=0, debug;
double robscale;

extern INT (*itype)();

void prefit() {}

void unitvec(x,k,p)
double *x;
INT k, p;
{ INT i;
  for (i=0; i<p; i++) x[i] = 0.0;
  x[k] = 1.0;
}

double median(x,n)
double *x;
INT n;
{ INT i, j, k, m;
  double lo, hi, s;
  lo = hi = x[0];
  k = n/2;
  for (i=0; i<n; i++)
  { lo = MIN(lo,x[i]);
    hi = MAX(hi,x[i]);
  }
  for (i=0; i<n; i++)
  { if ((x[i]>lo) & (x[i]<hi))
    { s = x[i]; m = 0;
      for (j=0; j<n; j++) m += (x[j]<=s);
      if (m==k) return(s);
      if (m<k) lo = s;
      if (m>k) hi = s;
    }
  }
  printf("median failed\n");
  return(0.0);
}

double robustscale(lf,des,c,tg)
struct tree *lf;
struct design *des;
INT c, tg;
{ INT i, ii, j, p;
  static double rs, os;
  if (tg==TROBT) return(1.0); /* not quasi - no scale */
  if (c) os = rs;
  p = des->p;
  for (i=0; i<des->n; i++)
  { des->th[i] = 0; ii = des->ind[i];
    for (j=0; j<p; j++) des->th[i] += des->cf[j]*des->X[i*p+j];
    des->res[i] = fabs(resp(lf,ii)-des->th[i])*sqrt(prwt(lf,ii));
  }
  rs = 3*median(des->res,des->n);

/* For stability, estimate must not change too much between calls. */
  if (c)
  { if (rs<0.8*os) rs = 0.8*os;
    if (rs>1.2*os) rs = 1.2*os;
  }
printf("robust scale: %8.5f\n",rs);
  return(rs);
}

double likereg(lf,des,x,zz,h,mi)
struct tree *lf;
struct design *des;
double *x, h;
INT *mi, *zz;
{ INT i, ii, j, p;
  double th, lk, ww, link[LLEN];
  lk = 0.0; p = des->p;
  for (i=0; i<p*p; i++) des->Z[i] = 0;
  for (i=0; i<p; i++) des->res[i] = 0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    des->th[i] = th = base(lf,ii)+innerprod(des->cf,&des->X[i*p],p);
    *zz = links(th,resp(lf,ii),mi[MTG],mi[MLINK],link,cens(lf,ii),prwt(lf,ii));
    if (*zz) return(0.0);
    ww = des->w[i];
    lk += ww*link[ZLIK];
    for (j=0; j<p; j++)
      des->res[j] += des->X[i*p+j]*ww*link[ZDLL];
    addouter(des->Z,&des->X[i*p],&des->X[i*p],p,ww*link[ZDDLL]);
  }
  switch(des->sm)
  { case 1:
      for (j=0; j<p; j++)
      { des->dg[j] = des->Z[j*p+j];
        if (des->dg[j]>0) des->dg[j] = 1/sqrt(des->dg[j]);
      }
      for (i=0; i<p; i++)
        for (j=0; j<p; j++)
          des->Z[i*p+j] *= des->dg[i]*des->dg[j];
      eigen(des->Z,des->Q,p,mi[MMXIT]);
      break;
    case 2: choldec(des->Z,p);
            break;
    default: ERROR(("likereg: no method for solving %d",des->sm))
  }
  return(lk);
}

void vxtwx(des,v,k) /* (X^T W X)^{-1} v */
struct design *des;
double *v;
INT k;
{ INT i;
  switch(des->sm)
  { case 1: /* eigenvalues on corr matrix */
      for (i=0; i<des->p; i++) v[i] *= des->dg[i];
      svdsolve(v,des->f2,des->Q,des->Z,des->Q,des->p,1.0e-8);
      for (i=0; i<des->p; i++) v[i] *= des->dg[i];
      return;
    case 2: /* chol decomposition of cov matrix */
      cholsolve(v,des->Z,k,des->p);
      return;
  }
  ERROR(("vxtwx: unknown method %d",des->sm))
}

void vmat(lf,des,M,tr)  /* M = X^T W^2 V X */
struct tree *lf;
struct design *des;
double *M, *tr;
{ INT i, ii, p, nk;
  double link[LLEN], h, ww;
  p = des->p;
  for (i=0; i<p*p; i++) M[i] = 0;
  if ((lf->mi[MTG]<THAZ) && (lf->mi[MLINK]==LLOG))
  { nk = -1;
    switch(lf->mi[MKER])
    { case WGAUS: nk = WGAUS; h = des->h/SQRT2; break;
      case WRECT: nk = WRECT; h = des->h; break;
      case WEPCH: nk = WBISQ; h = des->h; break;
      case WBISQ: nk = WQUQU; h = des->h; break;
      case WTCUB: nk = W6CUB; h = des->h; break;
    }
    if (nk != -1)
    { (des->itype)(NULL,M,des->P,lf,des->cf,h,lf->mi,nk);
      if (lf->mi[MTG]==TDEN)
        for (i=0; i<p*p; i++) M[i] *= des->na;
      tr[0] = des->ss[0];
      tr[1] = M[0]; /* n int W e^<a,A> */
      return;
    }
  }
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    links(des->th[i],resp(lf,ii),lf->mi[MTG],lf->mi[MLINK],link,cens(lf,ii),prwt(lf,ii));
    ww = SQR(des->w[i])*link[ZDDLL];
    tr[0] += des->w[i];
    tr[1] += SQR(des->w[i]);
    addouter(M,&des->X[i*p],&des->X[i*p],p,ww);
  }
}

void ldf(lf,des,tr,z,mi,t0)
struct tree *lf;
struct design *des;
INT z, *mi;
double *tr, *t0;
{ INT i, ii, j, k, p;
  double *m2, *V, ww, link[LLEN];
  tr[0] = tr[1] = tr[2] = tr[3] = tr[4] = tr[5] = 0.0;
  m2 = des->V; V = des->P; p = des->p;
  vmat(lf,des,m2,tr);  /* M = X^T W^2 V X  tr0=sum(W) tr1=sum(W*W) */
  for (j=0; j<p*p; j++) V[j] = m2[j];
  for (i=0; i<p; i++)
  { vxtwx(des,&m2[i*p],p);
    tr[2] += m2[i*(p+1)];  /* tr (XTWVX)^{-1}(XTW^2VX) */
  }

  switch(z)
  { case 1:
      unitvec(des->f1,0,p);
      vxtwx(des,des->f1,p);
      for (i=0; i<p; i++)
        for (j=0; j<p; j++)
        { tr[4] += m2[i*p+j]*m2[j*p+i];  /* tr(M^2) */
          tr[5] += des->f1[i]*V[i*p+j]*des->f1[j]; /* var(thetahat) */
        }
      tr[5] = sqrt(tr[5]);
      for (i=0; i<p*p; i++) m2[i] = 0.0;
      for (i=0; i<des->n; i++)
      { ii = des->ind[i];
        links(des->th[i],resp(lf,ii),mi[MTG],mi[MLINK],
          link,cens(lf,ii),prwt(lf,ii));
        ww = SQR(des->w[i])*des->w[i]*link[ZDDLL];
        addouter(m2,&des->X[i*p],&des->X[i*p],p,ww);
      }
      for (i=0; i<p; i++)
      { vxtwx(des,&m2[i*p],p);
        tr[3] += m2[i*(p+1)];
      }
      return;
    case 0:
      unitvec(des->f1,0,p);
      vxtwx(des,des->f1,p);
      for (i=0; i<=mi[MDIM]; i++) t0[i] = des->f1[i]; /* influnce and deriv */
      choldec(V,p);                                /* LL^T = X^T W^2 X */
      for (i=0; i<p; i++) vxtwx(des,&V[i*p],p); /* V = (X^TWX)^{-1}L */
      vxtwx(des,des->f1,p);
      for (i=0; i<p; i++)
      { for (j=0; j<p; j++)
        { m2[i*p+j] = 0;
          for (k=0; k<p; k++)
            m2[i*p+j] += V[k*p+i]*V[k*p+j]; /* ith column of covariance */
        }
      }
      if ((mi[MTG]==TDEN) && (mi[MLINK]==LIDENT))
        for (i=0; i<p*p; i++) m2[i] /= SQR(des->na);
      return;
  }
}

INT locfit(lf,des,x,h,noit)
struct tree *lf;
struct design *des;
double *x, h;
INT noit;
{ double f, s0, s1, lk1, lk0, *X, *Z, *cf, u[MXDIM], (*like)(), link[LLEN];
  INT i, ii, it, j, m, p, d, tg, cv, oob, zz, pf, *mi, st;
  mi  = lf->mi;
  d = mi[MDIM]; p = des->p; tg = mi[MTG];
  X = des->X; Z = des->Z; cf = des->cf;
  des->sm = 1+(mi[MDEG0]<mi[MDEG]);
  des->h = h;
  for (i=0; i<p; i++) cf[i] = des->res[i] = 0;

  m = des->n;
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    for (j=0; j<d; j++) u[j] = lf->x[j][ii]-x[j];
    fitfun(u,&X[i*p],lf->sca,d,mi[MDEG],mi[MKT],NULL,0,lf->sty);
  }

  if (tg<=THAZ)
  { st = densinit(lf,des,x,h,cf,mi,m);
    if (st)
    { cf[0] = NOSLN;
      for (j=1; j<p; j++) cf[j] = 0.0;
      return(1);
    }
    like = likeden;
  }
  else
  { s0 = s1 = des->ss[0] = 0;
    for (i=0; i<m; i++)
    { ii = des->ind[i];
      links(base(lf,ii),resp(lf,ii),mi[MTG],LINIT,link,cens(lf,ii),prwt(lf,ii));
      if (mi[MTG]==TCIRC)
      { s0 += des->w[i]*link[ZDLL]; /* sin term */
        s1 += des->w[i]*link[ZLIK]; /* cos term */
      }
      else
      { s1 += des->w[i]*link[ZDLL];
        s0 += des->w[i]*prwt(lf,ii);
      }
    }
    if (s0==0) return(5); /* no observations with W>0 */
    st = 0;
    des->ss[0] = (mi[MTG]==TCIRC) ? atan2(s0,s1) : s1/s0;
    switch(mi[MLINK])
    { case LIDENT: cf[0] = des->ss[0]; break;
      case LLOG:   if (des->ss[0]<=0) { cf[0] = -1000; st = 1; }
                   if (st==0) cf[0] = log(des->ss[0]);
                   break;
      case LLOGIT: if (des->ss[0]<=0) { cf[0] = -1000; st = 1; }
                   if (des->ss[0]>=1) { cf[0] = 1000;  st = 1; }
                   if (st==0)  cf[0] = logit(des->ss[0]);
                   break;
      case LINVER: cf[0] = 1/des->ss[0]; break;
      case LSQRT:  cf[0] = sqrt(des->ss[0]); break;
      default: ERROR(("locfit: invalid link %d",mi[MLINK]))
    }
    if (st==1)
    { for (j=1; j<p; j++) cf[j] = 0.0;
      return(1);
    }
    like = likereg;
  }
  
  if ((mi[MTG]&63)==TROBT) robscale = robustscale(lf,des,0,mi[MTG]);
  lk1 = like(lf,des,x,&zz,h,mi);

  for (it=0; it<mi[MMXIT]; it++)
  { for (j=0; j<p; j++)          /* store old coeffs and X^TW\dot{l} */
    { des->oc[j] = cf[j];
      des->f1[j] = des->res[j];
    }
    vxtwx(des,des->f1,p);        /* (XTWVX)^-1 (XTW\dot{l} */

    f = 1; lk0 = lk1;
    do
    { for (j=0; j<p; j++)
        cf[j] = des->oc[j]+f*des->f1[j]; /* update coefficients */
      lk1 = like(lf,des,x,&zz,h,mi);
      { if (((mi[MTG]&63)==TGAUS) && (noit))
        { des->llk = lk1;
          /* quick return for Gaussian - but iteration improves numerics */
          return(0);
        }
      }
      f /= 2.0;
      if (f<1.0e-5) /* seems to be a bad direction! */
      { WARN(("locfit: f problem")); break; }
    } while (zz || (lk0>lk1+1.0e-3));

    oob = pf = cv = 0;
    switch (tg&63)  /* check for convergence */
    { case TGAUS: cv = (it>0); /* max of 1 iteration */
                  break;
      case TLOGT: oob = 0; /* ((mi[MLINK]==LLOGIT) && (fabs(cf[0])>100)); */
                  pf = lk1>-1.0e-5*s0;
                  cv = (fabs(lk1-lk0)<1.0e-6*s0);
                  break;
      case TPOIS:
      case TGEOM:
      case TWEIB:
      case TGAMM: oob = ((mi[MLINK]==LLOG) && (fabs(cf[0])>100));
                  pf = lk1>-1.0e-5*s0;
                  cv = (fabs(lk1-lk0)<1.0e-6*s0);
                  break;
      case TCIRC: cv = (fabs(lk1-lk0)<1.0e-6);
                  break;
      case TDEN:
      case TRAT:
      case THAZ: if (mi[MLINK]==LLOG)
                 { oob = fabs(cf[0])>100;
                   cv = (fabs(lk1-lk0)<1.0e-6);
                 } else cv = 1;
                 break;
      case TROBT:
        robscale = robustscale(lf,des,1,mi[MTG]);
        cv = (fabs(des->f1[0])<=1.0e-6*robscale);
        lk1 = like(lf,des,x,&zz,h,mi);
        break;
      default: ERROR(("locfit: unknown target %d",tg))
    }
    if (oob) for (j=1; j<p; j++) cf[j] = 0.0;
    if (cv | oob | pf) break;
  }
  if (mi[MTG]==TDEN)
  { if (mi[MLINK]==LLOG) des->cf[0] -= log(des->na);
    else for (i=0; i<mi[MP]; i++) des->cf[i] /= des->na;
  }
  des->llk = lk1;
  if (oob) return(2); /* out of bounds */
  if (pf)  return(3); /* perfect fit */
  if (it==mi[MMXIT])
  { WARN(("locfit: not converged"))
    return(4);
  }
  return(0);
}

/*
void dercor(lf,des,h)
struct tree *lf;
struct design *des;
double h;
{ double s1, dc[MXDIM], z, wd, link[LLEN];
  INT i, ii, j, m, p, d; 
  if (lf->mi[MTG]<=THAZ) return;
  d = lf->mi[MDIM];
  p = des->p; m = des->n;
  unitvec(des->f1,0,p);
  vxtwx(des,des->f1,p);
  for (i=0; i<lf->mi[MDIM]; i++) dc[i] = 0.0;
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    s1 = 0.0;
    for (j=0; j<p; j++) s1 += des->f1[j]*des->X[i*p+j];
    links(des->th[i],resp(lf,ii),lf->mi[MTG],lf->mi[MLINK],link,cens(lf,ii),prwt(lf,ii));
    switch(lf->mi[MKT])
    { case KSPH:
        wd = Wd(des->di[ii]/h,lf->mi[MKER]);
        z = s1*wd*link[ZDLL];
        for (j=0; j<d; j++) dc[j] += z*des->X[i*p+j+1];
        break;
      case KPROD:
        z = s1*link[ZDLL];
        for (j=0; j<d; j++)
        { wd = fabs(des->X[i*p+j+1]/(h*lf->sca[j]));
          wd = des->w[i]*Wd(wd,lf->mi[MKER])/W(wd,lf->mi[MKER]);
          for (j=0; j<d; j++) dc[j] += z*wd*des->X[i*p+j+1];
        }
        break;
      case KANG:
        wd = Wd(des->di[ii]/h,lf->mi[MKER]);
        z = s1*wd*link[ZDLL];
        dc[0] += z*des->X[i*p+1];
    }
  }
  for (j=0; j<d; j++) des->cf[j+1] -= dc[j]/SQR(h*lf->sca[j]);
} */

void dercor(lf,des,x,h)
struct tree *lf;
struct design *des;
double h, *x;
{ double s1, dc[MXDIM], z, wd, link[LLEN];
  INT i, ii, j, m, p, d, *mi; 
  mi = lf->mi;
  if (mi[MTG]<=THAZ) return;
  d = mi[MDIM];
  p = des->p; m = des->n;
  unitvec(des->f1,0,p);
  vxtwx(des,des->f1,p);
  for (i=0; i<d; i++) dc[i] = 0.0;

  /* correction term is e1^T (XTWVX)^{-1} XTW' ldot. */
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    s1 = 0.0;
    for (j=0; j<p; j++) s1 += des->f1[j]*des->X[i*p+j];
    links(des->th[i],resp(lf,ii),mi[MTG],mi[MLINK],link,cens(lf,ii),prwt(lf,ii));
    for (j=0; j<d; j++)
    { wd = des->w[i]*weightd(lf->x[j][ii]-x[j],lf->sca[j],d,mi[MKER],mi[MKT],h,lf->sty[j],des->di[ii]);
      dc[j] += s1*wd*link[ZDLL];
    }

  }
  for (j=0; j<d; j++) des->cf[j+1] += dc[j];
}

void nnresproj(lf,des,u,m,p,mi)
struct tree *lf;
struct design *des;
double *u;
INT m, p, *mi;
{ INT i, j;
  double s, link[LLEN];
  for (i=0; i<p; i++) des->f1[i] = 0;
  for (j=0; j<m; j++)
  { i = des->ind[j];
    links(des->th[j],resp(lf,i),mi[MTG],mi[MLINK],link,cens(lf,i),prwt(lf,i));
    for (i=0; i<p; i++) des->f1[i] += link[ZDDLL]*des->X[j*p+i]*u[j];
  }
  vxtwx(des,des->f1,p);
  for (i=0; i<m; i++)
  { s = 0;
    for (j=0; j<p; j++) s += des->f1[j]*des->X[i*p+j];
    u[i] -= s*des->w[i];
  }
}

void makelxd(lf,des,x,lx,deg,deriv,nd,ty)
struct tree *lf;
struct design *des;
double *lx, *x;
INT deg, nd, *deriv, ty;
/* deg=0: l(x) only.
   deg=1: l(x), l'(x) (approx/exact ? mi[MDC] );
   deg=2: l(x), l'(x), l''(x);
   ty = 1: e1 (X^T WVX)^{-1} X^T W        -- hat matrix
   ty = 2: e1 (X^T WVX)^{-1} X^T WV^{1/2} -- scb's
*/
{ double w, *X, *Z, *lxd, *lxdd, wdd, wdw, *ulx, link[LLEN], h;
  double dfx[MXDIM], ei[MXDIM], hs[MXDIM];
  INT i, ii, j, k, l, m, d, p, *mi;
  mi = lf->mi; h = des->h;
  wd = des->wd;
  d = mi[MDIM]; p = des->p; X = des->X; Z = des->Z;
  ulx = des->res;
  m = des->n;
  for (i=0; i<m; i++) lx[i] = 0;
  for (i=0; i<d; i++) hs[i] = h*lf->sca[i];
  if (deg>0)
  { lxd = &lx[m];
    for (i=0; i<m*d; i++) lxd[i] = 0.0;
    if (deg>1)
    { lxdd = &lxd[d*m];
      for (i=0; i<m*d*d; i++) lxdd[i] = 0.0;
  } }
  for (i=0; i<d; i++) ei[i] = 0.0;
  fitfun(ei,des->f1,lf->sca,d,mi[MDEG],mi[MKT],deriv,nd,lf->sty); /* c(0) */
  vxtwx(des,des->f1,p);   /* c(0) (X^TWX)^{-1} */
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    for (j=0; j<p; j++) lx[i] += des->f1[j]*X[i*p+j]; /* c(0)(XTWX)^{-1}X^T */
    if ((deg>0) && (mi[MDC]))
    { wd[i] = Wd(des->di[ii]/h,mi[MKER]);
      for (j=0; j<d; j++)
      { dfx[j] = lf->x[j][ii]-x[j];
        lxd[j*m+i] = lx[i]*des->w[i]*weightd(dfx[j],lf->sca[j],
          d,mi[MKER],mi[MKT],h,lf->sty[j],des->di[ii]);
             /* c(0) (XTWX)^{-1}XTW' */
      }
      if (deg>1)
      { wdd = Wdd(des->di[ii]/h,mi[MKER]);
        for (j=0; j<d; j++)
          for (k=0; k<d; k++)
          { w = (des->di[ii]==0) ? 0 : h/des->di[ii];
            w = wdd*(x[k]-lf->x[k][ii])*(x[j]-lf->x[j][ii])*w*w/(hs[k]*hs[k]*hs[j]*hs[j]);
            if (j==k) w += wd[i]/(hs[j]*hs[j]);
            lxdd[(j*d+k)*m+i] = lx[i]*w;
              /* c(0)(XTWX)^{-1}XTW'' */
          }
      }
    }
    lx[i] *= des->w[i];
  }
  if ((deg==2) && (mi[MDC]))
  { for (i=0; i<d; i++)
    { deriv[nd] = i;
      fitfun(ei,des->f1,lf->sca,d,mi[MDEG],mi[MKT],deriv,nd+1,lf->sty);
      for (k=0; k<m; k++)
      { ii = des->ind[k];
        links(des->th[k],resp(lf,ii),mi[MTG],mi[MLINK],link,cens(lf,ii),prwt(lf,ii));
        for (j=0; j<p; j++)
          des->f1[j] -= link[ZDDLL]*lxd[i*m+k]*X[k*p+j];
        /* c'(x)-c(x)(XTWX)^{-1}XTW'X */
      }
      vxtwx(des,des->f1,p); /* (...)(XTWX)^{-1} */
      for (j=0; j<m; j++)
      { ulx[j] = 0;
        for (k=0; k<p; k++) ulx[j] += des->f1[k]*X[j*p+k]; /* (...)XT */
      }
      for (j=0; j<d; j++)
        for (k=0; k<m; k++)
        { ii = des->ind[k];
          dfx[j] = lf->x[j][ii]-x[j];
          wdw = des->w[k]*weightd(dfx[j],lf->sca[j],d,mi[MKER],mi[MKT],h,
            lf->sty[j],des->di[ii]);
          lxdd[(i*d+j)*m+k] += ulx[k]*wdw;
          lxdd[(j*d+i)*m+k] += ulx[k]*wdw;
        }
        /* + 2(c'-c(XTWX)^{-1}XTW'X)(XTWX)^{-1}XTW' */
    }
    for (j=0; j<d*d; j++) nnresproj(lf,des,&lxdd[j*m],m,p,mi);
        /* * (I-X(XTWX)^{-1} XTW */
  }
  if (deg>0)
  { if (mi[MDC]) for (j=0; j<d; j++) nnresproj(lf,des,&lxd[j*m],m,p,mi);
             /* c(0)(XTWX)^{-1}XTW'(I-X(XTWX)^{-1}XTW) */
    for (i=0; i<d; i++)
    { deriv[nd]=i;
      fitfun(ei,des->f1,lf->sca,d,mi[MDEG],mi[MKT],deriv,nd+1,lf->sty);
      vxtwx(des,des->f1,p);
      for (k=0; k<m; k++)
        for (l=0; l<p; l++)
          lxd[i*m+k] += des->f1[l]*X[k*p+l]*des->w[k];
            /* add c'(0)(XTWX)^{-1}XTW */
    }
  }
  if (deg==2)
  { for (i=0; i<d; i++)
    { deriv[nd]=i;
      for (j=0; j<d; j++)
      { deriv[nd+1]=j;
        fitfun(ei,des->f1,lf->sca,d,mi[MDEG],mi[MKT],deriv,nd+2,lf->sty);
        vxtwx(des,des->f1,p);
        for (k=0; k<m; k++)
          for (l=0; l<p; l++)
            lxdd[(i*d+j)*m+k] += des->f1[l]*X[k*p+l]*des->w[k];
        /* + c''(x)(XTWX)^{-1}XTW */
      }
    }
  }
  k = 1+d*(deg>0)+d*d*(deg==2);
  if (ty==1) return;
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    links(des->th[i],resp(lf,ii),mi[MTG],mi[MLINK],link,cens(lf,ii),prwt(lf,ii));
    link[ZDDLL] = sqrt(fabs(link[ZDDLL]));
    for (j=0; j<k; j++) lx[j*m+i] *= link[ZDDLL];
  }
  return;
}
