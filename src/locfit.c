/*
 *   Copyright (c) 1998 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

static double s0, s1;
static INT lf_status;
INT ident=0, (*like)();
double robscale;

void prefit() {}

void unitvec(x,k,p)
double *x;
INT k, p;
{ setzero(x,p);
  x[k] = 1.0;
}

double median(x,n)
double *x;
INT n;
{ INT i, j, lt, eq, gt;
  double lo, hi, s;
  lo = hi = x[0];
  for (i=0; i<n; i++)
  { lo = MIN(lo,x[i]);
    hi = MAX(hi,x[i]);
  }
  if (lo==hi) return(lo);
  lo -= (hi-lo);
  hi += (hi-lo);
  for (i=0; i<n; i++)
  { if ((x[i]>lo) & (x[i]<hi))
    { s = x[i]; lt = eq = gt = 0;
      for (j=0; j<n; j++)
      { lt += (x[j]<s);
        eq += (x[j]==s);
        gt += (x[j]>s);
      }
      if ((2*(lt+eq)>n) && (2*(gt+eq)>n)) return(s);
      if (2*(lt+eq)<=n) lo = s;
      if (2*(gt+eq)<=n) hi = s;
    }
  }
  return((hi+lo)/2);
}

double nrobustscale(lf,des,tg,rs)
lfit *lf;
design *des;
INT tg;
double rs;
{ int i, ii, p;
  double link[LLEN], sc, sd, sw, e;
  if ((tg&64)==0) return(1.0); /* not quasi - no scale */
  if (((tg&128)==0) & (((tg&63)!=TROBT) & ((tg&63)!=TCAUC))) return(1.0);
  p = des->p; sc = sd = sw = 0.0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    des->th[i] = base(lf,ii)+innerprod(des->cf,&des->X[i*p],p);
    e = resp(lf,ii)-des->th[i];
    stdlinks(link,lf,ii,des->th[i],rs);
    sc += des->w[i]*e*link[ZDLL];
    sd += des->w[i]*e*e*link[ZDDLL];
    sw += des->w[i];
  }

  /* newton-raphson iteration for log(s)
     -psi(ei/s) - log(s); s = e^{-th}
  */
  rs *= exp((sc-sw)/(sd+sc));
printf("ns: %8.5f\n",rs);
  return(rs);
}

double robustscale(lf,des,tg)
lfit *lf;
design *des;
INT tg;
{ INT i, ii, p;
  double rs;
  if ((tg&64)==0) return(1.0); /* not quasi - no scale */
  if (((tg&128)==0) & (((tg&63)!=TROBT) & ((tg&63)!=TCAUC))) return(1.0);
  p = des->p;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    des->th[i] = base(lf,ii) + innerprod(des->cf,&des->X[i*p],p);
    des->res[i] = fabs(resp(lf,ii)-des->th[i])*sqrt(prwt(lf,ii));
  }
  rs = median(des->res,des->n);
  if (rs==0.0) rs = 1.0;
printf("rs: %8.5f\n",rs);
  return(rs);
}

double lrobustscale(lf,des,tg)
lfit *lf;
design *des;
INT tg;
{ INT i, ii, p;
  double rs, link[LLEN];
  if ((tg&64)==0) return(1.0); /* not quasi - no scale */
  if (((tg&128)==0) & (((tg&63)!=TROBT) & ((tg&63)!=TCAUC))) return(1.0);
  p = des->p;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    des->th[i] = base(lf,ii) + innerprod(des->cf,&des->X[i*p],p);
    links(des->th[i],resp(lf,ii),tg&127,lf->mi[MLINK],link,cens(lf,ii),prwt(lf,ii),1.0);
    des->res[i] = -2*link[ZLIK];
  }
  rs = sqrt(median(des->res,des->n));
  if (rs==0.0) rs = 1.0;
/* printf("rs: %8.5f\n",rs); */
  return(rs);
}

INT likereg(lf,des,nop)
lfit *lf;
design *des;
INT nop;
{ INT i, ii, j, p, st;
  double lk, ww, link[LLEN];
  xtwxstruc *xtwx;
  lk = 0.0; p = des->p;
  xtwx = &des->xtwx;
  setzero(xtwx->Z,p*p);
  setzero(des->res,p);
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    des->th[i] = base(lf,ii)+innerprod(des->cf,&des->X[i*p],p);
    st = stdlinks(link,lf,ii,des->th[i],robscale);
    if (st>0) return(st);
    if (lf_error) return(LF_ERR);
    ww = des->w[i];
    lk += ww*link[ZLIK];
    for (j=0; j<p; j++)
      des->res[j] += des->X[i*p+j]*ww*link[ZDLL];
    addouter(xtwx->Z,&des->X[i*p],&des->X[i*p],p,ww*link[ZDDLL]);
  }
  if (lf->mi[MDEB]>2) prresp(des->cf,xtwx->Z,p);
  if (lf->mi[MDEB]>1) printf("  likelihood: %8.5f\n",lk);
  xtwx->p = p;
  xtwxdec(xtwx,xtwx->sm,des->res,nop);
  des->llk = lk;
  return(LF_OK);
}

INT vxtwx(xtwx,v,k) /* (X^T W X)^{-1} v */
xtwxstruc *xtwx;
double *v;
INT k;
{ INT i, rank;
  switch(xtwx->sm)
  { case 1: /* eigenvalues on corr matrix */
      for (i=0; i<xtwx->p; i++) v[i] *= xtwx->dg[i];
      rank = svdsolve(v,xtwx->f2,xtwx->Q,xtwx->Z,xtwx->Q,xtwx->p,1.0e-8);
      for (i=0; i<xtwx->p; i++) v[i] *= xtwx->dg[i];
      return(rank);
    case 2: /* chol decomposition of cov matrix */
      cholsolve(v,xtwx->Z,k,xtwx->p);
      return(xtwx->p);
  }
  ERROR(("vxtwx: unknown method %d",xtwx->sm));
  return(LF_ERR);
}

double vxtwxv(xtwx,v)  /* vT (xtwx)^{-1} v */
xtwxstruc *xtwx;
double *v;
{ INT i, j, p;
  double sum;
  if (xtwx->sm!=1) ERROR(("sm problem: vxtwxv"));
  p = xtwx->p;
  sum = 0.0;
  for (i=0; i<p; i++) v[i] *= xtwx->dg[i];
  for (i=0; i<p; i++)
  { xtwx->f2[i] = 0.0;
    for (j=0; j<p; j++) xtwx->f2[i] += xtwx->Q[j*p+i]*v[j];
    if (xtwx->Z[i*p+i]>1.0e-8) sum += xtwx->f2[i]*xtwx->f2[i]/xtwx->Z[i*p+i];
  }
  return(sum);
}

/*
  vmat() computes (after the local fit..) the matrix 
  Z = X^T W^2 V X.
  M = (X^T W V X)^{-1} Z
  Also, for convenience, tr[0] = sum(wi) tr[1] = sum(wi^2).
*/

void vmat(lf,des,M,Z,tr)
lfit *lf;
design *des;
double *M, *Z, *tr;
{ INT i, p, nk, ok;
  double link[LLEN], h, ww;
  p = des->p;
  setzero(M,p*p);

  nk = -1;

  /* for density estimation, use integral rather than
     sum form, if W^2 is programmed...
  */
  if ((lf->mi[MTG]<=THAZ) && (lf->mi[MLINK]==LLOG))
  { switch(lf->mi[MKER])
    { case WGAUS: nk = WGAUS; h = des->h/SQRT2; break;
      case WRECT: nk = WRECT; h = des->h; break;
      case WEPAN: nk = WBISQ; h = des->h; break;
      case WBISQ: nk = WQUQU; h = des->h; break;
      case WTCUB: nk = W6CUB; h = des->h; break;
      case WEXPL: nk = WEXPL; h = des->h/2; break;
    }
  }

  if (nk != -1)
  { ok = lf->mi[MKER]; lf->mi[MKER] = nk;
    (des->itype)(des->xev,M,des->P,lf,des->cf,h);
    lf->mi[MKER] = ok;
    if (lf->mi[MTG]==TDEN) multmatscal(M,lf->dp[DSWT],p*p);
    tr[0] = des->ss[0];
    tr[1] = M[0]; /* n int W e^<a,A> */
  }
  else
  { for (i=0; i<des->n; i++)
    { stdlinks(link,lf,des->ind[i],des->th[i],robscale);
      ww = SQR(des->w[i])*link[ZDDLL];
      tr[0] += des->w[i];
      tr[1] += SQR(des->w[i]);
      addouter(M,&des->X[i*p],&des->X[i*p],p,ww);
    }
  }

  memcpy(Z,M,p*p*sizeof(double));
  for (i=0; i<p; i++)
    vxtwx(&des->xtwx,&M[i*p],p);
}

/* compute influence, variance e.t.c. */
void ldf(lf,des,tr,z,mi,t0)
lfit *lf;
design *des;
INT z, *mi;
double *tr, *t0;
{ INT i, j, k, p;
  double *m2, *V, ww, link[LLEN];
  if (mi[MDEB]>1) printf("  in ldf\n");
  tr[0] = tr[1] = tr[2] = tr[3] = tr[4] = tr[5] = 0.0;
  m2 = des->V; V = des->P; p = des->p;
  vmat(lf,des,m2,V,tr);  /* M = X^T W^2 V X  tr0=sum(W) tr1=sum(W*W) */
  for (i=0; i<p; i++)
    tr[2] += m2[i*(p+1)];  /* tr (XTWVX)^{-1}(XTW^2VX) */
  if (mi[MDEB]>2) printf("    vmat ok\n");

  switch(z)
  { case 1:
      unitvec(des->f1,0,p);
      vxtwx(&des->xtwx,des->f1,p);
      for (i=0; i<p; i++)
        for (j=0; j<p; j++)
        { tr[4] += m2[i*p+j]*m2[j*p+i];  /* tr(M^2) */
          tr[5] += des->f1[i]*V[i*p+j]*des->f1[j]; /* var(thetahat) */
        }
      tr[5] = sqrt(tr[5]);
      setzero(m2,p*p);
      for (i=0; i<des->n; i++)
      { stdlinks(link,lf,des->ind[i],des->th[i],robscale);
        ww = SQR(des->w[i])*des->w[i]*link[ZDDLL];
        addouter(m2,&des->X[i*p],&des->X[i*p],p,ww);
      }
      for (i=0; i<p; i++)
      { vxtwx(&des->xtwx,&m2[i*p],p);
        tr[3] += m2[i*(p+1)];
      }
      return;
    case 0:
      unitvec(des->f1,0,p);
      vxtwx(&des->xtwx,des->f1,p);
      for (i=0; i<=mi[MDIM]; i++) t0[i] = des->f1[i]; /* influnce and deriv */
      choldec(V,p);                                /* LL^T = X^T W^2 X */
      for (i=0; i<p; i++) vxtwx(&des->xtwx,&V[i*p],p); /* V = (X^TWX)^{-1}L */
      vxtwx(&des->xtwx,des->f1,p);
      for (i=0; i<p; i++)
      { for (j=0; j<p; j++)
        { m2[i*p+j] = 0;
          for (k=0; k<p; k++)
            m2[i*p+j] += V[k*p+i]*V[k*p+j]; /* ith column of covariance */
        }
      }
      if ((mi[MTG]==TDEN) && (mi[MLINK]==LIDENT))
        multmatscal(m2,1/SQR(lf->dp[DSWT]),p*p);
      return;
  }
}

INT robustinit(lf,des)
lfit *lf;
design *des;
{ INT i;
  for (i=0; i<des->n; i++)
    des->res[i] = resp(lf,des->ind[i])-base(lf,des->ind[i]);
  des->cf[0] = median(des->res,des->n);
  for (i=1; i<des->p; i++) des->cf[i] = 0.0;
  return(LF_OK);
}

INT circinit(lf,des)
lfit *lf;
design *des;
{ INT i, ii;
  double s0, s1;
  s0 = s1 = 0.0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    s0 += des->w[i]*prwt(lf,ii)*sin(resp(lf,ii)-base(lf,ii));
    s1 += des->w[i]*prwt(lf,ii)*cos(resp(lf,ii)-base(lf,ii));
  }
  des->cf[0] = atan2(s0,s1);
  for (i=1; i<des->p; i++) des->cf[i] = 0.0;
  return(LF_OK);
}

INT reginit(lf,des)
lfit *lf;
design *des;
{ INT i, ii;
  double sb, link[LLEN];
  s0 = s1 = sb = 0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    links(base(lf,ii),resp(lf,ii),lf->mi[MTG],LINIT,link,cens(lf,ii),prwt(lf,ii),1.0);
    s1 += des->w[i]*link[ZDLL];
    s0 += des->w[i]*prwt(lf,ii);
    sb += des->w[i]*prwt(lf,ii)*base(lf,ii);
  }
  if (s0==0) return(LF_NOPT); /* no observations with W>0 */
  setzero(des->cf,des->p);
  switch(lf->mi[MLINK])
  { case LIDENT:
      des->cf[0] = (s1-sb)/s0;
      return(LF_OK);
    case LLOG:
      if (s1<=0.0)
      { des->cf[0] = -1000;
        return(LF_INFA);
      }
      des->cf[0] = log(s1/s0) - sb/s0;
      return(LF_OK);
    case LLOGIT:
      if (s1<=0.0)
      { des->cf[0] = -1000;
        return(LF_INFA);
      }
      if (s1>=s0)
      { des->cf[0] = -1000;
        return(LF_INFA);
      }
      des->cf[0] = logit(s1/s0)-sb/s0;
      return(LF_OK);
    case LINVER:
      if (s1<=0.0)
      { des->cf[0] = 1000;
        return(LF_INFA);
      }
      des->cf[0] = s0/s1-sb/s0;
      return(LF_OK);
    case LSQRT:
      des->cf[0] = sqrt(s1/s0)-sb/s0;
      return(LF_OK);
    case LASIN:
      des->cf[0] = asin(sqrt(s1/s0))-sb/s0;
      return(LF_OK);
    default:
      ERROR(("reginit: invalid link %d",lf->mi[MLINK]));
      return(LF_ERR);
  }
}

INT lfinit(lf,des)
lfit *lf;
design *des;
{ double u[MXDIM], *X;
  INT i, ii, j, d, p, *mi;

  mi = lf->mi;
  p = des->p; d = mi[MDIM];
  X = des->X;
  des->xtwx.sm = 1+(mi[MDEG0]<mi[MDEG]);

  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    for (j=0; j<d; j++) u[j] = datum(lf,j,ii);
    fitfun(lf,u,des->xev,&X[i*p],NULL,0);
  }

  like = likereg;
  switch(mi[MTG]&63)
  { case TDEN:
    case TRAT:
    case THAZ:
      like = likeden;
      return(densinit(lf,des,des->h,des->cf,des->n));
    case TCAUC:
    case TROBT:
      return(robustinit(lf,des));
    case TCIRC:
      return(circinit(lf,des));
    default:
      return(reginit(lf,des));
  }
}

#define acceptpar(st,lk0,lk1) ((st==LF_OK) && (lk1>lk0-1.0e-3))

INT updatecoefs(des,lf,lk0,cf,dlt,nop,noit)
design *des;
lfit *lf;
double lk0, *cf, *dlt;
INT nop, noit;
{ INT i, st;
  double nc, nd, cut, f;

  f = 1.0;
  do
  { for (i=0; i<des->p; i++)
      des->cf[i] = cf[i]+f*dlt[i];
    st = like(lf,des,des->p);
    if (acceptpar(st,lk0,des->llk)) return(LF_OK);
    if (f==1)
    { nc = innerprod(cf,cf,des->p);
      nd = innerprod(dlt,dlt,des->p);
      cut = 0.0001*MIN(1.0,sqrt(nc/nd));
    }

    f = (st==LF_OK) ? f/10.0 : f/2.0;
  } while (f>cut);

  return(LF_FPROB);
}

INT lfiter(lf,des,noit)
lfit *lf;
design *des;
INT noit;
{ INT it, i, st, p, *mi, tg, oob, pf, cv, nop, rank;
  double lk0, lk1, *cf;
  
  mi = lf->mi; p = des->p;
  tg = mi[MTG];
  cf = des->cf;

  st = like(lf,des,p);
  if (lf_error) return(LF_ERR);

  for (it=0; it<mi[MMXIT]; it++)
  { for (i=0; i<p; i++)          /* store old coeffs and X^TW\dot{l} */
    { des->oc[i] = cf[i];
      des->f1[i] = des->res[i];
    }
    rank = vxtwx(&des->xtwx,des->f1,p);     /* (XTWVX)^-1 (XTW\dot{l}) */

    lk0 = des->llk; nop = p;
    if (rank==0) /* NR won't move! */
      des->f1[0] = -lk0/des->res[0];

    st = updatecoefs(des,lf,lk0,des->oc,des->f1,nop,noit);
    if (st==LF_QR) return(LF_OK);
    if (lf_error) return(LF_ERR);
    if (st!=LF_OK) return(st);
    
    oob = pf = cv = 0;
    lk1 = des->llk;
    switch (tg&63)  /* check for convergence */
    { case TGAUS:
        if ((mi[MLINK]==LIDENT)&((tg&128)==0))
          cv = (it>0); /* max of 1 iteration */
        else
          cv = (fabs(lk1-lk0)<1.0e-6*s0);
        break;
      case TRBIN:
      case TLOGT: oob = 0; /* ((mi[MLINK]==LLOGIT) && (fabs(cf[0])>100)); */
                  pf = lk1>-1.0e-5*s0;
                  cv = (fabs(lk1-lk0)<1.0e-6*s0);
                  break;
      case TPOIS:
      case TGEOM:
      case TWEIB:
      case TGAMM: oob = ((mi[MLINK]==LLOG) && (fabs(cf[0])>100));
printf("s0 %8.5f  lk1 %8.5f\n",s0,lk1);
                  pf = lk1>-1.0e-5*s0;
                  cv = (fabs(lk1-lk0)<1.0e-6*s0);
                  break;
      case TCIRC: cv = (fabs(lk1-lk0)<1.0e-6);
                  break;
      case TDEN:
      case TRAT:
      case THAZ:
        switch(mi[MLINK])
        { case LLOG:
            oob = fabs(cf[0])>100;
            cv = (fabs(lk1-lk0)<1.0e-6);
            break;
          case LIDENT:
            cv = 1;
            break;
          case LSQRT:
            cv = it==5;
            break;
        }
        break;
      case TROBT:
      case TCAUC:
        cv = (fabs(lk1-lk0)<1.0e-6);
        break;
      default: ERROR(("locfit: unknown target %d",tg));
    }
    cv &= nop==p;
    if (oob) for (i=1; i<p; i++) cf[i] = 0.0;
    if (cv | oob | pf) break;
  }
  if ((tg&63)==TDEN) /* convert from rate to density */
  { switch(mi[MLINK])
    { case LLOG:
        des->cf[0] -= log(lf->dp[DSWT]);
        break;
      case LIDENT:
        multmatscal(des->cf,1.0/lf->dp[DSWT],p);
        break;
      default: ERROR(("Density adjustment; invalid link"));
    }
  }
  des->llk = lk1;
  if (oob) return(LF_OOB); /* out of bounds */
  if (pf)  return(LF_PF); /* perfect fit */
  if (it==mi[MMXIT]) return(LF_NCON);
  return(LF_OK);
}

INT locfit(lf,des,h,noit)
lfit *lf;
design *des;
double h;
INT noit;
{ INT i, tg;
  double s0[3], sh[3], slo, shi;

  if (lf->mi[MDEB]>0)
  { printf("locfit: ");
    for (i=0; i<lf->mi[MDIM]; i++) printf(" %10.6f",des->xev[i]);
    printf("  h = %8.5f\n",h);
  }

  des->h = h;
  lf_status = lfinit(lf,des);
  if (lf_status>0) return(lf_status);

  for (i=0; i<3; i++) s0[i] = sh[i] = 0.0;
  slo = shi = 0.0;
  tg = lf->mi[MTG];
  robscale = lrobustscale(lf,des,tg);
  for (i=0; i<lf->mi[MMXIT]; i++)
  { lf_status = lfiter(lf,des,noit);
    if (lf_status>0) return(lf_status);

    s0[0] = s0[1]; sh[0] = sh[1];
    s0[1] = robscale;
    /* sh[1] = robustscale(lf,des,tg)-s0[1];
    if (fabs(sh[1])<1.0e-6) return(0);
    if (sh[1]>0) slo = s0[1];
    if (sh[1]<0) shi = s0[1];
    if (s0[0] == 0.0)
      robscale = s0[1] + sh[1];
    else
      robscale = s0[1] - sh[1]*(s0[1]-s0[0])/(sh[1]-sh[0]);
    if (robscale>shi) robscale = (slo+shi)/2;
    if (robscale<slo) robscale = (shi==0.0) ? 2*slo : (slo+shi)/2.0;
    if (robscale<shi/2) robscale = shi/2; */
    robscale = lrobustscale(lf,des,lf->mi[MTG]);
    sh[1] = robscale - s0[1];
    if (fabs(sh[1])<1.0e-6) return(0);
  }

  lf_status = LF_OK;
  return(lf_status);
}

void dercor(lf,des,h)
lfit *lf;
design *des;
double h;
{ double s1, dc[MXDIM], wd, link[LLEN];
  INT i, ii, j, m, p, d, *mi; 
  mi = lf->mi;
  if (mi[MTG]<=THAZ) return;
  d = mi[MDIM];
  p = des->p; m = des->n;
  if (mi[MDEB]>1) printf("  Corrcting derivatives\n");
  unitvec(des->f1,0,p);
  vxtwx(&des->xtwx,des->f1,p);
  setzero(dc,d);

  /* correction term is e1^T (XTWVX)^{-1} XTW' ldot. */
  for (i=0; i<m; i++)
  { s1 = innerprod(des->f1,&des->X[i*p],p);
    ii = des->ind[i];
    stdlinks(link,lf,ii,des->th[i],robscale);
    for (j=0; j<d; j++)
    { wd = des->w[i]*weightd(datum(lf,j,ii)-des->xev[j],lf->sca[j],d,mi[MKER],mi[MKT],h,lf->sty[j],des->di[ii]);
      dc[j] += s1*wd*link[ZDLL];
    }

  }
  for (j=0; j<d; j++) des->cf[j+1] += dc[j];
}
