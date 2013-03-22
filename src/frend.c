/*
 *   Copyright (c) 1998 Lucent Technologies.
 *   See README file for details.
 */

#include <math.h>
#include <stdio.h>
#include "local.h"

extern INT ident;
static double *f;
static INT llf;

INT calcp(deg,dim,kt)
INT deg, dim, kt;
{ INT i, k;
  if (kt==KPROD) return(dim*deg+1);
  k = 1;
  for (i=1; i<=deg; i++) k = k*(dim+i)/i;
  return(k);
}

double getxi(x,i)
double *x;
INT i;
{ if (x==NULL) return(0.0);
  return(x[i]);
}

double resp(lf,i)
struct tree *lf;
INT i;
{ if (lf->y==NULL) return(0.0);
  return(lf->y[i]);
}

extern INT cvi;

double prwt(lf,i)
struct tree *lf;
INT i;
{ if (i==cvi) return(0.0);
  if (lf->w==NULL) return(1.0);
  return(lf->w[i]);
}

double base(lf,i)
struct tree *lf;
INT i;
{ if (lf->base==NULL) return(0.0);
  return(lf->base[i]);
}

double cens(lf,i)
struct tree *lf;
INT i;
{ if (lf->c==NULL) return(0.0);
  return(lf->c[i]);
}

void fitfun(x,f,sca,d,deg,kt,der,nd,sty)
double *x, *f, *sca;
INT d, deg, kt, *der, nd, *sty;
{ INT m, i, j, k, cd[MXDIM];
  double ff[MXDIM][1+MXDEG];

  m = 0;
  f[m++] = (nd==0);
  if (deg==0) return;

  for (i=0; i<d; i++) cd[i] = 0;
  for (i=0; i<nd; i++) cd[der[i]]++;

  for (i=0; i<d; i++)
  { if (sty[i]==KANG)
    { switch(cd[i])
      { case 0:
          ff[i][0] = 1;
          ff[i][1] = sin(x[i]/sca[i])*sca[i];
          ff[i][2] = (1-cos(x[i]/sca[i]))*sca[i]*sca[i];
          if (deg>=3)
            WARN(("Can't handle angular model with deg>=3"));
          break;
        case 1:
          ff[i][0] = 0;
          ff[i][1] = cos(x[i]/sca[i]);
          ff[i][2] = sin(x[i]/sca[i])*sca[i];
          break;
        case 2:
          ff[i][0] = 0;
          ff[i][1] = -sin(x[i]/sca[i])/sca[i];
          ff[i][2] = cos(x[i]/sca[i]);
          break;
        default: WARN(("Can't handle angular model with >2 derivs"));
      }
    }
    else
    { for (j=0; j<cd[i]; j++) ff[i][j] = 0.0;
      ff[i][cd[i]] = 1.0;
      for (j=cd[i]+1; j<=deg; j++)
        ff[i][j] = ff[i][j-1]*x[i]/(j-cd[i]);
    }
  }

  for (i=0; i<d; i++) f[m++] = (cd[i]==nd) ? ff[i][1] : 0.0;
  if (deg==1) return;

  for (i=0; i<d; i++)
  { f[m++] = (cd[i]==nd) ? ff[i][2] : 0.0;
    if (kt!=KPROD)
      for (j=i+1; j<d; j++)
        f[m++] = (cd[i]+cd[j]==nd) ? ff[i][1]*ff[j][1] : 0.0;
  }
  if (deg==2) return;

  for (i=0; i<d; i++)
  { f[m++] = (cd[i]==nd) ? ff[i][3] : 0.0;
    if (kt!=KPROD)
    { for (k=i+1; k<d; k++)
        f[m++] = (cd[i]+cd[k]==nd) ? ff[i][2]*ff[k][1] : 0.0;
      for (j=i+1; j<d; j++)
      { f[m++] = (cd[i]+cd[j]==nd) ? ff[i][1]*ff[j][2] : 0.0;
        for (k=j+1; k<d; k++)
          f[m++] = (cd[i]+cd[j]+cd[k]==nd) ? ff[i][1]*ff[j][1]*ff[k][1] : 0.0;
  } } }
  if (deg==3) return;

  if (d==1)
  { for (i=4; i<=deg; i++) f[m++] = ff[0][i];
    return;
  }
  ERROR(("fitfun: can't handle deg=%d",deg))
}

double vocri(lk,t0,t2,pen)
double lk, t0, t2, pen;
{ if (pen==0) return(-2*t0*lk/((t0-t2)*(t0-t2)));
  return((-2*lk+pen*t2)/t0);
}

INT procvraw(des,lf,v)
struct design *des;
struct tree *lf;
INT v;
{ INT i, k;
  double *xev;
  xev = &lf->xev[v*lf->mi[MDIM]];
  if (lf->dp[DADP]>0)
    lf->h[v] = afit(des,lf,xev);
  else
    lf->h[v] = nbhd(xev,lf,des,lf->dp[DALP],lf->dp[DFXH],0);
  if (lf->h[v]<=0)
  { WARN(("zero bandwidth in procvraw"))
    return(1);
  }
  k = locfit(lf,des,xev,lf->h[v],lf->mi[MDEG0]<lf->mi[MDEG]);
  if (lf->mi[MDC]) dercor(lf,des,xev,lf->h[v]);
  for (i=0; i<lf->mi[MP]; i++) lf->coef[i*lf->nvm+v] = des->cf[i];
  lf->deg[v] = lf->mi[MDEG];
  return(k);
}

INT procv(des,lf,v)
struct design *des;
struct tree *lf;
INT v;
{ INT d, p, nvm, *mi, i, k;
  double trc[6], t0[1+MXDIM];

  k = procvraw(des,lf,v);

  mi = lf->mi; d = mi[MDIM]; p = mi[MP];
  nvm = lf->nvm;

  if (k==0)
    ldf(lf,des,trc,0,mi,t0);
  else
  { trc[0] = trc[2] = 0.0;
    for (i=0; i<p*p; i++) des->V[i] = 0.0;
  }
  lf->lik[v] = des->llk;
  lf->lik[nvm+v] = trc[2];
  lf->lik[2*nvm+v] = trc[0]-trc[2];
  lf->nlx[v] = sqrt(des->V[0]);
  lf->t0[v] = sqrt(t0[0]);
  if (p>1)
  { if (lf->nlx[v]==0)
      for (i=1; i<=d; i++) lf->nlx[i*nvm+v] = 0.0;
    else
      for (i=1; i<=d; i++) lf->nlx[i*nvm+v] = des->V[i]/lf->nlx[v];
    for (i=1; i<p; i++)
      lf->nlx[(i+d)*nvm+v] = sqrt(des->V[i*(p+1)]);
    if (lf->t0[v]==0)
      for (i=1; i<=d; i++) lf->t0[i*nvm+v] = 0.0;
    else
      for (i=1; i<=d; i++) lf->t0[i*nvm+v] = t0[i]/lf->t0[v];
  }
  return(k);
}

void expand(l,n,ind,m)
double *l;
INT n, m, *ind;
{ INT i, j;
  j = m-1;
  for (i=n-1; i>=0; i--)
  { l[i] = ((j>=0) && (ind[j]==i)) ? l[j--] : 0.0; }
}

INT procvhatm(des,lf,v)
struct design *des;
struct tree *lf;
INT v;
{ INT k, n;
  n = (ident==0) ? lf->mi[MN] : des->p;
  k = procvraw(des,lf,v);
  makelxd(lf,des,&lf->xev[v*lf->mi[MDIM]],&lf->L[v*n],(INT)0,lf->deriv,
    lf->nd,(INT)1);
  if (ident==0) expand(&lf->L[v*n],n,des->ind,des->n);
  return(k);
}

double intvo(des,lf,c0,c1,a,p,t0,t20,t21)
struct design *des;
struct tree *lf;
double *c0, *c1, a, t0, t20, t21;
INT p;
{ double th, lk, link[LLEN];
  INT i, ii;
  lk = 0;
  for (i=0; i<des->n; i++)
  { th = (1-a)*innerprod(c0,&des->X[i*p],p) + a*innerprod(c1,&des->X[i*p],p);
    ii = des->ind[i];
    links(th,resp(lf,ii),lf->mi[MTG],lf->mi[MLINK],link,cens(lf,ii),prwt(lf,ii));
    lk += des->w[i]*link[ZLIK];
  }
  des->llk = lk;
  return(vocri(des->llk,t0,(1-a)*t20+a*t21,lf->dp[DADP]));
}

INT procvvord(des,lf,v)
struct design *des;
struct tree *lf;
INT v;
{ double *xev, tr[6], gcv, g0, ap, coef[4][10], t2[4], th, md;
  INT i, j, k, d1, *mi, i0, p1, ip;
  mi = lf->mi;
  xev = &lf->xev[v*mi[MDIM]];

  lf->h[v] = nbhd(xev,lf,des,lf->dp[DALP],lf->dp[DFXH],0);
  if (lf->h[v]<=0) WARN(("zero bandwidth in procvvord"))

  ap = lf->dp[DADP];
  if ((ap==0) & ((mi[MTG]&63)!=TGAUS)) ap = 2.0;
  d1 = mi[MDEG]; p1 = mi[MP];
  for (i=0; i<p1; i++) coef[0][i] = coef[1][i] = coef[2][i] = coef[3][i] = 0.0;
  i0 = 0; g0 = 0;
  ip = 1;
  for (i=mi[MDEG0]; i<=d1; i++)
  { mi[MDEG] = i; des->p = mi[MP] = calcp(i,mi[MDIM],mi[MKT]);
    k = locfit(lf,des,xev,lf->h[v],0);

    ldf(lf,des,tr,1,lf->mi,NULL);
    gcv = vocri(des->llk,tr[0],tr[2],ap);
    if ((i==mi[MDEG0]) || (gcv<g0)) { i0 = i; g0 = gcv; md = i; }

    for (j=0; j<des->p; j++) coef[i][j] = des->cf[j];
    t2[i] = tr[2];

    if ((ip) && (i>mi[MDEG0]))
    { for (j=1; j<10; j++)
      { gcv = intvo(des,lf,coef[i-1],coef[i],j/10.0,des->p,tr[0],t2[i-1],t2[i]);
        if (gcv<g0) { g0 = gcv; md = i-1+j/10.0; }
      }
    }

  }

  if (i0<d1) /* recompute the best fit */
  { mi[MDEG] = i0; des->p = mi[MP] = calcp(i0,mi[MDIM],mi[MKT]);
 /*   k = locfit(lf,des,xev,lf->h[v],0);
    for (i=mi[MP]; i<p1; i++) des->cf[i] = 0.0; */
    i0 = (INT)md; if (i0==d1) i0--;
    th = md-i0;
    for (i=0; i<p1; i++) des->cf[i] = (1-th)*coef[i0][i]+th*coef[i0+1][i];
    mi[MDEG] = d1; mi[MP] = p1;
  }

  for (i=0; i<p1; i++) lf->coef[i*lf->nvm+v] = des->cf[i];
  lf->deg[v] = md;
  return(k);
}

void trace(tr,des,tc,m)
struct tree *tr;
struct design *des;
double *tc;
INT m; /* m=1: tr(L), tr(L'L);  m=2: + tr(L^2), tr(L'LL'), tr(L'L)^2 */
{ double *L, s, q;
  INT i, j, k, nl, n;
  if (ident==1) /* parametric fit - all = p assuming nonsing. */
  { tc[0] = tc[1] = tr->mi[MP];
    if (m==2) tc[2] = tc[3] = tc[4] = tr->mi[MP];
    return;
  }
  nl = tr->nl; s = 0;
  tc[0] = tc[1] = 0.0;
  if (m==2) tc[2] = tc[3] = tc[4] = 0.0;
  n = tr->mi[MN];
  if ((tr->mi[MEV]==EDATA) && (nl>0))
  { L = tr->L;
WARN(("traces won't work right if non-unit weights"))
    for (i=0; i<n; i++)
    { tc[0] += L[n*i*nl+i];
      for (j=i; j<n; j++)
      { tc[1] += SQR(L[n*i*nl+j]);
        if (j>i) tc[1] += SQR(L[n*j*nl+i]);
        if (m==2)
        { tc[2] += (2-(i==j))*L[n*i*nl+j]*L[n*j*nl+i];
          s = 0;
          for (k=0; k<n; k++)
            s += L[n*i*nl+k]*L[n*j*nl+k];
          tc[3] += s*L[n*i*nl+j];
          if (j>i) tc[3] += s*L[n*j*nl+i];
          tc[4] += (2-(i==j))*s*s;
    } } }
    return;
  }

  checkvl(&f,&llf,n);
  tc[0] = tc[1] = 0;
  intd(tr,des,f,PT0,NULL,0);
  if (lf_error) return;
  for (i=0; i<n; i++) tc[0] += tr->w[i]*f[i];
  intd(tr,des,f,PNLX,NULL,0);
  if (lf_error) return;
  for (i=0; i<n; i++) tc[1] += tr->w[i]*SQR(f[i]);
  /* Assume p sing values=1, then cubic decay to p+q */
  q = 35*(tc[1]-tr->mi[MP])/13; 
  if (m==2)
  { tc[2] = tc[1];
    tc[3] = tr->mi[MP] + 43*q/140;
    tc[4] = tr->mi[MP] + 191*q/715;
  }
  return;
}

void ressumm(lf,des)
struct tree *lf;
struct design *des;
{ INT i, j;
  double *dp, pw, r1, r2, th, t0, t1, u[MXDIM], link[LLEN];
  dp = lf->dp;
  dp[DLK] = dp[DT0] = dp[DT1] = 0;
  r1 = r2 = 0.0;
  for (i=0; i<lf->mi[MN]; i++)
  { if ((lf->mi[MEV]==EDATA)|(lf->mi[MEV]==ECROS))
    { th = base(lf,i)+lf->coef[i];
      t0 = lf->t0[i]*lf->t0[i];
      t1 = lf->nlx[i];
    } else
    { for (j=0; j<lf->mi[MDIM]; j++) u[j] = lf->x[j][i];
      th = base(lf,i)+intp(lf,des,u,PCOEF,NULL,0,1);
      t0 = intp(lf,des,u,PT0,NULL,0,1);
      t1 = intp(lf,des,u,PNLX,NULL,0,1);
    }
    pw = prwt(lf,i);
    links(th,resp(lf,i),lf->mi[MTG],lf->mi[MLINK],link,cens(lf,i),pw);
    t1 = t1*t1*link[ZDDLL];
    t0 = t0*link[ZDDLL];
    if (t1>1) t1 = 1;
    if (t0>1) t0 = 1; /* no observation gives >1 deg.free */
    dp[DLK] += link[ZLIK];
    dp[DT0] += t0;
    dp[DT1] += t1;
    if (pw>0)
    { r1 += link[ZDLL]*link[ZDLL]/pw;
      r2 += link[ZDDLL]/pw;
    }
    if (lf->mi[MGETH]==4) /* for gam */
    { lf->y[i] = resp(lf,i)-th;
      des->di[i]  = t1;
      des->w[i] = 1.0;
      des->ind[i] = i;
    }
  }

  if (lf->mi[MGETH]==4) /* orthog. residuals */
  { des->n = lf->mi[MN];
    lf->mi[MDEG] = 1;
    lf->mi[MP] = des->p = 1+lf->mi[MDIM];
    for (i=0; i<lf->mi[MDIM]; i++) u[i] = lf->x[i][0]; /* can do better!! */
    locfit(lf,des,u,0.0,1);
    for (i=0; i<lf->mi[MN]; i++) lf->y[i] -= des->th[i];
    return;
  }

  if ((lf->mi[MTG]&64)==64) /* quasi family */
    dp[DRV] = r1/r2 * lf->mi[MN] / (lf->mi[MN]-2*dp[DT0]+dp[DT1]);
  else
    dp[DRV] = 1.0;

  /* try to ensure consistency for family="circ"! */
  if (((lf->mi[MTG]&63)==TCIRC) & (lf->mi[MDIM]==1))
  { INT *ind, nv;
    double dlt, th0, th1;
    ind = des->ind;
    nv = lf->nv;
    for (i=0; i<nv; i++) ind[i] = i;
    lforder(ind,lf->xev,0,nv-1);
    for (i=1; i<nv; i++)
    { dlt = lf->xev[ind[i]]-lf->xev[ind[i]-1];
      th0 = lf->coef[ind[i]]-dlt*lf->coef[ind[i]+nv]-lf->coef[ind[i-1]];
      th1 = lf->coef[ind[i]]-dlt*lf->coef[ind[i]-1+nv]-lf->coef[ind[i-1]];
      if ((th0>PI)&(th1>PI))
      { for (j=0; j<i; j++)
          lf->coef[ind[j]] += 2*PI;
        i--;
      }
      if ((th0<(-PI))&(th1<(-PI)))
      { for (j=0; j<i; j++)
          lf->coef[ind[j]] -= 2*PI;
        i--;
      }
    }
  }
}

double rss(lf,des,df)
struct tree *lf;
struct design *des;
double *df;
{ double ss;
  INT i;
  ss = 0;
  if (ident==1)
  { for (i=0; i<lf->mi[MN]; i++)
      ss += SQR(resp(lf,i)-lf->coef[i]);
    *df = lf->mi[MN]-lf->mi[MP];
    return(ss);
  }
  ressumm(lf,des);
  *df = lf->mi[MN] - 2*lf->dp[DT0] + lf->dp[DT1];
  return(-2*lf->dp[DLK]);
}
