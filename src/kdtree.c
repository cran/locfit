/*
 *   Copyright (c) 1998 Lucent Technologies.
 *   See README file for details.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "local.h"
extern double atan2();

void ksmall(), newcell();
INT cvi;

void fitdefault(lf,mi,dp,n,d)
struct tree *lf;
INT *mi, n, d;
double *dp;
{ INT i;
  if (lf!=NULL)
  { if (dp==NULL)
    { if (lf->dp==NULL) lf->dp = (double *)calloc(LEND,sizeof(double));
      dp = lf->dp;
    }
    if (mi==NULL)
    { if (lf->mi==NULL) lf->mi = (INT *)calloc(LENM,sizeof(INT));
      mi = lf->mi;
    }
  }
  if (mi!=NULL)
  { mi[MTG] = TNUL;
    if (lf!=NULL)
    { mi[MTG] = (lf->y==NULL) ? TDEN : (64+TGAUS);
    }
    mi[MLINK] = LDEFAU;
    mi[MACRI] = ACP;
    mi[MDEG] = mi[MDEG0] = 2;    mi[MEV] = (ident==1) ? EDATA : ETREE;
    mi[MKT] = KSPH; mi[MKER] = WTCUB;
    mi[MIT] = IDEFA; mi[MDC] = mi[MREN] = 0;
    mi[MK] = 50; mi[MMINT] = 20;
    mi[MMXIT] = 20;
    mi[MN] = n; mi[MDIM] = d;
  }

  if (dp!=NULL)
  { dp[DALP] = 0.7; dp[DFXH] = dp[DADP] = 0.0;
    dp[DCUT] = 0.8;
  }

  if (lf==NULL) return;
  if (d<=0) ERROR(("must set MDIM before calling fitdefault"))
  for (i=0; i<d; i++)
  { if (lf->sca==NULL) lf->sca= (double *)calloc(MXDIM,sizeof(double));
    if (lf->xl ==NULL) lf->xl = (double *)calloc(MXDIM,sizeof(double));
    if (lf->fl ==NULL) lf->fl = (double *)calloc(MXDIM,sizeof(double));
    lf->sca[i] = 1.0;
    lf->xl[i] = lf->xl[i+d] = 0.0;
    lf->fl[i] = lf->fl[i+d] = 0.0;
  }
}

void checkvl(z,n0,n1)
double **z;
INT *n0, n1;
{ if (n1==0) WARN(("checkvl for zero space - really??"))
  if ((z==NULL) || (n1>*n0))
  { *n0 = n1;
    *z = (double *)calloc(n1,sizeof(double));
  }
}

void trchck(tr,nvm,ncm,d,p,vc) /* check space for nvm vertices, ncm cells */
struct tree *tr;
INT nvm, ncm, d, p, vc;
{ INT rw, *k;
  double *z;
  rw = nvm*(2*p+3*d+6)+ncm;
  checkvl(&tr->tw,&tr->ltw,rw);
  z = tr->tw;
  tr->xev = z; z += nvm*d;
  tr->coef= z; z += nvm*p;
  tr->nlx = z; z += nvm*(p+d);
  tr->t0  = z; z += nvm*(d+1);
  tr->lik = z; z += 3*nvm;
  tr->h   = z; z += nvm;
  tr->deg = z; z += nvm;
  tr->sv  = z; z += ncm;
  if (z != tr->tw+rw) WARN(("trchck: double assign problem"))

  rw = ncm*vc+3*MAX(ncm,nvm);
  if (tr->liw<rw)
  { tr->iw = (INT *)calloc(rw,sizeof(INT));
    tr->liw = rw;
  }
  k = tr->iw;
  tr->ce = k; k += vc*ncm;
  tr->s  = k; k += MAX(ncm,nvm);
  tr->lo = k; k += MAX(ncm,nvm);
  tr->hi = k; k += MAX(ncm,nvm);
  if (k != tr->iw+rw) WARN(("trchck: int assign problem"))

  tr->nvm = nvm; tr->ncm = ncm; tr->mi[MDIM] = d; tr->mi[MP] = p; tr->vc = vc;
}

void deschk(des,n,p)
struct design *des;
INT n, p;
{ INT rw;
  double *z;
  rw = n*(p+5)+4*p*p+6*p;
  checkvl(&des->dw,&des->lw,rw);
  z = des->dw;
  des->X = z; z += n*p;
  des->Z = z; z += p*p;
  des->w = z; z += n;
  des->res=z; z += n;
  des->di =z; z += n;
  des->th =z; z += n;
  des->wd =z; z += n;
  des->V  =z; z += p*p;
  des->P  =z; z += p*p;
  des->Q  =z; z += p*p;
  des->f1 =z; z += p;
  des->f2 =z; z += p;
  des->ss =z; z += p;
  des->oc =z; z += p;
  des->cf =z; z += p;
  des->dg =z; z += p;

  if (des->li<n)
  { des->ind = (INT *)calloc(n,sizeof(INT));
    des->li = n;
  }
  des->n = n; des->p = p;
}

void bbox(tr,bx)
struct tree *tr;
double *bx;
{ INT i, j, d, n;
  double z, mx, mn;
  d = tr->mi[MDIM]; n = tr->mi[MN];
  if (tr->mi[MKT]==KANG)
  { bx[0] = 0.0; bx[1] = 2*PI*tr->sca[0];
    return;
  }
  for (i=0; i<d; i++)
    if (bx[i]==bx[i+d])
    { if (tr->sty[i]==KANG)
      { bx[i] = 0.0; bx[i+d] = 2*PI*tr->sca[i];
      }
      else
      { mx = mn = tr->x[i][0];
        for (j=1; j<n; j++)
        { mx = MAX(mx,tr->x[i][j]);
          mn = MIN(mn,tr->x[i][j]);
        }
        if (tr->xl[i]<tr->xl[i+d]) /* user set xlim; maybe use them. */
        { z = mx-mn;
          if (mn-0.2*z < tr->xl[i]) mn = tr->xl[i];
          if (mx+0.2*z > tr->xl[i+d]) mx = tr->xl[i+d];
        }
        bx[i] = mn;
        bx[i+d] = mx;
      }
    }
}

INT defaultlink(link,targ)
INT link, targ;
{ if (link != LDEFAU) return(link);
  switch(targ&63)
  { case TDEN:
    case TRAT:
    case THAZ:
    case TGAMM:
    case TGEOM:
    case TPOIS: return(LLOG);
    case TCIRC:
    case TGAUS:
    case TROBT: return(LIDENT);
    case TLOGT: return(LLOGIT);
  }
  ERROR(("Unknown Family in defaultlink"));
  return(0);
}

void preproc(tr)
struct tree *tr;
{ INT d, i, j, n;
  double xb;
  d = tr->mi[MDIM]; n = tr->mi[MN];
  for (i=0; i<d; i++)
    if (tr->sca[i]<=0) /* set automatic scales */
    { if (tr->sty[i]==KANG) tr->sca[i] = 1.0;
      else
      { xb = tr->sca[i] = 0.0;
        for (j=0; j<n; j++) xb += tr->x[i][j];
        xb /= n;
        for (j=0; j<n; j++) tr->sca[i] += SQR(tr->x[i][j]-xb);
        tr->sca[i] = sqrt(tr->sca[i]/(n-1));
      }
    }
  bbox(tr,tr->fl);
  tr->mi[MLINK] = defaultlink(tr->mi[MLINK],tr->mi[MTG]);
}

void dataf(des,tr)
struct design *des;
struct tree *tr;
{ INT d, i, j, nv;
  nv = tr->mi[MN]; d = tr->mi[MDIM];
  trchck(tr,nv,0,d,des->p,0);
  for (i=0; i<nv; i++)
    for (j=0; j<d; j++) tr->xev[i*d+j] = tr->x[j][i];
  for (i=0; i<nv; i++)
  { des->vfun(des,tr,i);
    tr->s[i] = 0;
  }
  tr->nv = nv; tr->nce = 0;
}

void preset(des,tr)
struct design *des;
struct tree *tr;
{ INT i, nv;
  double *tmp;
  nv = tr->nvm;
  tmp = tr->xev;
  trchck(tr,nv,0,tr->mi[MDIM],des->p,0);
  tr->xev = tmp;
  for (i=0; i<nv; i++)
  { des->vfun(des,tr,i);
    tr->s[i] = 0;
  }
  tr->nv = nv; tr->nce = 0;
}

void crossf(des,tr)
struct design *des;
struct tree *tr;
{ INT d, i, j, n;
  n = tr->mi[MN]; d = tr->mi[MDIM];
  trchck(tr,n,0,d,des->p,0);
  for (i=0; i<n; i++)
    for (j=0; j<d; j++) tr->xev[i*d+j] = tr->x[j][i];
  for (cvi=0; cvi<n; cvi++)
  { tr->s[cvi] = 0;
    des->vfun(des,tr,cvi);
  }
  tr->nv = n; tr->nce = 0; tr->mi[MN] = n;
}

void gridf(des,tr)
struct design *des;
struct tree *tr;
{ INT d, i, j, nv, u, z;
  nv = 1; d = tr->mi[MDIM];
  for (i=0; i<d; i++)
  { if (tr->mg[i]==0)
      tr->mg[i] = 2+(INT)((tr->fl[i+d]-tr->fl[i])/(tr->sca[i]*tr->dp[DCUT]));
    nv *= tr->mg[i];
  }
  trchck(tr,nv,0,d,des->p,1<<d);
  for (i=0; i<nv; i++)
  { z = i;
    for (j=0; j<d; j++)
    { u = z%tr->mg[j];
      tr->xev[i*d+j] = (tr->mg[j]==1) ? tr->fl[j] :
                       tr->fl[j]+(tr->fl[j+d]-tr->fl[j])*u/(tr->mg[j]-1);
      z = z/tr->mg[j];
    }
    tr->s[i] = 0;
    des->vfun(des,tr,i);
  }
  tr->nv = nv; tr->nce = 0;
}

double sprat(h0,h1) /* return proportion of h0 */
double h0, h1;
{ double rat;
  rat = 1/(1+sqrt(h0/h1));
  if (rat<0.2) rat = 0.2;
  if (rat>0.8) rat = 0.8;
  return(rat);
}

INT inrr(x,ll,ur,d)
double *x, *ll, *ur;
INT d;
{ INT i;
  for (i=0; i<d; i++)
  { if ((x[i]<ll[i]) && (x[i]<ll[i]-0.0001*(ur[i]-ll[i]))) return(0);
    if ((x[i]>ur[i]) && (x[i]>ur[i]+0.0001*(ur[i]-ll[i]))) return(0);
  }
  return(1);
}

void ksmall(l, r, k, x, pi)
INT l, r, k, *pi;
double *x;
{   static INT i, j;
    static double t;
/*     find the $k$-th smallest of elements in x[pi[l]..pi[r]] */
/*     Floyd+Rivest, CACM Mar '75, Algorithm 489 */
/*     Input, x is data and unchanged. pi = integers (0:n-1),
       any permutation seems to be OK.
       At output, x[pi[k]] is (k+1) smallest.
*/ 
    while(l<r)
    { t = x[pi[k]];
      i = l; j = r;
      ISWAP(pi[l],pi[k]);
      if (t < x[pi[r]]) ISWAP(pi[l],pi[r]);
      while (i<j)
      { ISWAP(pi[i],pi[j]);
        ++i; --j;
        while(x[pi[i]]<t) ++i;
        while(t<x[pi[j]]) --j;
      }
      if (x[pi[l]] == t) { ISWAP(pi[l],pi[j]); }
               else { ++j; ISWAP(pi[r],pi[j]); }
      if (j <= k) l = j + 1;
      if (k <= j) r = j - 1;
    }
    return;
}

INT terminal(tr,p,pi,fc,d,m,sv)
struct tree *tr;
INT p, *pi, d, fc, *m;
double *sv;
{ INT i, k, l, u, im;
  double al, be, mx, t;
  mx = 0.0; im = 0;
  l = tr->lo[p]; u = tr->hi[p];
  if (u-l < fc) return(-1);
  if (d>1)
    for (k=0; k<d; k++)
    { al = be = tr->x[k][pi[l]];
      for (i=l+1; i<=u; i++)
      { t = tr->x[k][pi[i]];
        if (t<al) al = t;
        if (t>be) be = t;
      }
      be = (be-al)/tr->sca[k];
      if (be>mx)
      { mx = be;
        im = k;
    } }
  *m = (u+l)/2;
  ksmall(l,u,*m, tr->x[im], pi);
  *sv = tr->x[im][pi[*m]];
  while ((*m<u-1) && (tr->x[im][pi[*m+1]]==*sv)) (*m)++;
  if (*m==u) im = -1;
  return(im);
}

void kdtree(des,tr)
struct design *des;
struct tree *tr;
{ INT i, j, vc, d, nc, nv, ncm, nvm, k, m, n, p, *pi, fc;
  double sv;
  d = tr->mi[MDIM]; n = tr->mi[MN]; pi = des->ind;
  switch(tr->mi[MEV])
  { case EKDTR:
      fc = (INT) (tr->dp[DCUT]/4*n*MIN(tr->dp[DALP],1.0));
      k = 2*n/fc; vc = 1<<d;
      ncm = 2*k+1; nvm = (k+2)*vc/2;
      break;
    case EKDCE:
      fc = (INT) (n*tr->dp[DALP]); vc = 1;
      nvm = 1+(INT)(2*n/fc);
      ncm = 2*nvm+1;
      break;
    case ETREE:
      vc = 1<<d; nvm = (tr->mi[MK]+2)*vc/2;
      ncm = 2*tr->mi[MK]+1;
      break;
    default: ERROR(("kdtree: invalid ev"))
  }
  trchck(tr,nvm,ncm,d,des->p,vc);
  nv = 0;
  if (tr->mi[MEV] != EKDCE)
  { for (i=0; i<vc; i++)
    { j = i;
      for (k=0; k<d; ++k)
      { tr->xev[d*i+k] = tr->fl[d*(j%2)+k];
        j >>= 1;
      }
    }
    nv = vc;
    for (j=0; j<vc; j++) tr->ce[j] = j;
  }
  if (tr->mi[MEV]==ETREE)
  { for (j=0; j<vc; j++)
    { des->vfun(des,tr,j);
      tr->s[j] = 0;
    }
    tr->nv = nv;
    growquad(des,tr,tr->ce,nvm,NULL,NULL,tr->fl,&tr->fl[d]);
    tr->nce = 1;
    return;
  }

  for (i=0; i<n; i++) pi[i] = i;
  p = 0; nc = 1;
  tr->lo[p] = 0; tr->hi[p] = n-1;
  tr->s[p] = -1;
  while (p<nc)
  { k = terminal(tr,p,pi,fc,d,&m,&sv);
    if (k>=0)
    { if ((ncm<nc+2) | (2*nvm<2*nv+vc))
      { WARN(("Insufficient space for full tree"))
        tr->nce = nc; tr->nv = nv;
        return;
      }
      tr->s[p] = k;
      tr->sv[p] = sv;
      tr->lo[nc] = tr->lo[p];
      tr->hi[nc] = m;
      tr->lo[nc+1] = m+1;
      tr->hi[nc+1] = tr->hi[p];
      tr->s[nc] = tr->s[nc+1] = -1;
      tr->lo[p] = nc; tr->hi[p] = nc+1;
      nc=nc+2; i = nv;
      if (tr->mi[MEV] != EKDCE)
        newcell(&nv,vc,tr->xev, d, k, sv,
             &tr->ce[p*vc], &tr->ce[(nc-2)*vc], &tr->ce[(nc-1)*vc]);
    }
    else if (tr->mi[MEV]==EKDCE) /* new vertex at cell center */
    { sv = 0;
      for (i=0; i<d; i++) tr->xev[d*nv+i] = 0;
      for (j=tr->lo[p]; j<=tr->hi[p]; j++)
      { sv += tr->w[pi[j]];
        for (i=0; i<d; i++)
          tr->xev[d*nv+i] += tr->x[i][pi[j]]*tr->w[pi[j]];
      }
      for (i=0; i<d; i++) tr->xev[d*nv+i] /= sv;
      tr->mi[MN] = tr->hi[p]-tr->lo[p]+1;
      des->ind = &pi[tr->lo[p]];
      des->vfun(des,tr,nv);
      tr->mi[MN] = n; des->ind = pi;
      nv++;
    }
    p++;
  }
  if (tr->mi[MEV]==EKDTR) for (i=0; i<nv; i++) des->vfun(des,tr,i);
  tr->nce = nc; tr->nv = nv;
  return;
}

void newcell(nv,vc,xev, d, k, t, cpar, clef, crig)
double *xev, t;
INT *nv, vc, d, k, *cpar, *clef, *crig;
{ INT i, ii, j, j2, tk, match;
  tk = 1<<k;
  for (i=0; i<vc; i++)
    if (xev[d*cpar[i]+k] < t)
    { for (j=0; j<d; j++) xev[*nv*d+j] = xev[d*cpar[i]+j];
      xev[*nv*d+k] = t;
      match = 0; j = vc; /* no matches in first vc points */
      while ((j<*nv) && (!match))
      { j2 = 0;
        while ((j2<d) && (xev[*nv*d+j2] == xev[j*d+j2])) j2++;
        match = (j2==d);
        if (!match) j++;
      }
      ii = i+tk;
      clef[i] = cpar[i];
      clef[ii]= crig[i] = j;
      crig[ii]= cpar[ii];
      if (!match) (*nv)++;
    }
  return;
}

INT newsplit(des,tr,i0,i1,nvm,pv)
struct design *des;
struct tree *tr;
INT i0, i1, nvm, pv;
{ double z[MXDIM];
  INT d, i, nv;
  d = tr->mi[MDIM];
  if (i0>i1) ISWAP(i0,i1);
  for (i=0; i<d; i++) z[i] = (tr->xev[i0*d+i]+tr->xev[i1*d+i])/2;
  for (i=i1+1; i<tr->nv; i++)
    if ((tr->lo[i]==i0) && (tr->hi[i]==i1)) return(i); /* duplicate */
  if (tr->nv==nvm) ERROR(("newsplit: out of vertex space"));
  if (lf_error) return(-1);
  nv = tr->nv;
  tr->lo[nv] = MIN(i0,i1); tr->hi[nv] = MAX(i0,i1);
  for (i=0; i<d; i++) tr->xev[nv*d+i] = (tr->xev[i0*d+i]+tr->xev[i1*d+i])/2;
  if (pv) /* pseudo vertex */
  { tr->h[nv] = (tr->h[i0]+tr->h[i1])/2;
    tr->s[nv] = 1; /* pseudo-vertex */
  }
  else /* real vertex */
  { des->vfun(des,tr,nv);
    tr->s[nv] = 0;
  }
  tr->nv++;
  return(nv);
}

INT needtosplitq(tr,ce,le,ll,ur)
struct tree *tr;
INT *ce;
double *le, *ll, *ur;
{ INT d, vc, i, is;
  double h;
  d = tr->mi[MDIM]; vc = 1<<d;
  is = 0;
  for (i=0; i<d; i++)
  { le[i] = (ur[i]-ll[i])/tr->sca[i];
    if (le[i]>le[is]) is = i;
  }
  for (i=0; i<vc; i++)
  { h = tr->h[ce[i]];
    if ((h>0) && (h*tr->dp[DCUT]<le[is])) return(is);
  }
  return(-1);
}

void growquad(des,tr,ce,nvm,ct,term,ll,ur)
struct design *des;
struct tree *tr;
INT nvm, *ce, *ct, *term;
double *ll, *ur;
{ INT i, i0, i1, d, vc, ns, tk, nce[1024], pv;
  double le[MXDIM], z;
  d = tr->mi[MDIM]; vc = 1<<d;
  ns = needtosplitq(tr,ce,le,ll,ur);
  if (ns==-1)
  { if (ct != NULL)
    { for (i=0; i<vc; i++) term[*ct*vc+i] = ce[i];
      (*ct)++;
    }
    return;
  }
  tk = 1<<ns;
  for (i=0; i<vc; i++)
  { if ((i&tk)==0) nce[i] = ce[i];
    else
    { i0 = ce[i];
      i1 = ce[i-tk];
      pv = le[ns] < (tr->dp[DCUT]*MIN(tr->h[i0],tr->h[i1]));
      nce[i] = newsplit(des,tr,i0,i1,nvm,pv);
      if (lf_error) return;
    }
  }
  z = ur[ns]; ur[ns] = (z+ll[ns])/2;
  growquad(des,tr,nce,nvm,ct,term,ll,ur);
  if (lf_error) return;
  ur[ns] = z;
  for (i=0; i<vc; i++)
    nce[i] = ((i&tk)== 0) ? nce[i+tk] : ce[i];
  z = ll[ns]; ll[ns] = (z+ur[ns])/2;
  growquad(des,tr,nce,nvm,ct,term,ll,ur);
  ll[ns] = z;
}

INT needtosplit(tr,ce,le)
struct tree *tr;
double *le;
INT *ce;
{ INT d, i, j, nts, vc;
  double di;
  nts = 0; d = tr->mi[MDIM]; vc = d+1;
  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
    { di = rho(&tr->xev[ce[i]*d],&tr->xev[ce[j]*d],tr->sca,d,KSPH,NULL);
      le[i*vc+j] = le[j*vc+i] = di/MIN(tr->h[ce[i]],tr->h[ce[j]]);
      nts = nts || le[i*vc+j]>tr->dp[DCUT];
    }
  return(nts);
}

void resort(pv,xev,dig)
double *xev;
INT *pv, *dig;
{ double d0, d1, d2;
  INT i;
  d0 = d1 = d2 = 0;
  for (i=0; i<3; i++)
  { d0 += (xev[3*pv[11]+i]-xev[3*pv[1]+i])*(xev[3*pv[11]+i]-xev[3*pv[1]+i]);
    d1 += (xev[3*pv[ 7]+i]-xev[3*pv[2]+i])*(xev[3*pv[ 7]+i]-xev[3*pv[2]+i]);
    d2 += (xev[3*pv[ 6]+i]-xev[3*pv[3]+i])*(xev[3*pv[ 6]+i]-xev[3*pv[3]+i]);
  }
  if ((d0<=d1) & (d0<=d2))
  { dig[0] = pv[1]; dig[1] = pv[11];
    dig[2] = pv[2]; dig[3] = pv[7];
    dig[4] = pv[3]; dig[5] = pv[6];
  }
  else if (d1<=d2)
  { dig[0] = pv[2]; dig[1] = pv[7];
    dig[2] = pv[1]; dig[3] = pv[11];
    dig[4] = pv[3]; dig[5] = pv[6];
  }
  else
  { dig[0] = pv[3]; dig[1] = pv[6];
    dig[2] = pv[2]; dig[3] = pv[7];
    dig[4] = pv[1]; dig[5] = pv[11];
  }
}

void growtri(des,tr,ce,nvm,ct,term)
struct design *des;
struct tree *tr;
INT nvm, *ce, *ct, *term;
{ double le[(1+MXDIM)*(1+MXDIM)], ml;
  INT pv[(1+MXDIM)*(1+MXDIM)], nce[1+MXDIM], d, i, j, im, jm, vc, dig[6];
  if (lf_error) return;
  d = tr->mi[MDIM]; vc = d+1;
  if (!needtosplit(tr,ce,le))
  { if (ct != NULL)
    { for (i=0; i<vc; i++) term[*ct*vc+i] = ce[i];
      (*ct)++;
    }
    return;
  }
  if (d>3)
  { ml = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<vc; j++)
        if (le[i*vc+j]>ml) { ml = le[i*vc+j]; im = i; jm = j; }
    pv[0] = newsplit(des,tr,ce[im],ce[jm],nvm,0);
    for (i=0; i<vc; i++) nce[i] = ce[i];
    nce[im] = pv[0]; growtri(des,tr,nce,nvm,ct,term); nce[im] = ce[im];
    nce[jm] = pv[0]; growtri(des,tr,nce,nvm,ct,term);
    return;
  }

  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
      pv[i*vc+j] = pv[j*vc+i]
        = newsplit(des,tr,ce[i],ce[j],nvm,le[i*vc+j]<=tr->dp[DCUT]);
  for (i=0; i<=d; i++) /* corners */
  { for (j=0; j<=d; j++) nce[j] = (j==i) ? ce[i] : pv[i*vc+j];
    growtri(des,tr,nce,nvm,ct,term);
  }
  
  if (d==2) /* center for d=2 */
  { nce[0] = pv[5]; nce[1] = pv[2]; nce[2] = pv[1];
    growtri(des,tr,nce,nvm,ct,term);
  }
  if (d==3) /* center for d=3 */
  { resort(pv,tr->xev,dig);
    nce[0] = dig[0]; nce[1] = dig[1];
    nce[2] = dig[2]; nce[3] = dig[4]; growtri(des,tr,nce,nvm,ct,term);
    nce[2] = dig[5]; nce[3] = dig[3]; growtri(des,tr,nce,nvm,ct,term);
    nce[2] = dig[2]; nce[3] = dig[5]; growtri(des,tr,nce,nvm,ct,term);
    nce[2] = dig[4]; nce[3] = dig[3]; growtri(des,tr,nce,nvm,ct,term);
  }
  if (d==1) return;
}

void descend(tr,xa,ce)
struct tree *tr;
double *xa;
INT *ce;
{ double le[(1+MXDIM)*(1+MXDIM)], ml;
  INT d, vc, i, j, pv[(1+MXDIM)*(1+MXDIM)], im, jm;
  struct design *des;
  des = NULL;
  if (!needtosplit(tr,ce,le)) return;
  d = tr->mi[MDIM]; vc = d+1;

  if (d>3) /* split longest edge */
  { ml = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<vc; j++)
        if (le[i*vc+j]>ml) { ml = le[i*vc+j]; im = i; jm = j; }
    pv[0] = newsplit(des,tr,ce[im],ce[jm],0,0);
    if (xa[im]>xa[jm])
    { xa[im] -= xa[jm]; xa[jm] *= 2; ce[jm] = pv[0]; }
    else
    { xa[jm] -= xa[im]; xa[im] *= 2; ce[im] = pv[0]; }
    descend(tr,xa,ce);
    return;
  }

  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
      pv[i*vc+j] = pv[j*vc+i]
        = newsplit(des,tr,ce[i],ce[j],0,le[i*d+j]<=tr->dp[DCUT]);
  for (i=0; i<=d; i++) if (xa[i]>=0.5) /* in corner */
  { for (j=0; j<=d; j++)
    { if (i!=j) ce[j] = pv[i*vc+j];
      xa[j] = 2*xa[j];
    }
    xa[i] -= 1;
    descend(tr,xa,ce);
    return;
  }
  if (d==1) { ERROR(("weights sum to < 1")); }
  if (d==2) /* center */
  { ce[0] = pv[5]; xa[0] = 1-2*xa[0];
    ce[1] = pv[2]; xa[1] = 1-2*xa[1];
    ce[2] = pv[1]; xa[2] = 1-2*xa[2];
    descend(tr,xa,ce);
  }
  if (d==3) /* center */
  { double z; INT dig[6];
    resort(pv,tr->xev,dig);
    ce[0] = dig[0]; ce[1] = dig[1];
    xa[0] *= 2; xa[1] *= 2; xa[2] *= 2; xa[3] *= 2;
    if (xa[0]+xa[2]>=1)
    { if (xa[0]+xa[3]>=1)
      { ce[2] = dig[2]; ce[3] = dig[4];
        z = xa[0];
        xa[3] += z-1; xa[2] += z-1; xa[0] = xa[1]; xa[1] = 1-z;
      }
      else
      { ce[2] = dig[2]; ce[3] = dig[5];
        z = xa[3]; xa[3] = xa[1]+xa[2]-1; xa[1] = z;
        z = xa[2]; xa[2] += xa[0]-1; xa[0] = 1-z;
    } }
    else
    { if (xa[1]+xa[2]>=1)
      { ce[2] = dig[5]; ce[3] = dig[3];
        xa[1] = 1-xa[1]; xa[2] -= xa[1]; xa[3] -= xa[1];
      }
      else
      { ce[2] = dig[4]; ce[3] = dig[3];
        z = xa[3]; xa[3] += xa[1]-1; xa[1] = xa[2];
        xa[2] = z+xa[0]-1; xa[0] = 1-z;
    } }
    descend(tr,xa,ce);
} }

void meanofdata(lf,z)
struct tree *lf;
double *z;
{ INT d, i, j;
  double s;
  s = 0; d = lf->mi[MDIM];
  for (i=0; i<d; i++) z[i] = 0;
  for (i=0; i<lf->mi[MN]; i++)
  { s += prwt(lf,i);
    for (j=0; j<d; j++) z[j] += prwt(lf,i)*lf->x[j][i];
  }
  for (i=0; i<d; i++) z[i] /= s;
}

void covrofdata(lf,V,mn) /* covar of data; mean in mn */
struct tree *lf;
double *V, *mn;
{ INT d, i, j, k;
  double s;
  s = 0; d = lf->mi[MDIM];
  for (i=0; i<d*d; i++) V[i] = 0;
  for (i=0; i<lf->mi[MN]; i++)
  { s += prwt(lf,i);
    for (j=0; j<d; j++)
      for (k=0; k<d; k++)
        V[j*d+k] += prwt(lf,i)*(lf->x[j][i]-mn[j])*(lf->x[k][i]-mn[k]);
  }
  for (i=0; i<d*d; i++) V[i] /= s;
}

INT intri(x,w,xev,xa,d) /* is x in triangle bounded by xd[0..d-1]? */
double *x, *xev, *xa;
INT d, *w;
{ INT i, j;
  double eps, *r, xd[MXDIM*MXDIM];
  eps = 1.0e-10;
  r = &xev[w[d]*d];
  for (i=0; i<d; i++)
  { xa[i] = x[i]-r[i];
    for (j=0; j<d; j++) xd[i*d+j] = xev[w[i]*d+j]-r[j];
  }
  solve(xd,xa,d);
  xa[d] = 1.0;
  for (i=0; i<d; i++) xa[d] -= xa[i];
  for (i=0; i<=d; i++) if ((xa[i]<-eps) | (xa[i]>1+eps)) return(0);
  return(1);
}

void phull(des,tr) /* Triangulation with polyhedral start */
struct design *des;
struct tree *tr;
{ INT i, j, k, n, d, nc, nvm, ncm, vc, *ce, ed[1+MXDIM];
  double V[MXDIM*MXDIM], P[MXDIM*MXDIM], sigma, z[MXDIM], xa[1+MXDIM];
  d = tr->mi[MDIM]; n = tr->mi[MN]; tr->nv = nc = 0;
  vc = d+1; nvm = tr->mi[MK]*d; ncm = nvm;
  trchck(tr,nvm,ncm,d,des->p,vc);
  ce = tr->ce;
  meanofdata(tr,des->xb); tr->nv++; /* construct centerpoint */
  covrofdata(tr,V,des->xb); /* fix this with scaling */
  eigen(V,P,d,tr->mi[MMXIT]);

  for (i=0; i<d; i++) /* add vertices +- 2sigma*eigenvect */
  { sigma = sqrt(V[i*(d+1)]);
    for (j=0; j<d; j++)
      tr->xev[tr->nv*d+j] = tr->xev[j]-2*sigma*P[j*d+i];
    tr->nv++;
    for (j=0; j<d; j++)
      tr->xev[tr->nv*d+j] = tr->xev[j]+2*sigma*P[j*d+i];
    tr->nv++;
  }

  for (i=0; i<n; i++) /* is point i inside? */
  { ed[0] = 0;
    for (j=0; j<d; j++)
    { z[j] = 0;
      for (k=0; k<d; k++) z[j] += P[k*d+j]*(tr->x[k][i]-tr->xev[k]);
      ed[j+1] = 2*j+1+(z[j]>0);
      for (k=0; k<d; k++) z[j] = tr->x[j][i];
    }
    k = intri(z,ed,tr->xev,xa,d);
    if (xa[0]<0)
    { for (j=1; j<=d; j++)
        for (k=0; k<d; k++)
          tr->xev[ed[j]*d+k] = xa[0]*tr->xev[k]+(1-xa[0])*tr->xev[ed[j]*d+k];
    }
  }

  nc = 1<<d; /* create initial cells */
  for (i=0; i<nc; i++)
  { ce[i*vc] = 0; k = i;
    for (j=0; j<d; j++)
    { ce[i*vc+j+1] = 2*j+(k%2)+1;
      k>>=1;
    }
  }

  for (i=0; i<tr->nv; i++)
  { des->vfun(des,tr,i);
    tr->s[i] = 0;
  }
  for (i=0; i<nc; i++)
  { for (j=0; j<=d; j++) ed[j] = tr->ce[i*vc+j];
    growtri(des,tr,&tr->ce[i*vc],nvm,(INT *)NULL,(INT *)NULL);
  }
  tr->nce = nc;
}

void evaluator(des,lf,vfun)
struct design *des;
struct tree *lf;
INT (*vfun)();
{ INT i, *mi;
  des->vfun = vfun;
  mi = lf->mi;
  mi[MP] = calcp(mi[MDEG],mi[MDIM],mi[MKT]);
  lf->nnl = mi[MDIM]+mi[MP];
  des->na = 0; des->pref = 0;
  cvi = -1; /* inhibit cross validation */
  preproc(lf);
  deschk(des,mi[MN],mi[MP]);
  if (mi[MTG]==TDEN)
  { des->na = 0.0;
    for (i=0; i<mi[MN]; i++)
      des->na += prwt(lf,i);
  }
  lf->ord = 0;
  if ((mi[MDIM]==1) && (lf->sty[0]!=KANG) && (mi[MKER]!=WGAUS))
  { i = 1;
    while ((i<mi[MN]) && (lf->x[0][i]>=lf->x[0][i-1])) i++;
    lf->ord = (i==mi[MN]);
  }

  switch(mi[MEV])
  { case EPHULL: phull(des,lf); break;
    case EDATA:  dataf(des,lf); break;
    case ECROS:  crossf(des,lf); break;
    case EGRID:  gridf(des,lf); break;
    case EKDCE:  mi[MKT] = KCE;
    case ETREE: 
    case EKDTR:  kdtree(des,lf); break;
    case EPRES:  preset(des,lf); break;
    default: ERROR(("evaluator: Invalid evaluation structure"));
  }

  /* renormalize for family=density */
  if ((mi[MREN]) && (mi[MTG]==TDEN)) densrenorm(lf,des);
}

void hermite1(x,z,phi)
double x, z, *phi;
{ if (z==0)
  { phi[0] = 1.0; phi[1];
    return;
  }
  phi[1] = x/z;
  phi[0] = 1-phi[1];
  return;
}

void hermite2(x,z,phi)
double x, z, *phi;
{ double h;
  if (z==0)
  { phi[0] = 1.0; phi[1] = phi[2] = phi[3] = 0.0;
    return;
  }
  h = x/z;
  if (h<0)
  { phi[0] = 1; phi[1] = 0;
    phi[2] = h; phi[3] = 0;
    return;
  }
  if (h>1)
  { phi[0] = 0; phi[1] = 1;
    phi[2] = 0; phi[3] = h-1;
    return;
  }
  phi[1] = h*h*(3-2*h);
  phi[0] = 1-phi[1];
  phi[2] = h*(1-h)*(1-h);
  phi[3] = h*h*(h - 1);
}

double cubint(h,f0,f1,d0,d1)
double h, f0, f1, d0, d1;
{ double phi[4];
  hermite2(h,1.0,phi);
  return(phi[0]*f0+phi[1]*f1+phi[2]*d0+phi[3]*d1);
}

double cubintd(h,f0,f1,d0,d1)
double h, f0, f1, d0, d1;
{ double phi[4];
  phi[1] = 6*h*(1-h);
  phi[0] = -phi[1];
  phi[2] = (1-h)*(1-3*h);
  phi[3] = h*(3*h-2);
  return(phi[0]*f0+phi[1]*f1+phi[2]*d0+phi[3]*d1);
}

double cubeint(h,v0,v1,z0,z1,dim)
double h, *v0, *v1, *z0, *z1;
INT dim;
{ double g0, g1;
  INT i;
  g0 = g1 = 0;
  for (i=0; i<dim; i++)
  { g0 += (v1[i]-v0[i])*z0[i+1];
    g1 += (v1[i]-v0[i])*z1[i+1];
  }
  return(cubint(h,z0[0],z1[0],g0,g1));
}

double interpcuk(xev,coef,ce,x,d,vc,ncc)
double *xev, **coef, *x;
INT *ce, d, vc, ncc;
{ double phi[4], g[2304], *ll, *ur;
  INT i, j, k, nc;
  nc = (ncc==1) ? 1 : d+1;
  for (i=0; i<nc; i++)
    for (j=0; j<vc; j++)
      g[j*nc+i] = coef[i][ce[j]];
  ll = &xev[d*ce[0]]; ur = &xev[d*ce[vc-1]];
  for (i=d-1; i>=0; i--)
  { vc >>= 1;
    if (nc==1)
      hermite1(x[i]-ll[i],ur[i]-ll[i],phi);
    else
      hermite2(x[i]-ll[i],ur[i]-ll[i],phi);
    for (j=0; j<vc; j++)
    { g[j*nc] = phi[0]*g[j*nc] + phi[1]*g[(j+vc)*nc];
      if (nc>1)
      { g[j*nc] += (phi[2]*g[j*nc+i+1]+phi[3]*g[(j+vc)*nc+i+1])*(ur[i]-ll[i]);
        for (k=1; k<=i; k++)
          g[j*nc+i] = phi[0]*g[j*nc+i] + phi[1]*g[(j+vc)*nc+i];
    } }
  }
  return(g[0]);
}

double blend(tr,s,x,z,ll,ur,j,nt,t,p,vc)
struct tree *tr;
double s, *x, *ll, *ur, **z;
INT j, nt, *t, p, vc;
{ INT i, k, k1, m, j0, j1, *ce;
  double v0, v1, xibar, g0[3], g1[3], gg[4], gp[4], h, phi[4];
  ce = tr->ce;
  for (k=0; k<4; k++)  /* North South East West */
  { k1 = (k>1);
    v0 = ll[k1]; v1 = ur[k1];
    j0 = ce[j+2*(k==0)+(k==2)];
    j1 = ce[j+3-2*(k==1)-(k==3)];
    xibar = (k%2==0) ? ur[k<2] : ll[k<2];
    m = nt;
    while ((m>=0) && ((tr->s[t[m]] != (k<=1)) | (tr->sv[t[m]] != xibar))) m--;
    if (m >= 0)
    { m = (k%2==1) ? tr->lo[t[m]] : tr->hi[t[m]];
      while (tr->s[m] != -1)
        m = (x[tr->s[m]] < tr->sv[m]) ? tr->lo[m] : tr->hi[m];
      if (v0 < tr->xev[2*ce[m*vc+2*(k==1)+(k==3)]+k1])
      { j0 = ce[m*vc+2*(k==1)+(k==3)];
        v0 = tr->xev[2*j0+k1];
      }
      if (tr->xev[2*ce[m*vc+3-2*(k==0)-(k==2)]+k1] < v1)
      { j1 = ce[m*vc+3-2*(k==0)-(k==2)];
        v1 = tr->xev[2*j1+k1];
      }
    }
    for (i=0; i<3-2*(p==1); i++)
    { g0[i] = z[i][j0];
      g1[i] = z[i][j1];
    }
    if (p==1)
    { hermite1(x[(k>1)]-v0,v1-v0,phi);
      gg[k] = phi[0]*g0[0] + phi[1]*g1[0];
    }
    else
    { hermite2(x[(k>1)]-v0,v1-v0,phi);
      gg[k] = phi[0]*g0[0]+phi[1]*g1[0]+(phi[2]*g0[1+k1]+phi[3]*g1[1+k1])*(v1-v0);
      gp[k] = phi[0]*g0[2-k1] + phi[1]*g1[2-k1];
    }
  }
  s = -s;
  if (p==1)
    for (k=0; k<2; k++)
    { hermite1(x[k]-ll[k],ur[k]-ll[k]);
      s += phi[0]*gg[3-2*k] + phi[1]*gg[2-2*k];
    }
    else
    for (k=0; k<2; k++) /* EW NS */
    { hermite2(x[k]-ll[k],ur[k]-ll[k],phi);
      s += phi[0]*gg[3-2*k] + phi[1]*gg[2-2*k]
          +(phi[2]*gp[3-2*k] + phi[3]*gp[2-2*k]) * (ur[k]-ll[k]);
    }
  return(s);
}

double kdint(tr,x,coef,nc)
struct tree *tr;
double *x, **coef;
INT nc;
{ INT d, vc, k, t[20], nt, *ce;
  double *ll, *ur, ff;
  d = tr->mi[MDIM];
  vc = tr->vc;
  if (d > 8) ERROR(("d too large in kdint"));
  nt = 0; t[nt] = 0; k = 0;
  while (tr->s[k] != -1)
  { nt++;
    if (nt>=20) { ERROR(("Too many levels in kdint")); return(NOSLN); }
    k = t[nt] = (x[tr->s[k]] < tr->sv[k]) ? tr->lo[k] : tr->hi[k];
  }
  ce = &tr->ce[k*vc];
  ll = &tr->xev[ce[0]*d];
  ur = &tr->xev[ce[vc-1]*d];
  ff = interpcuk(tr->xev,coef,ce,x,d,vc,nc);
  if (d==2) ff = blend(tr,ff,x,coef,ll,ur,k*vc,nt,t,nc,vc);
  return(ff);
}

double interptr(v,coef,w,d,nc,xxa)
double *v, **coef, *xxa;
INT d, *w, nc;
{ double sa, lb, cg[(MXDIM+1)*(MXDIM+1)];
  INT vc, i, j, k, ww[1+MXDIM];
  vc = d+1;
  if (nc==1) /* linear interpolate */
  { sa = 0;
    for (i=0; i<vc; i++) sa += xxa[i]*coef[0][w[i]];
    return(sa);
  }
  for (i=0; i<=d; i++) ww[i] = w[i];
  do
  { k=0;
    for (i=0; i<d; i++)
      if (ww[i]>ww[i+1])
      { j=ww[i]; ww[i]=ww[i+1]; ww[i+1]=j; k=1;
        lb = xxa[i]; xxa[i] = xxa[i+1]; xxa[i+1] = lb;
      }
  } while(k);
  for (j=0; j<=d; j++)
    for (k=0; k<=d; k++)
      cg[k*vc+j] = coef[j][ww[k]];
  sa = 1.0;
  for (j=d; j>0; j--)  /* eliminate v[ww[j]] */
  { lb = xxa[j]/sa;
    for (k=0; k<j; k++) /* Interpolate edge v[ww[k]],v[ww[j]] */
    { cg[k*vc] = cubeint(lb,&v[ww[k]*d],&v[ww[j]*d],&cg[k*vc],&cg[j*vc],d);
      for (i=1; i<=d; i++)
        cg[k*vc+i] = (1-lb)*((1-lb)*cg[k*vc+i]+lb*cg[j*vc+i]);
    }
    sa -= xxa[j];
    if (sa<=0) j = 0;
  }
  return(cg[0]);
}

double clotoch(xev,coef,ce,p,xxa)
double *xev, **coef, *xxa;
INT *ce, p;
{ double cfo[3], cfe[3], cg[9], *va, *vb, *vc,
    l0, nm[3], na, nb, nc, *xl, *xr, *xz, d0, d1, lb, dlt, gam;
  INT i, w[3], cfl, cfr;
  if (p==1)
    return(xxa[0]*coef[0][ce[0]]+xxa[1]*coef[0][ce[1]]+xxa[2]*coef[0][ce[2]]);
  if (xxa[2]<=MIN(xxa[0],xxa[1]))
  { w[0] = ce[0]; w[1] = ce[1]; w[2] = ce[2]; }
  else
  if (xxa[1]<xxa[0])
  { w[0] = ce[0]; w[1] = ce[2]; w[2] = ce[1];
    lb = xxa[1]; xxa[1] = xxa[2]; xxa[2] = lb;
  }
  else
  { w[0] = ce[2]; w[1] = ce[1]; w[2] = ce[0];
    lb = xxa[0]; xxa[0] = xxa[2]; xxa[2] = lb;
  }
  
/* set cg to values and derivatives on standard triangle */
  va = &xev[2*w[0]]; vb = &xev[2*w[1]]; vc = &xev[2*w[2]];
  for (i=0; i<3; i++)
  { cg[3*i] = coef[0][w[i]];
    cg[3*i+1] = ((vb[0]-va[0])*coef[1][w[i]]     
                +(vb[1]-va[1])*coef[2][w[i]])/2;  /* df/dx */
    cg[3*i+2] = ((2*vc[0]-vb[0]-va[0])*coef[1][w[i]]
                +(2*vc[1]-vb[1]-va[1])*coef[2][w[i]])/2.0; /* sqrt{3} df/dy */
  }
  dlt = (vb[0]-va[0])*(vc[1]-va[1])-(vc[0]-va[0])*(vb[1]-va[1]);
  /* Twice area; +ve if abc antic.wise  -ve is abc c.wise */
  cfo[0] = (cg[0]+cg[3]+cg[6])/3;
  cfo[1] = (2*cg[0]-cg[3]-cg[6])/4;
  cfo[2] = (2*cg[3]-cg[0]-cg[6])/4;
  na = -cg[1]+cg[2];  /* perp. deriv, rel. length 2 */
  nb = -cg[4]-cg[5];
  nc = 2*cg[7];
  cfo[1] += (nb-nc)/16;
  cfo[2] += (nc-na)/16;
  na = -cg[1]-cg[2]/3.0;  /* derivatives back to origin */
  nb =  cg[4]-cg[5]/3.0;
  nc =        cg[8]/1.5;
  cfo[0] -= (na+nb+nc)*7/54;
  cfo[1] += 13*(nb+nc-2*na)/144;
  cfo[2] += 13*(na+nc-2*nb)/144;
  for (i=0; i<3; i++)
  { /* Outward normals by linear interpolation on original triangle.
       Convert to outward normals on standard triangle.
       Actually, computed to opposite corner */
    switch(i)
    { case 0: xl = vc; xr = vb; xz = va; cfl = w[2]; cfr = w[1];
              break;
      case 1: xl = va; xr = vc; xz = vb; cfl = w[0]; cfr = w[2];
              break;
      case 2: xl = vb; xr = va; xz = vc; cfl = w[1]; cfr = w[0];
              break;
    }
    na = xr[0]-xl[0]; nb = xr[1]-xl[1];
    lb = na*na+nb*nb;
    d0 = 1.5*(coef[0][cfr]-coef[0][cfl]) - 0.25*(na*(coef[1][cfl]+coef[2][cfr])
         +nb*(coef[2][cfl]+coef[2][cfr]));
    d1 = 0.5*( na*(coef[2][cfl]+coef[2][cfr])-nb*(coef[1][cfl]+coef[1][cfr]) );
    l0 = (xz[0]-xl[0])*na+(xz[1]-xl[1])*nb-lb/2;
    nm[i] = (d1*dlt-l0*d0)/lb;
  }
  cfo[0] -= (nm[0]+nm[1]+nm[2])*4/81;
  cfo[1] += (2*nm[0]-nm[1]-nm[2])/27;
  cfo[2] += (2*nm[1]-nm[0]-nm[2])/27;

  gam = xxa[0]+xxa[1]-2*xxa[2];
  if (gam==0) return(cfo[0]);
  lb = (xxa[0]-xxa[2])/gam;
  d0 = -2*cg[4]; d1 = -2*cg[1];
  cfe[0] = cubint(lb,cg[3],cg[0],d0,d1);
  cfe[1] = cubintd(lb,cg[3],cg[0],d0,d1);
  cfe[2] = -(1-lb)*(1-2*lb)*cg[5] + 4*lb*(1-lb)*nm[2] - lb*(2*lb-1)*cg[2];
  d0 = 2*(lb*cfo[1]+(1-lb)*cfo[2]);
  d1 = (lb-0.5)*cfe[1]+cfe[2]/3.0;
  return(cubint(gam,cfo[0],cfe[0],d0,d1));
}

void getvertval(tr,coef,i,nc)
struct tree *tr;
double **coef;
INT i, nc;
{ double dx, P, le;
  INT d, il, ih, j;
  if (tr->s[i]==0) return;
  il = tr->lo[i]; getvertval(tr,coef,il,nc);
  ih = tr->hi[i]; getvertval(tr,coef,ih,nc);
  coef[0][i] = (coef[0][il]+coef[0][ih])/2;
  if (nc==1) return;
  d = tr->mi[MDIM];
  P = 1.5*(coef[0][ih]-coef[0][il]);
  le = 0.0;
  for (j=0; j<d; j++)
  { dx = tr->xev[ih*d+j]-tr->xev[il*d+j];
    coef[0][i] += dx*(coef[j+1][il]-coef[j+1][ih])/8;
    coef[j+1][i] = (coef[j+1][il]+coef[j+1][ih])/2;
    P -= 1.5*dx*coef[j+1][i];
    le += dx*dx;
  }
  for (j=0; j<d; j++)
    coef[j+1][i] += P*(tr->xev[ih*d+j]-tr->xev[il*d+j])/le;
}

void exvval(cf,vv,nv,nc,d)
double **cf, *vv;
INT nv, nc, d;
{ INT i, j, k;
  vv[0] = cf[0][nv];
  if (nc==1) return;
  for (j=0; j<d; j++)
    vv[1<<j] = cf[j+1][nv];
  if (nc==2) return;
  k = d+1;
  for (i=0; i<d; i++)
    for (j=i+1; j<d; j++)
      vv[(1<<i)+(1<<j)] = cf[k++][nv];
}

void exvvalpv(vv,vl,vr,d,k,dl)
double *vv, *vl, *vr, dl;
INT d, k;
{ INT i, tk, td;
  double f0, f1;
  tk = 1<<k;
  td = 1<<d;
  for (i=0; i<td; i++) if ((i&tk)==0)
  { f0 = (vl[i]+vr[i])/2 + dl*(vl[i+tk]-vr[i+tk])/8;
    f1 = 1.5*(vr[i]-vl[i])/dl - (vl[i+tk]+vr[i+tk])/4;
    vv[i] = f0;
    vv[i+tk] = f1;
  }
}

double intqgr(x,xev,vv,ll,ur,d,nc)
double *x, *xev, vv[64][64], *ll, *ur;
INT d, nc;
{ double xx[MXDIM], phi[4];
  INT i, j, k, tk;
  if (nc==1) /* no derivatives - use multilinear interpolation */
  { for (i=d-1; i>=0; i--)
    { tk = 1<<i;
      for (j=0; j<tk; j++)
      { hermite1(x[i]-ll[i],ur[i]-ll[i],phi);
        vv[j][0] = phi[0]*vv[j][0]+phi[1]*vv[j+tk][0];
      }
    }
    return(vv[0][0]);
  }
  for (i=d-1; i>=0; i--)
  { hermite2(x[i]-ll[i],ur[i]-ll[i],phi);
    tk = 1<<i;
    phi[2] *= ur[i]-ll[i];
    phi[3] *= ur[i]-ll[i];
    for (j=0; j<tk; j++)
      for (k=0; k<tk; k++)
        vv[j][k] = phi[0]*vv[j][k] + phi[1]*vv[j+tk][k]
                 + phi[2]*vv[j][k+tk] + phi[3]*vv[j+tk][k+tk];
  }
  return(vv[0][0]);
}

double htreint(tr,x,coef,nc)
struct tree *tr;
double *x, **coef;
INT nc;
{ double vv[64][64], *ll, *ur, xx[MXDIM];
  INT d, i, j, lo, tk, ns, nv, vc, ce[64];
  d = tr->mi[MDIM];
  vc = 1<<tr->mi[MDIM];
  for (i=0; i<vc; i++)
  { for (j=0; j<vc; j++)
      vv[i][j] = 0.0;
    exvval(coef,vv[i],i,nc,d);
    ce[i] = tr->ce[i];
  }
  ns = 0;
  while(ns!=-1)
  { ll = &tr->xev[ce[0]*d]; ur = &tr->xev[ce[vc-1]*d];
    ns = needtosplitq(tr,ce,xx,ll,ur);
    if (ns!=-1)
    { tk = 1<<ns;
      lo = (2*(x[ns]-ll[ns])) < (ur[ns]-ll[ns]);
      for (i=0; i<vc; i++) if ((tk&i)==0)
      { nv = newsplit((struct design *)NULL,tr,ce[i],ce[i+tk],0,0);
        if (lf_error) return(0.0);
        if (lo)
        { ce[i+tk] = nv;
          if (tr->s[nv]) exvvalpv(vv[i+tk],vv[i],vv[i+tk],d,ns,ur[ns]-ll[ns]);
                    else exvval(coef,vv[i+tk],nv,nc,d);
        }
        else
        { ce[i] = nv;
          if (tr->s[nv]) exvvalpv(vv[i],vv[i],vv[i+tk],d,ns,ur[ns]-ll[ns]);
                    else exvval(coef,vv[i],nv,nc,d);
      } }
  } }
  ll = &tr->xev[ce[0]*d]; ur = &tr->xev[ce[vc-1]*d];
  return(intqgr(x,tr->xev,vv,ll,ur,d,nc));
}

double gridint(tr,x,coef,nc)
struct tree *tr;
double *x, **coef;
INT nc;
{ INT d, i, j, jj, v[MXDIM], vc, z0, sk, nce[1024];
  double *ll, *ur, vv[64][64];
  d = tr->mi[MDIM];
  ll = tr->xev; ur = &tr->xev[(tr->nv-1)*d];
  z0 = 0; vc = 1<<d;
  for (j=d-1; j>=0; j--)
  { v[j] = (INT)((tr->mg[j]-1)*(x[j]-ll[j])/(ur[j]-ll[j]));
    if (v[j]<0) v[j]=0;
    if (v[j]>=tr->mg[j]-1) v[j] = tr->mg[j]-2;
    z0 = z0*tr->mg[j]+v[j];
  }
  nce[0] = z0; nce[1] = z0+1; sk = jj = 1; 
  for (i=1; i<d; i++)
  { sk *= tr->mg[i-1];
    jj<<=1;
    for (j=0; j<jj; j++)
      nce[j+jj] = nce[j]+sk;
  }
  for (i=0; i<vc; i++)
    exvval(coef,vv[i],nce[i],nc,d);
  ll = &tr->xev[nce[0]*d];
  ur = &tr->xev[nce[vc-1]*d];
  return(intqgr(x,tr->xev,vv,ll,ur,d,nc));
}

double triint(tr,x,coef,nc)
struct tree *tr;
double *x, **coef;
INT nc;
{ INT d, vc, i, j, *ce, nce[1+MXDIM];
  double xa[1+MXDIM];
  d = tr->mi[MDIM]; vc = d+1;
  ce = tr->ce;
  i = 0;
  while ((i<tr->nce) && (!intri(x,&ce[i*vc],tr->xev,xa,d))) i++;
  if (i==tr->nce) return(NOSLN);
  i *= vc;
  for (j=0; j<vc; j++) nce[j] = ce[i+j];
  descend(tr,xa,nce);
  for (i=0; i<vc; i++) getvertval(tr,coef,nce[i],nc);
  return((d==2) ? clotoch(tr->xev,coef,nce,nc,xa) :
                 interptr(tr->xev,coef,nce,d,nc,xa));
}

INT dcoef(dv,nd,d)
INT *dv, nd, d;
{ INT i, j, c;
  if (nd==0) return(0);
  if (nd==1) return(dv[0]);
  for (i=0; i<nd-1; i++)
    for (j=nd-1; j>i; j--)
      if (dv[j]>dv[j-1])
      { c = dv[j]; dv[j] = dv[j-1]; dv[j-1] = c; }
  if (nd==2) return(dv[1]+dv[0]*(2*d-dv[0]+1)/2);
  c = 0; dv[nd] = 0;
  for (i=nd; i>0; i--)
    c += calcp(i,d,KSPH)-calcp(i,d-dv[i-1]+dv[i],KSPH);
  return(c);
}

INT dvect(lf,z,dv,nd,what)
struct tree *lf;
double **z;
INT *dv, nd, what;
{ INT d, i, j, k, dw[2+MXDEG];
  d = lf->mi[MDIM];
  switch(what)
  { case PCOEF:
      for (i=0; i<nd; i++) dw[i] = dv[i];
      z[0] = &lf->coef[dcoef(dw,nd,d)*lf->nvm];
      if (nd==lf->mi[MDEG]) return(1);
      for (i=1; i<=d; i++)
      { for (j=0; j<nd; j++) dw[j] = dv[j];
        dw[nd] = i;
        z[i] = &lf->coef[dcoef(dw,nd+1,d)*lf->nvm];
      }
      if ((nd==0) & (lf->mi[MDEG]>=nd+2) & (lf->mi[MKT]==KSPH)) /* add mixed deriv's */
      { k = d+1;
        for (i=1; i<=d; i++)
          for (j=i+1; j<=d; j++)
          { dw[nd] = i; dw[nd+1] = j;
            z[k] = &lf->coef[dcoef(dw,nd+2,d)*lf->nvm];
            k++;
          }
        return(3);
      }
      return(2);
    case PT0:
      z[0] = lf->t0;
      if (lf->mi[MDEG]==0) return(1);
      for (i=1; i<=d; i++) z[i] = &lf->t0[i*lf->nvm];
      return(2);
    case PNLX:
      if (nd==0)
      { z[0] = lf->nlx;
        if (lf->mi[MDEG]==0) return(1);
        for (i=1; i<=d; i++) z[i] = &lf->nlx[i*lf->nvm];
        return(2);
      }
      z[0] = &lf->nlx[(d+dcoef(dv,nd,d))*lf->nvm];
      return(1);
    case PBAND:
      z[0] = lf->h;
      return(1);
    case PDEGR:
      z[0] = lf->deg;
      return(1);
  }
  ERROR(("dvect: what????")) return(0);
}

double dointpointpf(lf,des,x,what)
struct tree *lf;
struct design *des;
double *x;
INT what;
{ double trc[6], t0;
  locfit(lf,des,x,0.0,0,0,0);
  if (what==PCOEF) return(des->cf[0]);
  ldf(lf,des,trc,0,0,&t0);
  if (what==PNLX) return(sqrt(t0));
  if (what==PT0)  return(t0);
  ERROR(("dointpointpf: invalid what"));
  return(0.0);
}

double dointpoint(lf,des,x,z,nc,what)
struct tree *lf;
struct design *des;
double *x, **z;
INT nc, what;
{ double xf, f;
  INT i;
  for (i=0; i<lf->mi[MDIM]; i++) if (lf->sty[i]==KANG)
  { xf = floor(x[i]/(2*PI*lf->sca[i]));
    x[i] -= xf*2*PI*lf->sca[i];
  }
  if (ident==1) return(dointpointpf(lf,des,x,what));
  switch(lf->mi[MEV])
  { case EGRID: f = gridint(lf,x,z,nc); break;
    case EKDTR: f = kdint(lf,x,z,nc); break;
    case ETREE: f = htreint(lf,x,z,nc); break;
    case EPHULL: f = triint(lf,x,z,nc); break;
    default: ERROR(("dointpoint: cannot interpolate this structure"));
  }
  if (what==PT0) f = f*f;
  return(f);
}

double intp(lf,des,x,what,deriv,nd,fc) /* interpolate at point */
struct tree *lf;
struct design *des;
double *x;
INT what, *deriv, nd, fc;
{ static INT nc;
  static double *z[1+MXDIM];
  if (fc) nc = dvect(lf,z,deriv,nd,what);
  return(dointpoint(lf,des,x,z,nc,what));
}

void intv(lf,des,x,f,n,what,deriv,nd) /* interpolate a vector */
struct tree *lf;
struct design *des;
double **x, *f;
INT n, what, *deriv, nd;
{ INT i, j, d, k, nc;
  double *z[1+MXDIM], xx[MXDIM];
  nc = dvect(lf,z,deriv,nd,what);
  d = lf->mi[MDIM];
  if ((lf->mi[MEV]==EDATA)|(lf->mi[MEV]==ECROS))
  { k = (n==lf->nv);
    if (k)
      for (i=0; i<d; i++)
        for (j=0; j<n; j++) k = k & (lf->xev[j*d+i]==x[i][j]);
    if (!k)
      ERROR(("Cannot interpolate when struct = data"))
    else
      for (i=0; i<n; i++) f[i] = z[0][i];
  }
  else
    for (i=0; i<n; i++)
    { for (j=0; j<d; j++) xx[j] = x[j][i];
      f[i] = dointpoint(lf,des,xx,z,nc,what);
    }
}

void intg(lf,des,x,f,mg,what,deriv,nd) /* interpolate a grid given margins */
struct design *des;
struct tree *lf;
double **x, *f;
INT *mg, what, *deriv, nd;
{ INT i, ii, j, d, k, nc;
  double *z[1+MXDIM], xv[MXDIM];
  nc = dvect(lf,z,deriv,nd,what);
  d = lf->mi[MDIM];
  k = 1;
  for (i=0; i<d; i++) k *= mg[i];
  for (i=0; i<k; i++)
  { ii = i;
    for (j=d-1; j>=0; j--)
    { xv[j] = x[j][ii%mg[j]];
      ii /= mg[j];
    }
    f[i] = dointpoint(lf,des,xv,z,nc,what);
    if (lf_error) return;
  }
}

void intf(lf,des,f,what,deriv,nd) /* `interpolate' at fitted values */
struct tree *lf;
struct design *des;
double *f;
INT what, *deriv, nd;
{ INT i;
  double *z[1+MXDIM];
  dvect(lf,z,deriv,nd,what);
  for (i=0; i<lf->nv; i++) f[i] = (what==PT0) ? SQR(z[0][i]) : z[0][i];
}

void intd(lf,des,f,what,deriv,nd) /* interpolate at data points */
struct tree *lf;
struct design *des;
double *f;
INT what, *deriv, nd;
{ INT i, j, nc;
  double *z[1+MXDIM], xx[MXDIM];
  nc = dvect(lf,z,deriv,nd,what);
  if ((lf->mi[MEV]==EDATA)|(lf->mi[MEV]==ECROS))
    for (i=0; i<lf->mi[MN]; i++) f[i] = z[0][i];
  else
    for (i=0; i<lf->mi[MN]; i++)
    { for (j=0; j<lf->mi[MDIM]; j++) xx[j] = lf->x[j][i];
      f[i] = dointpoint(lf,des,xx,z,nc,what);
    }
}
