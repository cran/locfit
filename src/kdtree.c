/*
 *   Copyright (c) 1998-1999 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"
extern double atan2();

void ksmall(), newcell();
INT cvi=-1, exvval();

void trchck(tr,nvm,ncm,d,p,vc) /* check space for nvm vertices, ncm cells */
lfit *tr;
INT nvm, ncm, d, p, vc;
{ INT rw, *k;
  double *z;

  tr->xxev = checkvarlen(tr->xxev,d*nvm,"_lfxev",VDOUBLE);

  rw = nvm*(3*d+8)+ncm;
  tr->tw = checkvarlen(tr->tw,rw,"_lfwork",VDOUBLE);
  z = (double *)vdptr(tr->tw);
  tr->coef= z; z += nvm*(d+1);
  tr->nlx = z; z += nvm*(d+1);
  tr->t0  = z; z += nvm*(d+1);
  tr->lik = z; z += 3*nvm;
  tr->h   = z; z += nvm;
  tr->deg = z; z += nvm;
  tr->sv  = z; z += ncm;
  if (z != (double *)vdptr(tr->tw)+rw)
    WARN(("trchck: double assign problem"));

  rw = ncm*vc+3*MAX(ncm,nvm);
  tr->iw = checkvarlen(tr->iw,rw,"_lfiwork",VINT);
  k = (INT *)vdptr(tr->iw);
  tr->ce = k; k += vc*ncm;
  tr->s  = k; k += MAX(ncm,nvm);
  tr->lo = k; k += MAX(ncm,nvm);
  tr->hi = k; k += MAX(ncm,nvm);
  if (k != (INT *)vdptr(tr->iw)+rw)
    WARN(("trchck: int assign problem"));

  tr->nvm = nvm; tr->ncm = ncm; tr->mi[MDIM] = d; tr->mi[MP] = p; tr->vc = vc;
}

#ifdef CVERSION
void reassign(lf)
lfit *lf;
{ INT i, nvm, ncm, vc, d, k, p, *iw;
  double *tw, *ntw;
  setvarname(lf->tw,"__lfwork"); /* prevent overwrite */
  setvarname(lf->iw,"__lfiwork");
  nvm = lf->nvm; ncm = lf->ncm; vc = lf->vc;
  tw = (double *)vdptr(lf->tw);
  iw = (INT *)vdptr(lf->iw);
  d = lf->mi[MDIM];
  p = lf->mi[MP];
  trchck(lf,2*nvm,ncm,d,p,vc);
  ntw = vdptr(lf->tw);
/*
  xev is stored in blocks of d. other matrices by blocks on nvm
*/
  k = nvm*d;
  memcpy(vdptr(lf->xxev),tw,k*sizeof(double));
  tw += k; ntw += 2*k;
  for (i=0; i<2*p+2*d+6; i++)
  { memcpy(ntw,tw,nvm*sizeof(double));
    tw += nvm; ntw += 2*nvm;
  }
  k = ncm;       memcpy(lf->sv,tw,k*sizeof(double));   tw += k;

  k = vc*ncm;       memcpy(lf->ce,iw,k*sizeof(INT)); iw += k;
  k = MAX(ncm,nvm); memcpy(lf->s,iw,k*sizeof(INT));  iw += k;
  k = MAX(ncm,nvm); memcpy(lf->lo,iw,k*sizeof(INT)); iw += k;
  k = MAX(ncm,nvm); memcpy(lf->hi,iw,k*sizeof(INT)); iw += k;
  deletename("__lfwork");
  deletename("__lfiwork");
}
#endif

void dataf(des,lf)
design *des;
lfit *lf;
{ INT d, i, j, nv;
  nv = lf->mi[MN]; d = lf->mi[MDIM];
  trchck(lf,nv,0,d,des->p,0);
  for (i=0; i<nv; i++)
    for (j=0; j<d; j++) evptx(lf,i,j) = datum(lf,j,i);
  for (i=0; i<nv; i++)
  { des->vfun(des,lf,i);
    lf->s[i] = 0;
  }
  lf->nv = nv; lf->nce = 0;
}

void xbarf(des,lf)
design *des;
lfit *lf;
{ int i, d;
  d = lf->mi[MDIM];
  trchck(lf,1,0,lf->mi[MDIM],des->p,0);
  for (i=0; i<d; i++) evptx(lf,0,i) = lf->pc.xbar[i];
  des->vfun(des,lf,0);
  lf->s[0] = 0;
  lf->nv = 1; lf->nce = 0;
}

void preset(des,tr)
design *des;
lfit *tr;
{ INT i, nv;
  double *tmp;
  nv = tr->nvm;
  tmp = vdptr(tr->xxev);
  trchck(tr,nv,0,tr->mi[MDIM],des->p,0);
  tr->xxev->dpr = tmp;
  for (i=0; i<nv; i++)
  { des->vfun(des,tr,i);
    tr->s[i] = 0;
  }
  tr->nv = nv; tr->nce = 0;
}

void crossf(des,lf)
design *des;
lfit *lf;
{ INT d, i, j, n;
  n = lf->mi[MN]; d = lf->mi[MDIM];
  trchck(lf,n,0,d,des->p,0);
  for (i=0; i<n; i++)
    for (j=0; j<d; j++) evptx(lf,i,j) = datum(lf,j,i);
  for (cvi=0; cvi<n; cvi++)
  { lf->s[cvi] = 0;
    des->vfun(des,lf,cvi);
  }
  cvi = -1;
  lf->nv = n; lf->nce = 0; lf->mi[MN] = n;
}

void gridf(des,tr)
design *des;
lfit *tr;
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
      evptx(tr,i,j) = (tr->mg[j]==1) ? tr->fl[j] :
                      tr->fl[j]+(tr->fl[j+d]-tr->fl[j])*u/(tr->mg[j]-1);
      z = z/tr->mg[j];
    }
    tr->s[i] = 0;
    des->vfun(des,tr,i);
  }
  tr->nv = nv; tr->nce = 0;
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
lfit *tr;
INT p, *pi, d, fc, *m;
double *sv;
{ INT i, k, l, u, im;
  double al, be, mx, t;
  mx = 0.0; im = 0;
  l = tr->lo[p]; u = tr->hi[p];
  if (u-l < fc) return(-1);
  if (d>1)
    for (k=0; k<d; k++)
    { al = be = datum(tr,k,pi[l]);
      for (i=l+1; i<=u; i++)
      { t = datum(tr,k,pi[i]);
        if (t<al) al = t;
        if (t>be) be = t;
      }
      be = (be-al)/tr->sca[k];
      if (be>mx)
      { mx = be;
        im = k;
    } }
  *m = (u+l)/2;
  ksmall(l,u,*m, dvari(tr,im), pi);
  *sv = datum(tr,im,pi[*m]);
  while ((*m<u-1) && (datum(tr,im,pi[*m+1])==*sv)) (*m)++;
  if (*m==u) im = -1;
  return(im);
}

void kdtree(des,tr)
design *des;
lfit *tr;
{ INT i, j, vc, d, nc, nv, ncm, nvm, k, m, n, p, *pi, fc;
  double sv, ll[MXDIM], ur[MXDIM];
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
    default: ERROR(("kdtree: invalid ev"));
  }
  trchck(tr,nvm,ncm,d,des->p,vc);
  nv = 0;
  if (tr->mi[MEV] != EKDCE)
  { for (i=0; i<vc; i++)
    { 
      j = i;
      for (k=0; k<d; ++k)
      { evptx(tr,i,k) = tr->fl[d*(j%2)+k];
        j >>= 1;
      }
    }
    nv = vc;
    for (j=0; j<vc; j++) tr->ce[j] = j;
  }
  if (tr->mi[MEV]==ETREE)
  { for (j=0; j<vc; j++)
    { des->vfun(des,tr,j);
      if (lf_error) return;
      tr->s[j] = 0;
    }
    for (j=0; j<d; j++)
    { ll[j] = tr->fl[j];
      ur[j] = tr->fl[j+d];
    }
    tr->nv = nv;
    growquad(des,tr,tr->ce,NULL,NULL,ll,ur,tr->fl);
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
      { WARN(("Insufficient space for full tree"));
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
        newcell(&nv,vc,vdptr(tr->xxev), d, k, sv,
             &tr->ce[p*vc], &tr->ce[(nc-2)*vc], &tr->ce[(nc-1)*vc]);
    }
    else if (tr->mi[MEV]==EKDCE) /* new vertex at cell center */
    { sv = 0;
      for (i=0; i<d; i++) evptx(tr,nv,i) = 0;
      for (j=tr->lo[p]; j<=tr->hi[p]; j++)
      { sv += prwt(tr,pi[j]);
        for (i=0; i<d; i++)
          evptx(tr,nv,i) += datum(tr,i,pi[j])*prwt(tr,pi[j]);
      }
      for (i=0; i<d; i++) evptx(tr,nv,i) /= sv;
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

INT newsplit(des,lf,i0,i1,pv)
design *des;
lfit *lf;
INT i0, i1, pv;
{ double z[MXDIM];
  INT d, i, nv;
  d = lf->mi[MDIM];
  if (i0>i1) ISWAP(i0,i1);
  for (i=0; i<d; i++) z[i] = (evptx(lf,i0,i)+evptx(lf,i1,i))/2;
  nv = lf->nv;
  for (i=i1+1; i<nv; i++)
    if ((lf->lo[i]==i0) && (lf->hi[i]==i1)) return(i); /* duplicate */
  if (nv==lf->nvm)
  {
#ifdef CVERSION
    reassign(lf);
#else
    ERROR(("newsplit: out of vertex space"));
#endif
  }
  if (lf_error) return(-1);
  lf->lo[nv] = MIN(i0,i1); lf->hi[nv] = MAX(i0,i1);
  for (i=0; i<d; i++) evptx(lf,nv,i) = (evptx(lf,i0,i)+evptx(lf,i1,i))/2;
  if (pv) /* pseudo vertex */
  { lf->h[nv] = (lf->h[i0]+lf->h[i1])/2;
    lf->s[nv] = 1; /* pseudo-vertex */
  }
  else /* real vertex */
  { des->vfun(des,lf,nv);
    lf->s[nv] = 0;
  }
  lf->nv++;

  return(nv);
}

INT needtosplitq(lf,ce,le,ll,ur,flim)
lfit *lf;
INT *ce;
double *le, *ll, *ur, *flim;
{ INT d, vc, i, is;
  double h, hmin, score[MXDIM];
  d = lf->mi[MDIM]; vc = 1<<d;

  hmin = 0.0;
  for (i=0; i<vc; i++)
  { h = lf->h[ce[i]];
    if ((h>0) && ((hmin==0)|(h<hmin))) hmin = h;
  }

  is = 0;
  for (i=0; i<d; i++)
  { le[i] = (ur[i]-ll[i])/lf->sca[i];
    if ((lf->sty[i]==STCPAR) || (hmin==0))
      score[i] = 2*(ur[i]-ll[i])/(flim[i+d]-flim[i]);
    else
      score[i] = le[i]/hmin;
    if (score[i]>score[is]) is = i;
  }
  if (lf->dp[DCUT]<score[is]) return(is);
  return(-1);
}

void growquad(des,lf,ce,ct,term,ll,ur,flim)
design *des;
lfit *lf;
INT *ce, *ct, *term;
double *ll, *ur, *flim;
{ INT i, i0, i1, d, vc, ns, tk, nce[1024], pv;
  double le[MXDIM], z;
  d = lf->mi[MDIM]; vc = 1<<d;
  ns = needtosplitq(lf,ce,le,ll,ur,flim);
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
      pv = (lf->sty[i]!=STCPAR) &&
           (le[ns] < (lf->dp[DCUT]*MIN(lf->h[i0],lf->h[i1])));
      nce[i] = newsplit(des,lf,i0,i1,pv);
      if (lf_error) return;
    }
  }
  z = ur[ns]; ur[ns] = (z+ll[ns])/2;
  growquad(des,lf,nce,ct,term,ll,ur,flim);
  if (lf_error) return;
  ur[ns] = z;
  for (i=0; i<vc; i++)
    nce[i] = ((i&tk)== 0) ? nce[i+tk] : ce[i];
  z = ll[ns]; ll[ns] = (z+ur[ns])/2;
  growquad(des,lf,nce,ct,term,ll,ur,flim);
  ll[ns] = z;
}

INT needtosplit(tr,ce,le)
lfit *tr;
double *le;
INT *ce;
{ INT d, i, j, k, nts, vc;
  double di, dfx[MXDIM];
  nts = 0; d = tr->mi[MDIM]; vc = d+1;
  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
    { for (k=0; k<d; k++)
        dfx[k] = evptx(tr,ce[i],k)-evptx(tr,ce[j],k);
      di = rho(dfx,tr->sca,d,KSPH,NULL);
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

void growtri(des,tr,ce,ct,term)
design *des;
lfit *tr;
INT *ce, *ct, *term;
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
    pv[0] = newsplit(des,tr,ce[im],ce[jm],0);
    for (i=0; i<vc; i++) nce[i] = ce[i];
    nce[im] = pv[0]; growtri(des,tr,nce,ct,term); nce[im] = ce[im];
    nce[jm] = pv[0]; growtri(des,tr,nce,ct,term);
    return;
  }

  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
      pv[i*vc+j] = pv[j*vc+i]
        = newsplit(des,tr,ce[i],ce[j],le[i*vc+j]<=tr->dp[DCUT]);
  for (i=0; i<=d; i++) /* corners */
  { for (j=0; j<=d; j++) nce[j] = (j==i) ? ce[i] : pv[i*vc+j];
    growtri(des,tr,nce,ct,term);
  }
  
  if (d==2) /* center for d=2 */
  { nce[0] = pv[5]; nce[1] = pv[2]; nce[2] = pv[1];
    growtri(des,tr,nce,ct,term);
  }
  if (d==3) /* center for d=3 */
  { resort(pv,vdptr(tr->xxev),dig);
    nce[0] = dig[0]; nce[1] = dig[1];
    nce[2] = dig[2]; nce[3] = dig[4]; growtri(des,tr,nce,ct,term);
    nce[2] = dig[5]; nce[3] = dig[3]; growtri(des,tr,nce,ct,term);
    nce[2] = dig[2]; nce[3] = dig[5]; growtri(des,tr,nce,ct,term);
    nce[2] = dig[4]; nce[3] = dig[3]; growtri(des,tr,nce,ct,term);
  }
  if (d==1) return;
}

void descend(tr,xa,ce)
lfit *tr;
double *xa;
INT *ce;
{ double le[(1+MXDIM)*(1+MXDIM)], ml;
  INT d, vc, i, j, pv[(1+MXDIM)*(1+MXDIM)], im, jm;
  design *des;
  des = NULL;
  if (!needtosplit(tr,ce,le)) return;
  d = tr->mi[MDIM]; vc = d+1;

  if (d>3) /* split longest edge */
  { ml = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<vc; j++)
        if (le[i*vc+j]>ml) { ml = le[i*vc+j]; im = i; jm = j; }
    pv[0] = newsplit(des,tr,ce[im],ce[jm],0);
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
        = newsplit(des,tr,ce[i],ce[j],le[i*d+j]<=tr->dp[DCUT]);
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
    resort(pv,vdptr(tr->xxev),dig);
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

void covrofdata(lf,V,mn) /* covar of data; mean in mn */
lfit *lf;
double *V, *mn;
{ INT d, i, j, k;
  double s;
  s = 0; d = lf->mi[MDIM];
  for (i=0; i<d*d; i++) V[i] = 0;
  for (i=0; i<lf->mi[MN]; i++)
  { s += prwt(lf,i);
    for (j=0; j<d; j++)
      for (k=0; k<d; k++)
        V[j*d+k] += prwt(lf,i)*(datum(lf,j,i)-mn[j])*(datum(lf,k,i)-mn[k]);
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
design *des;
lfit *tr;
{ INT i, j, k, n, d, nc, nvm, ncm, vc, *ce, ed[1+MXDIM];
  double V[MXDIM*MXDIM], P[MXDIM*MXDIM], sigma, z[MXDIM], xa[1+MXDIM], *xev;
  xev = vdptr(tr->xxev);
  d = tr->mi[MDIM]; n = tr->mi[MN]; tr->nv = nc = 0;
  vc = d+1; nvm = tr->mi[MK]*d; ncm = nvm;
  trchck(tr,nvm,ncm,d,des->p,vc);
  ce = tr->ce;
  for (j=0; j<d; j++) xev[j] = tr->pc.xbar[j];
  tr->nv = 1;
  covrofdata(tr,V,tr->pc.xbar); /* fix this with scaling */
  eigen(V,P,d,tr->mi[MMXIT]);

  for (i=0; i<d; i++) /* add vertices +- 2sigma*eigenvect */
  { sigma = sqrt(V[i*(d+1)]);
    for (j=0; j<d; j++)
      xev[tr->nv*d+j] = xev[j]-2*sigma*P[j*d+i];
    tr->nv++;
    for (j=0; j<d; j++)
      xev[tr->nv*d+j] = xev[j]+2*sigma*P[j*d+i];
    tr->nv++;
  }

  for (i=0; i<n; i++) /* is point i inside? */
  { ed[0] = 0;
    for (j=0; j<d; j++)
    { z[j] = 0;
      for (k=0; k<d; k++) z[j] += P[k*d+j]*(datum(tr,k,i)-xev[k]);
      ed[j+1] = 2*j+1+(z[j]>0);
      for (k=0; k<d; k++) z[j] = datum(tr,j,i);
    }
    k = intri(z,ed,xev,xa,d);
    if (xa[0]<0)
    { for (j=1; j<=d; j++)
        for (k=0; k<d; k++)
          xev[ed[j]*d+k] = xa[0]*xev[k]+(1-xa[0])*xev[ed[j]*d+k];
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
    if (lf_error) return;
    tr->s[i] = 0;
  }
  for (i=0; i<nc; i++)
  { for (j=0; j<=d; j++) ed[j] = tr->ce[i*vc+j];
    growtri(des,tr,&tr->ce[i*vc],(INT *)NULL,(INT *)NULL);
  }
  tr->nce = nc;
}

void hermite1(x,z,phi)
double x, z, *phi;
{ if (z==0)
  { phi[0] = 1.0; phi[1] = 0.0;
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

double interpcuk(lf,ce,x,d,what)
lfit *lf;
double *x;
INT *ce, d, what;
{ double phi[4], g[2304], *ll, *ur;
  INT i, j, k, nc, vc;
  vc = lf->vc;
  for (j=0; j<vc; j++) nc = exvval(lf,&g[j*nc],ce[j],d,what,0);
  ll = evpt(lf,ce[0]); ur = evpt(lf,ce[vc-1]);
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

double blend(lf,s,x,ll,ur,j,nt,t,what)
lfit *lf;
double s, *x, *ll, *ur;
INT j, nt, *t, what;
{ INT k, k1, m, nc, j0, j1, *ce;
  double v0, v1, xibar, g0[3], g1[3], gg[4], gp[4], phi[4];
  ce = lf->ce;
  for (k=0; k<4; k++)  /* North South East West */
  { k1 = (k>1);
    v0 = ll[k1]; v1 = ur[k1];
    j0 = ce[j+2*(k==0)+(k==2)];
    j1 = ce[j+3-2*(k==1)-(k==3)];
    xibar = (k%2==0) ? ur[k<2] : ll[k<2];
    m = nt;
    while ((m>=0) && ((lf->s[t[m]] != (k<=1)) | (lf->sv[t[m]] != xibar))) m--;
    if (m >= 0)
    { m = (k%2==1) ? lf->lo[t[m]] : lf->hi[t[m]];
      while (lf->s[m] != -1)
        m = (x[lf->s[m]] < lf->sv[m]) ? lf->lo[m] : lf->hi[m];
      if (v0 < evptx(lf,ce[4*m+2*(k==1)+(k==3)],k1))
      { j0 = ce[4*m+2*(k==1)+(k==3)];
        v0 = evptx(lf,j0,k1);
      }
      if (evptx(lf,ce[4*m+3-2*(k==0)-(k==2)],k1) < v1)
      { j1 = ce[4*m+3-2*(k==0)-(k==2)];
        v1 = evptx(lf,j1,k1);
      }
    }
    nc = exvval(lf,g0,j0,2,what,0);
    nc = exvval(lf,g1,j1,2,what,0);
    if (nc==1)
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
  if (nc==1)
    for (k=0; k<2; k++)
    { hermite1(x[k]-ll[k],ur[k]-ll[k],phi);
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

double kdint(tr,x,what)
lfit *tr;
double *x;
INT what;
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
  ll = evpt(tr,ce[0]);
  ur = evpt(tr,ce[vc-1]);
  ff = interpcuk(tr,ce,x,d,what);
  if (d==2) ff = blend(tr,ff,x,ll,ur,k*vc,nt,t,what);
  return(ff);
}

double interptr(v,vv,w,d,nc,xxa)
double *v, *vv, *xxa;
INT d, *w, nc;
{ double sa, lb;
  INT i, j, k;
  if (nc==1) /* linear interpolate */
  { sa = 0;
    for (i=0; i<=d; i++) sa += xxa[i]*vv[i];
    return(sa);
  }
  sa = 1.0;
  for (j=d; j>0; j--)  /* eliminate v[w[j]] */
  { lb = xxa[j]/sa;
    for (k=0; k<j; k++) /* Interpolate edge v[w[k]],v[w[j]] */
    { vv[k*nc] = cubeint(lb,&v[w[k]*d],&v[w[j]*d],&vv[k*nc],&vv[j*nc],d);
      for (i=1; i<=d; i++)
        vv[k*nc+i] = (1-lb)*((1-lb)*vv[k*nc+i]+lb*vv[j*nc+i]);
    }
    sa -= xxa[j];
    if (sa<=0) j = 0;
  }
  return(vv[0]);
}

double clotoch(xev,vv,ce,p,xxa)
double *xev, *vv, *xxa;
INT *ce, p;
{ double cfo[3], cfe[3], cg[9], *va, *vb, *vc,
    l0, nm[3], na, nb, nc, *xl, *xr, *xz, d0, d1, lb, dlt, gam;
  INT i, w[3], cfl, cfr;
  if (p==1)
    return(xxa[0]*vv[0]+xxa[1]*vv[1]+xxa[2]*vv[2]);
  if (xxa[2]<=MIN(xxa[0],xxa[1]))
  { va = &xev[2*ce[0]]; vb = &xev[2*ce[1]]; vc = &xev[2*ce[2]];
    w[0] = 0; w[1] = 3; w[2] = 6;
  }
  else
  if (xxa[1]<xxa[0])
  { w[0] = 0; w[1] = 6; w[2] = 3;
    va = &xev[2*ce[0]]; vb = &xev[2*ce[2]]; vc = &xev[2*ce[1]];
    lb = xxa[1]; xxa[1] = xxa[2]; xxa[2] = lb;
  }
  else
  { w[0] = 6; w[1] = 3; w[2] = 0;
    va = &xev[2*ce[2]]; vb = &xev[2*ce[1]]; vc = &xev[2*ce[0]];
    lb = xxa[0]; xxa[0] = xxa[2]; xxa[2] = lb;
  }
  
/* set cg to values and derivatives on standard triangle */
  for (i=0; i<3; i++)
  { cg[3*i] = vv[w[i]];
    cg[3*i+1] = ((vb[0]-va[0])*vv[w[i]+1]
                +(vb[1]-va[1])*vv[w[i]+2])/2;  /* df/dx */
    cg[3*i+2] = ((2*vc[0]-vb[0]-va[0])*vv[w[i]+1]
                +(2*vc[1]-vb[1]-va[1])*vv[w[i]+2])/2.0; /* sqrt{3} df/dy */
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
    d0 = 1.5*(vv[cfr]-vv[cfl]) - 0.25*(na*(vv[cfl+1]+vv[cfr+1])
         +nb*(vv[cfl+2]+vv[cfr+2]));
    d1 = 0.5*( na*(vv[cfl+2]+vv[cfr+2])-nb*(vv[cfl+1]+vv[cfr+1]) );
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

INT getvertval(lf,vv,i,what)
lfit *lf;
double *vv;
INT i, what;
{ double dx, P, le, vl[1+MXDIM], vh[1+MXDIM];
  INT d, il, ih, j, nc;
  d = lf->mi[MDIM];
  if (lf->s[i]==0) return(exvval(lf,vv,i,d,what,0));

  il = lf->lo[i]; nc = getvertval(lf,vl,il,what);
  ih = lf->hi[i]; nc = getvertval(lf,vh,ih,what);
  vv[0] = (vl[0]+vh[0])/2;
  if (nc==1) return(nc);
  P = 1.5*(vh[0]-vl[0]);
  le = 0.0;
  for (j=0; j<d; j++)
  { dx = evptx(lf,ih,j)-evptx(lf,il,j);
    vv[0] += dx*(vl[j+1]-vh[j+1])/8;
    vv[j+1] = (vl[j+1]+vh[j+1])/2;
    P -= 1.5*dx*vv[j+1];
    le += dx*dx;
  }
  for (j=0; j<d; j++)
    vv[j+1] += P*(evptx(lf,ih,j)-evptx(lf,il,j))/le;
  return(nc);
}

INT exvval(lf,vv,nv,d,what,z)
lfit *lf;
double *vv;
INT nv, d, what, z;
{ INT i, k;
  double *values;
  k = (z) ? 1<<d : d+1;
  for (i=1; i<k; i++) vv[i] = 0.0;
  switch(what)
  { case PCOEF:
      values = lf->coef;
      break;
    case PVARI:
    case PNLX:
      values = lf->nlx;
      break;
    case PT0:
      values = lf->t0;
      break;
    case PBAND:
      vv[0] = lf->h[nv];
      return(1);
    case PDEGR:
      vv[0] = lf->deg[nv];
      return(1);
    case PLIK:
      vv[0] = lf->lik[nv];
      return(1);
    case PRDF:
      vv[0] = lf->lik[2*lf->nvm+nv];
      return(1);
    default:
      ERROR(("Invalid what in exvval"));
      return(0);
  }
  vv[0] = values[nv];
  if ((lf->mi[MDEG]==0) && (lf->mi[MDC]==0))return(1);
  if (z)
  { for (i=0; i<d; i++) vv[1<<i] = values[(i+1)*lf->nvm+nv];
    return(1<<d);
  }
  else
  { for (i=1;i<=d; i++) vv[i] = values[i*lf->nvm+nv];
    return(d+1);
  }
}

void exvvalpv(vv,vl,vr,d,k,dl,nc)
double *vv, *vl, *vr, dl;
INT d, k, nc;
{ INT i, tk, td;
  double f0, f1;
  if (nc==1)
  { vv[0] = (vl[0]+vr[0])/2;
    return;
  }
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
{ double phi[4];
  INT i, j, k, tk;
  tk = 1<<d;
  for (i=0; i<tk; i++) if (vv[i][0]==NOSLN) return(NOSLN);
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

double htreint(tr,x,what)
lfit *tr;
double *x;
INT what;
{ double vv[64][64], *ll, *ur, h, xx[MXDIM];
  INT d, i, lo, tk, ns, nv, nc, vc, ce[64];
  d = tr->mi[MDIM];
  vc = 1<<tr->mi[MDIM];
  for (i=0; i<vc; i++)
  { setzero(vv[i],vc);
    nc = exvval(tr,vv[i],i,d,what,1);
    ce[i] = tr->ce[i];
  }
  ns = 0;
  while(ns!=-1)
  { ll = evpt(tr,ce[0]); ur = evpt(tr,ce[vc-1]);
    ns = needtosplitq(tr,ce,xx,ll,ur,tr->fl);
    if (ns!=-1)
    { tk = 1<<ns;
      h = ur[ns]-ll[ns];
      lo = (2*(x[ns]-ll[ns])) < h;
      for (i=0; i<vc; i++) if ((tk&i)==0)
      { nv = newsplit((design *)NULL,tr,ce[i],ce[i+tk],0);
        if (lf_error) return(0.0);
        if (lo)
        { ce[i+tk] = nv;
          if (tr->s[nv]) exvvalpv(vv[i+tk],vv[i],vv[i+tk],d,ns,h,nc);
                    else exvval(tr,vv[i+tk],nv,d,what,1);
        }
        else
        { ce[i] = nv;
          if (tr->s[nv]) exvvalpv(vv[i],vv[i],vv[i+tk],d,ns,h,nc);
                    else exvval(tr,vv[i],nv,d,what,1);
      } }
  } }
  ll = evpt(tr,ce[0]); ur = evpt(tr,ce[vc-1]);
  return(intqgr(x,vdptr(tr->xxev),vv,ll,ur,d,nc));
}

double gridint(tr,x,what)
lfit *tr;
double *x;
INT what;
{ INT d, i, j, jj, v[MXDIM], vc, z0, sk, nce[1024], nc;
  double *ll, *ur, vv[64][64];
  d = tr->mi[MDIM];
  ll = evpt(tr,0); ur = evpt(tr,tr->nv-1);
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
    nc = exvval(tr,vv[i],nce[i],d,what,1);
  ll = evpt(tr,nce[0]);
  ur = evpt(tr,nce[vc-1]);
  return(intqgr(x,vdptr(tr->xxev),vv,ll,ur,d,nc));
}

double triint(tr,x,what)
lfit *tr;
double *x;
INT what;
{ INT d, vc, i, j, k, *ce, nc, nce[1+MXDIM];
  double xa[1+MXDIM], vv[(1+MXDIM)*(1+MXDIM)], lb;
  d = tr->mi[MDIM]; vc = d+1;
  ce = tr->ce;
  i = 0;
  while ((i<tr->nce) && (!intri(x,&ce[i*vc],vdptr(tr->xxev),xa,d))) i++;
  if (i==tr->nce) return(NOSLN);
  i *= vc;
  for (j=0; j<vc; j++) nce[j] = ce[i+j];
  descend(tr,xa,nce);

  /* order the vertices -- needed for asymmetric interptr */
  do
  { k=0;
    for (i=0; i<d; i++)
      if (nce[i]>nce[i+1])
      { j=nce[i]; nce[i]=nce[i+1]; nce[i+1]=j; k=1;
        lb = xa[i]; xa[i] = xa[i+1]; xa[i+1] = lb;
      }
  } while(k);
  for (i=0; i<vc; i++)
    nc =  getvertval(tr,&vv[i*nc],nce[i],what);
  return((d==2) ? clotoch(vdptr(tr->xxev),vv,nce,nc,xa) :
                 interptr(vdptr(tr->xxev),vv,nce,d,nc,xa));
}

double fitpint(lf,x,what,i)
lfit *lf;
double *x;
INT what, i;
{ double vv[1+MXDIM];
  exvval(lf,vv,i,lf->mi[MDIM],what,0);
  return(vv[0]);
}

double dointpointpf(lf,des,x,what)
lfit *lf;
design *des;
double *x;
INT what;
{ double trc[6], t0;
  locfit(lf,des,0.0,0);
  if (what==PCOEF) return(des->cf[0]);
  ldf(lf,des,trc,0,0,&t0);
  if ((what==PNLX)|(what==PT0)) return(sqrt(t0));
  ERROR(("dointpointpf: invalid what"));
  return(0.0);
}

double xbarint(lf,x,what)
lfit *lf;
double *x;
INT what;
{ INT i, nc;
  double vv[1+MXDIM], f;
  nc = exvval(lf,vv,0,lf->mi[MDIM],what,0);
  f = vv[0];
  if (nc>1)
    for (i=0; i<lf->mi[MDIM]; i++)
      f += vv[i+1]*(x[i]-evptx(lf,0,i));
  return(f);
}

double dointpoint(lf,des,x,what,ev,j)
lfit *lf;
design *des;
double *x;
INT what, ev, j;
{ double xf, f;
  INT i;
  for (i=0; i<lf->mi[MDIM]; i++) if (lf->sty[i]==STANGL)
  { xf = floor(x[i]/(2*PI*lf->sca[i]));
    x[i] -= xf*2*PI*lf->sca[i];
  }
  if (ident==1) return(dointpointpf(lf,des,x,what));
  switch(ev)
  { case EGRID: f = gridint(lf,x,what); break;
    case EKDTR: f = kdint(lf,x,what); break;
    case ETREE: f = htreint(lf,x,what); break;
    case EPHULL: f = triint(lf,x,what); break;
    case EFITP: f = fitpint(lf,x,what,j); break;
    case EXBAR: f = xbarint(lf,x,what); break;
    case ENONE: f = 0; break;
    default: ERROR(("dointpoint: cannot interpolate this structure"));
  }
  if ((what==PT0)|(what==PNLX))
  { if (f<0) f =0.0;
  }
  f += addparcomp(lf,x,what);
  return(f);
}
