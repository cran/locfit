/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

INT cvi=-1;

/*
 * the C version has its own checkvarlen, file src-c/vari.c.
 * For other versions, use the simple function here.
 */
#ifndef CVERSION
vari *createvar(name,type,n,mode)
varname name;
INT type, n, mode;
{ vari *v;
  v = (vari *)malloc(sizeof(vari));
  vdptr(v) = (double *)calloc( vlength(v)=n ,
              (mode==VDOUBLE) ? sizeof(double) : sizeof(INT));
  return(v);
}
#endif

/*
 * trchck checks the working space on the lfit structure 
 * has space for nvm vertices and ncm cells.
 */
void trchck(tr,nvm,ncm,d,p,vc)
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

/*
  add a new vertex at the midpoint of (x[i0],x[i1]).
  return the vertex number.
*/
INT newsplit(des,lf,i0,i1,pv)
design *des;
lfit *lf;
INT i0, i1, pv;
{ INT i, nv;

  /* first, check to see if the new point already exists */
  if (i0>i1) ISWAP(i0,i1);
  nv = lf->nv;
  for (i=i1+1; i<nv; i++)
    if ((lf->lo[i]==i0) && (lf->hi[i]==i1)) return(i);
  
  /* the point is new. Now check we have space for the new point. */
  if (nv==lf->nvm)
  {
#ifdef CVERSION
    reassign(lf);
#else
    ERROR(("newsplit: out of vertex space"));
    return(-1);
#endif
  }

  /* compute the new point, and evaluate the fit */
  lf->lo[nv] = i0;
  lf->hi[nv] = i1;
  for (i=0; i<lf->mi[MDIM]; i++)
    evptx(lf,nv,i) = (evptx(lf,i0,i)+evptx(lf,i1,i))/2;
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
