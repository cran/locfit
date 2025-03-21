/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

void solve(double *A, double *b, int d);
void triang_guessnv(int *nvm, int *ncm, int *vc, int d, int mk);
int triang_split(lfit *lf, Sint *ce, double *le);
void resort(int *pv, double *xev, int *dig);
void triang_grow(design *des, lfit *lf, Sint *ce, Sint *ct, Sint *term);
void triang_descend(lfit *tr, double *xa, Sint *ce);
void covrofdata(lfdata *lfd, double *V, double *mn);
int intri(double *x, Sint *w, double *xev, double *xa, int d);
void triang_start(design *des, lfit *lf);
double triang_cubicint(double *v, double *vv, Sint *w, int d, int nc, double *xxa);
double triang_clotoch(double *xev, double *vv, Sint *ce, int p, double *xxa);
int triang_getvertexvals(fitpt *fp, evstruc *evs, double *vv, int i, int what);
double triang_int(lfit *lf, double *x, int what);

void solve(double *A, double *b, int d) /* this is crude! A organized by column. */
/* solve(A,b,d) double *A, *b; int d; */
{ int i, j, k;
  double piv;
  for (i=0; i<d; i++)
  { piv = A[(d+1)*i];
    for (j=i; j<d; j++) A[j*d+i] /= piv;
    b[i] /= piv;
    for (j=0; j<d; j++) if (j != i)
    { piv = A[i*d+j];
      A[i*d+j] = 0;
      for (k=i+1; k<d; k++)
        A[k*d+j] -= piv*A[k*d+i];
      b[j] -= piv*b[i];
    }
  }
}

void triang_guessnv(int *nvm, int *ncm, int *vc, int d, int mk)
/* triang_guessnv(nvm,ncm,vc,d,mk) int *nvm, *ncm, *vc, d, mk; */
{
  *nvm = *ncm = mk*d;
  *vc = d+1;
  return;         
}

int triang_split(lfit *lf, Sint *ce, double *le)
/* triang_split(lf,ce,le) lfit *lf; Sint *ce; double *le; */
{ int d, i, j, k, nts, vc;
  double di, dfx[MXDIM];
  nts = 0; d = lf->fp.d; vc = d+1;
  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
    { for (k=0; k<d; k++)
        dfx[k] = evptx(&lf->fp,ce[i],k)-evptx(&lf->fp,ce[j],k);
      di = rho(dfx,lf->lfd.sca,d,KSPH,NULL);
      le[i*vc+j] = le[j*vc+i] = di/MIN(lf->fp.h[ce[i]],lf->fp.h[ce[j]]);
      nts = nts || le[i*vc+j]>cut(&lf->evs);
    }
  return(nts);
}

void resort(int *pv, double *xev, int *dig)
/* resort(pv,xev,dig) double *xev; int *pv, *dig; */
{ double d0, d1, d2;
  int i;
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

void triang_grow(design *des, lfit *lf, Sint *ce, Sint *ct, Sint *term)
/* triang_grow(des,lf,ce,ct,term) design *des; lfit *lf; Sint *ce, *ct, *term; */
{ double le[(1+MXDIM)*(1+MXDIM)], ml;
  int d, i, j, im=0, jm=0, vc, pv[(1+MXDIM)*(1+MXDIM)], dig[6];
  Sint nce[1+MXDIM];
  if (lf_error) return;
  d = lf->fp.d; vc = d+1;
  if (!triang_split(lf,ce,le))
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
    pv[0] = newsplit(des,lf,(int)ce[im],(int)ce[jm],0);
    for (i=0; i<vc; i++) nce[i] = ce[i];
    nce[im] = pv[0]; triang_grow(des,lf,nce,ct,term); nce[im] = ce[im];
    nce[jm] = pv[0]; triang_grow(des,lf,nce,ct,term);
    return;
  }

  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
      pv[i*vc+j] = pv[j*vc+i]
        = newsplit(des,lf,(int)ce[i],(int)ce[j],le[i*vc+j]<=cut(&lf->evs));
  for (i=0; i<=d; i++) /* corners */
  { for (j=0; j<=d; j++) nce[j] = (j==i) ? ce[i] : pv[i*vc+j];
    triang_grow(des,lf,nce,ct,term);
  }
  
  if (d==2) /* center for d=2 */
  { nce[0] = pv[5]; nce[1] = pv[2]; nce[2] = pv[1];
    triang_grow(des,lf,nce,ct,term);
  }
  if (d==3) /* center for d=3 */
  { resort(pv,evp(&lf->fp),dig);
    nce[0] = dig[0]; nce[1] = dig[1];
    nce[2] = dig[2]; nce[3] = dig[4]; triang_grow(des,lf,nce,ct,term);
    nce[2] = dig[5]; nce[3] = dig[3]; triang_grow(des,lf,nce,ct,term);
    nce[2] = dig[2]; nce[3] = dig[5]; triang_grow(des,lf,nce,ct,term);
    nce[2] = dig[4]; nce[3] = dig[3]; triang_grow(des,lf,nce,ct,term);
  }
}

void triang_descend(lfit *tr, double *xa, Sint *ce)
/* triang_descend(tr,xa,ce) lfit *tr; double *xa; Sint *ce; */
{ double le[(1+MXDIM)*(1+MXDIM)], ml;
  int d, vc, i, j, im=0, jm=0, pv[(1+MXDIM)*(1+MXDIM)];
  design *des;
  des = NULL;
  if (!triang_split(tr,ce,le)) return;
  d = tr->fp.d; vc = d+1;

  if (d>3) /* split longest edge */
  { ml = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<vc; j++)
        if (le[i*vc+j]>ml) { ml = le[i*vc+j]; im = i; jm = j; }
    pv[0] = newsplit(des,tr,(int)ce[im],(int)ce[jm],0);
    if (xa[im]>xa[jm])
    { xa[im] -= xa[jm]; xa[jm] *= 2; ce[jm] = pv[0]; }
    else
    { xa[jm] -= xa[im]; xa[im] *= 2; ce[im] = pv[0]; }
    triang_descend(tr,xa,ce);
    return;
  }

  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
      pv[i*vc+j] = pv[j*vc+i]
        = newsplit(des,tr,(int)ce[i],(int)ce[j],le[i*d+j]<=cut(&tr->evs));
  for (i=0; i<=d; i++) if (xa[i]>=0.5) /* in corner */
  { for (j=0; j<=d; j++)
    { if (i!=j) ce[j] = pv[i*vc+j];
      xa[j] = 2*xa[j];
    }
    xa[i] -= 1;
    triang_descend(tr,xa,ce);
    return;
  }
  if (d==1) { ERROR(("weights sum to < 1")); }
  if (d==2) /* center */
  { ce[0] = pv[5]; xa[0] = 1-2*xa[0];
    ce[1] = pv[2]; xa[1] = 1-2*xa[1];
    ce[2] = pv[1]; xa[2] = 1-2*xa[2];
    triang_descend(tr,xa,ce);
  }
  if (d==3) /* center */
  { double z; int dig[6];
    resort(pv,evp(&tr->fp),dig);
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
    triang_descend(tr,xa,ce);
} }

void covrofdata(lfdata *lfd, double *V, double *mn) /* covar of data; mean in mn */
/* covrofdata(lfd,V,mn) lfdata *lfd; double *V, *mn; */
{ int d, i, j, k;
  double s;
  s = 0; d = lfd->d;
  for (i=0; i<d*d; i++) V[i] = 0;
  for (i=0; i<lfd->n; i++)
  { s += prwt(lfd,i);
    for (j=0; j<d; j++)
      for (k=0; k<d; k++)
        V[j*d+k] += prwt(lfd,i)*(datum(lfd,j,i)-mn[j])*(datum(lfd,k,i)-mn[k]);
  }
  for (i=0; i<d*d; i++) V[i] /= s;
}

int intri(double *x, Sint *w, double *xev, double *xa, int d) /* is x in triangle bounded by xd[0..d-1]? */
/* intri(x,w,xev,xa,d) double *x, *xev, *xa; Sint *w; int d; */
/* is x in triangle bounded by xd[0..d-1]? */
{ int i, j;
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

void triang_start(design *des, lfit *lf)
/* triang_start(des,lf) design *des; lfit *lf; */
/* Triangulation with polyhedral start */
{
  int i, j, k, n, d, nc, nvm, ncm, vc;
  Sint *ce, ed[1+MXDIM];
  double V[MXDIM*MXDIM], P[MXDIM*MXDIM], sigma, z[MXDIM], xa[1+MXDIM], *xev;
  xev = evp(&lf->fp);
  d = lf->lfd.d; n = lf->lfd.n;
  lf->fp.nv = nc = 0;

  triang_guessnv(&nvm,&ncm,&vc,d,mk(&lf->evs));
  trchck(lf,nvm,ncm,vc);

  ce = lf->evs.ce;
  for (j=0; j<d; j++) xev[j] = lf->pc.xbar[j];
  lf->fp.nv = 1;
  covrofdata(&lf->lfd,V,xev); /* fix this with scaling */
  eig_dec(V,P,d);

  for (i=0; i<d; i++) /* add vertices +- 2sigma*eigenvect */
  { sigma = sqrt(V[i*(d+1)]);
    for (j=0; j<d; j++)
      xev[lf->fp.nv*d+j] = xev[j]-2*sigma*P[j*d+i];
    lf->fp.nv++;
    for (j=0; j<d; j++)
      xev[lf->fp.nv*d+j] = xev[j]+2*sigma*P[j*d+i];
    lf->fp.nv++;
  }

  for (i=0; i<n; i++) /* is point i inside? */
  { ed[0] = 0;
    for (j=0; j<d; j++)
    { z[j] = 0;
      for (k=0; k<d; k++) z[j] += P[k*d+j]*(datum(&lf->lfd,k,i)-xev[k]);
      ed[j+1] = 2*j+1+(z[j]>0);
      for (k=0; k<d; k++) z[j] = datum(&lf->lfd,j,i);
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

  for (i=0; i<lf->fp.nv; i++)
  { des->vfun(des,lf,i);
    if (lf_error) return;
    lf->evs.s[i] = 0;
  }
  for (i=0; i<nc; i++)
    triang_grow(des,lf,&ce[i*vc],NULL,NULL);
  lf->evs.nce = nc;
}

double triang_cubicint(double *v, double *vv, Sint *w, int d, int nc, double *xxa)
/* triang_cubicint(v,vv,w,d,nc,xxa) double *v, *vv, *xxa; int d, nc; Sint *w; */
{ double sa, lb, *vert0, *vert1, *vals0=NULL, *vals1, deriv0, deriv1;
  int i, j, k;
  if (nc==1) /* linear interpolate */
  { sa = 0;
    for (i=0; i<=d; i++) sa += xxa[i]*vv[i];
    return(sa);
  }
  sa = 1.0;
  for (j=d; j>0; j--)  /* eliminate v[w[j]] */
  { lb = xxa[j]/sa;
    for (k=0; k<j; k++) /* Interpolate edge v[w[k]],v[w[j]] */
    { vert0 = &v[w[k]*d];
      vert1 = &v[w[j]*d];
      vals0 = &vv[k*nc];
      vals1 = &vv[j*nc];
      deriv0 = deriv1 = 0;
      for (i=0; i<d; i++)
      { deriv0 += (vert1[i]-vert0[i])*vals0[i+1];
        deriv1 += (vert1[i]-vert0[i])*vals1[i+1];
      }
      vals0[0] = cubic_interp(lb,vals0[0],vals1[0],deriv0,deriv1);
      for (i=1; i<=d; i++)
        vals0[i] = (1-lb)*((1-lb)*vals0[i]+lb*vals1[i]);
    }
    sa -= xxa[j];
    if (sa<=0) j = 0;
  }
  return(vals0[0]);
}

double triang_clotoch(double *xev, double *vv, Sint *ce, int p, double *xxa)
/* triang_clotoch(xev,vv,ce,p,xxa) double *xev, *vv, *xxa; int p; Sint *ce; */
{ double cfo[3], cfe[3], cg[9], *va, *vb, *vc,
    l0, nm[3], na, nb, nc, *xl, *xr, *xz, d0, d1, lb, dlt, gam;
  int i, w[3], cfl, cfr;
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
  cfe[0] = cubic_interp(lb,cg[3],cg[0],d0,d1);
  cfe[1] = cubintd(lb,cg[3],cg[0],d0,d1);
  cfe[2] = -(1-lb)*(1-2*lb)*cg[5] + 4*lb*(1-lb)*nm[2] - lb*(2*lb-1)*cg[2];
  d0 = 2*(lb*cfo[1]+(1-lb)*cfo[2]);
  d1 = (lb-0.5)*cfe[1]+cfe[2]/3.0;
  return(cubic_interp(gam,cfo[0],cfe[0],d0,d1));
}

int triang_getvertexvals(fitpt *fp, evstruc *evs, double *vv, int i, int what)
/* triang_getvertexvals(fp,evs,vv,i,what) fitpt *fp; evstruc *evs; double *vv; int i, what; */
{ double dx, P, le, vl[1+MXDIM], vh[1+MXDIM];
  int d, il, ih, j, nc;
  d = fp->d;
  if (evs->s[i]==0) return(exvval(fp,vv,i,d,what,0));

  il = evs->lo[i]; nc = triang_getvertexvals(fp,evs,vl,il,what);
  ih = evs->hi[i]; nc = triang_getvertexvals(fp,evs,vh,ih,what);
  vv[0] = (vl[0]+vh[0])/2;
  if (nc==1) return(nc);
  P = 1.5*(vh[0]-vl[0]);
  le = 0.0;
  for (j=0; j<d; j++)
  { dx = evptx(fp,ih,j)-evptx(fp,il,j);
    vv[0] += dx*(vl[j+1]-vh[j+1])/8;
    vv[j+1] = (vl[j+1]+vh[j+1])/2;
    P -= 1.5*dx*vv[j+1];
    le += dx*dx;
  }
  for (j=0; j<d; j++)
    vv[j+1] += P*(evptx(fp,ih,j)-evptx(fp,il,j))/le;
  return(nc);
}

double triang_int(lfit *lf, double *x, int what)
/* triang_int(lf,x,what) lfit *lf; double *x; int what; */
{
  int d, i, j, k, vc, nc;
  Sint *ce, nce[1+MXDIM];
  double xa[1+MXDIM], vv[(1+MXDIM)*(1+MXDIM)], lb;
fitpt *fp;
evstruc *evs;
fp = &lf->fp;
evs= &lf->evs;

  d = fp->d; vc = d+1;
  ce = evs->ce;
  i = 0;
  while ((i<evs->nce) && (!intri(x,&ce[i*vc],evp(fp),xa,d))) i++;
  if (i==evs->nce) return(NOSLN);
  i *= vc;
  for (j=0; j<vc; j++) nce[j] = ce[i+j];
  triang_descend(lf,xa,nce);

  /* order the vertices -- needed for asymmetric interptr */
  do
  { k=0;
    for (i=0; i<d; i++)
      if (nce[i]>nce[i+1])
      { j=nce[i]; nce[i]=nce[i+1]; nce[i+1]=j; k=1;
        lb = xa[i]; xa[i] = xa[i+1]; xa[i+1] = lb;
      }
  } while(k);
  nc = 0;
  for (i=0; i<vc; i++)
    nc =  triang_getvertexvals(fp,evs,&vv[i*nc],nce[i],what);
  return((d==2) ? triang_clotoch(evp(fp),vv,nce,nc,xa) :
                 triang_cubicint(evp(fp),vv,nce,d,nc,xa));
}
