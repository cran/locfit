/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

extern INT ident;
extern double robscale;

INT calcp(mi,deg)
INT *mi, deg;
{ INT i, k;
  if (mi[MUBAS]) return(mi[MDEG]);

  if (mi[MKT]==KPROD) return(mi[MDIM]*deg+1);
  k = 1;
  for (i=1; i<=deg; i++) k = k*(mi[MDIM]+i)/i;
  return(k);
}

double resp(lf,i)
lfit *lf;
INT i;
{ if (lf->y==NULL) return(0.0);
  return(lf->y[i]);
}

extern INT cvi;

double prwt(lf,i)
lfit *lf;
INT i;
{ if (i==cvi) return(0.0);
  if (lf->w==NULL) return(1.0);
  return(lf->w[i]);
}

double base(lf,i)
lfit *lf;
INT i;
{ if (lf->base==NULL) return(0.0);
  return(lf->base[i]);
}

double cens(lf,i)
lfit *lf;
INT i;
{ if (lf->c==NULL) return(0.0);
  return(lf->c[i]);
}

INT coefnumber(deriv,nd,kt,d,deg)
INT *deriv, nd, kt, d, deg;
{ INT d0, d1, t;
  if (d==1)
  { if (nd<=deg) return(nd);
    return(-1);
  }
  if (nd==0) return(0);
  if (deg==0) return(-1);
  if (nd==1) return(1+deriv[0]);
  if (deg==1) return(-1);
  if (kt==KPROD) return(-1);
  if (nd==2)
  { d0 = deriv[0]; d1 = deriv[1];
    if (d0<d1) { t = d0; d0 = d1; d1 = t; }
    return((d+1)*(d0+1)-d0*(d0+3)/2+d1);
  }
  if (deg==2) return(-1);
  ERROR(("coefnumber not programmed for nd>=3"));
  return(-1);
}

void fitfunangl(dx,ff,sca,cd,deg)
double dx, *ff, sca;
INT deg, cd;
{ switch(cd)
  { case 0:
      ff[0] = 1;
      ff[1] = sin(dx/sca)*sca;
      ff[2] = (1-cos(dx/sca))*sca*sca;
      if (deg>=3)
        WARN(("Can't handle angular model with deg>=3"));
      return;
    case 1:
      ff[0] = 0;
      ff[1] = cos(dx/sca);
      ff[2] = sin(dx/sca)*sca;
      return;
    case 2:
      ff[0] = 0;
      ff[1] = -sin(dx/sca)/sca;
      ff[2] = cos(dx/sca);
      return;
    default: WARN(("Can't handle angular model with >2 derivs"));
  }
}

/* x is the data point.
   t is the centre point.
*/
void fitfun(lf,x,t,f,der,nd)
lfit *lf;
double *x, *t, *f;
INT *der, nd;
{ INT d, deg, m, i, j, k, cd[MXDIM];
  double ff[MXDIM][1+MXDEG], dx[MXDIM];

#ifdef SVERSION
  if (lf->mi[MUBAS])
  { basis(x,t,f,lf->mi[MDIM],lf->mi[MP]);
    return;
  }
#endif

  d = lf->mi[MDIM];
  deg = lf->mi[MDEG];
  m = 0;
  f[m++] = (nd==0);
  if (deg==0) return;

  for (i=0; i<d; i++)
  { cd[i] = 0;
    dx[i] = (t==NULL) ? x[i] : x[i]-t[i];
  }
  for (i=0; i<nd; i++) cd[der[i]]++;

  for (i=0; i<d; i++)
  { switch(lf->sty[i])
    {
      case STANGL:
        fitfunangl(dx[i],ff[i],lf->sca[i],cd[i],lf->mi[MDEG]);
        break;
      default:
        for (j=0; j<cd[i]; j++) ff[i][j] = 0.0;
        ff[i][cd[i]] = 1.0;
        for (j=cd[i]+1; j<=deg; j++)
          ff[i][j] = ff[i][j-1]*dx[i]/(j-cd[i]);
    }
  }

  for (i=0; i<d; i++) f[m++] = (cd[i]==nd) ? ff[i][1] : 0.0;
  if (deg==1) return;

  for (i=0; i<d; i++)
  { f[m++] = (cd[i]==nd) ? ff[i][2] : 0.0;
    if (lf->mi[MKT]!=KPROD)
      for (j=i+1; j<d; j++)
        f[m++] = (cd[i]+cd[j]==nd) ? ff[i][1]*ff[j][1] : 0.0;
  }
  if (deg==2) return;

  for (i=0; i<d; i++)
  { f[m++] = (cd[i]==nd) ? ff[i][3] : 0.0;
    if (lf->mi[MKT]!=KPROD)
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
  ERROR(("fitfun: can't handle deg=%d",deg));
}

double vocri(lk,t0,t2,pen)
double lk, t0, t2, pen;
{ if (pen==0) return(-2*t0*lk/((t0-t2)*(t0-t2)));
  return((-2*lk+pen*t2)/t0);
}

INT procvraw(des,lf,v)
design *des;
lfit *lf;
INT v;
{ INT lf_status;
  double h;
  des->xev = evpt(lf,v);

  lf_status = ainitband(des,lf);
  if (lf_error) return(lf_status);

  switch(lf->mi[MACRI])
  { case AKAT:
    case ACP:
    case AMDI:
      h = aband2(des,lf,des->h);
      h = aband3(des,lf,h);
      h = nbhd(lf,des,0,h,1);
      lf_status = locfit(lf,des,h,0);
      break;
    case ANONE:
    case AOK:
      break;
  }

  lf->h[v] = des->h;

  if (lf_error) return(lf_status);
  if (lf->mi[MDC]) dercor(lf,des,lf->h[v]);
  subparcomp(des,lf,v);
  lf->deg[v] = lf->mi[MDEG];
  return(lf_status);
}

void prvsetz(lf,nvm,v,d)
lfit *lf;
INT nvm, v;
{ INT i;
  lf->lik[v] = lf->lik[nvm+v] = 0;
  lf->lik[2*nvm+v] = 0; /* should use sum of weights here? */
  for (i=0; i<=d; i++)
    lf->t0[i*nvm+v] = lf->nlx[i*nvm+v] = 0.0;
}

INT procv(des,lf,v)
design *des;
lfit *lf;
INT v;
{ INT d, p, nvm, i, k;
  double trc[6], t0[1+MXDIM];
  k = procvraw(des,lf,v);
  if (lf_error) return(k);
   
  d = lf->mi[MDIM]; p = lf->mi[MP];
  nvm = lf->nvm;

  switch(k)
  { case LF_OK: break;
    case LF_NCON:
      WARN(("procv: locfit did not converge"));
      break;
    case LF_OOB:
      WARN(("procv: parameters out of bounds"));
      break;
    case LF_PF:
      WARN(("procv: perfect fit"));
      prvsetz(lf,nvm,v,d);
      return(k);
    case LF_NOPT:
      WARN(("procv: no points with non-zero weight"));
      prvsetz(lf,nvm,v,d);
      return(k);
    case LF_INFA:
      WARN(("procv: initial value problem"));
      prvsetz(lf,nvm,v,d);
      return(k);
    case LF_DEMP:
      WARN(("procv: density estimate, empty integration region"));
      prvsetz(lf,nvm,v,d);
      return(k);
    case LF_XOOR:
      WARN(("procv: fit point outside xlim region"));
      prvsetz(lf,nvm,v,d);
      return(k);
    case LF_DNOP:
      WARN(("density estimation -- insufficient points in smoothing window"));
      prvsetz(lf,nvm,v,d);
      return(k);
    case LF_FPROB:
      WARN(("procv: f problem; likelihood failure"));
      prvsetz(lf,nvm,v,d);
      return(k);
    default:
      WARN(("procv: unknown return code %d",k));
      return(k);
  }

  ldf(lf,des,trc,0,lf->mi,t0);
  lf->lik[v] = des->llk;
  lf->lik[nvm+v] = trc[2];
  lf->lik[2*nvm+v] = trc[0]-trc[2];
  lf->nlx[v] = sqrt(des->V[0]);
  lf->t0[v] = sqrt(t0[0]);
  if (p>1)
  { if (lf->t0[v]==0)
      for (i=1; i<=d; i++) lf->t0[i*nvm+v] = 0.0;
    else
      for (i=1; i<=d; i++) lf->t0[i*nvm+v] = t0[i]/lf->t0[v];
  }
  subparcomp2(des,lf,v);
  return(k);
}

double intvo(des,lf,c0,c1,a,p,t0,t20,t21)
design *des;
lfit *lf;
double *c0, *c1, a, t0, t20, t21;
INT p;
{ double th, lk, link[LLEN];
  INT i;
  lk = 0;
  for (i=0; i<des->n; i++)
  { th = (1-a)*innerprod(c0,&des->X[i*p],p) + a*innerprod(c1,&des->X[i*p],p);
    stdlinks(link,lf,des->ind[i],th,robscale);
    lk += des->w[i]*link[ZLIK];
  }
  des->llk = lk;
  return(vocri(des->llk,t0,(1-a)*t20+a*t21,lf->dp[DADP]));
}

INT procvvord(des,lf,v)
design *des;
lfit *lf;
INT v;
{ double tr[6], gcv, g0, ap, coef[4][10], t2[4], th, md;
  INT i, j, k, d1, *mi, i0, p1, ip;
  mi = lf->mi;
  des->xev = evpt(lf,v);

  lf->h[v] = nbhd(lf,des,(INT)(mi[MN]*lf->dp[DALP]),lf->dp[DFXH],0);
  if (lf->h[v]<=0) WARN(("zero bandwidth in procvvord"));

  ap = lf->dp[DADP];
  if ((ap==0) & ((mi[MTG]&63)!=TGAUS)) ap = 2.0;
  d1 = mi[MDEG]; p1 = mi[MP];
  for (i=0; i<p1; i++) coef[0][i] = coef[1][i] = coef[2][i] = coef[3][i] = 0.0;
  i0 = 0; g0 = 0;
  ip = 1;
  for (i=mi[MDEG0]; i<=d1; i++)
  { mi[MDEG] = i; des->p = mi[MP] = calcp(mi,i);
    k = locfit(lf,des,lf->h[v],0);

    ldf(lf,des,tr,1,lf->mi,NULL);
    gcv = vocri(des->llk,tr[0],tr[2],ap);
    if ((i==mi[MDEG0]) || (gcv<g0)) { i0 = i; g0 = gcv; md = i; }

    for (j=0; j<des->p; j++) coef[i][j] = des->cf[j];
    t2[i] = tr[2];

#ifdef RESEARCH
    if ((ip) && (i>mi[MDEG0]))
    { for (j=1; j<10; j++)
      { gcv = intvo(des,lf,coef[i-1],coef[i],j/10.0,des->p,tr[0],t2[i-1],t2[i]);
        if (gcv<g0) { g0 = gcv; md = i-1+j/10.0; }
      }
    }
#endif
  }

  if (i0<d1) /* recompute the best fit */
  { mi[MDEG] = i0; des->p = mi[MP] = calcp(mi,i0);
    k = locfit(lf,des,lf->h[v],0);
    for (i=mi[MP]; i<p1; i++) des->cf[i] = 0.0;
    i0 = (INT)md; if (i0==d1) i0--;
    th = md-i0;
    for (i=0; i<p1; i++) des->cf[i] = (1-th)*coef[i0][i]+th*coef[i0+1][i];
    mi[MDEG] = d1; mi[MP] = p1;
  }

  for (i=0; i<p1; i++) lf->coef[i*lf->nvm+v] = des->cf[i];
  lf->deg[v] = md;
  return(k);
}

/* special version of ressumm to estimate sigma^2, with derivative estimation */
void ressummd(lf,des)
lfit *lf;
design *des;
{ INT i;
  double s0, s1;
  s0 = s1 = 0.0;
  if ((lf->mi[MTG]&64)==0)
  { lf->dp[DRV] = 1.0;
    return;
  }
  for (i=0; i<lf->nv; i++)
  { s0 += lf->lik[2*lf->nvm+i];
    s1 += lf->lik[i];
  }
  if (s0==0.0)
    lf->dp[DRV] = 0.0;
  else
    lf->dp[DRV] = -2*s1/s0;
}

void ressumm(lf,des)
lfit *lf;
design *des;
{ INT i, j, ev, tg, orth;
  double *dp, *oy, pw, r1, r2, rdf, t0, t1, u[MXDIM], link[LLEN];
  dp = lf->dp;
  dp[DLK] = dp[DT0] = dp[DT1] = 0;
  if ((lf->mi[MEV]==EKDCE) | (lf->mi[MEV]==EPRES))
  { dp[DRV] = 1.0;
    return;
  }
  if (lf->nd>0)
  { ressummd(lf,des);
    return;
  }
  r1 = r2 = 0.0;
  ev = lf->mi[MEV];
  if ((ev==EDATA) | (ev==ECROS)) ev = EFITP;
  orth = (lf->mi[MGETH]==4) | (lf->mi[MGETH]==5);
  for (i=0; i<lf->mi[MN]; i++)
  { for (j=0; j<lf->mi[MDIM]; j++) u[j] = datum(lf,j,i);
    des->th[i] = base(lf,i)+dointpoint(lf,des,u,PCOEF,ev,i);
    des->wd[i] = resp(lf,i) - des->th[i];
    des->w[i] = 1.0;
    des->ind[i] = i;
  }

  tg = lf->mi[MTG];
  lf->dp[DRSC] = 1.0;
  if ((tg==TROBT+64) | (tg==TCAUC+64)) /* global robust scale */
  { oy = lf->y; lf->y = des->wd;
    des->xev = lf->pc.xbar;
    locfit(lf,des,0.0,1);
    lf->y = oy;
    lf->dp[DRSC] = robscale;
  }

  if (orth) /* orthog. residuals */
  { int od, op;
    des->n = lf->mi[MN];
    od = lf->mi[MDEG]; op = lf->mi[MP];
    lf->mi[MDEG] = 1;
    lf->mi[MP] = des->p = 1+lf->mi[MDIM];
    oy = lf->y; lf->y = des->wd;
    des->xev = lf->pc.xbar;
    locfit(lf,des,0.0,1);
    for (i=0; i<lf->mi[MN]; i++) oy[i] = resp(lf,i) - des->th[i];
    lf->y = oy;
    lf->mi[MDEG] = od; lf->mi[MP] = op;
  }

  for (i=0; i<lf->mi[MN]; i++)
  { for (j=0; j<lf->mi[MDIM]; j++) u[j] = datum(lf,j,i);
    t0 = dointpoint(lf,des,u,PT0,ev,i);
    t1 = dointpoint(lf,des,u,PNLX,ev,i);
    stdlinks(link,lf,i,des->th[i],lf->dp[DRSC]);
    t1 = t1*t1*link[ZDDLL];
    t0 = t0*t0*link[ZDDLL];
    if (t1>1) t1 = 1;
    if (t0>1) t0 = 1; /* no observation gives >1 deg.free */
    dp[DLK] += link[ZLIK];
    dp[DT0] += t0;
    dp[DT1] += t1;
    pw = prwt(lf,i);
    if (pw>0)
    { r1 += link[ZDLL]*link[ZDLL]/pw;
      r2 += link[ZDDLL]/pw;
    }
    if (orth) des->di[i]  = t1;
  }

  if (orth) return;

  dp[DRV] = 1.0;
  if ((lf->mi[MTG]&64)==64) /* quasi family */
  { rdf = lf->mi[MN]-2*dp[DT0]+dp[DT1];
    if (rdf<1.0)
    { WARN(("Estimated rdf < 1.0; not estimating variance"));
    }
    else
      dp[DRV] = r1/r2 * lf->mi[MN] / rdf;
  }

  /* try to ensure consistency for family="circ"! */
  if (((lf->mi[MTG]&63)==TCIRC) & (lf->mi[MDIM]==1))
  { INT *ind, nv;
    double dlt, th0, th1;
    ind = des->ind;
    nv = lf->nv;
    for (i=0; i<nv; i++) ind[i] = i;
    lforder(ind,vdptr(lf->xxev),0,nv-1);
    for (i=1; i<nv; i++)
    { dlt = evptx(lf,ind[i],0)-evptx(lf,ind[i-1],0);
      th0 = lf->coef[ind[i]]-dlt*lf->coef[ind[i]+nv]-lf->coef[ind[i-1]];
      th1 = lf->coef[ind[i]]-dlt*lf->coef[ind[i-1]+nv]-lf->coef[ind[i-1]];
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
lfit *lf;
design *des;
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
