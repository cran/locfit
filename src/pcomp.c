/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

/* functions for computing and subtracting, adding the
   parametric component
*/

#include "local.h"

INT noparcomp(lf)
lfit *lf;
{ INT tg;
  if (lf->sty[0]==STUSER) return(1);
  if (lf->mi[MDEG0]<lf->mi[MDEG]) return(1);
  tg = lf->mi[MTG] & 63;
  if (tg<=THAZ) return(1);
  if (tg==TROBT) return(1);
  if (tg==TCAUC) return(1);
  return(0);
}

INT hasparcomp(lf)
lfit *lf;
{ return(lf->mi[MPC]);
}

INT lenpc(d,p,lc)
INT d, p, lc;
{ return((lc) ? d + 4*p + 2*p*p : 0);
}

void pcchk(pc,d,p,lc)
paramcomp *pc;
INT d, p, lc;
{ INT k;
  double *z;
  pc->wk = checkvarlen(pc->wk,lenpc(d,p,lc),"_pcwork",VDOUBLE);
  z = vdptr(pc->wk);
  k = 0;
  pc->xbar = &z[k]; k = k+d;
  pc->coef = &z[k]; k = k+p;
  pc->f    = &z[k]; k = k+p;

  pc->xtwx.Z = &z[k]; k += p*p;
  pc->xtwx.Q = &z[k]; k += p*p;
  pc->xtwx.dg = &z[k]; k += p;
  pc->xtwx.f2 = &z[k]; k += p;
  pc->xtwx.p = p;
}

void compparcomp(des,lf,nopc)
design *des;
lfit *lf;
INT nopc;
{ INT i, j, k;
  double wt, sw;
  paramcomp *pc;
  pc = &lf->pc;
  pcchk(pc,lf->mi[MDIM],lf->mi[MP],1);
  for (i=0; i<lf->mi[MDIM]; i++) pc->xbar[i] = 0.0;
  sw = 0.0;
  for (i=0; i<lf->mi[MN]; i++)
  { wt = prwt(lf,i);
    sw += wt;
    for (j=0; j<lf->mi[MDIM]; j++)
      pc->xbar[j] += datum(lf,j,i)*wt;
    des->ind[i] = i;
    des->w[i] = 1.0;
  }
  for (i=0; i<lf->mi[MDIM]; i++) pc->xbar[i] /= sw;
  if ((nopc) || noparcomp(lf))
  { lf->mi[MPC] = 0;
    return;
  }
  lf->mi[MPC] = 1;
  des->xev = pc->xbar;
  k = locfit(lf,des,0.0,0);
  if (lf_error) return;
  switch(k)
  { case LF_NOPT:
      ERROR(("compparcomp: no points in dataset?"));
      return;
    case LF_INFA:
      ERROR(("compparcomp: infinite parameters in param. component"));
      return;
    case LF_NCON:
      ERROR(("compparcom: not converged"));
      return;
    case LF_OOB:
      ERROR(("compparcomp: parameters out of bounds"));
      return;
    case LF_PF:
      WARN(("compparcomp: perfect fit"));
    case LF_OK:
      for (i=0; i<lf->mi[MP]; i++)
      { pc->coef[i] = des->cf[i];
        pc->xtwx.dg[i] = des->xtwx.dg[i];
        pc->xtwx.f2[i] = des->xtwx.f2[i];
      }
      for (i=0; i<lf->mi[MP]*lf->mi[MP]; i++)
      { pc->xtwx.Z[i] = des->xtwx.Z[i];
        pc->xtwx.Q[i] = des->xtwx.Q[i];
      }
      pc->xtwx.sm = des->xtwx.sm;
      return;
    default:
      ERROR(("compparcomp: locfit unknown return status %d",k));
      return;
  }
}

void subparcomp(des,lf,v)
design *des;
lfit *lf;
{ double coef[1+MXDIM];
  INT d, deg, i, *deriv, nd, ncoef, *mi;

  mi = lf->mi;
  deg = mi[MDEG];
  d = mi[MDIM];
  deriv = lf->deriv;
  nd = lf->nd;
  ncoef = 1+d*(deg>nd);

  coef[0] = des->cf[coefnumber(deriv,nd,mi[MKT],d,deg)];
  for (i=0; i<d; i++)
  { deriv[nd] = i;
    coef[i+1] = des->cf[coefnumber(deriv,nd+1,mi[MKT],d,deg)];
  }

  if (hasparcomp(lf)) /* subtract parametric component */
  { fitfun(lf,evpt(lf,v),lf->pc.xbar,des->f1,deriv,nd);
    coef[0] -= innerprod(lf->pc.coef,des->f1,lf->mi[MP]);
    if (deg>0)
      for (i=0; i<d; i++)
      { deriv[nd] = i;
        fitfun(lf,evpt(lf,v),lf->pc.xbar,des->f1,deriv,nd+1);
        coef[i+1] -= innerprod(lf->pc.coef,des->f1,lf->mi[MP]);
      }
  }

  /* store coefficients in lf->coef */
  for (i=0; i<ncoef; i++)
    lf->coef[i*lf->nvm+v] = coef[i];
}

void subparcomp2(des,lf,v)
design *des;
lfit *lf;
{ double *coef, vcoef[1+MXDIM], t0, t1;
  INT d, deg, i, *deriv, nd, ncoef, *mi;

  mi = lf->mi;
  deg = mi[MDEG];
  d = mi[MDIM];
  deriv = lf->deriv;
  nd = lf->nd;
  ncoef = 1+d*(deg>nd);

  coef = &des->V[mi[MP]*coefnumber(deriv,nd,mi[MKT],d,deg)];
  vcoef[0] = sqrt(coef[coefnumber(deriv,nd,mi[MKT],d,deg)]);
  if (vcoef[0]==0.0)
    for (i=0; i<d; i++) vcoef[i] = 0.0;
  else
    for (i=0; i<d; i++)
    { deriv[nd] = i;
      vcoef[i+1] = coef[coefnumber(deriv,nd+1,mi[MKT],d,deg)]/vcoef[0];
    }

  if (hasparcomp(lf))
  { fitfun(lf,evpt(lf,v),lf->pc.xbar,des->f1,deriv,nd);
    for (i=0; i<mi[MP]; i++) lf->pc.f[i] = des->f1[i];
    vxtwx(&lf->pc.xtwx,des->f1,mi[MP]);
    t0 = sqrt(innerprod(lf->pc.f,des->f1,mi[MP]));
    vcoef[0] -= t0;
    lf->t0[v]  -= t0;
    if ((deg>0) && (t0>0))
      for (i=0; i<d; i++)
      { deriv[nd] = i;
        fitfun(lf,evpt(lf,v),lf->pc.xbar,lf->pc.f,deriv,nd+1);
        t1 = innerprod(lf->pc.f,des->f1,mi[MP])/t0;
        vcoef[i+1] -= t1;
        lf->t0[(i+1)*lf->nvm+v]  -= t1;
      }
  }

  /* store coefficients in lf->nlx */
  for (i=0; i<ncoef; i++)
    lf->nlx[i*lf->nvm+v] = vcoef[i];
}

double addparcomp(lf,x,c)
lfit *lf;
double *x;
int c;
{ if (!hasparcomp(lf)) return(0.0);
  fitfun(lf,x,lf->pc.xbar,lf->pc.f,lf->deriv,lf->nd);
  if (c==PCOEF) return(innerprod(lf->pc.coef,lf->pc.f,lf->mi[MP]));
  if ((c==PNLX)|(c==PT0)|(c==PVARI)) return(sqrt(vxtwxv(&lf->pc.xtwx,lf->pc.f)));
  return(0.0);
}
