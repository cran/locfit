#include "local.h"
#include <math.h>
#include <stdio.h>

FILE *ofile;

struct design des;
struct tree lf;
struct arc aru, vb;

INT lf_error;

void Sfit(x,y,c,w,b,lim,mi,dp,str,sca,wdes,wtre,nvc,
  iwork,lw,mg,L,kap,deriv,nd,sty)
double *x, *y, *c, *w, *b, *lim, *dp, *sca, *wdes, *wtre, *L, *kap;
INT *mi, *nvc, *iwork, *lw, *mg, *deriv, *nd, *sty;
char **str;
{ INT n, d, dv[10], i, kk;
  lf_error = 0;
  n = mi[MN]; d = mi[MDIM];
  for (i=0; i<d; i++)
  { lf.x[i] = &x[i*n];
    lf.sty[i] = sty[i];
  }
  lf.y = y; lf.w = w;
  lf.base = b;
  lf.c = c; lf.xl = lim;

  des.dw = wdes; des.lw = lw[0];
  des.ind= iwork;des.li = n; iwork += n;

  lf.tw = wtre; lf.ltw = lw[1];
  lf.fl = &lim[2*d];
  lf.iw = iwork; lf.liw = lw[2]-n;
  lf.mg = mg; lf.sca = sca;
  lf.mi = mi; lf.dp = dp;
  lf.L = L;
  lf.nvm = nvc[0];
  lf.xev = wtre; /* must do that for EPRESET */
  setstrval(mi,MKER ,str[0]);
  setstrval(mi,MTG  ,str[1]);
  setstrval(mi,MLINK,str[2]);
  setstrval(mi,MIT  ,str[3]);
  setstrval(mi,MACRI,str[4]);

  if (lf_error) return;

  switch(mi[MGETH])
  { case 0: /* the standard fit */
    case 4: /* for gam.lf, return residuals */
      if (mi[MDEG0]==mi[MDEG])
      { evaluator(&des,&lf,procv);
        if ((mi[MEV]!=EKDCE) & (mi[MEV]!=EPRES)) ressumm(&lf,&des);
      }
      else evaluator(&des,&lf,procvvord);
      break;
    case 1: /* hat martrix */
      lf.deriv = deriv; lf.nd = *nd;
      evaluator(&des,&lf,procvhatm);
      break;
    case 2: /* compute kappa0 */
      for (i=0; i< *nd; i++) dv[i] = deriv[i]-1;
      constants(&des,&lf,kap,dv,*nd);
      return;
    case 3: /* regression bandwidth selection */
      rband(&des,&lf,kap,deriv,nd,&kk);
      return;
  }

  nvc[0] = lf.nvm;
  nvc[1] = lf.ncm;
  nvc[3] = lf.nv;
  nvc[4] = lf.nce;
}

void Spred(xev,coef,sv,ce,x,res,cse,sca,m,nvc,mi,dp,mg,deriv,nd,sty)
double *xev, *coef, *sv, *x, *res, *sca, *dp;
INT *ce, *cse, *m, *nvc, *mi, *mg, *deriv, *nd, *sty;
{ INT p, i, vc;
  double *xx[MXDIM];
  lf_error = 0; p = mi[MP];
  lf.ncm = nvc[1]; lf.nnl = nvc[2]; lf.nv = lf.nvm = nvc[3]; lf.nce = nvc[4];
  lf.xev = xev;
  lf.coef = coef; coef += lf.nv*mi[MP];
  lf.nlx  = coef; coef += lf.nv*(mi[MDIM]+mi[MP]);
  lf.t0   = coef; coef += lf.nv*(mi[MDIM]+1);
  lf.lik  = coef; coef += lf.nv*3;
  lf.h    = coef; coef += lf.nv;
  lf.deg  = coef; coef += lf.nv;
  lf.sv = sv;   lf.sca = sca;
  lf.mg = mg;
  switch(mi[MEV])
  { case ETREE:
    case EKDTR:
    case EGRID: vc = 1<<mi[MDIM]; break;
    case EPHULL: vc = mi[MDIM]+1; break;
    case ECROS:
    case EDATA: break;
    default: ERROR(("Spred: Invalid ev"))
  }
  for (i=0; i<mi[MDIM]; i++)
  { xx[i] = &x[i**m];
    lf.sty[i] = sty[i];
  }

  lf.vc = vc;
  lf.ce = ce; ce += nvc[4]*vc;
  lf.s  = ce; ce += MAX(nvc[3],nvc[4]);
  lf.lo = ce; ce += MAX(nvc[3],nvc[4]);
  lf.hi = ce; ce += MAX(nvc[3],nvc[4]);
  lf.dp = dp; lf.mi = mi;
  if (lf_error) return;
  if ((lf.mi[MEV]==EDATA) | (lf.mi[MEV]==ECROS))
  { intf(&lf,&des,res,mi[MWH],deriv,*nd);
    if (*cse) intf(&lf,&des,&res[*m],PNLX,deriv,*nd);
  }
  else
  { intv(&lf,&des,xx,res,*m,mi[MWH],deriv,*nd);
    if (*cse) intv(&lf,&des,xx,&res[*m],*m,PNLX,deriv,*nd);
  }
}

void triterm(xev,h,ce,lo,hi,sca,nvc,mi,dp,nt,term)
double *xev, *h, *sca, *dp;
INT *ce, *lo, *hi, *nvc, *mi, *nt, *term;
{ INT d, i, t;
  lf_error = 0;
  lf.nv = lf.nvm = nvc[3];
  lf.xev = xev; lf.sca = sca;
  lf.h = h;
  lf.ce = ce;
  lf.lo = lo;
  lf.hi = hi;
  lf.dp = dp; lf.mi = mi; *nt = 0;
  d = mi[MDIM];
  if (mi[MEV]==ETREE)
    growquad(NULL,&lf,lf.ce,0,nt,term,xev,&xev[d*((1<<d)-1)]);
  else
    for (i=0; i<nvc[4]; i++)
      growtri(NULL,&lf,&lf.ce[i*lf.vc],0,nt,term);
}

void kdeb(x,mi,band,ind,h0,h1,meth,nmeth,ker)
double *x, *band, *h0, *h1;
INT *mi, *ind, *meth, *nmeth, *ker;
{ kdeselect(band,x,ind,*h0,*h1,meth,*nmeth,*ker,mi[MN]);
}

void compcv(k0,d,cov,z)
double *k0, *z, *cov;
INT *d;
{ *z = cv(k0,2+(*d>1),*d,1-*cov,10,2,0.0);
}
