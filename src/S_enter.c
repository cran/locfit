/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include "S.h"
#undef WARN
#undef ERROR
#undef INT

#include "local.h"

design des;
lfit lf;

INT lf_error;
extern INT cvi;

#ifdef RVERSION
#define CALL_S_FUNC char
#define CALL_S_ARGS void
#else
#define CALL_S_FUNC void
#define CALL_S_ARGS char
#endif
#define CALL_S_LEN  long int
#define CALL_S_NARG long int
#define CALL_S_MODE char
#define CALL_S_NRES long int
#define CALL_S_VALS char

static CALL_S_FUNC *bsfunc, *bsf2;

void basis(x,t,f,dim,p)
double *x, *t, *f;
INT dim, p;
{
  CALL_S_ARGS *args[2];
  CALL_S_LEN  length[2];
  CALL_S_NARG nargs;
  CALL_S_MODE *mode[2];
  CALL_S_NRES nres;
  CALL_S_VALS *values[1];
  double z0[1], z1[1], *vptr;
  int i;

  args[0] = (char *)x;
  mode[0] = "double";
  length[0] = dim;

  args[1] = (char *)t;
  mode[1] = "double";
  length[1] = dim;

  nargs = 2;
  nres = 1;

  call_S(bsfunc,nargs,args,mode,length,(char **)NULL,nres,values);

  vptr = (double *)values[0];
  for (i=0; i<p; i++) f[i] = vptr[i];
}

double vbasis(x,t,n,d,ind,m,p,X)
double **x, *t, *X;
int n, d, m, p, *ind;
{
  CALL_S_ARGS *args[3];
  CALL_S_LEN  length[3];
  CALL_S_NARG nargs;
  CALL_S_MODE *mode[3];
  CALL_S_NRES nres;
  CALL_S_VALS *values[1];
  double *vptr;
  int i;

  args[0] = (char *)x[0];
  mode[0] = "double";
  length[0] = n;

  args[1] = (char *)t;
  mode[1] = "double";
  length[1] = d;

  args[2] = (char *)ind;
  mode[2] = "integer";
  length[2] = m;

  nargs = 3;
  nres = 1;

  call_S(bsf2,nargs,args,mode,length,0,nres,values);
  vptr = (double *)values[0];
  for (i=0; i<m*p; i++) X[i] = vptr[i];
}

vari *setsvar(v,x,n)
vari *v;
double *x;
INT n;
{ v->dpr = x;
  v->n = n;
  return(v);
}

void slocfit(x,y,c,w,b,lim,mi,dp,str,sca,xev,wdes,wtre,wpc,nvc,
  iwork,lw,mg,L,kap,deriv,nd,sty,bs)
double *x, *y, *c, *w, *b, *lim, *dp, *sca, *xev, *L, *kap, *wdes, *wtre, *wpc;
INT *mi, *nvc, *iwork, *lw, *mg, *deriv, *nd, *sty;
char **str;
CALL_S_FUNC **bs;
{ INT n, d, i, kk;
  vari v1, v2, v3, vL, vi, vx, vxev;
  if (mi[MUBAS])
  { bsfunc = bs[0];
    bsf2 = bs[1];
  }
  lf_error = 0;
  n = mi[MN]; d = mi[MDIM];
  for (i=0; i<d; i++)
  { dvari(&lf,i) = &x[i*n];
    lf.sty[i] = sty[i];
  }
  lf.y = y; lf.w = w;
  lf.base = b;
  lf.c = c; lf.xl = lim;

  des.dw = setsvar(&v2,wdes,lw[0]);
  vx.dpr = (double *)iwork; vx.n = n; des.index = &vx;
  iwork += n;

  lf.xxev = setsvar(&vxev,xev,d*nvc[0]);
  lf.tw   = setsvar(&v1,wtre,lw[1]);
  lf.pc.wk= setsvar(&v3,wpc,lw[3]);
  vi.dpr = (double *)iwork; vi.n = lw[2]-n; lf.iw = &vi;
  lf.fl = &lim[2*d];
  lf.mg = mg; lf.sca = sca;
  lf.mi = mi; lf.dp = dp;
  vL.dpr = L; vL.n = lw[4]; lf.L = &vL;
  lf.nvm = nvc[0];
  setstrval(mi,MKER ,str[0]);
  setstrval(mi,MTG  ,str[1]);
  setstrval(mi,MLINK,str[2]);
  setstrval(mi,MIT  ,str[3]);
  setstrval(mi,MACRI,str[4]);

  lf.nd = *nd;
  for (i=0; i<*nd; i++) lf.deriv[i] = deriv[i]-1;
  if (lf_error) return;

  switch(mi[MGETH])
  { case 0: /* the standard fit */
    case 4: /* for gam.lf, return residuals */
    case 5: /* for gam.lf prediction */
      if (mi[MDEG0]==mi[MDEG])
      { startlf(&des,&lf,procv,0);
        if (!lf_error) ressumm(&lf,&des);
      }
      else startlf(&des,&lf,procvvord,0);
      break;
    case 1: /* hat matrix */
      startlf(&des,&lf,procvhatm,mi[MKER]!=WPARM);
      break;
    case 2: /* compute kappa0 */
      constants(&des,&lf,kap);
      return;
    case 3: /* regression bandwidth selection */
      rband(&des,&lf,kap,deriv,nd,&kk);
      return;
    case 6: /* lscv */
      startlf(&des,&lf,procv,1);
      dlscv(&des,&lf);
      return;
  }

  nvc[0] = lf.nvm;
  nvc[1] = lf.ncm;
  nvc[3] = lf.nv;
  nvc[4] = lf.nce;
}

void recoef(coef,ce,d,ev,nvc)
double *coef;
INT *ce, d, ev, *nvc;
{ INT vc;
  lf.coef = coef; coef += lf.nv*(d+1);
  lf.nlx  = coef; coef += lf.nv*(d+1);
  lf.t0   = coef; coef += lf.nv*(d+1);
  lf.lik  = coef; coef += lf.nv*3;
  lf.h    = coef; coef += lf.nv;
  lf.deg  = coef; coef += lf.nv;

  switch(ev)
  { case ETREE:
    case EKDTR:
    case EGRID: vc = 1<<d; break;
    case EPHULL: vc = d+1; break;
    case EXBAR:
    case ECROS:
    case EDATA:
    case EPRES:
    case ENONE: break;
    default: ERROR(("spreplot: Invalid ev"));
  }
  lf.vc = vc;
  lf.ce = ce; ce += nvc[4]*vc;
  lf.s  = ce; ce += MAX(nvc[3],nvc[4]);
  lf.lo = ce; ce += MAX(nvc[3],nvc[4]);
  lf.hi = ce; ce += MAX(nvc[3],nvc[4]);
}

void spreplot(xev,coef,sv,ce,x,res,se,wpc,sca,m,nvc,mi,dp,mg,deriv,nd,sty,where,what,bs)
double *xev, *coef, *sv, *x, *res, *se, *wpc, *sca, *dp;
INT *ce, *m, *nvc, *mi, *mg, *deriv, *nd, *sty, *where;
char **what;
void **bs;
{ INT p, i, vc;
  double *xx[MXDIM];
  vari v3, vxev;
  if (mi[MUBAS]) bsfunc = bs[0];
  lf_error = 0; p = mi[MP];
  lf.ncm = nvc[1]; lf.nv = lf.nvm = nvc[3]; lf.nce = nvc[4];

  vxev.dpr = xev; vxev.n = mi[MDIM]*lf.nv; lf.xxev = &vxev;

  recoef(coef,ce,mi[MDIM],mi[MEV],nvc);
  lf.sv = sv;   lf.sca = sca;
  lf.mg = mg;
  v3.dpr = wpc; v3.n = 0; lf.pc.wk = &v3;
  pcchk(&lf.pc,mi[MDIM],p,0);
  lf.pc.xtwx.sm = 1;
  for (i=0; i<mi[MDIM]; i++) lf.sty[i] = sty[i];

  /* set up the data structures right */
  switch (*where)
  { case 2: /* grid */
      for (i=0; i<mi[MDIM]; i++)
      { xx[i] = x;
        x += m[i];
      }
      break;
    case 1: /* vector */
    case 3: /* data */
      for (i=0; i<mi[MDIM]; i++) dvari(&lf,i) = xx[i] = &x[i**m];
      break;
    case 4: /* fit points, need nothing! */
      break;
    default:
      ERROR(("unknown where in spreplot"));
  }

  lf.dp = dp; lf.mi = mi;
  lf.nd = *nd;
  for (i=0; i<*nd; i++) lf.deriv[i] = deriv[i]-1;

  if (lf_error) return;
  preplot(&lf,&des,xx,res,se,what[1][0],m,*where,ppwhat(what[0]));
}

void sfitted(x,y,w,c,ba,fit,cv,st,xev,coef,sv,ce,wpc,sca,nvc,mi,dp,mg,deriv,nd,sty,what,bs)
double *x, *y, *w, *c, *ba, *fit, *xev, *coef, *sv, *wpc, *sca, *dp;
INT *cv, *st, *ce, *nvc, *mi, *mg, *deriv, *nd, *sty;
char **what;
void **bs;
{ INT i, n;
  vari v3, vxev;
  if (mi[MUBAS]) bsfunc = bs[0];
  n = mi[MN];
  lf_error = 0; cvi = -1;
  lf.ncm = nvc[1]; lf.nv = lf.nvm = nvc[3]; lf.nce = nvc[4];
  lf.mg = mg;
  for (i=0; i<mi[MDIM]; i++)
  { dvari(&lf,i) = &x[i*n];
    lf.sty[i] = sty[i];
  }
  lf.y = y; lf.w = w; lf.c = c; lf.base = ba;
  vxev.dpr = xev; vxev.n = mi[MDIM]*lf.nv; lf.xxev = &vxev;
  recoef(coef,ce,mi[MDIM],mi[MEV],nvc);
  lf.sv = sv;   lf.sca = sca;
  v3.dpr = wpc; v3.n = 0; lf.pc.wk = &v3;
  pcchk(&lf.pc,mi[MDIM],mi[MP],0);
  lf.pc.xtwx.sm = 1;
  lf.dp = dp; lf.mi = mi;
  lf.nd = *nd;
  for (i=0; i<*nd; i++) lf.deriv[i] = deriv[i]-1;
  fitted(&lf,&des,fit,ppwhat(what[0]),*cv,*st,restyp(what[1]));
}

void triterm(xev,h,ce,lo,hi,sca,nvc,mi,dp,nt,term)
double *xev, *h, *sca, *dp;
INT *ce, *lo, *hi, *nvc, *mi, *nt, *term;
{ INT d, i, t;
  vari vxev;
  lf_error = 0;
  lf.nv = lf.nvm = nvc[3];
  vxev.dpr = xev; vxev.n = mi[MDIM]*lf.nv; lf.xxev = &vxev;
  lf.sca = sca;
  lf.h = h;
  lf.ce = ce;
  lf.lo = lo;
  lf.hi = hi;
  lf.dp = dp; lf.mi = mi; *nt = 0;
  d = mi[MDIM];
  if (mi[MEV]==ETREE)
    atree_grow(NULL,&lf,lf.ce,nt,term,xev,&xev[d*((1<<d)-1)],lf.fl);
  else
    for (i=0; i<nvc[4]; i++)
      triang_grow(NULL,&lf,&lf.ce[i*lf.vc],nt,term);
}

void kdeb(x,mi,band,ind,h0,h1,meth,nmeth,ker)
double *x, *band, *h0, *h1;
INT *mi, *ind, *meth, *nmeth, *ker;
{ kdeselect(band,x,ind,*h0,*h1,meth,*nmeth,*ker,mi[MN]);
}
