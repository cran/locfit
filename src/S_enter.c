/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include "S.h"
#undef WARN
#undef ERROR

#include "local.h"
int deitype(char *);  /* in lfstr.c */


design des;
lfit lf;

int lf_error;

#ifdef RVERSION
typedef char * CALL_S_FUNC;
typedef void * CALL_S_ARGS;
#else
typedef void * CALL_S_FUNC;
typedef char * CALL_S_ARGS;
#endif
typedef long int CALL_S_LEN;
typedef long int CALL_S_NARG;
typedef char * CALL_S_MODE;
typedef long int CALL_S_NRES;
typedef char * CALL_S_VALS;

static CALL_S_FUNC bsfunc, bsf2;

void basis(x,t,f,dim,p)
double *x, *t, *f;
Sint dim, p;
{
  CALL_S_ARGS args[2];
  CALL_S_LEN  length[2];
  CALL_S_NARG nargs;
  CALL_S_MODE mode[2];
  CALL_S_NRES nres;
  CALL_S_VALS values[1];
  /*  double z0[1], z1[1], *vptr; */
  double *vptr;
  int i;

  args[0] = (CALL_S_ARGS)x;
  mode[0] = "double";
  length[0] = dim;

  args[1] = (CALL_S_ARGS)t;
  mode[1] = "double";
  length[1] = dim;

  nargs = 2;
  nres = 1;

  call_S(bsfunc,nargs,args,mode,length,(char **)NULL,nres,values);

  vptr = (double *)values[0];
  for (i=0; i<p; i++) f[i] = vptr[i];
}

void vbasis(x,t,n,d,ind,m,p,X)
double **x, *t, *X;
int n, d, m, p, *ind;
{
  CALL_S_ARGS args[MXDIM+3];
  CALL_S_LEN  length[MXDIM+3];
  CALL_S_NARG nargs;
  CALL_S_MODE mode[MXDIM+3];
  CALL_S_NRES nres;
  CALL_S_VALS values[1];
  double *vptr;
  int i;

  args[0] = (CALL_S_ARGS)(&d);
  mode[0] = "integer";
  length[0] = 1;

  args[1] = (CALL_S_ARGS)ind;
  mode[1] = "integer";
  length[1] = m;

  args[2] = (CALL_S_ARGS)t;
  mode[2] = "double";
  length[2] = d;

  for (i=0; i<d; i++)
  { args[3+i] = (CALL_S_ARGS)x[i];
    mode[3+i] = "double";
    length[3+i] = n;
  }

  nargs = d+3;
  nres = 1;

  call_S(bsf2,nargs,args,mode,length,0,nres,values);
  vptr = (double *)values[0];
  for (i=0; i<m*p; i++) X[i] = vptr[i];
}

void setevs(evs,mi,cut,mg,flim)
evstruc *evs;
int *mg;
Sint *mi;
double cut, *flim;
{ double *ll, *ur;
  int i, d;

  ev(evs) = mi[MEV];
  mk(evs) = mi[MK];
  d = mi[MDIM];

  if (flim != NULL)
  { ll = flim;
    ur = &flim[d];
    memcpy(evs->fl,ll,d*sizeof(double));
    memcpy(&evs->fl[d],ur,d*sizeof(double));
  }

  switch(ev(evs))
  { case ETREE:
    case EKDTR:
    case EKDCE:
    case EPHULL:
      cut(evs) = cut;
      return;
    case EGRID:
      for (i=0; i<d; i++)
        evs->mg[i] = mg[i];
      return;
    case ESPHR:
      for (i=0; i<2; i++) evs->mg[i] = mg[i];
      return;
    case EDATA:
    case ECROS:
    case EPRES:
    case EXBAR:
    case ENONE:
      return;
    default:
      printf("setevs: %2d not defined.\n",ev(evs));
  }
}

void setdata(lfd,x,y,c,w,b,n,d,sca,sty)
lfdata *lfd;
double *x, *y, *c, *w, *b, *sca;
Sint n, d, *sty;
{ int i;
  for (i=0; i<d; i++)
  { dvari(lfd,i) = &x[i*n];
    lfd->sca[i] = sca[i];
    lfd->sty[i] = sty[i];
  }
  lfd->y = y;
  lfd->w = w;
  lfd->b = b;
  lfd->c = c;
  lfd->n = n;
  lfd->d = d;
  lfd->ord = 0;
}

void setsmpar(sp,dp,mi)
smpar *sp;
double *dp;
Sint *mi;
{ nn(sp)  = dp[DALP];
  fixh(sp)= dp[DFXH];
  pen(sp) = dp[DADP];
  ker(sp) = mi[MKER];
  kt(sp)  = mi[MKT];
  acri(sp)= mi[MACRI];
  deg(sp) = mi[MDEG];
  deg0(sp) = mi[MDEG0];
  fam(sp)  = mi[MTG];
  link(sp) = mi[MLINK];
  ubas(sp) = mi[MUBAS];
  npar(sp) = mi[MP];
  lf.sp.vbasis = vbasis;
}

void slocfit(x,y,c,w,b,lim,mi,dp,str,sca,xev,wdes,wtre,wpc,nvc,
  iwk1, iwk2,lw,mg,L,kap,dv,nd,sty,bs)
double *x, *y, *c, *w, *b, *lim, *dp, *sca, *xev, *L, *kap, *wdes, *wtre, *wpc;
Sint *mi, *nvc, *iwk1, *iwk2, *lw, *mg, *dv, *nd, *sty;
char **str;
CALL_S_FUNC *bs;
{ Sint n, d, i;

  mi[MKER] = lfkernel(str[0]);
  mi[MTG]  = lffamily(str[1]);
  mi[MLINK]= lflink(str[2]);
  mi[MIT]  = deitype(str[3]);
  mi[MACRI]= lfacri(str[4]);
  mi[MKT]  = lfketype(str[5]);

  if (mi[MUBAS])
  { bsfunc = bs[0];
    bsf2 = bs[1];
  }
  lf_error = 0;
  n = mi[MN]; d = mi[MDIM];

  lfit_alloc(&lf);
  setdata(&lf.lfd,x,y,c,w,b,n,d,sca,sty);

  setsmpar(&lf.sp,dp,mi);
  setevs(&lf.evs,mi,dp[DCUT],mg,&lim[2*d]);

  lf_maxit = mi[MMXIT];
  lf_debug = mi[MDEB];
  de_mint  = mi[MMINT];
  de_itype = mi[MIT];
  de_renorm= mi[MREN];
  dc(&lf.fp) = mi[MDC];
  geth(&lf.fp)=mi[MGETH];

  des.wk = wdes;  des.lwk = lw[0];
  des.ind= iwk2; des.lind = lw[6];
  des.des_init_id = DES_INIT_ID;

  lf.fp.xev = xev;   lf.fp.lev = d*nvc[0];
  lf.fp.coef= wtre;  lf.fp.lwk = lw[1];
  lf.pc.wk  = wpc;   lf.pc.lwk = lw[3];

  lf.evs.iwk = iwk1; lf.evs.liw = lw[2];

  lf.fp.L = L; lf.fp.ll = lw[4];

  lf.fp.nvm = nvc[0];

  lf.dv.nd = *nd;
  for (i=0; i<lf.dv.nd; i++) lf.dv.deriv[i] = dv[i]-1;
  if (lf_error) return;

  memcpy(lf.lfd.xl ,lim,2*d*sizeof(double));

 if (mi[MGETH] >= 70)
    scb(&des,&lf);
  else switch(mi[MGETH])
  { case GSTD: /* the standard fit */
    case GAMF: /* for gam.lf, return residuals */
    case GAMP: /* for gam.lf prediction */
      if (mi[MDEG0]==mi[MDEG])
      { startlf(&des,&lf,procv,0);
        if (!lf_error) ressumm(&lf,&des);
      }
      else startlf(&des,&lf,procvvord,0);
      break;
    case GSMP:
      startlf(&des,&lf,procvraw,0);
      break;
    case GHAT:
      startlf(&des,&lf,procvhatm,(int)mi[MKER]!=WPARM);
      break;
    case GKAP:
      constants(&des,&lf);
      for(i=0; i<lw[5]; i++) kap[i] = lf.fp.kap[i];
      return;
    case GRBD:
      rband(&des,&lf,kap,lf.dv.deriv,lf.dv.nd);
      return;
    case GLSC:
      startlf(&des,&lf,procv,1);
      dens_lscv(&des,&lf);
      return;
  }

  nvc[0] = lf.fp.nvm;
  nvc[1] = lf.evs.ncm;
  nvc[3] = lf.fp.nv;
  nvc[4] = lf.evs.nce;
  mi[MEV]= ev(&lf.evs);
  mi[MP] = npar(&lf.sp);
  mi[MLINK] = link(&lf.sp);
  mi[MPC] = haspc(&lf.pc);
  dp[DLK] = llk(&lf.fp);
  dp[DT0] = df0(&lf.fp);
  dp[DT1] = df1(&lf.fp);
  dp[DRV] = rv(&lf.fp);
  dp[DRSC]= rsc(&lf.fp);
  memcpy(sca,lf.lfd.sca,d*sizeof(double));
  memcpy(&lim[2*d],lf.evs.fl,2*d*sizeof(double));
  for(i=0; i<lw[5]; i++) kap[i] = lf.fp.kap[i];
}

void recoef(xev,coef,cell,nvc,mi,dp)
double *xev, *coef, *dp;
Sint *cell, *nvc, *mi;
{ int d, vc=0;

  d = mi[MDIM];
  lf.fp.nv = lf.fp.nvm = nvc[3];
  lf.fp.xev = xev;
  lf.fp.d   = d;
  lf.fp.coef = coef; coef += lf.fp.nv*(d+1);
  lf.fp.nlx  = coef; coef += lf.fp.nv*(d+1);
  lf.fp.t0   = coef; coef += lf.fp.nv*(d+1);
  lf.fp.lik  = coef; coef += lf.fp.nv*3;
  lf.fp.h    = coef; coef += lf.fp.nv;
  lf.fp.deg  = coef; coef += lf.fp.nv;
  rv(&lf.fp) = dp[DRV];
  rsc(&lf.fp)= dp[DRSC];
  dc(&lf.fp) = mi[MDC];
  lf.fp.hasd = (mi[MDEG]>0) | dc(&lf.fp);

  switch(mi[MEV])
  { case ETREE:
    case EKDTR:
    case EGRID:
    case ESPHR: vc = 1<<d; break;
    case EPHULL: vc = d+1; break;
    case EXBAR:
    case ECROS:
    case EDATA:
    case EPRES:
    case ENONE: vc=0; break;
    default: ERROR(("spreplot: Invalid ev"));
  }

  lf.evs.ce = cell; cell += nvc[4]*vc;
  lf.evs.s  = cell; cell += MAX(nvc[3],nvc[4]);
  lf.evs.lo = cell; cell += MAX(nvc[3],nvc[4]);
  lf.evs.hi = cell; cell += MAX(nvc[3],nvc[4]);
}

void spreplot(xev,coef,sv,cell,x,res,se,wpc,sca,m,nvc,mi,dp,
  mg,dv,nd,sty,where,what,bs)
double *xev, *coef, *sv, *x, *res, *se, *wpc, *sca, *dp;
Sint *cell, *m, *nvc, *mi, *mg, *dv, *nd, *sty, *where;
char **what;
void **bs;
{ Sint i, p;
  double *xx[MXDIM];

  for (i=0; i<mi[MDIM]; i++)
  { lf.lfd.sty[i] = sty[i];
    lf.lfd.sca[i] = sca[i];
  }
  lf.lfd.d = mi[MDIM];

  setsmpar(&lf.sp,dp,mi);
  setevs(&lf.evs,mi,dp[DCUT],mg,NULL);

  if (mi[MUBAS]) bsfunc = bs[0];

  lf_error = 0; p = mi[MP];
  lf.evs.ncm = nvc[1]; lf.evs.nce = nvc[4];

  recoef(xev,coef,cell,nvc,mi,dp);
  lf.evs.sv = sv;

  lf.pc.wk = wpc;
  lf.pc.lwk = pc_reqd(mi[MDIM],p);
  pcchk(&lf.pc,(int)mi[MDIM],p,0);
  haspc(&lf.pc) = mi[MPC];
  lf.pc.xtwx.st = JAC_EIGD;

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
      for (i=0; i<mi[MDIM]; i++) dvari(&(lf.lfd),i) = xx[i] = &x[i**m];
      break;
    case 4: /* fit points, need nothing! */
      break;
    default:
      ERROR(("unknown where in spreplot"));
  }

  lf.dv.nd = *nd;
  for (i=0; i<*nd; i++) lf.dv.deriv[i] = dv[i]-1;

  if (lf_error) return;
  preplot(&lf,xx,res,se,what[1][0],m,*where,ppwhat(what[0]));
}

void sfitted(x,y,w,c,ba,fit,cv,st,xev,coef,sv,cell,wpc,sca,nvc,mi,dp,mg,dv,nd,sty,what,bs)
double *x, *y, *w, *c, *ba, *fit, *xev, *coef, *sv, *wpc, *sca, *dp;
Sint *cv, *st, *cell, *nvc, *mi, *mg, *dv, *nd, *sty;
char **what;
void **bs;
{ Sint i, n;

  recoef(xev,coef,cell,nvc,mi,dp);
  setsmpar(&lf.sp,dp,mi);
  setevs(&lf.evs,mi,dp[DCUT],mg,NULL);

  if (mi[MUBAS]) bsfunc = bs[0];
  n = mi[MN];
  lf_error = 0;
  lf.evs.ncm = nvc[1]; lf.evs.nce = nvc[4];

  setdata(&lf.lfd,x,y,c,w,ba,mi[MN],mi[MDIM],sca,sty);

  lf.evs.sv = sv;

  lf.pc.wk = wpc;
  lf.pc.lwk= pc_reqd(mi[MDIM],mi[MP],0);
  pcchk(&lf.pc,mi[MDIM],mi[MP],0);
  haspc(&lf.pc) = mi[MPC];
  lf.pc.xtwx.st = JAC_EIGD;

  lf.dv.nd = *nd;
  for (i=0; i<*nd; i++) lf.dv.deriv[i] = dv[i]-1;

  fitted(&lf,fit,ppwhat(what[0]),*cv,*st,restyp(what[1]));
}

void triterm(xev,h,ce,lo,hi,sca,nvc,mi,dp,nt,term,box)
double *xev, *h, *sca, *dp, *box;
Sint *ce, *lo, *hi, *nvc, *mi, *nt, *term;
{ int i, d, vc;
  Sint mg;

  lf_error = 0;
  d = mi[MDIM];

  lf.fp.xev = xev;
  lf.fp.h = h;
  lf.fp.d = d;
  lf.fp.nv = lf.fp.nvm = nvc[3];

  memcpy(lf.lfd.sca,sca,d*sizeof(double));
  setevs(&lf.evs,mi,dp[DCUT],&mg,NULL);

  lf.evs.ce = ce;
  lf.evs.lo = lo;
  lf.evs.hi = hi;
  *nt = 0;

  if (mi[MEV]==ETREE)
    atree_grow(NULL, &lf, lf.evs.ce, nt, term, box, &box[d]);
  else
  { vc = d+1;
    for (i=0; i<nvc[4]; i++)
      triang_grow(NULL,&lf,&lf.evs.ce[i*vc],nt,term);
  }
}

void guessnv(lw,evt,dp,mi,nvc,mg)
double *dp;
char **evt;
int *lw, *mi, *nvc, *mg;
{ int n, d, i, nvm, ncm, vc;
  smpar sp;
  evstruc evs;

  mi[MEV] = ev(&evs) = lfevstr(evt[0]);
  mi[MKT] = lfketype(evt[1]);
  mk(&evs) = mi[MK];
  d = mi[MDIM];
  n = mi[MN];

  switch(mi[MEV])
  { case ETREE:
      cut(&evs) = dp[DCUT];
      atree_guessnv(&evs,&nvm,&ncm,&vc,d,dp[DALP]);
      break;
    case EPHULL:
      nvm = ncm = mi[MK]*mi[MDIM];
      vc = mi[MDIM]+1;
      break;
    case EDATA:
    case ECROS:
      nvm = mi[MN];
      ncm = vc = 0;
      break;
    case EKDTR:
    case EKDCE:
      cut(&evs) = dp[DCUT];
      kdtre_guessnv(&evs,&nvm,&ncm,&vc,n,d,dp[DALP]);
      break;
    case EGRID:
      nvm = 1;
      for (i=0; i<d; i++) nvm *= mg[i];
      ncm = 0;
      vc = 1<<d;
      break;
    case EXBAR:
    case ENONE:
      nvm = 1;
      ncm = vc = 0;
      break;
    case EPRES:
      nvm = mg[0];
      ncm = vc = 0;
      break;
    default:
      ERROR(("guessnv: I don't know this evaluation structure."));
  }

  ubas(&sp)= mi[MUBAS];
  deg(&sp) = mi[MDEG];
  kt(&sp)  = mi[MKT];
  npar(&sp)= mi[MDEG]; /* for user basis */
  mi[MP] = calcp(&sp,d);
  lw[0] = des_reqd(n,(int)mi[MP]);
  lw[1] = lfit_reqd(d,nvm,ncm,(int)mi[MGETH]);
  lw[2] = lfit_reqi(nvm,ncm,vc);
  lw[6] = des_reqi(n,(int)mi[MP]);
  lw[3] = pc_reqd(d,(int)mi[MP]);
  lw[5] = 1;

  if (mi[MGETH] >= 70)
  { lw[4] = k0_reqd(d,n,0);
    if (lw[4]<2*nvm) lw[4] = 2*nvm;
    lw[5] = d+1;
  }
  else switch(mi[MGETH])
  { case GSTD:  lw[4] = 1; break;                  /* standard fit */
    case GSMP:  lw[4] = 1; break;                  /* simple fit   */
    case GHAT:  lw[4] = nvm*n; break;              /* hat matrix   */
    case GKAP:  lw[4] = k0_reqd(d,n,0);            /* kappa0       */
             lw[5] = 1+d;
             break;
    case GRBD:  lw[5] = 10;                        /* regband      */
    case GAMF:                                     /* gam.lf fit   */
    case GAMP:  lw[4] = 1; break;                  /* gam.lf pred  */
    case GLSC:  lw[4] = 2; break;                  /* lscv         */
    default:
      printf("sguessnv: invalid geth\n");
      lw[4] = 0;
  }

  nvc[0] = nvm;
  nvc[1] = ncm;
  nvc[2] = vc;
  nvc[3] = nvc[4] = 0;

}
