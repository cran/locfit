/*
 *   Copyright (c) 1998-1999 Lucent Technologies.
 *   See README file for details.
 */

/*
  routines for performing scb max. simulations
*/

#include "local.h"

#ifdef CVERSION
extern lfit lf;
extern design des;
extern vari *aru;
static INT corr, side;

INT nqmax(x,f,z,xl,mbe)
double *x, *f, *xl;
INT z, mbe; /* z=0, first call; z=1, subsequent calls  mbe, may be end */
{ INT im;
  double c, h1, h2, d1, d2;
/* printf("%6.4f %8.5f  %6.4f %8.5f  %6.4f %8.5f  %6.4f %8.5f\n",x[0],f[0],x[1],f[1],x[2],f[2],x[3],f[3]); */
  if (z==1)
  { if (f[3]<f[2])
    { x[3] = (x[3]+x[0])/2;
      return(0);
    }
    im = 2;
    if (f[3]>f[1]) { im = 1; f[2] = f[1]; x[2] = x[1]; }
    if (f[3]>f[0]) { im = 0; f[1] = f[0]; x[1] = x[0]; }
    f[im] = f[3]; x[im] = x[3];
  }
  else /* sort on f */
  { if (f[0]<f[1])
    { x[3] = x[0]; x[0] = x[1]; x[1] = x[3];
      f[3] = f[0]; f[0] = f[1]; f[1] = f[3];
    }
    if (f[1]<f[2])
    { x[3] = x[2]; x[2] = x[1]; x[1] = x[3];
      f[3] = f[2]; f[2] = f[1]; f[1] = f[3];
    }
    if (f[0]<f[1])
    { x[3] = x[0]; x[0] = x[1]; x[1] = x[3];
      f[3] = f[0]; f[0] = f[1]; f[1] = f[3];
    }
  }

  d2 = f[2]-f[0]; d1 = f[1]-f[0];
  h2 = x[2]-x[0]; h1 = x[1]-x[0];
  if (d2*h1-d1*h2==0) /* linear, return x0 */
  { x[3] = x[0]; f[3] = f[0];
    return(0);
  }
  x[3] = x[0]+(d2*h1*h1-d1*h2*h2)/(d2*h1-d1*h2);
  if (mbe)
  { if (x[3]<xl[0]) { x[3] = xl[0]; return(1); }
    if (x[3]>xl[1]) { x[3] = xl[1]; return(1); }
    c = ((x[2]-x[1])*f[0]-h2*f[1]+h1*f[2])/((x[2]-x[1])*h2*h1);
    if (c>0.0) return(1); /* local minimum */
  }

  while ((x[3]<xl[0]) | (x[3]>xl[1])) x[3] = (x[3]+x[0])/2;
  if (fabs(x[3]-x[0])*10<fabs(h1)) x[3] = (x[3]+x[1])/2;
  if (fabs(x[3]-x[1])*10<fabs(h1)) x[3] = (x[3]+x[0])/2;
  if (fabs(x[3]-x[2])*10<fabs(h2)) x[3] = (x[3]+x[0])/2;
  return(0);
}

double b2(th,tg,w)
double th, w;
INT tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(w);
    case TPOIS: return(w*exp(th));
    case TLOGT:
      y = expit(th);
      return(w*y*(1-y));
  }
  ERROR(("b2: invalid family %d",tg));
  return(0.0);
}

double b3(th,tg,w)
double th, w;
INT tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(0.0);
    case TPOIS: return(w*exp(th));
    case TLOGT:
      y = expit(th);
      return(w*y*(1-y)*(1-2*y));
  }
  ERROR(("b3: invalid family %d",tg));
  return(0.0);
}

double b4(th,tg,w)
double th, w;
INT tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(0.0);
    case TPOIS: return(w*exp(th));
    case TLOGT:
      y = expit(th); y = y*(1-y);
      return(w*y*(1-6*y));
  }
  ERROR(("b4: invalid family %d",tg));
  return(0.0);
}

double cumulant(lf,des,w)
lfit *lf;
design *des;
double w;
{ double c1, c2, c3, c4, c5, c6, c7, c8, c9;
  double b2i, b3i, b3j, b4i, p2, s[10], *ui, *uj;
  double ss, si, sj, uii, uij, ujj, k1, k2, k4;
  INT i, j, *mi;
  c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = c9 = 0;
  k1 = 0;
  mi = lf->mi;

  locfit(lf,des,0.0,0);
  wdiag(lf,des,s,0,2,0);
  ss = innerprod(s,s,mi[MP]);

  for (i=0; i<mi[MN]; i++)
  { b2i = b2(des->th[i],mi[MTG],prwt(lf,i));
    b3i = b3(des->th[i],mi[MTG],prwt(lf,i));
    b4i = b4(des->th[i],mi[MTG],prwt(lf,i));
    ui = (double *)viptr(lf->L,i*mi[MP]);
    si = innerprod(s,ui,mi[MP]);
    uii= innerprod(ui,ui,mi[MP]);
    if (lf_error) return(0.0);

    c2 += b4i*si*si*uii;
    c6 += b4i*si*si*si*si;
    c7 += b3i*si*uii;
    c8 += b3i*si*si*si;
    c9 += b2i*b2i*si*si*si*si;
    k1 += b3i*si*(si*si/ss-uii);
/* printf("b3i %8.5f si %8.5f ss %8.5f uii %8.5f k1 %8.5f\n",b3i,si,ss,uii,k1);
*/

    /* i=j components */
    c1 += b3i*b3i*si*si*uii*uii;
    c3 += b3i*b3i*si*si*si*si*uii;
    c4 += b3i*b3i*si*si*uii*uii;

    for (j=i+1; j<mi[MN]; j++)
    { b3j = b3(des->th[j],mi[MTG],prwt(lf,j));
      uj = (double *)viptr(lf->L,j*mi[MP]);
      sj = innerprod(s,uj,mi[MP]);
      uij= innerprod(ui,uj,mi[MP]);
      ujj= innerprod(uj,uj,mi[MP]);

      c1 += 2*b3i*b3j*si*sj*uij*uij;
      c3 += 2*b3i*b3j*si*si*sj*sj*uij;
      c4 += b3i*b3j*uij*(si*si*ujj+sj*sj*uii);
      if (lf_error) return(0.0);
    }
  }
  c5 = c1;
  c7 = c7*c8;
  c8 = c8*c8;

  c1 /= ss; c2 /= ss; c3 /= ss*ss; c4 /= ss;
  c5 /= ss; c6 /= ss*ss; c7 /= ss*ss; c8 /= ss*ss*ss;
  c9 /= ss*ss;
/* printf("%7.5f %7.5f  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  ",x[0],w,c1,c2,c3,c4,c5,c6,c7,c8,c9); */

  k1 = k1/(2*sqrt(ss));
  k2 = 1+c1/2-c2/2-3*c3+c4/2+c6-c7/2+1.75*c8;
  k4 = -9*c3+3*c6+6*c8+3*c9;
/* printf("kappa: %8.5f %8.5f %8.5f %8.5f\n",k1,k2,sqrt(c8),k4); */
  if (k2<0) WARN(("cumul: k2<0 in cumul"));

  w = (w-k1)/sqrt(k2);
/* printf("%7.5f %7.5f %7.5f  %7.5f  ",k1,k2,k4,w); */
  p2 = -w*(k4*(w*w-3)/24 + c8*((w*w-10)*w*w+15)/72);
/* printf("%7.5f  %7.5f\n",p2,fabs(w)+p2); */
  switch(side)
  { case -1: return(-w*sqrt(k2));
    case  0: return(fabs(w)+p2);
    case  1: return(w*sqrt(k2));
  }
  return(0);
}



void procvscb(des,lf,v)
design *des;
lfit *lf;
INT v;
{ double mean, w, thhat, sd;
  des->xev = evpt(lf,v);
  if ((lf->mi[MKER]==WPARM) && (hasparcomp(lf)))
  { thhat = addparcomp(lf,des->xev,PCOEF);
    sd = addparcomp(lf,des->xev,PNLX);
  }
  else
  { procv(des,lf,v);
    thhat = lf->coef[v];
    sd = lf->nlx[v];
  }
  mean = dareval(aru,0,des->xev);
  switch(lf->mi[MLINK])
  { case LIDENT: break;
    case LLOG:   mean = log(mean); break;
    case LLOGIT: mean = logit(mean); break;
    default: ERROR(("procvscb: invalid link %d",lf->mi[MLINK]));
      return;
  }
  w = (thhat-mean)/sd;
  if (corr)
  { lf->coef[v] = cumulant(lf,des,w);
    return;
  }
  switch(side)
  { case -1:lf->coef[v] = -w; return;
    case 0: lf->coef[v] = fabs(w); return;
    case 1: lf->coef[v] = w; return;
  }
}

void scbmax(lf,des,co,si)
lfit *lf;
design *des;
INT co, si;
{ double max, xmx, x[4], f[4], kap[3];
  INT c, i, im, nv, v, dv[MXDIM], mbe;
  if ((co) & (ident==0))
  { WARN(("Correction doesn't work with locfit; ignoring."));
    co = 0;
  }
  corr = co;
  side = si;
  if (corr)
  { v = lf->mi[MEV]; lf->mi[MEV] = EDATA;
    i = calcp(lf->mi[MDEG],lf->mi[MDIM],lf->mi[MKT]);
    lf->L = checkvarlen(lf->L,lf->mi[MN]*i,"_hatmat",VDOUBLE);
    startlf(des,lf,procvhatm,lf->mi[MKER]!=WPARM);
    lf->mi[MEV] = v;
  }
  startlf(des,lf,procvscb,lf->mi[MKER]!=WPARM);
  if (lf_error) return;
  max = 0.0; im = 0;
  nv = lf->nv;
  for (i=0; i<nv; i++)
    if (lf->coef[i]>max) { max = lf->coef[i]; im = i; }
  if (nv>=4)
  { mbe = 0;
    xmx = evptx(lf,im,0);
    if (im==0) { mbe = 1; im = 1; }
    if (im==nv-1) { mbe = 1; im = nv-2; }
    if (im<=1) v = 3; else v = 0;
    x[0] = evptx(lf,im,0);   f[0] = lf->coef[im];
    x[1] = evptx(lf,im-1,0); f[1] = lf->coef[im-1];
    x[2] = evptx(lf,im+1,0); f[2] = lf->coef[im+1];
    i = 0;
    do
    { c = nqmax(x,f,(i>0),lf->fl,mbe);
      i = 1;
      if (!c)
      { evptx(lf,v,0) = x[3];
        procvscb(des,lf,v);
        f[3] = lf->coef[v];
        if (f[3]>max) { max = f[3]; xmx = x[3]; }
      }
    } while ((!c) && (f[0]-f[2]>0.000001));
  }

  /* now, compute kappa0 and tail probability */

  constants(des,lf,kap,dv,0);
  printf("xmx: %10.6f  max: %10.6f  k0 %10.6f %10.6f  pr %10.6f\n",xmx,max,kap[0],kap[1],tailp(max,kap,1+(side==0),lf->mi[MDIM],0.0));
}

void cscbmax(v)
vari *v;
{ INT i, corr, side;
  corr = side = 0;
  fitoptions(&lf,v,0);
 
  i = getarg(v,"mean",1);
  if (i==0) { ERROR(("cscbmax: no mean function given")); }
  else
  { aru = arbuild(argval(v,i),0,strlen(argval(v,i))-1,NULL,0,1);
    setvarname(aru,"_aru");
  }

  i = getarg(v,"corr",1);
  if (i>0) corr = getlogic(v,i);
  if (lf_error) return;

  i = getarg(v,"side",1);
  if (i>0) sscanf(argval(v,i),"%d",&side);
  if (lf_error) return;

  scbmax(&lf,&des,corr,side);
}
#endif
