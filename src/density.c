/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

static double u[MXDIM], ilim[2*MXDIM], *ff, tmax;

INT multint(), prodint(), gausint(), mlinint(), hazint(), harint();

void prresp(coef,resp,p)
double *coef, *resp;
INT p;
{ INT i, j;
  printf("Coefficients:\n");
  for (i=0; i<p; i++) printf("%8.5f ",coef[i]);
  printf("\n");
  printf("Response matrix:\n");
  for (i=0; i<p; i++)
  { for (j=0; j<p; j++) printf("%9.6f, ",resp[i+j*p]);
    printf("\n");
  }
}

INT multint(t,resp1,resp2,lf,cf,h)
lfit *lf;
double *t, *resp1, *resp2, *cf, h;
{ INT d, p, i, j, k, m, m1, w, z, z1;
  double th, wt, dj[MXDIM];
  d = lf->mi[MDIM]; p = lf->mi[MP];
  setzero(resp1,p*p);
  m = 1; m1 = lf->mi[MMINT]+1;
  for (i=0; i<d; i++)
  { m *= m1;
    dj[i] = (ilim[i+d]-ilim[i])/lf->mi[MMINT];
  }
  for (i=0; i<m; i++)
  { z = i; w = 1;
    for (j=d-1; j>=0; j--)
    { z1 = z%m1;
      u[j] = t[j]+ilim[j]+dj[j]*z1;
      w *= (4-2*(z1%2==0)-(z1==0)-(z1==lf->mi[MMINT]));
      z /= m1;
    }
    wt = w*weight(lf,u,t,h,0,0.0);
    if (wt>0)
    { fitfun(lf,u,t,ff,NULL,0);
      th = innerprod(ff,cf,p);
      switch(lf->mi[MLINK])
      { case LLOG:
          addouter(resp1,ff,ff,p,wt*exp(th));
          break;
        case LIDENT:
          addouter(resp1,ff,ff,p,wt);
          break;
        case LSQRT:
          if (th<0) return(LF_BADP);
          addouter(resp1,ff,ff,p,2*wt*th);
          break;
        default:
          ERROR(("multint: Invalid link"));
          return(LF_LNK);
      }
    }
  }
  wt = 1;
  for (j=0; j<d; j++) wt *= dj[j]/3;
  for (j=0; j<p; j++)
    for (k=j; k<p; k++)
      resp1[p*k+j] = resp1[p*j+k] = resp1[p*j+k]*wt;
  return(LF_OK);
}

INT mlinint(t,resp1,resp2,lf,cf,h)
lfit *lf;
double *t, *resp1, *resp2, *cf, h;
{ double P[MXDIM*MXDIM], C[MXDIM*MXDIM], pu[MXDIM];
  double hd, nb, c, wt, wu, g[4], w0, w1, v, *sca;
  INT d, p, i, j, jmax, k, l, w, z, z1, jj[2];
  d = lf->mi[MDIM]; p = lf->mi[MP]; sca = lf->sca;
  hd = 1; for (i=0; i<d; i++) hd *= h*sca[i];
  if (lf->mi[MLINK]==LIDENT)
  { setzero(resp1,p*p);
    resp1[0] = wint(d,NULL,0,lf->mi[MKER])*hd;
    if (lf->mi[MDEG]==0) return(LF_OK);
    jj[0] = 2; w0 = wint(d,jj,1,lf->mi[MKER])*hd*h*h;
    for (i=0; i<d; i++) resp1[(i+1)*p+i+1] = w0*sca[i]*sca[i];
    if (lf->mi[MDEG]==1) return(LF_OK);
    for (i=0; i<d; i++)
    { j = p-(d-i)*(d-i+1)/2;
      resp1[j] = resp1[p*j] = w0*sca[i]*sca[i]/2;
    }
    if (d>1) { jj[1] = 2; w0 = wint(d,jj,2,lf->mi[MKER])*hd*h*h*h*h; }
    jj[0] = 4; w1 = wint(d,jj,1,lf->mi[MKER])*hd*h*h*h*h/4;
    z = d+1;
    for (i=0; i<d; i++)
    { k = p-(d-i)*(d-i+1)/2;
      for (j=i; j<d; j++)
      { l = p-(d-j)*(d-j+1)/2;
        if (i==j) resp1[z*p+z] = w1*SQR(sca[i])*SQR(sca[i]);
        else
        { resp1[z*p+z] = w0*SQR(sca[i])*SQR(sca[j]);
          resp1[k*p+l] = resp1[k+p*l] = w0/4*SQR(sca[i])*SQR(sca[j]);
        }
        z++;
    } }
    return(LF_OK);
  }
  switch(lf->mi[MDEG])
  { case 0:
      resp1[0] = exp(cf[0])*wint(d,NULL,0,lf->mi[MKER]);
      for (i=0; i<d; i++) resp1[0] *= h*sca[i];
      return(LF_OK);
    case 1:
      nb = 0.0;
      for (i=1; i<=d; i++)
      { v = h*cf[i]*sca[i-1];
        nb += v*v;
      }
      if (lf->mi[MKER]==WGAUS)
      { w0 = 1/(GFACT*GFACT);
        g[0] = exp(cf[0]+w0*nb/2+d*log(S2PI/2.5));
        g[1] = g[3] = g[0]*w0;
        g[2] = g[0]*w0*w0;
      }
      else
      { wt = wu = exp(cf[0]);
        w0 = wint(d,NULL,0,lf->mi[MKER]); g[0] = wt*w0;
        g[1] = g[2] = g[3] = 0.0;
        j = 0; jmax = (d+2)*lf->mi[MMINT];
        while ((j<jmax) && (wt*w0/g[0]>1.0e-8))
        { j++;
          jj[0] = 2*j; w0 = wint(d,jj,1,lf->mi[MKER]);
          if (d==1) g[3] += wt * w0;
          else
          { jj[0] = 2; jj[1] = 2*j-2; w1 = wint(d,jj,2,lf->mi[MKER]);
            g[3] += wt*w1;
            g[2] += wu*(w0-w1);
          }
          wt /= (2*j-1.0); g[1] += wt*w0;
          wt *= nb/(2*j); g[0] += wt*w0;
          wu /= (2*j-1.0)*(2*j);
          if (j>1) wu *= nb;
        }
        if (j==jmax) WARN(("mlinint: series not converged"));
      }
      g[0] *= hd; g[1] *= hd;
      g[2] *= hd; g[3] *= hd;
      resp1[0] = g[0];
      for (i=1; i<=d; i++)
      { resp1[i] = resp1[(d+1)*i] = cf[i]*SQR(h*sca[i-1])*g[1];
        for (j=1; j<=d; j++)
        { resp1[(d+1)*i+j] = (i==j) ? g[3]*SQR(h*sca[i-1]) : 0;
          resp1[(d+1)*i+j] += g[2]*SQR(h*h*sca[i-1]*sca[j-1])*cf[i]*cf[j];
        }
      }
      return(LF_OK);
    case 2:
      if (lf->mi[MKER]==WGAUS) return(gausint(t,resp1,resp2,cf,h,lf->mi,sca));
      k = d+1;
      for (i=0; i<d; i++)
      { C[i*(d+1)] = cf[k++]*sca[i]*sca[i];
        for (j=i+1; j<d; j++) C[i*d+j] = C[j*d+i] = cf[k++]*sca[i]*sca[j];
      }
      eigen(C,P,d,20);
      k = 1;
      for (i=0; i<d; i++) k *= lf->mi[MMINT]+1;
      setzero(resp1,p*p);
      for (i=0; i<k; i++)
      { z = i; w = 1;
        for (j=d-1; j>=0; j--)
        { z1 = z%(lf->mi[MMINT]+1);
          pu[j] = h*sca[j]*(-1+2.0*z1/lf->mi[MMINT]);
          w *= (4-2*(z1%2==0)-(z1==0)-(z1==lf->mi[MMINT]));
          z /= (lf->mi[MMINT]+1);
        }
        wt = w*weight(lf,pu,NULL,h,0,0.0);
        if (wt>0)
        { for (j=0; j<d; j++)
          { u[j] = 0.0;
            for (l=0; l<d; l++) u[j] += sca[j]*P[j*d+l]*pu[l]/sca[l];
          }
          c = cf[0];
          for (j=0; j<d; j++)
            c += cf[j+1]*u[j]+SQR(pu[j]/sca[j])*C[j*(d+1)]/2;
          fitfun(lf,u,NULL,ff,NULL,0);
          addouter(resp1,ff,ff,p,wt*exp(c));
      } }
      for (i=0; i<p*p; i++)
        for (j=0; j<d; j++)
          resp1[i] *= 2.0*h*sca[j]/(3*lf->mi[MMINT]);
      return(LF_OK);
    default: ERROR(("mlinint: deg<=2 only"));
  }
  return(LF_ERR);
}

INT prodint(t,resp1,resp2,lf,coef,h)
lfit *lf;
double *t, *resp1, *resp2, *coef, h;
{ INT d, deg, p, i, j, j0, j1, k0, k1, f0, m, st;
  double ea, cf[MXDEG+1], hj, hs;
  d = lf->mi[MDIM]; deg = lf->mi[MDEG];
  p = lf->mi[MP];
  for (i=0; i<p*p; i++) resp1[i] = 1;
  cf[0] = coef[0];
  for (i=0; i<d; i++)
  { hj = 1; hs = h*lf->sca[i];
    for (j=0; j<deg; j++)
    { hj *= hs;
      cf[j+1] = hj*coef[j*d+i+1];
    }
    st = onedint(cf,lf->mi,ilim[i]/hs,ilim[i+d]/hs,resp2);
    if (st != LF_OK) printf("onedint return %d\n",st);
    hj = 1;
    for (j=0; j<=2*deg; j++)
    { hj *= hs;
      resp2[j] *= hj;
    }
    resp1[0] *= resp2[0];
    for (k0=1; k0<=deg; k0++)
    { m = 1+(k0-1)*d;
      for (j0=0; j0<d; j0++) resp1[m +j0] *= resp2[k0*(i==j0)];
      for (k1=1; k1<=k0; k1++) /* update (k0,k1) block */
      { m = p*(1+(k1-1)*d)+1+(k0-1)*d;
        for (j0=0; j0<d; j0++)
          for (j1=0; j1<d; j1++)
            resp1[m+p*j1+j0] *= resp2[k0*(i==j0)+k1*(i==j1)];
      }
  } }
  ea = exp(-(d-1)*coef[0]);
/* factorial divisions */

  ff[0] = f0 = 1.0;
  for (k0=0; k0<deg; k0++)
  { f0 *= k0+1;
    for (j0=1; j0<=d; j0++) ff[k0*d+j0] = f0;
  }
  for (k0=0; k0<p; k0++)
    for (k1=0; k1<=k0; k1++)
    { resp1[k1*p+k0] *= ea/(ff[k0]*ff[k1]);
      resp1[k0*p+k1] = resp1[k1*p+k0];
    }

  return(LF_OK);
}

INT gausint(t,resp,C,cf,h,mi,sca)
double *t, *resp, *C, *cf, h, *sca;
INT *mi;
{ double nb, det, z, *P;
  INT d, p, i, j, k, l, m1, m2, f;
  d = mi[MDIM]; p = mi[MP];
  m1 = d+1; nb = 0;
  P = &C[d*d];
  resp[0] = 1;
  for (i=0; i<d; i++)
  { C[i*d+i] = SQR(GFACT/(h*sca[i]))-cf[m1++];
    for (j=i+1; j<d; j++) C[i*d+j] = C[j*d+i] = -cf[m1++];
  }
  eigen(C,P,d,20);
  det = 1;
  for (i=1; i<=d; i++)
  { det *= C[(i-1)*(d+1)];
    if (det <= 0) return(LF_BADP);
    resp[i] = cf[i];
    for (j=1; j<=d; j++) resp[j+i*p] = 0;
    resp[i+i*p] = 1;
    svdsolve(&resp[i*p+1],u,P,C,P,d,0.0);
  }
  svdsolve(&resp[1],u,P,C,P,d,0.0);
  det = sqrt(det);
  for (i=1; i<=d; i++)
  { nb += cf[i]*resp[i];
    resp[i*p] = resp[i];
    for (j=1; j<=d; j++)
      resp[i+p*j] += resp[i]*resp[j];
  }
  m1 = d;
  for (i=1; i<=d; i++)
    for (j=i; j<=d; j++)
    { m1++; f = 1+(i==j);
      resp[m1] = resp[m1*p] = resp[i*p+j]/f;
      m2 = d;
      for (k=1; k<=d; k++)
      { resp[m1+k*p] = resp[k+m1*p] =
        ( resp[i]*resp[j*p+k] + resp[j]*resp[i*p+k]
        + resp[k]*resp[i*p+j] - 2*resp[i]*resp[j]*resp[k] )/f;
        for (l=k; l<=d; l++)
        { m2++; f = (1+(i==j))*(1+(k==l));
          resp[m1+m2*p] = resp[m2+m1*p] = ( resp[i+j*p]*resp[k+l*p]
            + resp[i+k*p]*resp[j+l*p] + resp[i+l*p]*resp[j+k*p]
            - 2*resp[i]*resp[j]*resp[k]*resp[l] )/f;
    } } }
  z = exp(d*0.918938533+cf[0]+nb/2)/det;
  multmatscal(resp,z,p*p);
  return(LF_OK);
}

INT hrao(lf,dfx,cf,h,r1)
lfit *lf;
double *dfx, *cf, h, *r1;
{ double s, t0, t1, wt, th;
  INT j, p, *mi;
  mi = lf->mi;
  s = 0; p = mi[MP];
  switch (mi[MKT])
  { case KPROD:
      for (j=1; j<mi[MDIM]; j++) s = MAX(s,fabs(dfx[j]/(h*lf->sca[j])));
      break;
    case KSPH:
      for (j=1; j<mi[MDIM]; j++) s += SQR(dfx[j]/(h*lf->sca[j]));
      break;
    default:
      ERROR(("harint not available for kt=%d",mi[MKT]));
      return(LF_ERR);
  }
  if (s>1) return(LF_OK);
  if (mi[MKT]==KPROD) t1 = h*lf->sca[0];
    else t1 = sqrt(1-s)*h*lf->sca[0];
  t0 = -t1;
  if (t0<ilim[0]) t0 = ilim[0];
  if (t1>ilim[mi[MDIM]]) t1 = ilim[mi[MDIM]];
  if (t1>dfx[0]) t1 = dfx[0];
  if (t1<t0) return(LF_OK);
  setzero(r1,p*p);
  for (j=0; j<=mi[MMINT]; j++)
  { dfx[0] = t0+(t1-t0)*j/mi[MMINT];
    wt = weight(lf,dfx,NULL,h,0,0.0);
    fitfun(lf,dfx,NULL,ff,NULL,0);
    th = innerprod(cf,ff,p);
    wt *= 2+2*(j&1)-(j==0)-(j==mi[MMINT]);
    addouter(r1,ff,ff,p,wt*exp(th));
  }
  multmatscal(r1,(t1-t0)/(3*mi[MMINT]),p*p);
  return(LF_OK);
}

INT harint(t,resp,r1,lf,cf,h)
lfit *lf;
double *t, *resp, *r1, *cf, h;
{ INT i, j, p;
  double dfx[MXDIM];
  p = lf->mi[MP];
  setzero(resp,p*p);
  for (i=0; i<=lf->mi[MN]; i++)
  { if (i==lf->mi[MN])
    { dfx[0] = tmax-t[0];
      for (j=1; j<lf->mi[MDIM]; j++) dfx[j] = 0.0;
    }
    else
      for (j=0; j<lf->mi[MDIM]; j++) dfx[j] = datum(lf,j,i)-t[j];
    if (hrao(lf,dfx,cf,h,r1))
      for (j=0; j<p*p; j++) resp[j] += r1[j];
  }
  return(LF_OK);
}

void hadd(resp,r,p,d,deg,dfx)
double *resp, *r, *dfx;
INT p, d, deg;
{ INT i, j;
  resp[0] += r[0];
  if (deg==0) return;
  resp[1] += r[1];     /* t */
  resp[p+1] += r[2];   /* t*t */
  for (i=1; i<d; i++)
  { resp[i+1] += r[0]*dfx[i];                     /* x */
    resp[p+i+1] += r[1]*dfx[i];                   /* x*t */
    for (j=i; j<d; j++)
      resp[(i+1)*p+j+1] += r[0]*dfx[i]*dfx[j];    /* x*x */
  }
  if (deg==1) return;
  resp[d+1] += r[2]/2;         /* t^2/2 */
  resp[p+d+1] += r[3]/2;       /* t^2*t/2 */
  resp[(d+1)*p+d+1] += r[4]/4; /* t^2*t^2/4 */
  for (i=1; i<d; i++)
  { resp[d+1+i] += dfx[i]*dfx[i]*r[0]/2;    /* x^2/2 */
    resp[p+d+1+i] += dfx[i]*dfx[i]*r[1]/2;  /* x^2*t/2 */
    resp[p*(i+1)+d+1] += dfx[i]*r[2]/2;     /* x*t^2/2 */
    resp[p*(d+1)+d+1+i] += r[2]*dfx[i]*dfx[i]/4; /* x^2*t^2/4 */
    for (j=1; j<d; j++)
      resp[p*(i+1)+d+1+j] += dfx[i]*dfx[j]*dfx[j]*r[0]/2; /* x^2*x/2 */
    for (j=i; j<d; j++)
      resp[(i+d+1)*p+j+d+1] += r[0]*dfx[i]*dfx[i]*dfx[j]*dfx[j]/4;
              /* x^2*x^2/4 */
  }
  if (deg==2) return;
  WARN(("hazint for deg<=2 only"));
}

INT hazint(t,resp,x,lf,cf,h)
lfit *lf;
double *t, *resp, *x, *cf, h;
{ INT d, p, i, j, k, st;
  double r[2*MXDEG+1], rr[2*MXDEG+1], dfx[MXDIM],
         hj, hs, ncf[MXDEG], wx, ef, lt;
  p = lf->mi[MP]; d = lf->mi[MDIM];
  setzero(resp,p*p);
  hj = hs = h*lf->sca[0];
  ncf[0] = cf[0];
  for (i=1; i<=lf->mi[MDEG]; i++)
  { ncf[i] = hj*cf[(i-1)*d+1]; hj *= hs;
  }
  k = 0; lt = -1.0;
  for (i=0; i<=lf->mi[MN]; i++)
  { if (i==lf->mi[MN])
    { dfx[0] = tmax-t[0];
      for (j=1; j<d; j++) dfx[j] = 0.0;
    }
    else
      for (j=0; j<d; j++) dfx[j] = datum(lf,j,i)-t[j];
    if ((d==1) & (dfx[0]>=hs)) k++;
    else if (dfx[0]>ilim[0])
    { wx = 1;
      ef = 0.0;
      for (j=1; j<d; j++)
      { wx *= W(dfx[j]/(h*lf->sca[j]),lf->mi[MKER]);
        if (lf->mi[MDEG]>=1) ef += cf[j+1]*dfx[j];
        if (lf->mi[MDEG]>=2) ef += cf[j+d+1]*dfx[j]*dfx[j]/2;
      }
      ef = exp(ef);
      if (wx>0)
      { if (dfx[0]>hs) dfx[0] = hs;
        if ((lt==-1.0) || (dfx[0]!=lt))
        { st = onedint(ncf,lf->mi,ilim[0]/hs,dfx[0]/hs,r);
          if (st>0) return(st);
          lt = dfx[0];
        }
        hj = 1;
        for (j=0; j<=2*lf->mi[MDEG]; j++) { hj *= hs; rr[j] = r[j]*wx*hj*ef; }
        hadd(resp,rr,p,d,lf->mi[MDEG],dfx);
      }
    }
  }
  if (k>0)
  { st = onedint(ncf,lf->mi,ilim[0]/hs,1.0,r);
    if (st>0) return(st);
    lt = 1.0;
    hj = 1;
    for (i=0; i<=2*lf->mi[MDEG]; i++) { hj *= hs; r[i] *= k*hj; }
    for (i=1; i<d; i++) dfx[i] = 0.0;
    hadd(resp,r,p,d,lf->mi[MDEG],dfx);
  }
  for (i=0; i<p; i++)
    for (j=0; j<i; j++)
      resp[i*p+j] = resp[j*p+i];
  return(LF_OK);
}

INT likeden(lf,des,nop)
lfit *lf;
design *des;
INT nop;
{ double lk, r, *A;
  INT i, j, p, st;
  p = des->p; A = des->xtwx.Z;
  if ((lf->mi[MLINK]==LIDENT) && (des->cf[0] != 0.0)) return(LF_OK);
  st = (des->itype)(des->xev,A,des->xtwx.Q,lf,des->cf,des->h);
  des->xtwx.p = p;
  if (lf->mi[MDEB]>2) prresp(des->cf,A,p);
  if (st!=LF_OK) return(st);
  switch(lf->mi[MLINK])
  { case LLOG:
      r = des->ss[0]/A[0];
      des->cf[0] += log(r);
      multmatscal(A,r,p*p);
      A[0] = des->ss[0];
      lk = -A[0];
      for (i=0; i<p; i++)
      { lk += des->cf[i]*des->ss[i];
        des->res[i] = des->ss[i]-A[i];
      }
      break;
    case LIDENT:
      lk = 0.0;
      for (i=0; i<p; i++)
      { des->res[i] = des->ss[i];
        for (j=0; j<p; j++)
          des->res[i] -= A[i*p+j]*des->cf[j];
      }
      break;
    case LSQRT:
      lk = 0.0;
      for (i=0; i<p; i++)
      { des->res[i] = des->ss[i];
        for (j=0; j<p; j++)
          des->res[i] -= A[i*p+j]*des->cf[j]/2;
      }
      break;
  }
  xtwxdec(&des->xtwx,des->xtwx.sm,des->res,nop);
  des->llk = lk;

  return(LF_OK);
}

INT inre(x,bound,d)
double *x, *bound;
INT d;
{ INT i, z;
  z = 1;
  for (i=0; i<d; i++)
    if (bound[i]<bound[i+d])
      z &= (x[i]>=bound[i]) & (x[i]<=bound[i+d]);
  return(z);
}

INT setintlimits(lf, x, h, ang, lset)
lfit *lf;
INT *ang, *lset;
double *x, h;
{ INT d, i;
  d = lf->mi[MDIM];
  *ang = *lset = 0;
  for (i=0; i<d; i++)
  { if (lf->sty[i]==STANGL)
    { ilim[i+d] = ((h<2) ? 2*asin(h/2) : PI)*lf->sca[i];
      ilim[i] = -ilim[i+d];
      *ang = 1;
    }
    else
    { ilim[i+d] = h*lf->sca[i];
      ilim[i] = -ilim[i+d];

      if (lf->sty[i]==STLEFT) { ilim[i+d] = 0; *lset = 1; }
      if (lf->sty[i]==STRIGH) { ilim[i] = 0;   *lset = 1; }
      if (lf->sty[i]==STUSER) { *ang = 1; } /* crude! to force multint */

      if (lf->xl[i]<lf->xl[i+d]) /* user limits for this variable */
      { if (lf->xl[i]-x[i]> ilim[i])
        { ilim[i] = lf->xl[i]-x[i]; *lset=1; }
        if (lf->xl[i+d]-x[i]< ilim[i+d])
        { ilim[i+d] = lf->xl[i+d]-x[i]; *lset=1; }
      }
    }
    if (ilim[i]==ilim[i+d]) return(LF_DEMP); /* empty integration */
  }
  return(LF_OK);
}

INT selectintmeth(mi,lset,ang)
INT *mi, lset, ang;
{
  if (mi[MIT]==IDEFA) /* select the default method */
  { if (mi[MTG]==THAZ)
    { if (ang) return(IDEFA);
      return( ((mi[MDIM]==1)|(mi[MKT]==KPROD)) ? IHAZD : IHARD );
    }

    if ((ang) | (mi[MLINK]==LSQRT)) return(IMULT);

    if (iscompact(mi[MKER]))
    { if (mi[MKT]==KPROD) return(IPROD);
      if (lset)
        return( (mi[MDIM]==1) ? IPROD : IMULT );
      if (mi[MDEG]<=1) return(IMLIN);
      if (mi[MDIM]==1) return(IPROD);
      return(IMULT);
    }

    if (mi[MKER]==WGAUS)
    { if (lset) WARN(("Integration for Gaussian weights ignores limits"));
      if ((mi[MDIM]==1)|(mi[MKT]==KPROD)) return(IPROD);
      return(IMLIN);
    }

    return(IDEFA);
  }

  /* user provided an integration method, check it is valid */

  if (mi[MTG]==THAZ)
  { if (ang) return(INVLD);
    if (!iscompact(mi[MKER])) return(INVLD);
    switch(mi[MIT])
    { case IHAZD: return( ((mi[MDIM]==1) | (mi[MKT]==KPROD)) ? IHAZD : INVLD );
      case IHARD: return( IHARD );
      default: return(INVLD);
    }
  }

  if ((ang) && (mi[MIT] != IMULT)) return(INVLD);

  switch(mi[MIT])
  { case IMULT: return( iscompact(mi[MKER]) ? IMULT : INVLD );
    case IPROD: return( ((mi[MDIM]==1) | (mi[MKT]==KPROD)) ? IPROD : INVLD );
    case IMLIN: return( ((mi[MKT]==KSPH) && (!lset) &&
      (mi[MDEG]<=2)) ? IMLIN : INVLD );
  }

  return(INVLD);
}

INT densinit(lf,des,h,cf,m)
lfit *lf;
design *des;
double h, *cf;
INT m;
{ INT deg, p, i, ii, j, nnz, rnz, lset, ang, status;
  double w;

  p = des->p; deg = lf->mi[MDEG];
  ff = des->xtwx.f2;
  cf[0] = NOSLN;
  for (i=1; i<p; i++) cf[i] = 0.0;

  if (!inre(des->xev,lf->xl,lf->mi[MDIM])) return(LF_XOOR);

  status = setintlimits(lf,des->xev,h,&ang,&lset);
  if (status != LF_OK) return(status);

  switch(selectintmeth(lf->mi,lset,ang))
  { case IMULT: des->itype = multint; break;
    case IPROD: des->itype = prodint; break;
    case IMLIN: des->itype = mlinint; break;
    case IHARD: des->itype = harint; break;
    case IHAZD: des->itype = hazint; break;
    case INVLD: ERROR(("Invalid integration method %d",lf->mi[MIT]));
                break;
    case IDEFA: ERROR(("No integration type available for this model"));
                break;
    default: ERROR(("densinit: unknown integral type"));
  }

  switch(deg)
  { case 0: rnz = 1; break;
    case 1: rnz = 1; break;
    case 2: rnz = lf->mi[MDIM]+1; break;
    case 3: rnz = lf->mi[MDIM]+2; break;
    default: ERROR(("densinit: invalid degree %d",deg));
  }
  if (lf_error) return(LF_ERR);

  setzero(des->ss,p);
  nnz = 0;
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    if (!cens(lf,ii))
    { w = des->w[i]*prwt(lf,ii);
      for (j=0; j<p; j++) des->ss[j] += des->X[i*p+j]*w;
      if (des->w[i]>0.00001) nnz++;
  } }

  if (lf->mi[MTG]==THAZ)
  { tmax = datum(lf,0,0);
    for (i=1; i<lf->mi[MN]; i++) tmax = MAX(tmax,datum(lf,0,i));
  }

  if (lf->mi[MDEB]>2)
  { printf("    LHS: ");
    for (i=0; i<p; i++) printf(" %8.5f",des->ss[i]);
    printf("\n");
  }

  switch(lf->mi[MLINK])
  { case LIDENT:
      cf[0] = 0.0;
      return(LF_OK);
    case LLOG:
      if (nnz<rnz) { cf[0] = -1000; return(LF_DNOP); }
      cf[0] = 0.0;
      return(LF_OK);
    case LSQRT:
      if (nnz<rnz) { cf[0] = 0.0; return(LF_DNOP); }
      cf[0] = 1.0;
      return(LF_OK);
    default:
      ERROR(("unknown link in densinit"));
      return(LF_ERR);
  }
}

void lforder(ind,x,l,r)
INT *ind, l, r;
double *x;
{ double piv;
  INT i, i0, i1;
  piv = (x[ind[l]]+x[ind[r]])/2;
  i0 = l; i1 = r;
  while (i0<=i1)
  { while ((i0<=i1) && (x[ind[i0]]<=piv)) i0++;
    while ((i0<=i1) && (x[ind[i1]]>piv))  i1--;
    if (i0<i1)
    { ISWAP(ind[i0],ind[i1]);
      i0++; i1--;
    }
  }
  /* now, x[ind[l..i1]] <= piv < x[ind[i0..r]].
     put the ties in the middle */
  while ((i1>=l) && (x[ind[i1]]==piv)) i1--;
  for (i=l; i<=i1; i++)
    if (x[ind[i]]==piv)
    { ISWAP(ind[i],ind[i1]);
      while (x[ind[i1]]==piv) i1--;
    }

  if (l<i1) lforder(ind,x,l,i1);
  if (i0<r) lforder(ind,x,i0,r);
}

double estdiv(coef,x,i0,i1,nvm,xlim,z)
double *coef, *x, *xlim, z;
INT i0, i1, nvm;
{ double cf[3], I[2], dlt, e0, e1;
  if (x[i0]==x[i1]) return(0.0);

  if (i0==-1) /* left boundary */
  { if (xlim[0]<xlim[1])
    { cf[0] = coef[i1]; cf[1] = -coef[i1+nvm]; cf[2] = 0.0;
      if (cf[0]==NOSLN) return(0.0);
      cf[0] *= z; cf[1] *= z; cf[2] *= z;
      recurint(0.0,x[i1]-xlim[0],cf,I,0,WRECT);
      return(I[0]);
    }
    if (coef[i1+nvm]<=0.0) return(0.0); /* neg slope, no boundary */
    return(exp(z*coef[i1])/(z*coef[i1+nvm]));
  }
  if (i1==-1) /* right boundary */
  { if (xlim[0]<xlim[1])
    { cf[0] = coef[i0]; cf[1] = coef[i0+nvm]; cf[2] = 0.0;
      if (cf[0]==NOSLN) return(0.0);
      cf[0] *= z; cf[1] *= z; cf[2] *= z;
      recurint(0.0,xlim[1]-x[i0],cf,I,0,WRECT);
      return(I[0]);
    }
    if (coef[i0+nvm]>=0.0) return(0.0);
    return(-exp(z*coef[i0])/(z*coef[i0+nvm]));
  }

  dlt =(x[i1]-x[i0])/2;
  cf[0] = z*coef[i0];
  cf[1] = z*coef[i0+nvm];
  cf[2] = z*(2*(coef[i1]-coef[i0])-dlt*(coef[i1+nvm]+3*coef[i0+nvm]))/(4*dlt*dlt);
  recurint(0.0,dlt,cf,I,0,WRECT);
  e0 = I[0];

  cf[0] = z*coef[i1];
  cf[1] = -z*coef[i1+nvm];
  cf[2] = z*(2*(coef[i0]-coef[i1])+dlt*(coef[i0+nvm]+3*coef[i1+nvm]))/(4*dlt*dlt);
  recurint(0.0,dlt,cf,I,0,WRECT);
  e1 = I[0];

  return(e0+e1);
}

double densintgl(lf,des,z)
lfit *lf;
design *des;
double z;
{ INT i, nv, *ind;
  double sum;
  if ((lf->mi[MDIM]>=2)|(lf->mi[MDEG]==0)|(lf->mi[MLINK]==LIDENT)) return(0.0);
  nv = lf->nv;
  if (lf->mi[MN]<nv) return(0.0);
  ind = des->ind;
  for (i=0; i<nv; i++) ind[i] = i;
  lforder(ind,vdptr(lf->xxev),0,nv-1);
  sum = estdiv(lf->coef,vdptr(lf->xxev),-1,ind[0],lf->nvm,lf->xl,z)
      + estdiv(lf->coef,vdptr(lf->xxev),ind[nv-1],-1,lf->nvm,lf->xl,z);
  for (i=1; i<nv; i++)
    sum += estdiv(lf->coef,vdptr(lf->xxev),ind[i-1],ind[i],lf->nvm,lf->xl,z);
  return(sum);
}

void densrenorm(lf,des)
lfit *lf;
design *des;
{ INT i;
  double sum;
  sum = densintgl(lf,des,1.0);
  if (sum==0.0) return;
  sum = log(sum);
  for (i=0; i<lf->nv; i++) lf->coef[i] -= sum;
}

void dlscv(des,lf)
lfit *lf;
design *des;
{ double df, fh, t0, z0, z1, x[MXDIM];
  int i, j, ev;
  z0 = densintgl(lf,des,2.0);
  z1 = df = 0.0;
  ev = lf->mi[MEV];
  if ((ev==EDATA) | (ev==ECROS)) ev = EFITP;
  for (i=0; i<lf->mi[MN]; i++)
  { for (j=0; j<lf->mi[MDIM]; j++) x[j] = datum(lf,j,i);
    fh = base(lf,i)+dointpoint(lf,des,x,PCOEF,ev,i);
    if (lf->mi[MLINK]==LLOG) fh = exp(fh);
    t0 = dointpoint(lf,des,x,PT0,ev,i);
    t0 = t0*t0;
    if (t0>1) t0 = 1;
    z1 += fh*(1-t0);
    df += t0;
  }
  vdptr(lf->L)[0] = z0-2*z1/(lf->mi[MN]-1);
  vdptr(lf->L)[1] = df;
}
