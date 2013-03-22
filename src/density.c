/*
 *   Copyright (c) 1998 Lucent Technologies.
 *   See README file for details.
 */

#include <math.h>
#include <stdio.h>
#include "local.h"

static double u[MXDIM], ilim[2*MXDIM], *ff, tmax;
static INT p, debug;

INT multint(), prodint(), gausint(), mlinint(), hazint(), harint();

void prresp(coef,resp)
double *coef, *resp;
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

INT multint(t,resp1,resp2,lf,cf,h,mi,ker)
struct tree *lf;
double *t, *resp1, *resp2, *cf, h;
INT *mi, ker;
{ INT d, i, j, k, m, m1, w, z, z1;
  double th, wt, dj[MXDIM];
  d = mi[MDIM];
  for (i=0; i<p*p; i++) resp1[i] = 0;
  m = 1; m1 = mi[MMINT]+1;
  for (i=0; i<d; i++)
  { m *= m1;
    dj[i] = (ilim[i+d]-ilim[i])/mi[MMINT];
  }
  for (i=0; i<m; i++)
  { z = i; w = 1;
    for (j=d-1; j>=0; j--)
    { z1 = z%m1;
      u[j] = ilim[j]+dj[j]*z1;
      w *= (4-2*(z1%2==0)-(z1==0)-(z1==mi[MMINT]));
      z /= m1;
    }
    wt = w*weight(u,lf->sca,d,ker,mi[MKT],h,lf->sty);
    if (wt>0)
    { fitfun(u,ff,lf->sca,d,mi[MDEG],mi[MKT],NULL,0,lf->sty);
      th = innerprod(ff,cf,p);
      addouter(resp1,ff,ff,p,wt*exp(th));
    }
  }
  wt = 1;
  for (j=0; j<d; j++) wt *= dj[j]/3;
  for (j=0; j<p; j++)
    for (k=j; k<p; k++)
      resp1[p*k+j] = resp1[p*j+k] = resp1[p*j+k]*wt;
  return(0);
}

INT mlinint(t,resp1,resp2,lf,cf,h,mi,ker)
struct tree *lf;
double *t, *resp1, *resp2, *cf, h;
INT *mi, ker;
{ double P[MXDIM*MXDIM], C[MXDIM*MXDIM], pu[MXDIM];
  double hd, nb, c, wt, wu, g[4], w0, w1, v, *sca;
  INT d, i, j, jmax, k, l, w, z, z1, jj[2];
  d = mi[MDIM]; sca = lf->sca;
  hd = 1; for (i=0; i<d; i++) hd *= h*sca[i];
  if (mi[MLINK]==LIDENT)
  { for (i=0; i<p*p; i++) resp1[i] = 0.0;
    resp1[0] = wint(d,NULL,0,ker)*hd;
    if (mi[MDEG]==0) return(0);
    jj[0] = 2; w0 = wint(d,jj,1,ker)*hd*h*h;
    for (i=0; i<d; i++) resp1[(i+1)*p+i+1] = w0*sca[i]*sca[i];
    if (mi[MDEG]==1) return(0);
    for (i=0; i<d; i++)
    { j = p-(d-i)*(d-i+1)/2;
      resp1[j] = resp1[p*j] = w0*sca[i]*sca[i]/2;
    }
    if (d>1) { jj[1] = 2; w0 = wint(d,jj,2,ker)*hd*h*h*h*h; }
    jj[0] = 4; w1 = wint(d,jj,1,ker)*hd*h*h*h*h/4;
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
    return(0);
  }
  switch(mi[MDEG])
  { case 0:
      resp1[0] = exp(cf[0])*wint(d,NULL,0,ker);
      for (i=0; i<d; i++) resp1[0] *= h*sca[i];
      return(0);
    case 1:
      nb = 0.0;
      for (i=1; i<=d; i++)
      { v = h*cf[i]*sca[i-1];
        nb += v*v;
      }
      if (ker==WGAUS)
      { w0 = 1/(GFACT*GFACT);
        g[0] = exp(cf[0]+w0*nb/2+d*log(S2PI/2.5));
        g[1] = g[3] = g[0]*w0;
        g[2] = g[0]*w0*w0;
      }
      else
      { wt = wu = exp(cf[0]);
        w0 = wint(d,NULL,0,ker); g[0] = wt*w0;
        g[1] = g[2] = g[3] = 0.0;
        j = 0; jmax = (d+2)*mi[MMINT];
        while ((j<jmax) && (wt*w0/g[0]>1.0e-8))
        { j++;
          jj[0] = 2*j; w0 = wint(d,jj,1,ker);
          if (d==1) g[3] += wt * w0;
          else
          { jj[0] = 2; jj[1] = 2*j-2; w1 = wint(d,jj,2,ker);
            g[3] += wt*w1;
            g[2] += wu*(w0-w1);
          }
          wt /= (2*j-1.0); g[1] += wt*w0;
          wt *= nb/(2*j); g[0] += wt*w0;
          wu /= (2*j-1.0)*(2*j);
          if (j>1) wu *= nb;
        }
        if (j==jmax) WARN(("mlinint: series not converged"))
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
      return(0);
    case 2:
      if (ker==WGAUS) return(gausint(t,resp1,resp2,cf,h,d,sca));
      k = d+1;
      for (i=0; i<d; i++)
      { C[i*(d+1)] = cf[k++]*sca[i]*sca[i];
        for (j=i+1; j<d; j++) C[i*d+j] = C[j*d+i] = cf[k++]*sca[i]*sca[j];
      }
      eigen(C,P,d,20);
      k = 1;
      for (i=0; i<d; i++) k *= (mi[MMINT]+1);
      for (i=0; i<p*p; i++) resp1[i] = 0;
      for (i=0; i<k; i++)
      { z = i; w = 1;
        for (j=d-1; j>=0; j--)
        { z1 = z%(mi[MMINT]+1);
          pu[j] = h*sca[j]*(-1+2.0*z1/mi[MMINT]);
          w *= (4-2*(z1%2==0)-(z1==0)-(z1==mi[MMINT]));
          z /= (mi[MMINT]+1);
        }
        wt = w*weight(pu,sca,d,ker,KSPH,h,lf->sty);
        if (wt>0)
        { for (j=0; j<d; j++)
          { u[j] = 0.0;
            for (l=0; l<d; l++) u[j] += sca[j]*P[j*d+l]*pu[l]/sca[l];
          }
          c = cf[0];
          for (j=0; j<d; j++)
            c += cf[j+1]*u[j]+SQR(pu[j]/sca[j])*C[j*(d+1)]/2;
          fitfun(u,ff,sca,d,mi[MDEG],KSPH,NULL,0,lf->sty);
          addouter(resp1,ff,ff,p,wt*exp(c));
      } }
      for (i=0; i<p*p; i++)
        for (j=0; j<d; j++)
          resp1[i] *= 2.0*h*sca[j]/(3*mi[MMINT]);
      return(0);
    default: ERROR(("mlinint: deg<=2 only"))
  }
  return(1);
}

INT exbctay(b,c,n,z) /* n-term taylor series of e^(bx+cx^2) */
double b, c, *z;
INT n;
{ double ec[20];
  INT i, j;
  z[0] = 1;
  for (i=1; i<=n; i++) z[i] = z[i-1]*b/i;
  if (c==0.0) return(n);
  if (n>=40)
  { WARN(("exbctay limit to n<40"));
    n = 39;
  }
  ec[0] = 1;
  for (i=1; 2*i<=n; i++) ec[i] = ec[i-1]*c/i;
  for (i=n; i>1; i--)
    for (j=1; 2*j<=i; j++)
      z[i] += ec[j]*z[i-2*j];
  return(n);
}

double explinjtay(l0,l1,j,cf)
/* int_l0^l1 x^j e^(a+bx+cx^2); exbctay aroud l1 */
double l0, l1, *cf;
INT j;
{ double tc[40], f, s;
  INT k, n;
  if ((l0!=0.0) | (l1!=1.0)) WARN(("explinjtay: invalid l0, l1"));
  n = exbctay(cf[1]+2*cf[2]*l1,cf[2],20,tc);
  s = tc[0]/(j+1);
  f = 1/(j+1);
  for (k=1; k<=n; k++)
  { f *= -k/(j+k+1.0);
    s += tc[k]*f;
  }
  return(f);
}

void explint1(l0,l1,cf,I,p) /* int x^j exp(a+bx); j=0..p-1 */
double l0, l1, *cf, *I;
INT p;
{ double y0, y1, f;
  INT j, k, k1;
  y0 = exp(cf[0]+l0*cf[1]);
  y1 = exp(cf[0]+l1*cf[1]);
  if (p<2*fabs(cf[1])) k = p; else k = (INT)fabs(cf[1]);

  if (k>0)
  { I[0] = (y1-y0)/cf[1];
    for (j=1; j<k; j++) /* forward steps for small j */
    { y1 *= l1; y0 *= l0;
      I[j] = (y1-y0-j*I[j-1])/cf[1];
    }
    if (k==p) return;
    y1 *= l1; y0 *= l0;
  }

  f = 1; k1 = k;
  while ((k<50) && (f>1.0e-8)) /* initially Ik = diff(x^{k+1}e^{a+bx}) */
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
    if (k>=p) f *= fabs(cf[1])/(k+1);
    k++;
  }
  if (k==50) WARN(("explint1: want k>50"))
  I[k] = 0.0;
  for (j=k-1; j>=k1; j--) /* now do back step recursion */
    I[j] = (I[j]-cf[1]*I[j+1])/(j+1);
}

void explintyl(l0,l1,cf,I,p) /* small c, use taylor series and explint1 */
double l0, l1, *cf, *I;
INT p;
{ INT i;
  double c;
  explint1(l0,l1,cf,I,p+8);
  c = cf[2];
  for (i=0; i<p; i++)
    I[i] = (((I[i+8]*c/4+I[i+6])*c/3+I[i+4])*c/2+I[i+2])*c+I[i];
}

void solvetrid(X,y,m)
double *X, *y;
INT m;
{ INT i;
  double s;
  for (i=1; i<m; i++)
  { s = X[3*i]/X[3*i-2];
    X[3*i] = 0; X[3*i+1] -= s*X[3*i-1];
    y[i] -= s*y[i-1];
  }
  for (i=m-2; i>=0; i--)
  { s = X[3*i+2]/X[3*i+4];
    X[3*i+2] = 0;
    y[i] -= s*y[i+1];
  }
  for (i=0; i<m; i++) y[i] /= X[3*i+1];
}

void initi0i1(I,cf,y0,y1,l0,l1)
double *I, *cf, y0, y1, l0, l1;
{ double a0, a1, c, d, bi;
  d = -cf[1]/(2*cf[2]); c = sqrt(2*fabs(cf[2]));
  a0 = c*(l0-d); a1 = c*(l1-d);
  if (cf[2]<0)
  { bi = exp(cf[0]+cf[1]*d+cf[2]*d*d)/c;
    if (a0>0)
    { if (a0>6) I[0] = (y0*ptail(-a0)-y1*ptail(-a1))/c;
      else I[0] = S2PI*(pnorm(-a0,0.0,1.0)-pnorm(-a1,0.0,1.0))*bi;
    }
    else
    { if (a1< -6) I[0] = (y1*ptail(a1)-y0*ptail(a0))/c;
      else I[0] = S2PI*(pnorm(a1,0.0,1.0)-pnorm(a0,0.0,1.0))*bi;
    }
  }
  else
    I[0] = (y1*daws(a1)-y0*daws(a0))/c;
  I[1] = (y1-y0)/(2*cf[2])+d*I[0];
}

void explinsid(l0,l1,cf,I,p) /* large b; don't use fwd recursion */
double l0, l1, *cf, *I;
INT p;
{ INT k, k0, k1, k2;
  double y0, y1, Z[150];
if (debug) printf("side: %8.5f %8.5f %8.5f    limt %8.5f %8.5f  p %2d\n",cf[0],cf[1],cf[2],l0,l1,p);
 
  k0 = 2;
  k1 = (INT)(fabs(cf[1])+fabs(2*cf[2]));
  if (k1<2) k1 = 2;
  if (k1>p+20) k1 = p+20;
  k2 = p+20;

  if (debug) printf("k0 %2d  k1 %2d  k2 %2d  p %2d\n",k0,k1,k2,p);

  y0 = exp(cf[0]+l0*(cf[1]+l0*cf[2]));
  y1 = exp(cf[0]+l1*(cf[1]+l1*cf[2]));
  initi0i1(I,cf,y0,y1,l0,l1);
if (debug) printf("i0 %8.5f  i1 %8.5f\n",I[0],I[1]);

  y1 *= l1; y0 *= l0; /* should be x^(k1)*exp(..) */
  if (k0<k1) /* center steps; initially x^k*exp(...) */
    for (k=k0; k<k1; k++)
    { y1 *= l1; y0 *= l0;
      I[k] = y1-y0;
      Z[3*k] = k; Z[3*k+1] = cf[1]; Z[3*k+2] = 2*cf[2];
    }
   
  y1 *= l1; y0 *= l0; /* should be x^(k1)*exp(..) */
if (debug) printf("k1 %2d  y0 %8.5f  y1 %8.5f\n",k1,y0,y1);
  for (k=k1; k<k2; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }
  I[k2] = I[k2+1] = 0.0;
  for (k=k2-1; k>=k1; k--)
    I[k] = (I[k]-cf[1]*I[k+1]-2*cf[2]*I[k+2])/(k+1);

  if (k0<k1)
  { I[k0] -= k0*I[k0-1];
    I[k1-1] -= 2*cf[2]*I[k1];
    Z[3*k0] = Z[3*k1-1] = 0;
    solvetrid(&Z[3*k0],&I[k0],k1-k0);
  }
if (debug)
{ printf("explinsid:\n");
  for (k=0; k<p; k++) printf("  %8.5f\n",I[k]);
}
}

void explinbkr(l0,l1,cf,I,p) /* small b,c; use back recursion */
double l0, l1, *cf, *I;
INT p;
{ INT k, km;
  double y0, y1;
  y0 = exp(cf[0]+l0*(cf[1]+cf[2]*l0));
  y1 = exp(cf[0]+l1*(cf[1]+cf[2]*l1));
  km = p+10;
  for (k=0; k<=km; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }
  I[km+1] = I[km+2] = 0;
  for (k=km; k>=0; k--)
    I[k] = (I[k]-cf[1]*I[k+1]-2*cf[2]*I[k+2])/(k+1);
}

void explinfbk0(l0,l1,cf,I,p) /* fwd and bac recur; b=0; c<0 */
double l0, l1, *cf, *I;
INT p;
{ double y0, y1, f1, f2, f, ml2;
  INT k, ks;

  y0 = exp(cf[0]+l0*l0*cf[2]);
  y1 = exp(cf[0]+l1*l1*cf[2]);
  initi0i1(I,cf,y0,y1,l0,l1);

  ml2 = MAX(l0*l0,l1*l1);
  ks = 1+(INT)(2*fabs(cf[2])*ml2);
  if (ks<2) ks = 2;
  if (ks>p-3) ks = p;

  /* forward recursion for k < ks */
  for (k=2; k<ks; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = (y1-y0-(k-1)*I[k-2])/(2*cf[2]);
  }
  if (ks==p) return;

  y1 *= l1*l1; y0 *= l0*l0;
  for (k=ks; k<p; k++) /* set I[k] = x^{k+1}e^(a+cx^2) | {l0,l1} */
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }

  /* initialize I[p-2] and I[p-1] */
  f1 = 1.0/p; f2 = 1.0/(p-1);
  I[p-1] *= f1; I[p-2] *= f2;
  k = p; f = 1.0;
  while (f>1.0e-8)
  { y1 *= l1; y0 *= l0;
    if ((k-p)%2==0) /* add to I[p-2] */
    { f2 *= -2*cf[2]/(k+1);
      I[p-2] += (y1-y0)*f2;
    }
    else /* add to I[p-1] */
    { f1 *= -2*cf[2]/(k+1);
      I[p-1] += (y1-y0)*f1;
      f *= 2*fabs(cf[2])*ml2/(k+1);
    }
    k++;
  }
  
  /* use back recursion for I[ks..(p-3)] */
  for (k=p-3; k>=ks; k--)
    I[k] = (I[k]-2*cf[2]*I[k+2])/(k+1);
}

void explinfbk(l0,l1,cf,I,p) /* fwd and bac recur; b not too large */
double l0, l1, *cf, *I;
INT p;
{ double y0, y1;
  INT k, ks, km;

  y0 = exp(cf[0]+l0*(cf[1]+l0*cf[2]));
  y1 = exp(cf[0]+l1*(cf[1]+l1*cf[2]));
  initi0i1(I,cf,y0,y1,l0,l1);

  ks = (INT)(3*fabs(cf[2]));
  if (ks<3) ks = 3;
  if (ks>0.75*p) ks = p; /* stretch the forward recurs as far as poss. */
  /* forward recursion for k < ks */
  for (k=2; k<ks; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = (y1-y0-cf[1]*I[k-1]-(k-1)*I[k-2])/(2*cf[2]);
  }
  if (ks==p) return;

  km = p+15;
  y1 *= l1*l1; y0 *= l0*l0;
  for (k=ks; k<=km; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }
  I[km+1] = I[km+2] = 0.0;
  for (k=km; k>=ks; k--)
    I[k] = (I[k]-cf[1]*I[k+1]-2*cf[2]*I[k+2])/(k+1);
}

void recent(I,resp,wt,p,s,x)
double *I, *resp, *wt, x;
INT p, s;
{ INT i, j;

  /* first, use W taylor series I -> resp */
  for (i=0; i<=p; i++)
  { resp[i] = 0.0;
    for (j=0; j<s; j++) resp[i] += wt[j]*I[i+j];
  }

  /* now, recenter x -> 0 */
  if (x==0) return;
  for (j=0; j<=p; j++) for (i=p; i>j; i--) resp[i] += x*resp[i-1];
}

void recurint(l0,l2,cf,resp,p,ker)
double l0, l2, *cf, *resp;
INT p, ker;
{ INT i, s;
  double l1, d0, d1, d2, dl, z0, z1, z2, wt[20], ncf[3], I[50], r1[5], r2[5];
if (debug) printf("\nrecurint: %8.5f %8.5f %8.5f   %8.5f %8.5f\n",cf[0],cf[1],cf[2],l0,l2);

  if (cf[2]==0) /* go straight to explint1 */
  { s = wtaylor(wt,0.0,ker);
if (debug) printf("case 1\n");
    explint1(l0,l2,cf,I,p+s);
    recent(I,resp,wt,p,s,0.0);
    return;
  }

  dl = l2-l0;
  d0 = cf[1]+2*l0*cf[2];
  d2 = cf[1]+2*l2*cf[2];
  z0 = cf[0]+l0*(cf[1]+l0*cf[2]);
  z2 = cf[0]+l2*(cf[1]+l2*cf[2]);

  if ((fabs(cf[1]*dl)<1) && (fabs(cf[2]*dl*dl)<1))
  { ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
if (debug) printf("case 2\n");
    s = wtaylor(wt,l0,ker);
    explinbkr(0.0,dl,ncf,I,p+s);
    recent(I,resp,wt,p,s,l0);
    return;
  }

  if (fabs(cf[2]*dl*dl)<0.001) /* small c, use explint1+tay.ser */
  { ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
if (debug) printf("case small c\n");
    s = wtaylor(wt,l0,ker);
    explintyl(0.0,l2-l0,ncf,I,p+s);
    recent(I,resp,wt,p,s,l0);
    return;
  }

  if (d0*d2<=0) /* max/min in [l0,l2] */
  { l1 = -cf[1]/(2*cf[2]);
    z1 = cf[0]+l1*(cf[1]+l1*cf[2]);
    d1 = 0.0;
    if (cf[2]<0) /* peak, integrate around l1 */
    { s = wtaylor(wt,l1,ker);
      ncf[0] = z1; ncf[1] = 0.0; ncf[2] = cf[2];
if (debug) printf("case peak  p %2d  s %2d\n",p,s);
      explinfbk0(l0-l1,l2-l1,ncf,I,p+s);
      recent(I,resp,wt,p,s,l1);
      return;
    }
  }

  if ((d0-2*cf[2]*dl)*(d2+2*cf[2]*dl)<0) /* max/min is close to [l0,l2] */
  { l1 = -cf[1]/(2*cf[2]);
    z1 = cf[0]+l1*(cf[1]+l1*cf[2]);
    if (l1<l0) { l1 = l0; z1 = z0; }
    if (l1>l2) { l1 = l2; z1 = z2; }

    if ((z1>=z0) & (z1>=z2)) /* peak; integrate around l1 */
    { s = wtaylor(wt,l1,ker);
if (debug) printf("case 4\n");
      d1 = cf[1]+2*l1*cf[2];
      ncf[0] = z1; ncf[1] = d1; ncf[2] = cf[2];
      explinfbk(l0-l1,l2-l1,ncf,I,p+s);
      recent(I,resp,wt,p,s,l1);
      return;
    }

    /* trough; integrate [l0,l1] and [l1,l2] */
    for (i=0; i<=p; i++) r1[i] = r2[i] = 0.0;
    if (l0<l1)
    { s = wtaylor(wt,l0,ker);
if (debug) printf("case 5\n");
      ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
      explinfbk(0.0,l1-l0,ncf,I,p+s);
      recent(I,r1,wt,p,s,l0);
    }
    if (l1<l2)
    { s = wtaylor(wt,l2,ker);
if (debug) printf("case 6\n");
      ncf[0] = z2; ncf[1] = d2; ncf[2] = cf[2];
      explinfbk(l1-l2,0.0,ncf,I,p+s);
      recent(I,r2,wt,p,s,l2);
    }
    for (i=0; i<=p; i++) resp[i] = r1[i]+r2[i];
    return;
  }

  /* Now, quadratic is monotone on [l0,l2]; big b; moderate c */
  if (z2>z0+3) /* steep increase, expand around l2 */
  { s = wtaylor(wt,l2,ker);
if (debug) printf("case 7\n");


    ncf[0] = z2; ncf[1] = d2; ncf[2] = cf[2];
    explinsid(l0-l2,0.0,ncf,I,p+s);
    recent(I,resp,wt,p,s,l2);
if (debug) printf("7 resp: %8.5f %8.5f %8.5f %8.5f\n",resp[0],resp[1],resp[2],resp[3]);
    return;
  }

  /* bias towards expansion around l0, because it's often 0 */
if (debug) printf("case 8\n");
  s = wtaylor(wt,l0,ker);
  ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
  explinsid(0.0,l2-l0,ncf,I,p+s);
  recent(I,resp,wt,p,s,l0);
  return;
}

void onedint(cf,mi,l0,l1,resp,ker) /* int W(u)u^j exp(..), j=0..2*deg */
double *cf, l0, l1, *resp;
INT *mi, ker;
{ double d, u, y, ncf[4], rr[5];
  INT deg, i, j;
if (debug) printf("onedint: %f %f %f   %f %f\n",cf[0],cf[1],cf[2],l0,l1);
  deg = mi[MDEG];

  if (deg<=2)
  { for (i=0; i<3; i++) ncf[i] = (i>deg) ? 0.0 : cf[i];
    ncf[2] /= 2;

    if (l1>0)
      recurint(MAX(l0,0.0),l1,ncf,resp,2*deg,ker);
    else for (i=0; i<=2*deg; i++) resp[i] = 0;

    if (l0<0)
    { ncf[1] = -ncf[1];
      l0 = -l0; l1 = -l1;
      recurint(MAX(l1,0.0),l0,ncf,rr,2*deg,ker);
    }
    else for (i=0; i<=2*deg; i++) rr[i] = 0.0;

    for (i=0; i<=2*deg; i++)
      resp[i] += (i%2==0) ? rr[i] : -rr[i];

    return;
  }

  for (i=0; i<=mi[MMINT]; i++)
  { u = l0+(l1-l0)*i/mi[MMINT];
    y = 0; d = 1;
    for (j=0; j<=deg; j++)
    { y += cf[j]*d;
      d *= u/(j+1);
    }
    y = (4-2*(i%2==0)-(i==0)-(i==mi[MMINT])) * W(fabs(u),ker)*exp(MIN(y,300));
    for (j=0; j<=2*deg; j++)
    { resp[j] += y;
      y *= u;
    }
  }
  for (j=0; j<=2*deg; j++) resp[j] = resp[j]*(l1-l0)/(3*mi[MMINT]);
}

INT prodint(t,resp1,resp2,lf,coef,h,mi,ker)
struct tree *lf;
double *t, *resp1, *resp2, *coef, h;
INT *mi, ker;
{ INT d, deg, i, j, j0, j1, k0, k1, f0, f1, m;
  double ea, cf[MXDEG+1], hj, hs, *sca;
  d = mi[MDIM]; deg = mi[MDEG]; sca = lf->sca;
  for (i=0; i<p*p; i++) resp1[i] = 1;
  cf[0] = coef[0];
  for (i=0; i<d; i++)
  { hj = 1; hs = h*sca[i];
    for (j=0; j<deg; j++)
    { hj *= hs;
      cf[j+1] = hj*coef[j*d+i+1];
    }
    onedint(cf,mi,ilim[i]/hs,ilim[i+d]/hs,resp2,ker);
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
    } }
  }
  ea = exp(-(d-1)*coef[0]);
/* factorial divisions */
  for (k0=0; k0<=deg; k0++)
  { f0 = (k0==0) ? 1 : f0*k0;
    for (k1=0; k1<=k0; k1++)
    { f1 = (k1==0) ? f0 : f1*k1;
      m = p*(1+(k1-1)*d)+1+(k0-1)*d;
      for (j0=(k0==0)*(d-1); j0<d; j0++)
        for (j1=(k1==0)*(d-1); j1<d; j1++)
          resp1[m+p*j1+j0] *= ea/f1;
  } }
  for (i=0; i<p; i++) /* symmetry */
    for (j=i; j<p; j++)
      resp1[p*j+i] = resp1[i*p+j];
  return(0);
}

INT gausint(t,resp,C,cf,h,d,sca)
double *t, *resp, *C, *cf, h, *sca;
INT d;
{ double nb, det, z, *P;
  INT i, j, k, l, m1, m2, f;
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
    if (det <= 0) return(1);
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
  for (i=0; i<p*p; i++) resp[i] *= z;
  return(0);
}

INT hrao(dfx,cf,mi,h,sca,ker,r1,sty)
double *dfx, *cf, h, *sca, *r1;
INT *mi, ker, *sty;
{ double s, t0, t1, wt, th;
  INT j;
  s = 0;
  switch (mi[MKT])
  { case KPROD:
      for (j=1; j<mi[MDIM]; j++) s = MAX(s,fabs(dfx[j]/(h*sca[j])));
      break;
    case KSPH:
      for (j=1; j<mi[MDIM]; j++) s += SQR(dfx[j]/(h*sca[j]));
      break;
    default: ERROR(("harint not available for kt=%d",mi[MKT]))
  }
  if ((lf_error) | (s>1)) return(0);
  if (mi[MKT]==KPROD) t1 = h*sca[0];
    else t1 = sqrt(1-s)*h*sca[0];
  t0 = -t1;
  if (t0<ilim[0]) t0 = ilim[0];
  if (t1>ilim[mi[MDIM]]) t1 = ilim[mi[MDIM]];
  if (t1>dfx[0]) t1 = dfx[0];
  if (t1<t0) return(0);
  for (j=0; j<p*p; j++) r1[j] = 0;
  for (j=0; j<=mi[MMINT]; j++)
  { dfx[0] = t0+(t1-t0)*j/mi[MMINT];
    wt = weight(dfx,sca,mi[MDIM],ker,mi[MKT],h,sty);
    fitfun(dfx,ff,sca,mi[MDIM],mi[MDEG],mi[MKT],NULL,0,sty);
    th = innerprod(cf,ff,p);
    wt *= 2+2*(j&1)-(j==0)-(j==mi[MMINT]);
    addouter(r1,ff,ff,p,wt*exp(th));
  }
  for (j=0; j<p*p; j++) r1[j] *= (t1-t0)/(3*mi[MMINT]);
  return(1);
}

INT harint(t,resp,r1,lf,cf,h,mi,ker)
struct tree *lf;
double *t, *resp, *r1, *cf, h;
INT *mi, ker;
{ INT i, j;
  double dfx[MXDIM];
  for (i=0; i<p*p; i++) resp[i] = 0;
  for (i=0; i<=mi[MN]; i++)
  { if (i==mi[MN])
    { dfx[0] = tmax-t[0];
      for (j=1; j<mi[MDIM]; j++) dfx[j] = 0.0;
    }
    else
      for (j=0; j<mi[MDIM]; j++) dfx[j] = lf->x[j][i]-t[j];
    if (hrao(dfx,cf,mi,h,lf->sca,ker,r1,lf->sty))
      for (j=0; j<p*p; j++) resp[j] += r1[j];
  }
  return(0);
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

INT hazint(t,resp,x,lf,cf,h,mi,ker)
struct tree *lf;
double *t, *resp, *x, *cf, h;
INT *mi, ker;
{ INT d, i, j, k;
  double r[2*MXDEG+1], rr[2*MXDEG+1], dfx[MXDIM],
         hj, hs, ncf[MXDEG], wx, ef, lt;
  for (i=0; i<p*p; i++) resp[i] = 0;
  hj = hs = h*lf->sca[0];
  d = mi[MDIM];
  ncf[0] = cf[0];
  for (i=1; i<=mi[MDEG]; i++)
  { ncf[i] = hj*cf[(i-1)*d+1]; hj *= hs;
  }
  k = 0; lt = -1.0;
  for (i=0; i<=mi[MN]; i++)
  { if (i==mi[MN])
    { dfx[0] = tmax-t[0];
      for (j=1; j<d; j++) dfx[j] = 0.0;
    }
    else
      for (j=0; j<d; j++) dfx[j] = lf->x[j][i]-t[j];
    if ((d==1) & (dfx[0]>=hs)) k++;
    else if (dfx[0]>ilim[0])
    { wx = 1;
      ef = 0.0;
      for (j=1; j<d; j++)
      { wx *= W(dfx[j]/(h*lf->sca[j]),ker);
        if (mi[MDEG]>=1) ef += cf[j+1]*dfx[j];
        if (mi[MDEG]>=2) ef += cf[j+d+1]*dfx[j]*dfx[j]/2;
      }
      ef = exp(ef);
      if (wx>0)
      { if (dfx[0]>hs) dfx[0] = hs;
        if ((lt==-1.0) || (dfx[0]!=lt))
        { onedint(ncf,mi,ilim[0]/hs,dfx[0]/hs,r,ker);
          lt = dfx[0];
        }
        hj = 1;
        for (j=0; j<=2*mi[MDEG]; j++) { hj *= hs; rr[j] = r[j]*wx*hj*ef; }
        hadd(resp,rr,p,d,mi[MDEG],dfx);
      }
    }
  }
  if (k>0)
  { onedint(ncf,mi,ilim[0]/hs,1.0,r,ker);
    lt = 1.0;
    hj = 1;
    for (i=0; i<=2*mi[MDEG]; i++) { hj *= hs; r[i] *= k*hj; }
    for (i=1; i<d; i++) dfx[i] = 0.0;
    hadd(resp,r,p,d,mi[MDEG],dfx);
  }
  for (i=0; i<p; i++)
    for (j=0; j<i; j++)
      resp[i*p+j] = resp[j*p+i];
  return(0);
}

double likeden(lf,des,x,zz,h,mi)
struct tree *lf;
struct design *des;
double *x, h;
INT *zz, *mi;
{ double lk, r, *A;
  INT i, j, p;
  p = des->p; A = des->Z;
  if ((mi[MLINK]==LIDENT) && (des->cf[0] != 0.0)) return(0.0);
  *zz = (des->itype)(x,A,des->Q,lf,des->cf,h,mi,mi[MKER]);
  if (debug) prresp(des->cf,A);
  if (*zz) return(0.0);
  if (mi[MLINK]==LLOG)
  { r = des->ss[0]/A[0];
    des->cf[0] += log(r);
    for (i=1; i<p*p; i++) A[i] *= r;
    A[0] = des->ss[0];
    lk = -A[0];
    for (i=0; i<p; i++) lk += des->cf[i]*des->ss[i];
  }
  else lk = 0.0;
  for (i=0; i<p; i++)
  { if (A[i*(p+1)]<=0)
    { WARN(("likeden: negative diagonal, zeroing"))
      des->dg[i] = 0;
    }
    else
      des->dg[i] = 1/sqrt(A[i*(p+1)]);
    if (mi[MLINK]==LLOG)
      des->res[i]= des->ss[i]-A[i];
    else
    { des->res[i] = des->ss[i];
      for (j=0; j<p; j++)
        des->res[i] -= A[i*p+j]*des->cf[j];
    }
  }
  switch(des->sm)
  { case 1:
      for (i=0; i<p; i++)
        for (j=0; j<p; j++)
          A[i*p+j] *= des->dg[i]*des->dg[j];
      eigen(A,des->Q,p,20);
      break;
    case 2:
      choldec(A,p);
      break;
    default: ERROR(("likeden: unknown solution method %d",des->sm))
  }
  return(lk);
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

INT densinit(lf,des,x,h,cf,mi,m)
struct tree *lf;
struct design *des;
double *x, h, *cf;
INT *mi, m;
{ INT d, deg, i, ii, j, nnz, rnz, va, it, lset, hzd, kt;
  double w;
debug = 0;
  hzd = mi[MTG]==THAZ;
  d = mi[MDIM]; p = des->p;
  deg = mi[MDEG];
  kt = mi[MKT];   it = mi[MIT]; 
  ff = des->f2;
  lset = 0;
  for (i=0; i<d; i++)
  { if (lf->sty[i]==KANG)
    { ilim[i+d] = ((h<2) ? 2*asin(h/2) : PI)*lf->sca[i];
      ilim[i] = -ilim[i+d];
      if (it==IDEFA) it = IMULT;
    }
    else
    { ilim[i+d] = h*lf->sca[i];
      ilim[i] = -ilim[i+d];
      if (lf->sty[i]==KLEF) { ilim[i+d] = 0; lset = 1; }
      if (lf->sty[i]==KRIG) { ilim[i] = 0; lset = 1; }
      if (lf->xl[i]<lf->xl[i+d]) /* user limits for this variable */
      { if (lf->xl[i]-x[i]> ilim[i])
        { ilim[i] = lf->xl[i]-x[i]; lset=1; }
        if (lf->xl[i+d]-x[i]< ilim[i+d])
        { ilim[i+d] = lf->xl[i+d]-x[i]; lset=1; }
      }
    }
    if (ilim[i]==ilim[i+d]) return(1); /* empty integration */
  }

  if (it==IDEFA) /* choose default integration method */
  { if (hzd) it = ((d==1)|(kt==KPROD)) ? IHAZD : IHARD;
    else
    { it = IMULT;
      if ((d==1) || (kt==KPROD)) it = IPROD;
      if ((!lset) & (deg<=1)) it = IMLIN;
      if ((!lset) & (deg==2) & (mi[MLINK]==LIDENT)) it = IMLIN;
      if (mi[MKER]==WGAUS) it = IMLIN;
    }
  }
  switch(it)
  { case IMULT: va = (mi[MKER] != WGAUS); break;
    case IPROD: va = (d==1) | (kt==KPROD); break;
    case IMLIN: va = (kt==KSPH) && (!lset) && (deg<=2); break;
    case IHAZD: va = (hzd) & ((d==1) | (kt==KPROD)); break;
    case IHARD: va = hzd; break;
    default: ERROR(("densinit: unknown integration method %d",it))
  }
  if (va==0) ERROR(("densinit: invalid integration method %d",it))
  switch(deg)
  { case 0: rnz = 1; break;
    case 1: rnz = 1; break;
    case 2: rnz = d+1; break;
    case 3: rnz = d+2; break;
    default: ERROR(("densinit: invalid degree %d",deg))
  }
  switch(it)
  { case IMULT: des->itype = multint; break;
    case IPROD: des->itype = prodint; break;
    case IMLIN: des->itype = mlinint; break;
    case IHARD: des->itype = harint; break;
    case IHAZD: des->itype = hazint; break;
    default: ERROR(("densinit: invalid integral type %d",it))
  }
  if ((lf_error) | (!inre(x,lf->xl,d))) return(1);
  for (i=0; i<p; i++) des->ss[i] = 0;
  nnz = 0;
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    if (!cens(lf,ii))
    { w = des->w[i]*prwt(lf,ii);
      for (j=0; j<p; j++) des->ss[j] += des->X[i*p+j]*w;
      if (des->w[i]>0.00001) nnz++;
  } }
  if (hzd)
  { tmax = lf->x[0][0];
    for (i=1; i<mi[MN]; i++) tmax = MAX(tmax,lf->x[0][i]);
  }
if (debug) printf("ss: %8.5f %8.5f %8.5f\n",des->ss[0],des->ss[1],des->ss[2]);
  if ((mi[MLINK]==LLOG) && (nnz<rnz)) return(1);
  return(0);
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

double estdiv(coef,x,i0,i1,nvm,xlim)
double *coef, *x, *xlim;
INT i0, i1, nvm;
{ double cf[3], I[2], dlt, e0, e1, fl;
  if (x[i0]==x[i1]) return(0.0);

  if (i0==-1) /* left boundary */
  { if (xlim[0]<xlim[1])
    { cf[0] = coef[i1]; cf[1] = -coef[i1+nvm]; cf[2] = 0.0;
      recurint(0.0,x[i1]-xlim[0],cf,I,0,WRECT);
      return(I[0]);
    }
    if (coef[i1+nvm]<=0.0) return(0.0); /* neg slope, no boundary */
    return(exp(coef[i1])/coef[i1+nvm]);
  }
  if (i1==-1) /* right boundary */
  { if (xlim[0]<xlim[1])
    { cf[0] = coef[i0]; cf[1] = coef[i0+nvm]; cf[2] = 0.0;
      recurint(0.0,xlim[1]-x[i0],cf,I,0,WRECT);
      return(I[0]);
    }
    if (coef[i0+nvm]>=0.0) return(0.0);
    return(-exp(coef[i0])/coef[i0+nvm]);
  }

  dlt =(x[i1]-x[i0])/2;
  cf[0] = coef[i0];
  cf[1] = coef[i0+nvm];
  cf[2] = (2*(coef[i1]-coef[i0])-dlt*(coef[i1+nvm]+3*coef[i0+nvm]))/(4*dlt*dlt);
  recurint(0.0,dlt,cf,I,0,WRECT);
  e0 = I[0];

  cf[0] = coef[i1];
  cf[1] = -coef[i1+nvm];
  cf[2] = (2*(coef[i0]-coef[i1])+dlt*(coef[i0+nvm]+3*coef[i1+nvm]))/(4*dlt*dlt);
  recurint(0.0,dlt,cf,I,0,WRECT);
  e1 = I[0];

  return(e0+e1);
}

void densrenorm(lf,des)
struct tree *lf;
struct design *des;
{ INT i, nv, *ind;
  double sum;
  if ((lf->mi[MDIM]>=2)|(lf->mi[MDEG]==0)|(lf->mi[MLINK]==LIDENT)) return;
  nv = lf->nv;
  if (lf->mi[MN]<nv) return;
  ind = des->ind;
  for (i=0; i<nv; i++) ind[i] = i;
  lforder(ind,lf->xev,0,nv-1);
  sum = estdiv(lf->coef,lf->xev,-1,ind[0],lf->nvm,lf->xl)
      + estdiv(lf->coef,lf->xev,ind[nv-1],-1,lf->nvm,lf->xl);
  for (i=1; i<nv; i++)
    sum += estdiv(lf->coef,lf->xev,ind[i-1],ind[i],lf->nvm,lf->xl);
  sum = log(sum);
  for (i=0; i<nv; i++) lf->coef[i] -= sum;
}
