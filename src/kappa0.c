#include <stdio.h>
#include <math.h>
#include "local.h"
#include <malloc.h>

extern void preproc();
extern int cvi;

static double *fd, *ft, *lij, *d1a;
static INT par;

void assignk0(z,d,n) /* z should be n*(2*d*d+2*d+2); */
double *z;
INT d, n;
{ d1a= z; z += d*d*n;
  ft = z; z += n*(d*(d+1)+1);
  fd = z; z += n*(d+1);
}

void christ(d,nn,nl)  /* lij[i][j] = res proj. of Tij to (T1...Td) */
double nl;
INT d, nn;
{ INT i, j, k, l;
  double p4, *ll, v[1+MXDIM];
  for (i=0; i<d; i++)
    for (j=i; j<d; j++)
    { ll = &lij[(i*d+j)*nn];
      for (k=0; k<=d; k++)
      { v[k] = 0;
        for (l=0; l<nn; l++)
          v[k] += ft[k*nn+l]*ll[l];
      }
      bacT(fd,v,d+1,0,d+1);
      for (k=0; k<nn; k++)
        for (l=0; l<=d; l++)
          ll[k] -= ft[l*nn+k]*v[l];
      p4 = 0;
      for (k=0; k<=i+1; k++)
        p4 += fd[k*(d+1)+i+1]*fd[k*(d+1)+j+1];
      p4 = (fd[i+1]*fd[j+1]-p4)/(nl*nl);
      for (k=0; k<nn; k++)
        ll[k] = lij[(j*d+i)*nn+k] = ll[k] + p4*ft[k];
    }
}

void d1(n,d)   /* d1[i][j] = e_i^T (A^T A)^{-1} B_j^T */
INT n, d;
{ INT a, b, i, j;
  double *dd, v[MXDIM];
  for (i=0; i<d; i++)
  { for (j=0; j<d; j++) v[j] = 0;
    v[i] = 1;
    bacT(fd,v,d+1,1,d+1);
    for (j=0; j<d; j++)
    { dd = &d1a[(i*d+j)*n];
      for (a=0; a<n; a++)
      { dd[a] = 0;
        for (b=0; b<d; b++)
          dd[a] += v[b] * lij[(j*d+b)*n+a]; 
} } } }

void k2x(z,lf,des,kap,dv,nd)
struct tree *lf;
struct design *des;
double *z, *kap;
INT *dv, nd;
{ double det, s;
  INT i, j, k, l, d, m;
  d = lf->mi[MDIM];
  m = des->n;
if (m==1)
{ kap[0] = kap[2] = 0.0;
  return;
}
  makelxd(lf,des,z,ft,1+(d>1),dv,nd,2);
  lij = &ft[(d+1)*m];
  for (i=0; i<m; i++)
    for (j=0; j<=d; j++)
      fd[i*(d+1)+j] = ft[j*m+i];
  QR1(fd,m,d+1,NULL);
  s = 0;
  if (d>1)
  { christ(d,m,fd[0]);
    d1(m,d);
    for (j=0; j<d; j++)
      for (k=0; k<j; k++)
        for (l=0; l<m; l++)
          s += d1a[(j*d+k)*m+l]*d1a[(k*d+j)*m+l]
             - d1a[(j*d+j)*m+l]*d1a[(k*d+k)*m+l];
  }
  det = 1;
  for (j=1; j<=d; j++)
    det *= fd[j*(d+2)]/fd[0];
  kap[0] = det;
  kap[2] = s*det*fd[0]*fd[0];
}

void l1x(z,lf,des,lap,dv,nd,re)
struct tree *lf;
struct design *des;
double *z, *lap;
INT nd, *dv, re;
{ double det, t, sumcj, nu, *u, v[MXDIM];
  INT i, j, j1, k, d, m;
  d = lf->mi[MDIM]; u = des->res;
  m = des->n;
  makelxd(lf,des,z,ft,2,dv,nd,2);
  lij = &ft[(d+1)*m];
  for (i=0; i<m; i++)
  { t = ft[(re+1)*m+i];
    ft[(re+1)*m+i] = ft[d*m+i];
    ft[d*m+i] = t;
    for (j=0; j<d; j++) /* don't copy last column */
      fd[i*d+j] = ft[j*m+i];
    u[i] = ft[d*m+i];
  }
  QR1(fd,m,d,&ft[d*m]);
  bacK(fd,&ft[d*m],d);
  nu = 0;
  for (i=0; i<m; i++)
  { for (j=0; j<d; j++)
      u[i] -= ft[j*m+i]*ft[d*m+j];
    nu += u[i]*u[i];
  }    /* now u is outward vector, nu = ||u|| */
  sumcj = 0;
  for (i=0; i<d; i++) /* copy l(d-1,i) to l(re,i) */
    for (j=0; j<m; j++)
      lij[(re*d+i)*m+j] = lij[((d-1)*d+i)*m+j];
  for (j=0; j<d; j++)
  { if (j != re)
    { j1 = (j==(d-1)) ? re : j;
      for (i=0; i<d-1; i++)
      { v[i] = 0;
        for (k=0; k<m; k++)
          v[i] += lij[(i*d+j)*m+k]*u[k];
      }
      bacT(fd,v,d,1,d);
      sumcj += -v[j1];
    }
  }                                   /* stage 3,4 now complete */
  det = 1;
  for (j=1; j<d; j++)
    det *= fd[j*(d+1)]/fd[0];
  lap[0] = det;
  lap[1] = sumcj*det*fd[0]/sqrt(nu);
}

void m0x(z,lf,des,m0,dv,nd,re,rg)
struct tree *lf;
struct design *des;
double *z, *m0;
INT nd, *dv, re, rg;
{ double det, t;
  INT d, m, i, j;
  d = lf->mi[MDIM];
  m = des->n;
  makelxd(lf,des,z,ft,1,dv,nd,2);
  for (i=0; i<m; i++)
  { t=ft[(rg+1)*m+i]; ft[(rg+1)*m+i]=ft[d*m+i]; ft[d*m+i]=t;
    t=ft[(re+1)*m+i]; ft[(re+1)*m+i]=ft[(d-1)*m+i]; ft[(d-1)*m+i]=t;
    for (j=0; j<=d; j++)
      fd[i*(d+1)+j] = ft[j*m+i];
  }
  det = 1;
  QR1(fd,m,d+1,NULL);
  for (j=1; j<d-1; j++)
    det *= fd[j*(d+2)]/fd[0];
  m0[0] = det*atan2(fd[d*(d+2)],-par*fd[d*(d+1)-1]);
}

INT constants(des,lf,kap,dv,nd)
struct design *des;
struct tree *lf;
double *kap;
INT nd, *dv;
{ double h, k0[3], k1[3], l0[2], l1[2], m0[1], m1[1];
  double z[MXDIM], delt[MXDIM], mk, ml, mm;
  INT d, i, j, nnn, wt, index[MXDIM], *mi, pe, re, rg;
  if (dv==NULL) ERROR(("constants: don't provide dv=NULL"))
  mi = lf->mi;
  d = mi[MDIM];
  if (lf_error) return(0);
  if ((ident!=1) && (lf->dp[DALP]>0))
    WARN(("constants only work right for fixed h"))
  preproc(lf);
  mi[MP] = calcp(mi[MDEG],mi[MDIM],mi[MKT]);
  nnn = (ident==1) ? lf->mi[MP] : lf->mi[MN];
  checkvl(&lf->L,&lf->ll,2*nnn*(d*d+d+1));
  assignk0(lf->L,d,nnn);
  deschk(des,mi[MN],mi[MP]);
  mi[MDC] = 1;

  mk = 1.0;
  for (i=0; i<d; i++)
  { index[i] = 0;
    z[i] = lf->fl[i];
    delt[i] = (lf->fl[i+d]-z[i])/(3*mi[MMINT]);
    mk *= delt[i];
  }
  i = 0;
  
  k0[0] = k0[1] = k0[2] = 0.0;
  l0[0] = l0[1] = 0.0;
  m0[0] = 0.0;
  cvi = -1; /* avoid cross valid */
  if (mi[MIT]==IMONT)
  { for (i=0; i<mi[MMINT]; i++)
    { for (j=0; j<d; j++) z[j] = lf->fl[j]+(lf->fl[j+d]-lf->fl[j])*runif();
      h = nbhd(z,lf,des,lf->dp[DALP],lf->dp[DFXH],0);
      locfit(lf,des,z,h,1);
      k2x(z,lf,des,k1,dv,nd);
      k0[0] += k1[0];
    }
    for (j=0; j<d; j++) k0[0] *= lf->fl[j+d]-lf->fl[j];
    kap[0] = k0[0]/mi[MMINT];
    return(1);
  }
  while(1)
  { wt = 1;
    for (i=0; i<d; i++)
      wt *= (4-2*(index[i]%2==0)-(index[i]==0)-(index[i]==mi[MMINT]));
    h = nbhd(z,lf,des,lf->dp[DALP],lf->dp[DFXH],0);
    locfit(lf,des,z,h,1);
    k2x(z,lf,des,k1,dv,nd);
    k0[0] += wt*mk*k1[0];
    k0[2] += wt*mk*k1[2];

    for (re=0; re<d; re++) if ((index[re]==0) | (index[re]==mi[MMINT]))
    { l1x(z,lf,des,l1,dv,nd,re);
      ml = 1;
      for (i=0; i<d; i++) if (i!=re) ml *= delt[i];
      pe = 1-2*(index[re]==0);
      l0[0] += wt*ml*l1[0];
      l0[1] += wt*ml*pe*l1[1];

      for (rg=re+1; rg<d; rg++) if ((index[rg]==0) | (index[rg]==mi[MMINT]))
      { par = pe*(1-2*(index[rg]==0));
        m0x(z,lf,des,m1,dv,nd,re,rg);
        mm = 1;
        for (i=0; i<d; i++) if ((i!=re) & (i!=rg)) mm *= delt[i];
        m0[0] += wt*mm*m1[0];
      }
    }

    /* compute next grid point */
    for (i=0; i<d; i++)
    { index[i]++;
      z[i] = lf->fl[i]+3*delt[i]*index[i];
      if (index[i]>mi[MMINT])
      { index[i] = 0;
        z[i] = lf->fl[i];
        if (i==d-1) /* done */
        { kap[0] = k0[0];
          kap[1] = l0[0]/2;
printf("%8.5f %8.5f\n",kap[0],kap[1]);
          if (d==1) return(2);
          k0[2] = -k0[2] - d*(d-1)*k0[0]/2;
printf("k0: %8.5f  k2: %8.5f\n",k0[0],k0[2]);
printf("l0: %8.5f  l1: %8.5f\n",l0[0],l1[1]);
printf("m0: %8.5f\n",m0[0]);
printf("check: %8.5f\n",(k0[0]+k0[2]+l0[1]+m0[0])/(2*PI));
          kap[2] = (k0[2]+l0[1]+m0[0])/(2*PI);
          return(3);
        }
      }
      else i = d;
    }
  }
}

double tailp(c,k0,m,d,nu)
double c, *k0, nu;
INT m, d;
{ INT i;
  double p;
  p = 0;
  if (nu==0)
    for (i=0; i<m; i++)
      p += k0[i]*exp(LGAMMA((d+1-i)/2.0)-(d+1-i)*0.5723649429)
          *(1-pchisq(c*c,(double) d+1-i));
  else
    for (i=0; i<m; i++)
      p += k0[i]*exp(LGAMMA((d+1-i)/2.0)-(d+1-i)*0.5723649429)
          *(1-pf(c*c/(d+1-i),(double) (d+1-i), nu));
  return(p);
}

double taild(c,k0,m,d,nu)
double c, *k0, nu;
INT m, d;
{ double p;
  INT i;
  p = 0;
  if (nu==0)
    for (i=0; i<m; i++)
      p += k0[i]*exp(LGAMMA((d+1-i)/2.0)-(d+1-i)*0.5723649429)
          *2*c*dchisq(c*c,(double) (d+1-i));
  else
    for (i=0; i<m; i++)
      p += k0[i]*exp(LGAMMA((d+1-i)/2.0)-(d+1-i)*0.5723649429)
          *2*c*df(c*c/(d+1-i),(double) (d+1-i), nu)/(d+1-i);
  return(-p);
}

double cv(k0,m,d,al,it,s,nu)
double *k0, al, nu;
INT m, d, it, s;
{ double c;
  INT j;
  if ((m<0) | (m>d+1)) ERROR(("cv: invalid no. of terms %d",m))
  if ((al<=0) | (al>=1)) ERROR(("cv: invalid alpha %8.5f",al))
  if (lf_error) return(0.0);
  if (al>0.5) WARN(("cv: A mighty large tail probability al=%8.5f",al))
  if (s==1) al = 2*al;
  if (m==0) { d = 0; k0[0] = 1; m = 1; }
  c = 2.0;
  for (j=0; j<it; j++)
    c -= (tailp(c,k0,m,d,nu)-al)/taild(c,k0,m,d,nu);
  return(c);
}

double cvc(k0,m,d,al,it,s,nu,d0,d1)
double *k0, al, nu, d0, d1;
INT m, d, it, s;
{ double c, c1, c2;
  INT j;
  if ((m<0) | (m>d+1)) ERROR(("cvc: invalid no. of terms %d",m))
  if ((al<=0) | (al>=1)) ERROR(("cvc: invalid alpha %8.5f",al))
  if (lf_error) return(0.0);
  if (al>0.5) WARN(("cvc: A mighty large tail probability al=%8.5f",al))
  if (m==0) { d = 1; k0[0] = 1; m = 1; }
  if (s==1) al = 2*al;
  c = sqrt(-2*log(al*3.14159265/k0[0]))+d0;
  for (j=0; j<it; j++)
  { c1 = c-d0;
    c2 = c+d0-d1*d1/c;
    if (c2<c1)
      c -= (tailp(c1,k0,m,d,nu)-al)/taild(c1,k0,m,d,nu);
    else
      c -= (tailp(c1,k0,m,d,nu)+tailp(c2,k0,m,d,nu)-2*al) /
           (taild(c1,k0,m,d,nu)+taild(c2,k0,m,d,nu)*(1+d1*d1/(2*c*c)));
  }
  return(c);
}
