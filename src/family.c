#include <math.h>
#include <stdio.h>
#include "local.h"

#ifdef RVERSION
extern double pgamma(), pbeta();
#define igamma(x,p) pgamma(x,p,1.0)
#define ibeta(x,a,b) pbeta(x,a,b)
#endif

extern double robscale;

INT famdens(mean,th,link,res,cens,w)
double mean, th, *res, w;
INT link, cens;
{ if (cens)
    res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
  else
  { res[ZLIK] = w*th;
    res[ZDLL] = res[ZDDLL] = w;
  }
  return(0);
}

INT famgaus(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double z, pz, dp;
  if (link==LINIT)
  { res[ZDLL] = w*y;
    return(0);
  }
  z = y-mean;
  if (cens)
  { pz = pnorm(-z,0.0,1.0);
    dp = ((z>6) ? ptail(-z) : exp(-z*z/2)/pz)/2.5066283;
    res[ZLIK] = w*log(pz);
    res[ZDLL] = w*dp;
    res[ZDDLL]= w*dp*(dp-z);
    return(0);
  }
  res[ZDLL] = w*z;           /*  T'(th)Y-psi'(th) */
  res[ZLIK] = -w*z*z/2;      /* log-likelihood */
  res[ZDDLL]= w;         /*-(T"(th)Y-psi"(th) */
  return(0);                /* 0 for success, 1 for fail */
}

INT famrobu(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double z, sw;
  if (link==LINIT)
  { res[ZDLL] = w*y;
    return(0);
  }
  sw = (w==1.0) ? 1.0 : sqrt(w); /* don't want unnecess. sqrt! */
  z = sw*(y-mean)/robscale;
  res[ZLIK] = robscale*robscale*((fabs(z)<1) ? -z*z/2 : 0.5-fabs(z));
  if (z< -1.0)
  { res[ZDLL] = -sw*robscale;
    res[ZDDLL]= 0.0;
    return(0);
  }
  if (z> 1.0)
  { res[ZDLL] = sw*robscale;
    res[ZDDLL]= 0.0;
    return(0);
  }
  res[ZDLL] =  sw*z*robscale;
  res[ZDDLL] = w;
  return(0);
}

INT fambino(y,p,th,link,res,cens,w)
double y, p, th, *res, w;
INT link, cens;
{ double wp;
  if (link==LINIT)
  { res[ZDLL] = y;
    return(0);
  }
  wp = w*p;
  if (link==LIDENT)
  { if ((p<=0) && (y>0)) return(1);
    if ((p>=1) && (y<w)) return(1);
    res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
    if (y>0)
    { res[ZLIK] += y*log(wp/y);
      res[ZDLL] += y/p;
      res[ZDDLL]+= y/(p*p);
    }
    if (y<w)
    { res[ZLIK] += (w-y)*log((w-wp)/(w-y));
      res[ZDLL] -= (w-y)/(1-p);
      res[ZDDLL]+= (w-y)/SQR(1-p);
    }
    return(0);
  }
  if (link==LLOGIT)
  { if ((y<0) | (y>w)) /* goon observation; delete it */
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
      return(0);
    }
    res[ZLIK] = (th<0) ? th*y-w*log(1+exp(th)) : th*(y-w)-w*log(1+exp(-th));
    if (y>0) res[ZLIK] -= y*log(y/w);
    if (y<w) res[ZLIK] -= (w-y)*log(1-y/w);
    res[ZDLL] = (y-wp);
    res[ZDDLL]= wp*expit(-th);
    return(0);
  }
  ERROR(("link %d invalid for binomial family",link))
  return(1);
}

INT fampois(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double wmu, pt, dp, dq;
  if (link==LINIT)
  { res[ZDLL] = y;
    return(0);
  }
  wmu = w*mean;
  if (cens)
  { if (y<=0)
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
      return(0);
    }
    pt = igamma(wmu,y);
    dp = exp((y-1)*log(wmu)-wmu-LGAMMA(y))/pt;
    dq = dp*((y-1)/wmu-1);
    res[ZLIK] = log(pt);
    if (link==LLOG)
    { res[ZDLL] = dp*wmu;
      res[ZDDLL]= -(dq-dp*dp)*wmu*wmu-dp*wmu;
      return(0);
    }
    if (link==LIDENT)
    { res[ZDLL] = dp*w;
      res[ZDDLL]= -(dq-dp*dp)*w*w;
      return(0);
    }
    if (link==LSQRT)
    { res[ZDLL] = dp*2*w*th;
      res[ZDDLL]= -(dq-dp*dp)*(4*w*w*mean)-2*dp*w;
      return(0);
  } }
  if (link==LLOG)
  { res[ZLIK] = res[ZDLL] = y-wmu;
    if (y>0) res[ZLIK] += y*(th-log(y/w));
    res[ZDDLL] = wmu;
    return(0);
  }
  if (link==LIDENT)
  { if ((mean<=0) && (y>0)) return(1);
    res[ZLIK] = y-wmu;
    res[ZDLL] = -w;
    res[ZDDLL] = 0;
    if (y>0)
    { res[ZLIK] += y*log(wmu/y);
      res[ZDLL] += y/mean;
      res[ZDDLL]= y/(mean*mean);
    }
    return(0);
  }
  if (link==LSQRT)
  { if ((mean<=0) && (y>0)) return(1);
    res[ZLIK] = y-wmu;
    res[ZDLL] = -2*w*th;
    res[ZDDLL]= 2*w;
    if (y>0)
    { res[ZLIK] += y*log(wmu/y);
      res[ZDLL] += 2*y/th;
      res[ZDDLL]+= 2*y/mean;
    }
    return(0);
  }
  ERROR(("link %d invalid for Poisson family",link))
  return(1);
}

INT famgamm(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double pt, dg;
  if (link==LINIT)
  { res[ZDLL] = y;
    return(0);
  }
  if ((mean<=0) & (y>0)) return(1);
  if (cens)
  { if (y<=0)
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
      return(0);
    }
    if (link==LLOG)
    { pt = 1-igamma(y/mean,w);
      dg = exp((w-1)*log(y/mean)-y/mean-LGAMMA(w));
      res[ZLIK] = log(pt);
      res[ZDLL] = y*dg/(mean*pt);
      res[ZDDLL]= dg*(w*y/mean-y*y/(mean*mean))/pt+SQR(res[ZDLL]);
      return(0);
    }
    if (link==LINVER)
    { pt = 1-igamma(th*y,w);
      dg = exp((w-1)*log(th*y)-th*y-LGAMMA(w));
      res[ZLIK] = log(pt);
      res[ZDLL] = -y*dg/pt;
      res[ZDDLL]= dg*y*((w-1)*mean-y)/pt+SQR(res[ZDLL]);
      return(0);
    }
  }
  else
  { if (y<0) WARN(("Negative Gamma observation"))
    if (link==LLOG)
    { res[ZLIK] = -y/mean+w*(1-th);
      if (y>0) res[ZLIK] += log(y/w);
      res[ZDLL] = y/mean-w;
      res[ZDDLL]= y/mean;
      return(0);
    }
    if (link==LINVER)
    { res[ZLIK] = -y/mean+w-w*log(mean);
      if (y>0) res[ZLIK] += log(y);
      res[ZDLL] = -y+w*mean;
      res[ZDDLL]= w*mean*mean;
      return(0);
    }
    if (link==LIDENT)
    { res[ZLIK] = -y/mean+w-w*log(mean);
      res[ZDLL] = (y-mean)/(mean*mean);
      res[ZDDLL]= w/(mean*mean);
      return(0);
    }
  }
  ERROR(("link %d invalid for Gamma family",link))
  return(1);
}

INT famgeom(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double p, pt, dp, dq;
  if (link==LINIT)
  { res[ZDLL] = y;
    return(0);
  }
  p = 1/(1+mean);
  if (cens) /* censored observation */
  { if (y<=0)
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0;
      return(0);
    }
    pt = 1-ibeta(p,w,y);
    dp = -exp(LGAMMA(w+y)-LGAMMA(w)-LGAMMA(y)+(y-1)*th+(w+y-2)*log(p))/pt;
    dq = ((w-1)/p-(y-1)/(1-p))*dp;
    res[ZLIK] = log(pt);
    res[ZDLL] = -dp*p*(1-p);
    res[ZDDLL]= (dq-dp*dp)*p*p*(1-p)*(1-p)+dp*(1-2*p)*p*(1-p);
    res[ZDDLL]= -res[ZDDLL];
    return(0);
  }
  else
  { res[ZLIK] = (y+w)*log((y/w+1)/(mean+1));
    if (y>0) res[ZLIK] += y*log(w*mean/y);
    if (link==LLOG)
    { res[ZDLL] = (y-w*mean)*p;
      res[ZDDLL]= (y+w)*p*(1-p);
      return(0);
    }
    if (link==LIDENT)
    { res[ZDLL] = (y-w*mean)/(mean*(1+mean));
      res[ZDDLL]= w/(mean*(1+mean));
      return(0);
    }
  }
  ERROR(("link %d invalid for geometric family",link))
  return(1);
}

INT famweib(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double yy;
  yy = pow(y,w);
  if (link==LINIT)
  { res[ZDLL] = yy;
    return(0);
  }
  if (cens)
  { res[ZLIK] = -yy/mean;
    res[ZDLL] = res[ZDDLL] = yy/mean;
    return(0);
  }
  res[ZLIK] = 1-yy/mean-th;
  if (yy>0) res[ZLIK] += log(w*yy);
  res[ZDLL] = -1+yy/mean;
  res[ZDDLL]= yy/mean;
  return(0);
}

INT famcirc(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ if (link==LINIT)
  { res[ZDLL] = w*sin(y);
    res[ZLIK] = w*cos(y);
    return(0);
  }
  res[ZDLL] = w*sin(y-mean);
  res[ZDDLL]= w*cos(y-mean);
  res[ZLIK] = res[ZDDLL]-w;
  return(0);
}

INT links(th,y,fam,lin,res,cd,w) /* the link and various related functions */
double th, y, *res, w, cd;
INT fam, lin;
{ double mean;
  INT c, link;
  c = (INT)cd; link = (INT)lin;
  switch(link)
  { case LIDENT: mean = res[ZMEAN] = th; break; /* mean */
    case LLOG:   mean = res[ZMEAN] = exp(MIN(th,300)); break;
    case LLOGIT: mean = res[ZMEAN] = expit(th); break;
    case LINVER: mean = res[ZMEAN] = 1/th; break;
    case LSQRT:  mean = res[ZMEAN] = th*fabs(th); break;
    case LINIT:  mean = 0; break;
    default: ERROR(("links: unknown link %d",link)) return(1);
  }
  switch(fam&63)
  { case THAZ:
    case TDEN:
    case TRAT: return(famdens(mean,th,link,res,c,w));
    case TGAUS: return(famgaus(y,mean,th,link,res,c,w));
    case TLOGT: return(fambino(y,mean,th,link,res,c,w));
    case TPOIS: return(fampois(y,mean,th,link,res,c,w));
    case TGAMM: return(famgamm(y,mean,th,link,res,c,w));
    case TGEOM: return(famgeom(y,mean,th,link,res,c,w));
    case TWEIB: return(famweib(y,mean,th,link,res,c,w));
    case TCIRC: return(famcirc(y,mean,th,link,res,c,w));
    case TROBT:
      return(famrobu(y,mean,th,link,res,c,w));
  }
  ERROR(("links: invalid family %d",fam)) return(1);
}
