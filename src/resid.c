/*
 *   Copyright (c) 1998 Lucent Technologies.
 *   See README file for details.
 */

/*
  Functions for computing residuals and fitted values from
  the locfit object.

  resid(y,c,w,th,mi,ty) computes residual for response y; fitted
    value th; censoring indicator c; prior weight w. ty controls
    the type of residual.
  compres(...) is SVERSION front end to resid; y etc are vectors.
  fitted(v,ty) is CVERSION front end, interpreting command
    line arguments and computing th.
  vfitted() and vresid() for use by arithmetic interpreter.
*/

#include <math.h>
#include <stdlib.h>
#include "local.h"

extern struct tree lf;
extern struct design des;

double resid(y,c,w,th,mi,ty)
INT *mi, ty;
double y, c, w, th;
{ double res[LLEN], raw;
  links(th,y,mi[MTG],mi[MLINK],res,c,w);

  if (((mi[MTG]&63)==TGAUS) || ((mi[MTG]&63)==TROBT))
    raw = y-res[ZMEAN];
  else
    raw = y-w*res[ZMEAN];
  switch(ty)
  { case 1: /* deviance */
      if (res[ZDLL]>0) return(sqrt(-2*res[ZLIK]));
            else return(-sqrt(-2*res[ZLIK]));
    case 2: /* pearson */
      if (res[ZDDLL]<=0)
      { if (res[ZDLL]==0) return(0);
        return(NOSLN);
      }
      return(res[ZDLL]/sqrt(res[ZDDLL]));
    case 3: /* raw */
      return(raw);
    case 4: /* ldot */
      return(res[ZDLL]);
    case 5: /* fitted value */
      return(res[ZMEAN]);
    case 6: /* dev^2 */
      return(-2*res[ZLIK]);
    case 7: /* -ddot{l} */
      return(res[ZDDLL]);
    default: ERROR(("resid: unknown residual type %d",ty))
  }
  return(0.0);
}

void compres(y,c,w,th,mi,ty,m)
double *y, *w, *th;
INT *c, *mi, *ty, *m;
{ INT i;
  for (i=0; i<*m; i++)
    th[i] = resid(y[i],(double)c[i],w[i],th[i],mi,*ty);
}

#ifdef CVERSION
vari *vfitted()
{ vari *v;
  INT i, n;
  n = lf.mi[MN];
  v = createvar("vfitted",STHIDDEN,n,VDOUBLE);
  recondat(1,&n);
  if (lf_error) return(NULL);

  intv(&lf,&des,lf.x,vdptr(v),n,PCOEF,NULL,0);
  return(v);
}

vari *vresid()
{ vari *vf, *vr;
  INT i, n, ty;
  n = lf.mi[MN];
  vr = createvar("vresid",STHIDDEN,n,VDOUBLE);
  vf = vfitted();
  if (lf_error) return(NULL);

  for (i=0; i<n; i++)
    vassn(vr,i,resid(resp(&lf,i),cens(&lf,i),prwt(&lf,i),vitem(vf,i),lf.mi,1));
  return(vr);
}

void fitted(v,ty)
vari *v;
INT ty;
{ double *f, *infl, ldt;
  vari *vr;
  INT i, j, n, cv;

  i = getarg(v,"type",1); if (i>0) switch(argval(v,i)[0])
  { case 'd': ty = (argval(v,i)[1]=='2') ? 6 : 1; break;
    case 'p': ty = 2; break;
    case 'r': ty = 3; break;
    case 'l': ty = 4; break;
    case 'f': ty = 5; break;
    default: ERROR(("Unknown residual type %s",argval(v,i)));
  }

  i = getarg(v,"cv",1); cv = (i>0) ? getlogic(v,i) : 0;

  recondat(ty==5,&n);
  if (lf_error) return;

  vr = createvar("fitted",STHIDDEN,n,VDOUBLE);
  f = vdptr(vr);
  intv(&lf,&des,lf.x,f,n,PCOEF,NULL,0);
  if (cv)
  { infl = vdptr(createvar("_duminfl",STHIDDEN,n,VDOUBLE));
    intv(&lf,&des,lf.x,infl,n,PT0,NULL,0);
    for (i=0; i<n; i++)
    { ldt = resid(resp(&lf,i),cens(&lf,i),prwt(&lf,i),f[i],lf.mi,4);
      f[i] -= infl[i]*ldt;
    }
  }

  for (i=0; i<n; i++)
    f[i] = resid(resp(&lf,i),cens(&lf,i),prwt(&lf,i),f[i],lf.mi,ty);
  saveresult(vr,argarg(v,0),STREGULAR);
}
#endif
