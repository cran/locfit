/*
 *   Copyright (c) 1998-1999 Lucent Technologies.
 *   See README file for details.
 */

/*
  functions for handling locfit variables. Note the definition
  of the vari type, and some functions, such as createvar(),
  differ for different versions of locfit.
*/

#include "local.h"

#define MAXV 1000
#define LF_WORK 1024

#ifdef CVERSION
static char *db;
static INT lfwptr, lf_work;
vari root;
#endif

double *vdptr(v)
vari *v;
{ if (v==NULL) return(NULL);
  return(v->dpr);
}

#ifdef CVERSION
void initdb() /* initialize locfit's work space */
{ char *z;
  z = getenv("LFWORK");
  if (z==NULL) lf_work = LF_WORK;
    else sscanf(z,"%d",&lf_work);
  lf_work <<= 10;
  db = (char *)malloc(lf_work);
  root.stat = STSYSTEM;
  root.mode = VVARI;
  root.dpr = (double *)db;
  lfwptr = root.len = MAXV*sizeof(vari);
  root.n = 0;
}

vari *growvar(vold,n)
vari *vold;
INT n;
{ vari *vnew;
  INT rlen;
  if (vold==NULL)
  { ERROR(("growvar: NULL old"));
    return(NULL);
  }
  rlen = vlen(n,vold->mode);
  if (rlen <= vold->len) return(vold);
  vnew = createvar("_grow",vold->stat,n,vold->mode);
  memcpy(vdptr(vnew),vdptr(vold),vlen(vold->n,vold->mode));
  setvarname(vnew,vold->name);
  vnew->n = vold->n;
  deletevar(vold);
  return(vnew);
}

INT vlen(n,mode)
INT n, mode;
{ switch(mode)
  { case VDOUBLE: return(n*sizeof(double));
    case VINT:    return(n*sizeof(INT));
    case VCHAR:   return(n);
    case VARGL:   return(n*sizeof(carg));
    case VPREP:   return(sizeof(pplot));
    case VARC:    return(n*sizeof(arstruct));
    case VVARI:   return(n*sizeof(vari));
    case VXYZ:    return(n*sizeof(plxyz));
  }
  ERROR(("unknown mode in vlen"));
  return(0);
}

void *viptr(v,i) /* return pointer to ith data item, take account of mode */
vari *v;
INT i;
{ switch(v->mode)
  { case VDOUBLE: return(&v->dpr[i]);
    case VCHAR: return(&((char *)v->dpr)[i]);
    case VARGL: return(&((carg *)v->dpr)[i]);
    case VARC:  return(&((arstruct *)v->dpr)[i]);
    case VVARI: return(&((vari *)v->dpr)[i]);
    case VXYZ:  return(&((plxyz *)v->dpr)[i]);
  }
  ERROR(("Unknown mode %d in viptr",v->mode));
  return(NULL);
}

void setvarname(v,name)
vari *v;
varname name;
{ if (strcmp(v->name,name)==0) return;
  deletename(name);
  strcpy(v->name,name);
}

/*
  findvar finds the variable name.
  err=0, keep quiet if not found; 1 produce error message.
  *n returns length of variable (if initially>0, checks length)
*/

vari *findvar(name,err,n)
varname name;
INT err, *n;
{ INT i, status;
  vari *v;
  if (strcmp(name,"_NuLl")==0) return(NULL);
  for (i=0; i<root.n; i++)
  { v = viptr(&root,i);
    status = v->stat;
    if ((strcmp(v->name,name)==0) &&
      ((status!=STHIDDEN)&(status!=STEMPTY)) )
    { if (n==NULL) return(v);
      if (*n==-1) *n = v->n;
      if ((*n==0) | (*n==v->n)) return(v);
      if (err) ERROR(("Variable %s has wrong length",name));
      return(NULL);
    }
  }
  if (err) ERROR(("Variable %s not found",name));
  return(NULL);
}

void deletevar(v) /* delete variable, or top variable if NULL */
vari *v;
{ if (root.n==0) return;
  if (v!=NULL) v->stat = STEMPTY;
  if ((v==NULL) || (v==viptr(&root,root.n-1))) /* top variable */
  { root.n--;
    lfwptr -= ((vari *)viptr(&root,root.n))->len;
  }
}

void deleteifhidden(v)
vari *v;
{ if (v==NULL) return;
  if (v->stat == STHIDDEN) deletevar(v);
}

void deletename(name) /* delete variable name, or top variable if NULL */
varname name;
{ vari *v;
  v = findvar(name,0,NULL);
  if (v!=NULL) deletevar(v);
}

vari *createvar(name,status,n,mode)
varname name;
INT status, n, mode;
{ INT i, len;
  vari *v;
  len = vlen(n,mode);
  while (len%8>0) len++;
  if (lf_error) return(NULL);
  if ((status==STSYSTEM)|(status==STREGULAR)|(status==STPLOTVAR))
    deletename(name);

  if (status!=STREADFI)
  { for (i=0; i<root.n; i++) /* does unused variable have space allocated? */
    { v = viptr(&root,i);
      if ((v->stat == STEMPTY) && (v->len >= len))
      { strcpy(v->name,name);
        v->n = n;
        v->stat = status;
        v->mode = mode;
        return(v);
      }
    }
  }

  /* must allocate next variable. First, is there space? */
  if (root.n==MAXV) ERROR(("Too many variables"));
  if ((status!=STSYSPEC) && (lfwptr+len>lf_work))
    ERROR(("Insufficient space for variable creation"));
  if (lf_error) return(NULL);

  v = viptr(&root,root.n);
  strcpy(v->name,name);
  v->n = n;
  v->stat = status;
  v->len = len;
  v->mode = mode;
  if (status!=STSYSPEC)
  { v->dpr = (double *)(&db[lfwptr]);
    lfwptr += len;
  }
  root.n++;
  return(v);
}
#else
void *viptr(v,i)
vari *v;
INT i;
{ return(&v->dpr[i]);
}

vari *createvar(name,type,n,mode)
varname name;
INT type, n, mode;
{ vari *v;
  v = (vari *)malloc(sizeof(vari));
  v->n = n;
  switch(mode)
  { case VDOUBLE:
      v->dpr = (double *)calloc(n,sizeof(double));
      break;
    case VINT:
      v->dpr = (double *)calloc(n,sizeof(INT));
      break;
  }
  return(v);
}
#endif

vari *checkvarlen(v,n,name,mode)
vari *v;
INT n, mode;
varname name;
{ 
#ifndef CVERSION
  if ((v!=NULL) && (v->n >= n)) return(v);
#endif
  return(createvar(name,STSYSTEM,n,mode));
}
