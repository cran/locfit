/*
  strval() function for converting string arguments to Locfit's
  numeric values.
  a typical call will be setstrval(lf.mi,MKER,"gauss")

  components that can be set in this manner are
    MKER  (weight function)
    MKT   (kernel type -- spherical or product)
    MTG   (local likelihood family)
    MLINK (link function)
    MIT   (integration type for density estimation)
    MEV   (evaluation structure)

  INT ppwhat(str) interprets the preplot what argument.
  INT restyp(str) interprets the residual type argument.

  Also includes other misc. string functions.
*/

#include "local.h"

INT stm(u,v,k)
char *u, *v;
INT k;
{ int i;
  for (i=0; i<k; i++) if (u[i]!=v[i]) return(0);
  return(1);
}

INT lffamily(z)
char *z;
{ INT quasi, robu, f;
  quasi = robu = 0;
  while ((z[0]=='q') | (z[0]=='r'))
  { quasi |= (z[0]=='q');
    robu  |= (z[0]=='r');
    z++;
  }
  switch(z[0])
  { /* 'a' and 'o' are for old 'rate' and 'robu'. */
    case 'a':
      f = TRAT; robu = 0;
      break;
    case 'o':
      f = TROBT; robu = 0;
      break;
    case 'd': f = TDEN; break;
    case 'h':
      f = THAZ;
      if (z[1]=='u') f = TROBT;
      break;
    case 'b': f = TLOGT; break;
    case 'c':
      f = TCIRC;
      if (z[1]=='a') f = TCAUC;
      break;
    case 'p':
      f = TPOIS;
      if (z[1]=='r') f = TPROB;
      break;
    case 'w': f = TWEIB; break;
    case 'g':
      if (z[1]=='e') f = TGEOM;
      else
        f = (stm(z,"gau",3)) ? TGAUS : TGAMM;
      break;
    default:
      WARN(("unknown family %s",z));
      f = TGAUS;
  }
  if (quasi) f += 64;
  if (robu)  f += 128;
  return(f);
}

void getlffam(z,x)
char **z;
INT *x;
{ *x = lffamily(z[0]);
}

void setstrval(mi,v,z)
INT *mi, v;
char *z;
{ switch(v)
  { case MKER:
      switch(z[0])
      { case 'r': mi[v] = WRECT; return;
        case 'e': mi[v] = (z[1]=='x') ? WEXPL : WEPAN; return;
        case 'b': mi[v] = WBISQ; return;
        case 'g': mi[v] = WGAUS; return;
        case 't': if (z[1]=='r')
                    mi[v] = (z[2]=='i') ? WTRIA : WTRWT;
                  else mi[v] = WTCUB;
                  return;
        case 'q': mi[v] = WQUQU; return;
        case '6': mi[v] = W6CUB; return;
        case 'm': mi[v] = (z[1]=='c') ? WMACL : WMINM;
                  return;
        case 'p': mi[v] = WPARM; return;
      }
      WARN(("Unknown weight function %s; using default",z));
      mi[v] = WTCUB;
      return;

    case MKT:
      switch(z[0])
      { case 's': mi[v] = KSPH;  return;
        case 'p': mi[v] = KPROD; return;
        case 'c': mi[v] = KCE;   return;
      }

    case MTG:
      mi[v] = lffamily(z);
      return;

    case MLINK:
      switch(z[0])
      { case 'd' : mi[v] = LDEFAU; return;
        case 'c' : mi[v] = LCANON; return;
        case 'i' : mi[v] = (z[1]=='n') ? LINVER : LIDENT; return;
        case 'l' : if (stm(z,"logi",4)) mi[v] = LLOGIT;
                                   else mi[v] = LLOG;
                   return;
        case 's' : mi[v] = LSQRT; return;
        case 'a' : mi[v] = LASIN; return;
      }
      WARN(("unknown link %s",z));
      mi[v] = LDEFAU;
      return;

    case MIT:
      switch(z[0])
      { case 'p': mi[v] = IPROD; return;
        case 'd': mi[v] = IDEFA; return;
        case 'h':
          mi[v] = IHAZD;
          if (stm(z,"har",3)) mi[v] = IHARD;
          return;
        case 'm':
          mi[v] = IMULT;
          if (z[1]=='l') mi[v] = IMLIN;
          if (z[1]=='o') mi[v] = IMONT;
          return;
      }
      WARN(("unknown integration type %s",z));
      mi[v] = IDEFA;
      return;

    case MEV:
      switch(z[0])
      { case 't': mi[v] = ETREE; return;
        case 'p': mi[v] = EPHULL;return;
        case 'd': mi[v] = EDATA; return;
        case 'g': mi[v] = EGRID; return;
        case 'k': mi[v] = (stm(z,"kdc",3)) ? EKDCE : EKDTR;
                  return;
        case 'c': mi[v] = ECROS; return;
        case 'x': mi[v] = EXBAR; return;
        case 'n': mi[v] = ENONE; return;
      }
      ERROR(("unknown evaluation structure %s",z));
      mi[v] = ETREE;
      return;

    case MACRI:
      switch(z[0])
      { case 'n': mi[v] = ANONE; return;
        case 'c': mi[v] = ACP; return;
        case 'i': mi[v] = AKAT;return;
        case 'm': mi[v] = AMDI;return;
        case 'o': mi[v] = AOK; return;
      }
      WARN(("unknown adaptive criterion %s",z));
      mi[v] = 0;
      return;
  }

  WARN(("setstrval: invalid value %d",v));
  return;
}

INT restyp(z)
char *z;
{
printf("restyp: %s\n",z); 
  switch(z[0])
  { case 'd':
      if (z[1]=='2') return(RDEV2);
                else return(RDEV);
    case 'p': return(RPEAR);
    case 'r': return(RRAW);
    case 'l':
      if ((z[1]=='d') && (z[2]=='d')) return(RLDDT);
         else return(RLDOT);
    case 'f': return(RFIT);
    case 'm': return(RMEAN);
  }
  ERROR(("Unknown type= argument"));
  return(0);
}

INT ppwhat(z)
char *z;
{
  switch(z[0])
  { case 'c': return(PCOEF);
    case 'n': return(PNLX);
    case 'i': return(PT0);
    case 'b': return(PBAND);
    case 'd': return(PDEGR);
    case 'l': return(PLIK);
    case 'r': return(PRDF);
    case 'v': return(PVARI);
  }
  ERROR(("Unknown what= argument"));
  return(0);
}
