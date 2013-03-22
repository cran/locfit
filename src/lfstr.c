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

void setstrval(mi,v,z)
INT *mi, v;
char *z;
{ INT quasi;
  char *zz;
  switch(v)
  { case MKER:
      switch(z[0])
      { case 'r': mi[v] = WRECT; return;
        case 'e': mi[v] = WEPCH; return;
        case 'b': mi[v] = WBISQ; return;
        case 'g': mi[v] = WGAUS; return;
        case 't': if (z[1]=='r')
                    mi[v] = (z[2]=='i') ? WTRIA : WTRWT;
                  else mi[v] = WTCUB;
                  return;
        case 'q': mi[v] = WQUQU; return;
        case '6': mi[v] = W6CUB; return;
        case 'm': mi[v] = WMINM; return;
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
      quasi = (z[0]=='q');
      zz = &z[quasi];
      switch(zz[0])
      { case 'd': mi[v] = TDEN; return;
        case 'r':
          mi[v] = (zz[1]=='o') ? TROBT : TRAT;
          break;
        case 'h': mi[v] = THAZ; return;
        case 'b': mi[v] = TLOGT; break;
        case 'c': mi[v] = TCIRC; break;
        case 'p': mi[v] = TPOIS; break;
        case 'w': mi[v] = TWEIB; break;
        case 'g':
          if (zz[1]=='e') mi[v] = TGEOM;
            else
              mi[v] = (stm(zz,"gau",3)) ? TGAUS : TGAMM;
        break;
        default:
          WARN(("unknown family %s",z));
          mi[v] = TGAUS;
      }
      if (quasi) mi[v] += 64;
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
      }
      WARN(("unknown evaluation structure %s",z));
      mi[v] = ETREE;
      return;

    case MACRI:
      switch(z[0])
      { case 'c': mi[v] = ACP; return;
        case 'k': mi[v] = AKAT;return;
        case 'm': mi[v] = AMDI;return;
      }
      WARN(("unknown adaptive criterion %s",z));
      mi[v] = ACP;
      return;
  }

  WARN(("setstrval: invalid value %d",v));
  return;
}
