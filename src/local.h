/*
   Which version of locfit is being built?
   One of the following should be uncommented.
   If you change this after building one version
   of Locfit, then you'll also need to
   touch *.c
   at the unix prompt.
   For R compilation, both RVERSION and SVERSION should be uncommented.
   INTERFACE (not CVERSION) should be used, for user-written interfaces.
*/

/* #define CVERSION */
/* #define RVERSION */
#define SVERSION
/* #define INTERFACE */


/*
   Some older math libraries have no lgamma() function, and gamma(arg)
   actually returns log(gamma(arg)). If so, you need to change
   LGAMMA macro below.

   If all else fails, you can also use lflgamma().

   uncomment the definitions for erf, erfc and daws only if your
   math libraries don't include these functions.
*/
#define LGAMMA(arg) lgamma(arg)
/* #define erf(x) lferf(x)   */
/* #define erfc(x) lferfc(x) */
#define daws(x) lfdaws(x)

/*
   Does your system support popen() and pclose()
   for pipes? For most flavours of unix, yes.
   For other OS's, usually not, and you'll need to
   uncomment the following line.

   (This is only relevant if making the C version).
*/
/* #define NOPIPES */


/*
   the INT type is used for all integers provided in .C() calls
   from S.
   For the S version on 64 bit systems, this should be long int.
   For the R version, should be int.
   For the C version, either is adequate.
*/

typedef long int INT;

/*
  DIRSEP: '/' for unix; '\\' for DOS
*/
#define DIRSEP '/'

/******** NOTHING BELOW HERE NEEDS CHANGING **********/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "lfcons.h"
#include "lfstruc.h"
#include "lffuns.h"

#ifdef CVERSION
#undef printf
#define printf lfprintf
extern int lfprintf(const char *format, ...);
extern int printe(const char *format, ...);
#endif

#ifdef SVERSION
#define printe printf
#endif

#ifdef INTERFACE
#define printe printf
#endif

#define ERROR(args) printe("Error: "), printe args , printe("\n"), lf_error=1
#define WARN(args)  printe("Warning: "),printe args, printe("\n")

#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SGN(x) (((x)>0) ? 1 : -1)
#define SQR(x) ((x)*(x))
#define FSWAP(a,b) { double zz; zz = a; a = b; b = zz; }
#define ISWAP(a,b) { INT zz; zz = a; a = b; b = zz; }
#define NOSLN 0.1278433
#define GFACT 2.5
#define EFACT 3.0

#define MAXCOLOR 20
#define MAXWIN 5

extern INT lf_error;
