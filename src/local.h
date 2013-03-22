/* Which version of locfit is being built?
   One of the following should be uncommented.
   If you change this after building one version
   of Locfit, then you'll also need to
   touch *.c
   at the unix prompt.
   For R compilation, both RVERSION and SVERSION should be uncommented.
   INTERFACE (not CVERSION) should be used, for user-written interfaces.
*/
/* #define CVERSION */
#define RVERSION
#define SVERSION
/* #define INTERFACE */


/* The definitions of gamma() and lgamma() functions vary by math library.
   The LGAMMA(arg) macro must be set to a legitimate log(gamma(arg)) function.
   If all else fails, change lgamma() to lfgamma() in the following.
*/
extern double lgamma();
#define LGAMMA(arg) lgamma(arg)


/* Does your system support popen() and pclose()
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

typedef int INT;

/******** NOTHING BELOW HERE NEEDS CHANGING **********/

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

#define ERROR(args) {printe("Error: "); printe args ; printe("\n"); lf_error=1; }
#define WARN(args) {printe("Warning: ");printe args ; printe("\n"); }

#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SGN(x) (((x)>0) ? 1 : -1)
#define SQR(x) ((x)*(x))
#define FSWAP(a,b) { double zz; zz = a; a = b; b = zz; }
#define ISWAP(a,b) { INT zz; zz = a; a = b; b = zz; }
#define MXDIM 15
#define MXDEG 7
#define NOSLN 0.1278433
#ifndef PI
#define PI    3.141592653589793238462643
#endif
#define S2PI  2.506628274631000502415765
#define SQRT2 1.4142135623730950488
#define SQRPI 1.77245385090552
#define LOGPI 1.144729885849400174143427
#define GOLDEN 0.61803398874989484820
#define GFACT 2.5

#define MAXCOLOR 20
#define MAXWIN 5

#define ACP  0
#define AKAT 1
#define AMDI 2

#define DALP 0
#define DFXH 1
#define DADP 2
#define DCUT 3
#define DLK  4
#define DT0  5
#define DT1  6
#define DRV  7
#define LEND 8

#define ENULL  0
#define ETREE  1
#define EPHULL 2
#define EDATA  3
#define EGRID  4
#define EKDTR  5
#define EKDCE  6
#define ECROS  7
#define EPRES  8

#define MN     0
#define MP     1
#define MDEG0  2
#define MDEG   3
#define MDIM   4
#define MACRI  5
#define MKER   6
#define MKT    7
#define MIT    8
#define MMINT  9
#define MMXIT 10
#define MREN  11
#define MEV   12
#define MTG   13
#define MLINK 14
#define MDC   15
#define MK    16
#define MWH   17
#define MGETH 18
#define LENM  19

#define LINIT  0
#define LDEFAU 1
#define LCANON 2
#define LIDENT 3
#define LLOG   4
#define LLOGIT 5
#define LINVER 6
#define LSQRT  7

#define LLEN  4
#define ZLIK  0
#define ZMEAN 1
#define ZDLL  2
#define ZDDLL 3

#define WRECT 1
#define WEPCH 2
#define WBISQ 3
#define WTCUB 4
#define WTRWT 5
#define WGAUS 6
#define WTRIA 7
#define WQUQU 8
#define W6CUB 9
#define WMINM 10
#define WPARM 11

#define KSPH  1
#define KPROD 2
#define KANG  3
#define KCE   4
#define KLEF  5
#define KRIG  6

#define TNUL 0
#define TDEN 1
#define TRAT 2
#define THAZ 3
#define TGAUS 4
#define TLOGT 5
#define TPOIS 6
#define TGAMM 7
#define TGEOM 8
#define TWEIB 9
#define TCIRC 10
#define TROBT 11

#define IDEFA 1
#define IMULT 2
#define IPROD 3
#define IMLIN 4
#define IHAZD 5
#define IHARD 6
#define IMONT 7

#define PCOEF 1
#define PT0   2
#define PNLX  3
#define PBAND 4
#define PDEGR 5

#define CBAK 0
#define CAXI 1
#define CTEX 2
#define CLIN 3
#define CPOI 4
#define CCON 5
#define CCLA 6
#define CSEG 7
#define CPA1 8
#define CPA2 9

#define VDOUBLE 0
#define VINT    1
#define VCHAR   2
#define VARGL   3

typedef char varname[15];
typedef struct {
  char *arg, *val;
  INT used; } carg;

typedef struct {
  void (*AddColor)(), (*SetColor)(), (*ClearScreen)(), (*TextDim)(), (*DoText)();
  void (*DrawPoint)(), (*DrawLine)(), (*DrawPatch)(), (*wrapup)();
  INT (*makewin)(), ticklength, defth, deftw;
} device;

struct design {
  double *dw, *X, *Z, *w, *di, *res, *th, *wd, h, xb[MXDIM];
  double *V, *P, *Q, *f1, *f2, *ss, *oc, *cf, *dg, llk, na;
  INT *ind, n, p, lw, li, sm, pref, (*itype)();
  INT (*vfun)(); };

struct tree {
  double *x[MXDIM], *y, *w, *base, *c, *xl;
  double *tw, *xev, *coef, *nlx, *t0, *lik, *h, *deg;
  double *sv, *L, *fl, *sca, *dp, kap[3];
  INT *iw, *ce, *s, *lo, *hi, sty[MXDIM];
  INT ltw, liw, ll, *mg, nvm, ncm, vc;
  INT nl, nv, nnl, nce, nk, nn, *mi, ord, *deriv, nd;
  varname yname, xname[MXDIM], wname, bname, cname; };

#define STEMPTY   0
#define STREGULAR 1
#define STRESULT  2
#define STHIDDEN  3
#define STPLOTVAR 4
#define STSYSTEM  5
#define STSYSPEC  6

typedef struct {
  varname name;
  INT n, len, mode, stat;
  double *dpr; } vari;

struct arc {
  char cmd;
  double x, *v, (*f)();
  INT m;
  vari *vv;
  struct arc *ar[3]; };

typedef struct arc arstruct;

typedef struct {
  double *wk, *x[5], *y[5], *z[5], theta, phi, xl[2], yl[2], zl[2], sl[10];
  INT d, r, lw, add, mx[5], my[5], mz[5], ty, nsl;
  char type[5], main[50], xlab[50], ylab[50], zlab[50];
  vari *track; } plots;
/* pl.ty for validation purposes.. */
#define PLNONE 0
#define PLDATA 1
#define PLFIT  2
#define PLTRK  4

struct lfcol {
  char name[10];
  INT n, r, g, b;
};

extern INT lf_error;

/* main.c */
extern void SetWinDev();

/* post.c */
extern void SetPSDev();

/* lfd.c */
extern void doreaddata(), dosavedata(), dosavefit();
extern INT  setfilename();

/* cmd.c */
extern INT setuplf(), dispatch(), getlogic();
extern INT getarg(), argused(), argvalis();
extern void deletevar();
extern vari *createvar(), *findvar(), *cmdsplit();
extern char *argarg(), *argval();

/* lfstr.c */
extern void setstrval();
extern INT stm();

/* arith.c */
extern void arbuild(), vassn(), *viptr();
extern INT arvect(), vlen(), intitem();
extern double areval(), arith(), darith(), dareval(), vitem(), *vdptr();
extern vari *varith(), *saveresult();

/* adap.c */
extern double afit();

/* density.c */
extern double likeden();
extern INT densinit();
extern void densrenorm();

/* frend.c */
extern void fitfun(), degfree(), ressumm(), trace();
extern INT procv(), procvraw(), procvhatm(), procvvord();
extern double base(), cens(), prwt(), resp(), getxi(), rss();
extern INT calcp();

/* kappa0.c */
extern double cv(), cvc(), tailp(), taild();
extern INT constants();

/* kdtree.c */
extern void fitdefault(), trchck(), deschk(), intv(), intg(), intf(), intd();
extern void bbox(), evaluator(), checkvl(), growtri(), growquad();
extern double intp();

/* family.c */
extern INT links();

/* locfit.c or parfit.c (most) */
extern void prefit(), makelxd(), dercor(), ldf(), vxtwx();
extern INT ident, locfit();

/* nbhd.c */
extern double kordstat(), nbhd(), rho();

/* linalg.c */
extern void eigen(), svd(), svdsolve(), addouter(), choldec(), cholsolve();
extern void QRupd(), QR1(), bacu1(), bacK(), bacT(), solve(), grsc();
extern double innerprod();
extern INT factorial();

/* simul.c */
extern void liksim(), scbsim(), scbmax(), kdeselect(), regband(), rband();

/* pout.c */
extern INT pretty();
extern void displayplot();
extern void plotmaple(), plotmathe(), plotmatlb(), plotgnup(), plotxwin();

/* random.c */
extern double igamma(), ibeta();
extern double rgamma(), rbeta(), rt(), rnorm(), rexp(), runif(), rpois();
extern double pt(), pf(), pchisq(), pnorm(), daws(), ptail();
extern double expit(), logit(), df(), dchisq();
extern void rseed();

/* resid.c */
extern double resid();
extern void fitted();
extern vari *vfitted(), *vresid();

/* weight.c */
extern double W(), weightmm(), weight(), weightd(), Wd(), Wdd(), wint();
extern double Wconv(), Wconv1(), Wconv4(), Wconv5(), Wconv6(), Wikk();
extern double minmax();
extern INT wtaylor();

/* help.c */
extern void example();
