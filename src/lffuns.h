


/* FILES IN THE src DIRECTORY */

/* adap.c */
extern double afit(), aband2(), aband3();
extern INT ainitband();

/* band.c */
extern void band(), kdeselect();

/* density.c */
extern INT densinit(), likeden();
extern void densrenorm(), dlscv(), lforder(), prresp();

/* dist.c */
extern double igamma(), ibeta();
extern double pf(), pchisq(), pnorm();
extern double df(), dchisq();

/* ev_atree.c */
extern void atree_start(), atree_grow(), guessnv();
extern double atree_int();

/* ev_interp.c */
extern double dointpoint(), cubintd();
extern double linear_interp(), cubic_interp(), rectcell_interp();
extern INT exvval();
extern void exvvalpv(), hermite2();

/* ev_kdtre.c */
extern void kdtre_start();
extern double kdtre_int();

/* ev_main.c */
extern void trchck();
extern void dataf(), gridf(), crossf(), xbarf(), preset();
extern INT newsplit();
#ifndef CVERSION
extern vari *createvar();
#endif

/* ev_trian.c */
extern void triang_start(), triang_grow();
extern double triang_int();

/* family.c */
extern INT links(), stdlinks(), defaultlink(), validlinks();

/* frend.c */
extern void fitfun(), degfree(), ressumm();
extern INT procv(), procvraw(), procvvord();
extern double base(), cens(), prwt(), resp(), getxi(), rss();
extern INT calcp(), coefnumber();

/* kappa0.c */
extern double critval(), critvalc(), tailp(), taild();
extern INT constants();

/* lfbasis.c */
extern void basis();

/* lfstr.c */
extern void setstrval();
extern INT stm(), ppwhat(), restyp();

/* linalg.c */
extern void eigen(), svd(), hsvdsolve();
extern void addouter(), multmatscal(), setzero(), choldec(), cholsolve();
extern void QRupd(), QR1(), bacu1(), bacK(), bacT(), solve(), grsc();
extern void xtwxdec();
extern double innerprod();
extern INT factorial(), svdsolve();

/* locfit.c or parfit.c (most) */
extern void prefit(), dercor(), ldf(), unitvec();
extern double vxtwxv();
extern INT ident, locfit(), vxtwx();

/* math.c */
extern double lflgamma(), lferf(), lferfc(), lfdaws();
extern double ptail(), logit(), expit();
extern double lgamma(), erf(), erfc();

/* minmax.c */
extern double ipower(), minmax(), weightmm();

/* nbhd.c */
extern double kordstat(), nbhd(), rho();

/* odint.c */
extern INT onedint();
extern void recurint();

/* pcomp.c */
extern double addparcomp();
extern void compparcomp(), subparcomp(), subparcomp2(), pcchk();
extern INT lenpc(), noparcomp(), hasparcomp();

/* preplot.c */
extern void preplot(), cpreplot();
extern INT setpppoints();

/* resid.c */
extern double resid();
extern void cfitted();
extern vari *vfitted(), *vresid();

/* scbmax.c */
extern void cscbmax();

/* simul.c */
extern void liksim(), scbsim(), scbmax(), regband(), rband();

/* startlf.c */
extern void bbox(), deschk(), startlf(), preproc(), fitdefault();
extern void fitoptions(), clocfit(), endfit();
extern INT nofit();

/* wdiag.c */
extern INT wdiag(), procvhatm();
extern void hvxtwx(), cwdiag();

/* weight.c */
extern double W(), weight(), weightd(), Wd(), Wdd(), wint();
extern double Wconv(), Wconv1(), Wconv4(), Wconv5(), Wconv6(), Wikk();
extern INT iscompact(), wtaylor();






/* FILES IN THE src-c DIRECTORY */

#ifdef CVERSION
/* arith.c */
extern void vassn();
extern INT arvect(), intitem();
extern double areval(), arith(), darith(), dareval(), vitem();
extern vari *varith(), *saveresult(), *arbuild();

/* cmd.c */
extern INT dispatch(), getlogic();
extern INT getarg(), argused(), argvalis(), readilist();
extern vari *cmdsplit(), *setnextline();
extern char *argarg(), *argval();
extern void setuplf(), setused(), recondat(), cmdint();

/* help.c */
extern void example();

/* lfd.c */
extern void doreaddata(), dosavedata(), dosavefit();
extern INT  setfilename();

/* main.c */
extern void SetWinDev();

/* post.c */
extern void SetPSDev();

/* pout.c */
extern INT pretty();
extern void displayplot();
extern void plotmaple(), plotmathe(), plotmatlb(), plotgnup(), plotxwin();

/* random.c */
extern double rnorm(), rexp(), runif(), rpois();
extern void rseed();

/* vari.c */
extern INT vbytes();
extern vari *createvar(), *findvar(), *growvar();
extern void initdb(), deletevar(), deletename(), deleteifhidden(), setvarname();
extern void *viptr();
extern double *vdptrn();

#endif
