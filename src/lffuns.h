/* band.c */
extern void band(), kdeselect();

/* main.c */
extern void SetWinDev();

/* post.c */
extern void SetPSDev();

/* lfd.c */
extern void doreaddata(), dosavedata(), dosavefit();
extern INT  setfilename();

/* cmd.c */
extern INT setuplf(), dispatch(), getlogic();
extern INT getarg(), argused(), argvalis(), readilist();
extern vari *cmdsplit(), *setnextline();
extern char *argarg(), *argval();
extern void setused(), recondat(), cmdint();

/* lfstr.c */
extern void setstrval();
extern INT stm(), ppwhat(), restyp();

/* arith.c */
extern void vassn();
extern INT arvect(), intitem();
extern double areval(), arith(), darith(), dareval(), vitem();
extern vari *varith(), *saveresult(), *arbuild();

/* adap.c */
extern double afit(), aband2(), aband3();
extern INT ainitband();

/* density.c */
extern INT densinit(), likeden();
extern void densrenorm(), dlscv(), lforder(), prresp();

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

/* kdtree.c */
extern void trchck(), growtri(), growquad();
extern void kdtree(), dataf(), phull(), gridf(), crossf(), xbarf(), preset();
extern double dointpoint();

/* lfbasis.c */
extern void basis();

/* linalg.c */
extern void eigen(), svd(), hsvdsolve();
extern void addouter(), multmatscal(), setzero(), choldec(), cholsolve();
extern void QRupd(), QR1(), bacu1(), bacK(), bacT(), solve(), grsc();
extern void xtwxdec();
extern double innerprod();
extern INT factorial(), svdsolve();

/* locfit.c or parfit.c (most) */
extern void prefit(), dercor(), ldf();
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

/* pout.c */
extern INT pretty();
extern void displayplot();
extern void plotmaple(), plotmathe(), plotmatlb(), plotgnup(), plotxwin();

/* preplot.c */
extern void preplot(), cpreplot();
extern INT setpppoints();

/* random.c */
extern double igamma(), ibeta();
extern double rgamma(), rbeta(), rt(), rnorm(), rexp(), runif(), rpois();
extern double pt(), pf(), pchisq(), pnorm();
extern double df(), dchisq();
extern void rseed();

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

/* vari.c */
extern double *vdptr();
extern INT vlen();
extern vari *createvar(), *checkvarlen(), *findvar(), *growvar();
extern void initdb(), deletevar(), deletename(), deleteifhidden(), *viptr(), setvarname();

/* wdiag.c */
extern INT wdiag(), procvhatm();
extern void hvxtwx(), cwdiag();

/* weight.c */
extern double W(), weight(), weightd(), Wd(), Wdd(), wint();
extern double Wconv(), Wconv1(), Wconv4(), Wconv5(), Wconv6(), Wikk();
extern INT iscompact(), wtaylor();

/* help.c */
extern void example();
