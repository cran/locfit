/*
 *   Copyright (c) 1998-1999 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

static unsigned long cc, tv, ss=1;

double runif()
{	if (ss)
	{ WARN(("runif: No seed set."));
          return(0.0);
	}
        cc = cc * 69069;    /* congruential part */
        tv ^= tv >> 15;       /* tausworthe part */
        tv ^= tv << 17;
	return(((tv ^ cc) >> 1) / 2147483648.0);
}

void rseed(seed)
/*
  Seed should be string of at least 8 characters.
*/
char *seed;
{	ss = 0;
	tv = seed[0];
        tv = (tv<<8) | seed[1];
        tv = (tv<<8) | seed[2];
        tv = (tv<<8) | seed[3];
 	cc = seed[4];
        cc = (cc<<8) | seed[5];
        cc = (cc<<8) | seed[6];
        cc = (cc<<8) | seed[7];
        if(cc % 2 == 0)
             cc++;
}

double dchisq(x, df)
double x, df;
{ return(exp(log(x/2)*(df/2-1) - x/2 - LGAMMA(df/2) - LOG_2));
}

double df(x, df1, df2)
double x, df1, df2;
{ double p;
  p = exp(LGAMMA((df1+df2)/2) + df1/2*log(df1/df2) + (df1/2-1)*log(x)
               - LGAMMA(df1/2) - LGAMMA(df2/2) - (df1+df2)/2*log(1+x*df1/df2));
  return(p);
}

double ibeta(x, a, b)
double x, a, b;
{ int flipped = 0, i, k, count;
  double I = 0, temp, pn[6], ak, bk, next, prev, factor, val;
  if (x <= 0) return(0);
  if (x >= 1) return(1);
/* use ibeta(x,a,b) = 1-ibeta(1-x,b,z) */
  if ((a+b+1)*x > (a+1))
  { flipped = 1;
    temp = a;
    a = b;
    b = temp;
    x = 1 - x;
  }
  pn[0] = 0.0;
  pn[2] = pn[3] = pn[1] = 1.0;
  count = 1;
  val = x/(1.0-x);
  bk = 1.0;
  next = 1.0;
  do
  { count++;
    k = count/2;
    prev = next;
    if (count%2 == 0)
      ak = -((a+k-1.0)*(b-k)*val)/((a+2.0*k-2.0)*(a+2.0*k-1.0));
    else
      ak = ((a+b+k-1.0)*k*val)/((a+2.0*k)*(a+2.0*k-1.0));
    pn[4] = bk*pn[2] + ak*pn[0];
    pn[5] = bk*pn[3] + ak*pn[1];
    next = pn[4] / pn[5];
    for (i=0; i<=3; i++)
      pn[i] = pn[i+2];
    if (fabs(pn[4]) >= IBETA_LARGE)
      for (i=0; i<=3; i++)
        pn[i] /= IBETA_LARGE;
    if (fabs(pn[4]) <= IBETA_SMALL)
      for (i=0; i<=3; i++)
        pn[i] /= IBETA_SMALL;
  } while (fabs(next-prev) > DOUBLE_EP*prev);
  factor = a*log(x) + (b-1)*log(1-x);
  factor -= LGAMMA(a+1) + LGAMMA(b) - LGAMMA(a+b);
  I = exp(factor) * next;
  return(flipped ? 1-I : I);
}

/*
 * Incomplete gamma function.
 * Reference:  Abramowitz and Stegun.
 * Assumptions: x >= 0; df > 0.
 */

double igamma(x, df)
double x, df;
{ double factor, term, gintegral, pn[6], rn, ak, bk;
  double increment, df1;
  int i, count, k;
  if (x <= 0.0) return(0.0);
  if (df < 1.0)
  { increment = exp(df*log(x) - x - LGAMMA(df + 1.0));
    df1 = df + 1.0;
  } else
  { increment = 0.0;
    df1 = df;
  }
  factor = exp(df1*log(x) - x - LGAMMA(df1));
  if (x > 1.0 && x >= df1)
  { pn[0] = 0.0;
    pn[2] = pn[1] = 1.0;
    pn[3] = x;
    count = 1;
    rn = 1.0 / x;
    do
    { count++;
      k = count / 2;
      gintegral = rn;
      if (count%2 == 0)
      { bk = 1.0;
        ak = (double)k - df1;
      } else
      { bk = x;
        ak = (double)k;
      }
      pn[4] = bk*pn[2] + ak*pn[0];
      pn[5] = bk*pn[3] + ak*pn[1];
      rn = pn[4] / pn[5];
      for (i=0; i<4; i++)
        pn[i] = pn[i+2];
      if (pn[4] > IGAMMA_LARGE)
        for (i=0; i<4; i++)
          pn[i] /= IGAMMA_LARGE;
    } while (fabs(gintegral-rn) > DOUBLE_EP*rn);
    gintegral = 1.0 - factor*rn;
  } else
  { gintegral = term = 1.0;
    rn = df1;
    do
    { rn += 1.0;
      term *= x/rn;
      gintegral += term;
    } while (term > DOUBLE_EP*gintegral);
    gintegral *= factor/df1;
  }
  return(increment + gintegral);
}

double pf(q, df1, df2)
double q, df1, df2;
{ return(ibeta(q*df1/(df2+q*df1), df1/2, df2/2));
}

double pchisq(q, df)
double q, df;
{ return(igamma(q/2, df/2));
}

double pnorm(x,mu,s)
double x, mu, s;
{ if(x == mu)
    return(0.5);
  x = (x-mu)/s;
  if(x > 0) return((1 + erf(x/SQRT2))/2);
  return(erfc(-x/SQRT2)/2);
}

/*
 * Gaussian random variable.
 * Reference: Kinderman & Monahan, Proceedings of
 * the ASA, Statistical Computing Section, 1975, 128-131.
 */
double rnorm(mu,s)
double mu, s;
{
	double rnormk, u, x2;

	do {
		u = runif();
		rnormk = 1.715527769 * (runif()-0.5) / u;
		x2 = rnormk * rnormk / 4;
		if(x2 <= 1-u)
			break;
	} while(x2 > -log(u));
	return(mu+s*rnormk);
}

double rexp(lb)
double lb;
{ return(-log(runif())/lb);
}

/*
 * Poisson random variable.
 * Simple algorithm for small lambda, else complex algorithm.
 * Crossover point must be at least 5 for the complex algorithm
 * to work correctly.
 * Reference: Devroye, pages 504, 511 and 516 (with corrections!)
 */
double rpois(lambda)
double lambda;
{
	static double olambda = -1, a, mu, delta, d, c1, c2, c3, c4, c5;
	double u, e, n, x, y, w, t, p, q;
	int new = lambda != olambda;

	olambda = lambda;
	if(lambda < 8) {
		if(new)
			a = exp(-lambda);
		q = 1;
		x = -1;
		do {
			q *= runif();
			x++;
		} while(q >= a);
		return(x);
	}

	if(new) {
		mu = floor(lambda);
		delta = sqrt(2 * mu * log(mu * PI128));
		delta = MAX(6.0, MIN(mu, floor(delta)));
		d = 2*mu + delta;
		c1 = sqrt(mu * PI_HALF);
		c2 = c1 + sqrt(d * PI_QUARTER) * exp(1/d);
		c3 = c2 + 1;
		c4 = c3 + EXP78;
		c5 = c4 + 2 * d * exp(-delta*(1+delta/2)/d) / delta;
	}
	while(1) {
		u = c5 * runif();
		e = -log(runif());
		if(u <= c1) {
			n = rnorm(0.0,1.0);
			x = floor(-fabs(n) * sqrt(mu));
			if(x < -mu)
				continue;
			w = n*n/2 + e + x*log(lambda/mu);
		} else if(u <= c2) {
			y = 1 + fabs(rnorm(0.0,1.0)) * sqrt(d/2);
			x = ceil(y);
			if(x > delta)
				continue;
			w = y*(y-2)/d + e + x*log(lambda/mu);
		} else if(u <= c3) {
			x = 0;
			w = e;
		} else if(u <= c4) {
			x = 1;
			w = e + log(lambda/mu);
		} else {
			y = delta - 2*d*log(runif())/delta;
			x = ceil(y);
			w = delta*(1+y/2)/d + e + x*log(lambda/mu);
		}
		w = -w;
		t = x*(x+1) / (2*mu);
		if(x >= 0 && w <= -t)
			return(x+mu);
		if(x < 0 && w > -t)
			continue;
		q = t * ((2*x+1)/(6*mu) - 1);
		if(w > q)
			continue;
		p = x+1 <= 0 ? x+1 : 0;
		p = q - t*t/(3*(mu+p));
		if(w <= p)
			return(x+mu);
		if(w <= x*log(mu) - LGAMMA(mu+x+1) + LGAMMA(mu+1))
			return(x+mu);
	}
}
