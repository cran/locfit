/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *  solve f(x)=c by various methods, with varying stability etc...
 *    xlo and xhi should be initial bounds for the solution.
 *    convergence criterion is |f(x)-c| < tol.
 *
 *  double solve_secant(f,c,xlo,xhi,tol,bd_flag,err)
 *    secant method solution of f(x)=c.
 *    xlo and xhi are starting values and bound for solution.
 *    tol = convergence criterion, |f(x)-c| < tol.
 *    bd_flag = if (xlo,xhi) doesn't bound a solution, what action to take?
 *      BDF_NONE returns error.
 *      BDF_EXPRIGHT increases xhi.
 *      BDF_EXPLEFT  decreases xlo.
 *    err = error flag.
 *    The (xlo,xhi) bound is not formally necessary for the secant method.
 *    But having such a bound vastly improves stability; the code performs
 *    a bisection step whenever the iterations run outside the bounds.
 *
 *  double solve_nr(f,f1,c,x0,tol,err)
 *    Newton-Raphson solution of f(x)=c.
 *    f1 = f'(x).
 *    x0 = starting value.
 *    tol = convergence criteria, |f(x)-c| < tol.
 *    err = error flag.
 *    No stability checks at present.
 *
 *  double solve_fp(f,x0,tol)
 *    fixed-point iteration to solve f(x)=x.
 *    x0 = starting value.
 *    tol = convergence criteria, stops when |f(x)-x| < tol.
 *    Convergence requires |f'(x)|<1 in neighborhood of true solution;
 *      f'(x) \approx 0 gives the fastest convergence.
 *    No stability checks at present.
 *
 *  TODO: additional error checking, non-convergence stop.
 */

#include <R.h>
#include <math.h>
#include <stdio.h>
#include "mutil.h"

double solve_secant(double (*f)(), double c, double xlo, double xhi, double tol,
                    int bd_flag, int *err)
/* double solve_secant(f,c,xlo,xhi,tol,bd_flag,err)
   double (*f)(), c, xhi, xlo, tol;
   int bd_flag, *err; */
{ double ylo, yhi, x1, x2, x, y1, y2, y;
  *err = 0;
  ylo = f(xlo)-c;
  yhi = f(xhi)-c;

  switch(bd_flag)
  { case BDF_EXPRIGHT:
      while (yhi*ylo > 0)
      { x1 = xhi + (xhi-xlo);
        y1 = f(x1) - c;
        xlo = xhi; xhi = x1;
        ylo = yhi; yhi = y1;
      }
      break;
    case BDF_EXPLEFT:
      while (yhi*ylo > 0)
      { x1 = xlo - (xhi-xlo);
        y1 = f(x1) - c;
        xhi = xlo; xlo = x1;
        yhi = ylo; ylo = y1;
      }
      break;
    case BDF_NONE:
    default:
      if (yhi*ylo > 0)
      { *err = 1;
        return((xlo+xhi)/2);
      }
      break;
  }

  x1 = xlo; y1 = ylo;
  x2 = xhi; y2 = yhi;

  while (1)
  { x = x2 + (x1-x2)*y2/(y2-y1);
    if ((x<=xlo) | (x>=xhi)) x = (xlo+xhi)/2;
    y = f(x)-c;
    if (fabs(y) < tol) return(x);
    if (y*ylo>0) { xlo = x; ylo = y; }
            else { xhi = x; yhi = y; }
if (y2==y)
{ Rprintf("secant: y2 %12.9f\n",y2);
  return(x);
}
    x1 = x2; y1 = y2;
    x2 = x;  y2 = y;
  }
}

double solve_nr(double (*f)(), double (*f1)(), double c, double x0, double tol,
                int *err)
/* double solve_nr(f,f1,c,x0,tol,err)
   double (*f)(), (*f1)(), c, x0, tol;
   int *err; */
{ double y;
  do
  { y = f(x0)-c;
    x0 -= y/f1(x0);
  } while (fabs(y)>tol);
  return(x0);
}

double solve_fp(double (*f)(), double x0, double tol, int maxit)
/* double solve_fp(f,x0,tol,maxit)
   double (*f)(), x0, tol;
   int maxit; */
{ double x1=0.0;
  int i;
  for (i=0; i<maxit; i++)
  { x1 = f(x0);
    if (fabs(x1-x0)<tol) return(x1);
    x0 = x1;
  }
  return(x1); /* although it hasn't converged */
}
