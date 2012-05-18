/*
 *      Copyright (C) 2005-2009 Jonathan Marshall
 *      http://ifs.massey.ac.nz
 *
 *  This Program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *
 *  This Program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XBMC; see the file COPYING.  If not, write to
 *  the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *  http://www.gnu.org/copyleft/gpl.html
 *
 */

#include <math.h>
#include "mathutils.h"

#ifndef HAS_R
double erf(double x)
//erf(x) = 2/sqrt(pi)*integral(exp(-t^2),t,0,x)
//       = 2/sqrt(pi)*[x - x^3/3 + x^5/5*2! - x^7/7*3! + ...]
//       = 1-erfc(x)
{
static const double two_sqrtpi=  1.128379167095512574;        // 2/sqrt(pi)
  if (fabs(x) > 2.2) {
  return 1.0 - erfc(x);        //use continued fraction when fabs(x) > 2.2
  }
  double sum= x, term= x, xsqr= x*x;
  int j= 1;
  do {
    term*= xsqr/j;
    sum-= term/(2*j+1);
    ++j;
    term*= xsqr/j;
    sum+= term/(2*j+1);
    ++j;
  } while (fabs(term/sum) > rel_error);   // CORRECTED LINE
  return two_sqrtpi*sum;
}

double erfc(double x)
//erfc(x) = 2/sqrt(pi)*integral(exp(-t^2),t,x,inf)
//        = exp(-x^2)/sqrt(pi) * [1/x+ (1/2)/x+ (2/2)/x+ (3/2)/x+ (4/2)/x+ ...]
//        = 1-erf(x)
//expression inside [] is a continued fraction so '+' means add to denominator only
{
  static const double one_sqrtpi=  0.564189583547756287;        // 1/sqrt(pi)
  if (fabs(x) < 2.2)
  {
    return 1.0 - erf(x);        //use series when fabs(x) < 2.2
  }
  if (x < 0)
  {               //continued fraction only valid for x>0
    return 2.0 - erfc(-x);
  }
  double a=1, b=x;                //last two convergent numerators
  double c=x, d=x*x+0.5;          //last two convergent denominators
  double q1, q2= b/d;             //last two convergents (a/c and b/d)
  double n= 1.0, t;
  do
  {
    t= a*n+b*x;
    a= b;
    b= t;
    t= c*n+d*x;
    c= d;
    d= t;
    n+= 0.5;
    q1= q2;
    q2= b/d;
  } while (fabs(q1-q2)/q2 > rel_error);

  return one_sqrtpi*exp(-x*x)*q2;
}

bool isnan(double x)
{
  return x != x;
}

#endif

double phi(double x)
{
  // normal cdf
  return 0.5*(1 + erf(x * ONE_ON_SQRT_2));
}

