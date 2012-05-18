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
#include <limits>
#include <algorithm>
#include <cstdio>
#include "mathutils.h"

#define NMATLAB

#ifndef NMATLAB
#include "mex.h"      // for matlab
#endif

static const double one_on_sqrt_pi = 0.39894228040143267793994605993438;

void integrate_gaussian(double x0, double y0, double x1, double y1, bool adaptive, double *result)
{
  // we don't bother integrating rectangles that are far from our origin, as they will contribute very little to the integral
  const double bound_for_inclusion = 5; // around 10^(-6) for A_{00}

  double dx = (x0 < 0 && x1 >= 0 || x0 >= 0 && x1 < 0) ? 0 : std::min(abs(x0), abs(x1));
  double dy = (y0 < 0 && y1 >= 0 || y0 >= 0 && y1 < 0) ? 0 : std::min(abs(y0), abs(y1));

  if (std::max(dx, dy) > bound_for_inclusion)
    return;

  // int_y0^y1 int_x0^x1 1/(2*pi*h^2)*((x-ox)/h)^m*((y-oy)/h)^n*exp(-((x - ox)^2 + (y - oy)^2)/(2*h^2)) dx dy

  // = 1/(2*pi) * int_u0^u1 u^m exp(-u^2/2) du * int_v0^v1 v^n exp(-v^2/2) dv

  // = [phi(u1) - phi(u0)] * [phi(v1) - phi(v0)]

  // where u = (x - ox)/h and v = (y - oy)/h, and phi is the normal cdf

  // u0 = (x0 - ox)/h, u1 = (x1 - ox)/h

  // h*du = dy, h*dv = dx

  // int x exp(-x^2/2)dx = -1/sqrt(2*pi)*exp(-x^2/2)
  // int x^2 exp(-x^2/2)dx = -1/sqrt(2*pi)*x*exp(-x^2/2) + int exp(-x^2/2)dx
  // int x^3 exp(-x^2/2)dx = 1/sqrt(2*pi)*(-x^2*exp(-x^2/2) - 2*exp(-x^2/2))
  // int x^4 exp(-x^2/2)dx = - x^3*exp(-x^2/2) + 3*int x^2 exp(-x^2/2)dx

  double exp_x0 = one_on_sqrt_pi * exp(-0.5*x0*x0);
  double exp_x1 = one_on_sqrt_pi * exp(-0.5*x1*x1);
  double exp_y0 = one_on_sqrt_pi * exp(-0.5*y0*y0);
  double exp_y1 = one_on_sqrt_pi * exp(-0.5*y1*y1);

  double Kx = phi(x1) - phi(x0);
  double Ky = phi(y1) - phi(y0);
  double xKx = exp_x0 - exp_x1;
  double yKy = exp_y0 - exp_y1;
  double xxKx = (x0*exp_x0 - x1*exp_x1) + Kx;
  double yyKy = (y0*exp_y0 - y1*exp_y1) + Ky;

  if (adaptive)
  {
    double xxxKx = (x0*x0*exp_x0 - x1*x1*exp_x1) + 2*xKx;
    double yyyKy = (y0*y0*exp_y0 - y1*y1*exp_y1) + 2*yKy;
    double xxxxKx = (x0*x0*x0*exp_x0 - x1*x1*x1*exp_x1) + 3*xxKx;
    double yyyyKy = (y0*y0*y0*exp_y0 - y1*y1*y1*exp_y1) + 3*yyKy;

    // f[0] = K
    result[0] += Kx * Ky;

    // f[1] = 2xK + 1/2x^2dK/dx + 1/2xy dK/dy
    result[1] += 2 * xKx*Ky + 0.5 * -xxxKx * Ky + 0.5 * xKx * -yyKy;

    // f[2] = 2yK + 1/2xy dK/dx + 1/2y^2dK/dy
    result[2] += 2 * Kx * yKy + 0.5 * -xxKx * yKy + 0.5 * Kx * -yyyKy;

    // f[3] = xK
    result[3] += xKx * Ky;

    // f[4] = 2x(xK) + 1/2x^2 d(xK)/dx + 1/2xy d(xK)/dy
    result[4] += 2 * xxKx * Ky + 0.5 * (xxKx - xxxxKx)*Ky + 0.5 * xxKx * -yyKy;

    // f[5] = 2x(yK) + 1/2xy d(xK)/dx + 1/2y^2 d(xK)/dy
    double result5 = 2 * xKx * yKy + 0.5 * (xKx - xxxKx)*yKy + 0.5 * xKx * -yyyKy;
    result[5] += result5;

    // f[6] = yK
    result[6] += Kx * yKy;

    // f[7] = 2x(yK) + 1/2x^2 d(yK)/dx + 1/2xy d(yK)/dy
    result[7] += result5;

    // f[8] = 2y(yK) + 1/2xy d(yK)/dx + 1/2y^2 d(yK)/dy
    result[8] += 2 * Kx * yyKy + 0.5 * -xxKx * yyKy + 0.5 * Kx * (yyKy - yyyyKy);
  }
  else
  {
    result[0] += Kx * Ky;
    double result1 = xKx * Ky;
    result[1] += result1;
    double result2 = Kx * yKy;
    result[2] += result2;
    result[3] += result1;
    result[4] += xxKx * Ky;
    double result5 = xKx * yKy;
    result[5] += result5;
    result[6] += result2;
    result[7] += result5;
    result[8] += Kx * yyKy;
  }
}

void biweight_segment(double x1, double y1, double x2, double y2, bool adaptive, double *result, bool flipX, bool flipY)
{
  double y22 = y2*y2;
  double y23 = y2*y22;
  double y24 = y2*y23;
  double y25 = y2*y24;
  double y26 = y2*y25;
  double y27 = y2*y26;

  double x22 = x2*x2;
  double x23 = x2*x22;
  double x24 = x2*x23;
  double x25 = x2*x24;
  double x26 = x2*x25;
  double x27 = x2*x26;
  double x28 = x2*x27;

  double y12 = y1*y1;
  double y13 = y1*y12;
  double y14 = y1*y13;
  double y15 = y1*y14;
  double y16 = y1*y15;
  double y17 = y1*y16;

  double x12 = x1*x1;
  double x13 = x1*x12;
  double x14 = x1*x13;
  double x15 = x1*x14;
  double x16 = x1*x15;
  double x17 = x1*x16;
  double x18 = x1*x17;

  double asinx1 = asin(x1);
  double asinx2 = asin(x2);

  double A00 = 4.0/45*y15*x2+1.0/9*y13*x2+1.0/6*x2*y1+1.0/6*asinx2
               -1.0/5*y15*x2-1.0/3*y13*(-2*x2+2.0/3*x23)-y1*(x2+1.0/5*x25-2.0/3*x23)
               -4.0/45*x1*y25-1.0/9*y23*x1-1.0/6*x1*y2-1.0/6*asinx1
               +1.0/5*y15*x1+1.0/3*y13*(-2*x1+2.0/3*x13)+y1*(x1+1.0/5*x15-2.0/3*x13);
  
  double A10 = -6.0/35*y17-2.0/21*x22*y15+2.0/21*y15-1.0/6*y1*x26
               -1.0/4*(2.0/3*y13-2*y1)*x24-1.0/2*(1.0/5*y15-2.0/3*y13+y1)*x22
               +6.0/35*y27+2.0/21*x12*y25-2.0/21*y25+1.0/6*y1*x16
               +1.0/4*(2.0/3*y13-2*y1)*x14+1.0/2*(1.0/5*y15-2.0/3*y13+y1)*x12;
      
  double A01 = 1.0/6*x2-1.0/42*x27+1.0/10*x25-1.0/6*x23
              -1.0/6*y16*x2-1.0/4*y14*(-2*x2+2.0/3*x23)-1.0/2*y12*(x2+1.0/5*x25-2.0/3*x23)
              -1.0/6*x1+1.0/42*x17-1.0/10*x15
              +1.0/6*x13+1.0/6*y16*x1+1.0/4*y14*(-2*x1+2.0/3*x13)+1.0/2*y12*(x1+1.0/5*x15-2.0/3*x13);

  double A11 = -1.0/48*x28+1.0/12*x26-1.0/8*x24+1.0/12*x22
               -1.0/12*y12*x26-1.0/4*(1.0/2*y14-y12)*x24-1.0/2*(1.0/6*y16-1.0/2*y14+1.0/2*y12)*x22
               +1.0/48*x18-1.0/12*x16+1.0/8*x14-1.0/12*x12
               +1.0/12*y12*x16+1.0/4*(1.0/2*y14-y12)*x14+1.0/2*(1.0/6*y16-1.0/2*y14+1.0/2*y12)*x12;

  double A20 = -3.0/20*x2*y17+17.0/180*x2*y15+1.0/72*y13*x2+1.0/48*x2*y1+1.0/48*asinx2
               -1.0/12*x23*y15-1.0/7*y1*x27-1.0/5*(2.0/3*y13-2*y1)*x25-1.0/3*(1.0/5*y15-2.0/3*y13+y1)*x23
               +3.0/20*x1*y27-17.0/180*x1*y25-1.0/72*y23*x1-1.0/48*x1*y2-1.0/48*asinx1
               +1.0/12*x13*y25+1.0/7*y1*x17+1.0/5*(2.0/3*y13-2*y1)*x15+1.0/3*(1.0/5*y15-2.0/3*y13+y1)*x13;

  double A02 =  1.0/105*x2*y17+1.0/90*x2*y15+1.0/72*y13*x2+1.0/48*x2*y1+1.0/48*asinx2
               -1.0/7*y17*x2-1.0/5*y15*(-2*x2+2.0/3*x23)-1.0/3*y13*(x2+1.0/5*x25-2.0/3*x23)
               -1.0/105*x1*y27-1.0/90*x1*y25-1.0/72*y23*x1-1.0/48*x1*y2-1.0/48*asinx1
               +1.0/7*y17*x1+1.0/5*y15*(-2*x1+2.0/3*x13)+1.0/3*y13*(x1+1.0/5*x15-2.0/3*x13);

  if (flipX)
  {
    A10 *= -1;
    A11 *= -1;
  }
  if (flipY)
  {
    A01 *= -1;
    A11 *= -1;
  }

  // speedup: Don't need the scaling
  /*
  const static double three_on_pi = 0.95492965855137201461330258023509;
  A00 *= three_on_pi;
  A01 *= three_on_pi;
  A10 *= three_on_pi;
  A20 *= three_on_pi;
  A11 *= three_on_pi;
  A02 *= three_on_pi;
  */
  if (adaptive)
  {
    // need a bunch more calcs here
    double A30_10 = 1.0/3*x23*y15+1.0/6*x2*y15-1.0/24*y13*x2-1.0/16*x2*y1-1.0/16*asinx2
                   -4.0/7*y1*x27+4.0/5*(y1-1.0/3*y13)*x25
                   -1.0/3*x13*y25-1.0/6*x1*y25+1.0/24*y23*x1+1.0/16*x1*y2+1.0/16*asinx1
                   +4.0/7*y1*x17-4.0/5*(y1-1.0/3*y13)*x15;

    double A21_01 = 1.0/15*x2*y17-1.0/90*x2*y15-1.0/72*y13*x2-1.0/48*x2*y1-1.0/48*asinx2
                   -4.0/15*y13*x25+4.0/3*(-1.0/5*y15+1.0/3*y13)*x23
                   -1.0/15*x1*y27+1.0/90*x1*y25+1.0/72*y23*x1+1.0/48*x1*y2+1.0/48*asinx1
                   +4.0/15*y13*x15-4.0/3*(-1.0/5*y15+1.0/3*y13)*x13;

    double A12_10 = A21_01;

    double A03_01 = -1.0/35*x2*y17-1.0/30*x2*y15-1.0/24*y13*x2-1.0/16*x2*y1-1.0/16*asinx2
                    -4.0/7*y17*x2+4.0/5*y15*(x2-1.0/3*x23)
                    +1.0/35*x1*y27+1.0/30*x1*y25+1.0/24*y23*x1+1.0/16*x1*y2+1.0/16*asinx1
                    +4.0/7*y17*x1-4.0/5*y15*(x1-1.0/3*x13);

    double A21_10 = -1.0/8*x28+1.0/3*x26-1.0/4*x24-1.0/3*y12*x26+(-1.0/4*y14+1.0/2*y12)*x24
                    +1.0/8*x18-1.0/3*x16+1.0/4*x14+1.0/3*y12*x16-(-1.0/4*y14+1.0/2*y12)*x14;

    double A12_01 = 1.0/24*x28-1.0/6*x26+1.0/4*x24-1.0/6*x22-1.0/4*y14*x24+2*(-1.0/6*y16+1.0/4*y14)*x22
                   -1.0/24*x18+1.0/6*x16-1.0/4*x14+1.0/6*x12+1.0/4*y14*x14-2*(-1.0/6*y16+1.0/4*y14)*x12;

    double A11_10 = -1.0/7*x27+2.0/5*x25-1.0/3*x23-2.0/5*y12*x25+4.0/3*(-1.0/4*y14+1.0/2*y12)*x23
                    +1.0/7*x17-2.0/5*x15+1.0/3*x13+2.0/5*y12*x15-4.0/3*(-1.0/4*y14+1.0/2*y12)*x13;

    double A02_01 = -1.0/3*x2+1.0/21*x27-1.0/5*x25+1.0/3*x23-2.0/3*y16*x2+y14*(x2-1.0/3*x23)
                    +1.0/3*x1-1.0/21*x17+1.0/5*x15-1.0/3*x13+2.0/3*y16*x1-y14*(x1-1.0/3*x13);

    double A20_10 = 8.0/21*x22*y15+16.0/105*y15-2.0/3*y1*x26+(y1-1.0/3*y13)*x24
                   -8.0/21*x12*y25-16.0/105*y25+2.0/3*y1*x16-(y1-1.0/3*y13)*x14;

    double A11_01 = 8.0/105*y17-1.0/3*y13*x24+2*(-1.0/5*y15+1.0/3*y13)*x22
                   -8.0/105*y27+1.0/3*y13*x14-2*(-1.0/5*y15+1.0/3*y13)*x12;

    if (flipX)
    {
      A30_10 *= -1;
      A21_01 *= -1;
      A12_10 *= -1;
      A03_01 *= -1;
      A11_10 *= -1;
      A02_01 *= -1;
    }
    if (flipY)
    {
      A30_10 *= -1;
      A21_01 *= -1;
      A12_10 *= -1;
      A03_01 *= -1;
      A20_10 *= -1;
      A11_01 *= -1;
    }
    result[0] += A00;
    result[1] += A20_10 + A11_01 + 4*A10;
    result[2] += A11_10 + A02_01 + 4*A01;
    result[3] += A10;
    result[4] += A30_10 + A21_01 + 5*A20;
    result[5] += A21_10 + A12_01 + 5*A11;
    result[6] += A01;
    result[7] += A21_10 + A12_01 + 5*A11;
    result[8] += A12_10 + A03_01 + 5*A02;
  }
  else
  {
    result[0] += A00;
    result[1] += A10;
    result[2] += A01;
    result[3] += A10;
    result[4] += A20;
    result[5] += A11;
    result[6] += A01;
    result[7] += A11;
    result[8] += A02;
    printf("results from biweight_segment: %f, %f, %f, %f, %f, %f\n", A00, A01, A10, A11, A20, A02);
  }
}

void biweight_rectangle(double x1, double y1, double x2, double y2, bool adaptive, double *result, bool flipX, bool flipY)
{
  double y22 = y2*y2;
  double y23 = y2*y22;
  double y24 = y2*y23;
  double y25 = y2*y24;
  double y26 = y2*y25;
  double y27 = y2*y26;

  double x22 = x2*x2;
  double x23 = x2*x22;
  double x24 = x2*x23;
  double x25 = x2*x24;
  double x26 = x2*x25;
  double x27 = x2*x26;

  double y12 = y1*y1;
  double y13 = y1*y12;
  double y14 = y1*y13;
  double y15 = y1*y14;
  double y16 = y1*y15;
  double y17 = y1*y16;

  double x12 = x1*x1;
  double x13 = x1*x12;
  double x14 = x1*x13;
  double x15 = x1*x14;
  double x16 = x1*x15;
  double x17 = x1*x16;

  double A00 = 1.0/5*x2*y25+1.0/3*y23*(-2*x2+2.0/3*x23)+y2*(x2+1.0/5*x25-2.0/3*x23)
              -1.0/5*x2*y15-1.0/3*y13*(-2*x2+2.0/3*x23)-y1*(x2+1.0/5*x25-2.0/3*x23)
              -1.0/5*y25*x1-1.0/3*y23*(-2*x1+2.0/3*x13)-y2*(x1+1.0/5*x15-2.0/3*x13)
              +1.0/5*x1*y15+1.0/3*y13*(-2*x1+2.0/3*x13)+y1*(x1+1.0/5*x15-2.0/3*x13);

  double A01 = 1.0/6*y26*x2+1.0/4*y24*(-2*x2+2.0/3*x23)+1.0/2*y22*(x2+1.0/5*x25-2.0/3*x23)
              -1.0/6*y16*x2-1.0/4*y14*(-2*x2+2.0/3*x23)-1.0/2*y12*(x2+1.0/5*x25-2.0/3*x23)
              -1.0/6*y26*x1-1.0/4*y24*(-2*x1+2.0/3*x13)-1.0/2*y22*(x1+1.0/5*x15-2.0/3*x13)
              +1.0/6*y16*x1+1.0/4*y14*(-2*x1+2.0/3*x13)+1.0/2*y12*(x1+1.0/5*x15-2.0/3*x13);

  double A10 = 1.0/6*y2*x26+1.0/4*(2.0/3*y23-2*y2)*x24+1.0/2*(1.0/5*y25-2.0/3*y23+y2)*x22
              -1.0/6*y1*x26-1.0/4*(2.0/3*y13-2*y1)*x24-1.0/2*(1.0/5*y15-2.0/3*y13+y1)*x22
              -1.0/6*y2*x16-1.0/4*(2.0/3*y23-2*y2)*x14-1.0/2*(1.0/5*y25-2.0/3*y23+y2)*x12
              +1.0/6*y1*x16+1.0/4*(2.0/3*y13-2*y1)*x14+1.0/2*(1.0/5*y15-2.0/3*y13+y1)*x12;

  double A11 = 1.0/12*y22*x26+1.0/4*(1.0/2*y24-y22)*x24+1.0/2*(1.0/6*y26-1.0/2*y24+1.0/2*y22)*x22
              -1.0/12*y12*x26-1.0/4*(1.0/2*y14-y12)*x24-1.0/2*(1.0/6*y16-1.0/2*y14+1.0/2*y12)*x22
              -1.0/12*y22*x16-1.0/4*(1.0/2*y24-y22)*x14-1.0/2*(1.0/6*y26-1.0/2*y24+1.0/2*y22)*x12
              +1.0/12*y12*x16+1.0/4*(1.0/2*y14-y12)*x14+1.0/2*(1.0/6*y16-1.0/2*y14+1.0/2*y12)*x12;

  double A20 = 1.0/7*y2*x27+1.0/5*(2.0/3*y23-2*y2)*x25+1.0/3*(1.0/5*y25-2.0/3*y23+y2)*x23
              -1.0/7*y1*x27-1.0/5*(2.0/3*y13-2*y1)*x25-1.0/3*(1.0/5*y15-2.0/3*y13+y1)*x23
              -1.0/7*y2*x17-1.0/5*(2.0/3*y23-2*y2)*x15-1.0/3*(1.0/5*y25-2.0/3*y23+y2)*x13
              +1.0/7*y1*x17+1.0/5*(2.0/3*y13-2*y1)*x15+1.0/3*(1.0/5*y15-2.0/3*y13+y1)*x13;

  double A02 = 1.0/7*x2*y27+1.0/5*y25*(-2*x2+2.0/3*x23)+1.0/3*y23*(x2+1.0/5*x25-2.0/3*x23)
              -1.0/7*y17*x2-1.0/5*y15*(-2*x2+2.0/3*x23)-1.0/3*y13*(x2+1.0/5*x25-2.0/3*x23)
              -1.0/7*y27*x1-1.0/5*y25*(-2*x1+2.0/3*x13)-1.0/3*y23*(x1+1.0/5*x15-2.0/3*x13)
              +1.0/7*x1*y17+1.0/5*y15*(-2*x1+2.0/3*x13)+1.0/3*y13*(x1+1.0/5*x15-2.0/3*x13);

  if (flipX)
  {
    A10 *= -1;
    A11 *= -1;
  }
  if (flipY)
  {
    A01 *= -1;
    A11 *= -1;
  }

  // speedup: Don't need the scaling
  /*
  const static double three_on_pi = 0.95492965855137201461330258023509;
  A00 *= three_on_pi;
  A01 *= three_on_pi;
  A10 *= three_on_pi;
  A20 *= three_on_pi;
  A11 *= three_on_pi;
  A02 *= three_on_pi;
  */
  if (adaptive)
  {
    // need a bunch more calcs here
    double A30_10 = 4.0/7*y2*x27-4.0/5*(y2-1.0/3*y23)*x25
                   -4.0/7*y1*x27+4.0/5*(y1-1.0/3*y13)*x25
                   -4.0/7*y2*x17+4.0/5*(y2-1.0/3*y23)*x15
                   +4.0/7*y1*x17-4.0/5*(y1-1.0/3*y13)*x15;

    double A21_01 = 4.0/15*y23*x25-4.0/3*(-1.0/5*y25+1.0/3*y23)*x23
                   -4.0/15*y13*x25+4.0/3*(-1.0/5*y15+1.0/3*y13)*x23
                   -4.0/15*y23*x15+4.0/3*(-1.0/5*y25+1.0/3*y23)*x13
                   +4.0/15*y13*x15-4.0/3*(-1.0/5*y15+1.0/3*y13)*x13;

    double A12_10 = A21_01;

    double A03_01 = 4.0/7*y27*x2-4.0/5*y25*(x2-1.0/3*x23)
                   -4.0/7*x2*y17+4.0/5*y15*(x2-1.0/3*x23)
                   -4.0/7*x1*y27+4.0/5*y25*(x1-1.0/3*x13)
                   +4.0/7*y17*x1-4.0/5*y15*(x1-1.0/3*x13);

    double A21_10 = 1.0/3*y22*x26-(-1.0/4*y24+1.0/2*y22)*x24
                   -1.0/3*y12*x26+(-1.0/4*y14+1.0/2*y12)*x24
                   -1.0/3*y22*x16+(-1.0/4*y24+1.0/2*y22)*x14
                   +1.0/3*y12*x16-(-1.0/4*y14+1.0/2*y12)*x14;

    double A12_01 = 1.0/4*y24*x24-2*(-1.0/6*y26+1.0/4*y24)*x22
                   -1.0/4*y14*x24+2*(-1.0/6*y16+1.0/4*y14)*x22
                   -1.0/4*y24*x14+2*(-1.0/6*y26+1.0/4*y24)*x12
                   +1.0/4*y14*x14-2*(-1.0/6*y16+1.0/4*y14)*x12;

    double A11_10 = 2.0/5*y22*x25-4.0/3*(-1.0/4*y24+1.0/2*y22)*x23
                   -2.0/5*y12*x25+4.0/3*(-1.0/4*y14+1.0/2*y12)*x23
                   -2.0/5*y22*x15+4.0/3*(-1.0/4*y24+1.0/2*y22)*x13
                   +2.0/5*y12*x15-4.0/3*(-1.0/4*y14+1.0/2*y12)*x13;

    double A02_01 = 2.0/3*y26*x2-y24*(x2-1.0/3*x23)
                   -2.0/3*y16*x2+y14*(x2-1.0/3*x23)
                   -2.0/3*y26*x1+y24*(x1-1.0/3*x13)
                   +2.0/3*y16*x1-y14*(x1-1.0/3*x13);

    double A20_10 = 2.0/3*y2*x26-(y2-1.0/3*y23)*x24
                   -2.0/3*y1*x26+(y1-1.0/3*y13)*x24
                   -2.0/3*y2*x16+(y2-1.0/3*y23)*x14
                   +2.0/3*y1*x16-(y1-1.0/3*y13)*x14;

    double A11_01 = 1.0/3*y23*x24-2*(-1.0/5*y25+1.0/3*y23)*x22
                   -1.0/3*y13*x24+2*(-1.0/5*y15+1.0/3*y13)*x22
                   -1.0/3*y23*x14+2*(-1.0/5*y25+1.0/3*y23)*x12
                   +1.0/3*y13*x14-2*(-1.0/5*y15+1.0/3*y13)*x12;

    if (flipX)
    {
      A30_10 *= -1;
      A21_01 *= -1;
      A12_10 *= -1;
      A03_01 *= -1;
      A11_10 *= -1;
      A02_01 *= -1;
    }
    if (flipY)
    {
      A30_10 *= -1;
      A21_01 *= -1;
      A12_10 *= -1;
      A03_01 *= -1;
      A20_10 *= -1;
      A11_01 *= -1;
    }

    result[0] += A00;
    result[1] += A20_10 + A11_01 + 4*A10;
    result[2] += A11_10 + A02_01 + 4*A01;
    result[3] += A10;
    result[4] += A30_10 + A21_01 + 5*A20;
    result[5] += A21_10 + A12_01 + 5*A11;
    result[6] += A01;
    result[7] += A21_10 + A12_01 + 5*A11;
    result[8] += A12_10 + A03_01 + 5*A02;
  }
  else
  {
    result[0] += A00;
    result[1] += A10;
    result[2] += A01;
    result[3] += A10;
    result[4] += A20;
    result[5] += A11;
    result[6] += A01;
    result[7] += A11;
    result[8] += A02;
    printf("results from biweight_rectangle: %f, %f, %f, %f, %f, %f\n", A00, A01, A10, A11, A20, A02);
  }
}

void integrate_biweight(double x1, double y1, double x2, double y2, bool adaptive, double *result, bool flipX = false, bool flipY = false)
{
  // at this point we know that x1 < x2 and y1 < y2
  if (x1 < 0 && x2 < 0)
  { // rectangle is to the left of the origin, so integrate by reflecting
    printf("reflecting in x\n");
    integrate_biweight(-x2, y1, -x1, y2, adaptive, result, !flipX, flipY);
    return;
  }
  if (y1 < 0 && y2 < 0)
  { // rectangle is to the bottom of the origin, so integrate by reflecting
    printf("reflecting in y\n");
    integrate_biweight(x1, -y2, x2, -y1, adaptive, result, flipX, !flipY);
    return;
  }
  if (x1 < 0)
  { // split into two rectangles
    printf("splitting in x (%f,%f)->(%f,%f)\n", x1, y1, x2, y2);
    integrate_biweight(0, y1, x2, y2, adaptive, result, flipX, flipY);
    integrate_biweight(0, y1, -x1, y2, adaptive, result, !flipX, flipY);
    return;
  }
  if (y1 < 0)
  { // split into two rectangles
    printf("splitting in y (%f,%f)->(%f,%f)\n", x1, y1, x2, y2);
    integrate_biweight(x1, 0, x2, y2, adaptive, result, flipX, flipY);
    integrate_biweight(x1, 0, x2, -y1, adaptive, result, flipX, !flipY);
    return;
  }
  // finally, we're assured that our rectangle is the appropriate way around
  if (x1 < 0 || y1 < 0)
    printf("x1 < 0 or y1 < 0!!\n");

  // check for intersection - no intersection == 0 area
  if (x1*x1 + y1*y1 > 1)
    return; // zero area

  x2 = std::min(x2, 1.0);
  y2 = std::min(y2, 1.0);

  // integrating
  printf("integrating over (%f,%f)->(%f,%f) ", x1, y1, x2, y2);
  // we have an intersection as the bottom left corner is in the circle

  bool tl = x1*x1 + y2*y2 < 1;
  bool tr = x2*x2 + y2*y2 < 1;
  bool br = x2*x2 + y1*y1 < 1;

  if (tl && tr && br)
  { // all inside -> rectangular integration
    printf("using simple rect integration\n");
    biweight_rectangle(x1, y1, x2, y2, adaptive, result, flipX, flipY);
    return;
  }
  else if (tl && !tr && br)
  { // top right outside only - intersection points on top and right.
    printf("using 2 rect + segment integration\n");
    double x3 = sqrt(1 - y2*y2);
    double y3 = sqrt(1 - x2*x2);
    biweight_segment(x3,y3,x2,y2,adaptive,result, flipX, flipY);
    biweight_rectangle(x1,y1,x3,y2,adaptive,result, flipX, flipY);
    biweight_rectangle(x3,y1,x2,y3,adaptive,result, flipX, flipY);
    return;
  }
  else if (!tl && !tr && br)
  { // top is outside only - intersection points on left and right.
    printf("using rect + segment integration\n");
    double y3 = sqrt(1 - x1*x1);
    double y4 = sqrt(1 - x2*x2);
    biweight_segment(x1,y4,x2,y3,adaptive,result, flipX, flipY);
    biweight_rectangle(x1,y1,x2,y4,adaptive,result, flipX, flipY);
    return;
  }
  else if (tl && !tr && !br)
  { // right is outside only - intersection point top and bottom.
    printf("using rect + segment integration\n");
    double x3 = sqrt(1 - y1*y1);
    double x4 = sqrt(1 - y2*y2);
    biweight_segment(x4,y1,x3,y2,adaptive,result, flipX, flipY);
    biweight_rectangle(x1,y1,x4,y2,adaptive,result, flipX, flipY);
    return;
  }
  else if (!tl && !tr && !br)
  { // all but bottom left is outside - intersection points bottom and left.
    printf("using segment integration\n");
    double y3 = sqrt(1 - x1*x1);
    double x3 = sqrt(1 - y1*y1);
    biweight_segment(x1,y1,x3,y3,adaptive,result, flipX, flipY);
    return;
  }
  else
  {
    printf("unable to integrate!\n");
  }
}

void integrate_rectangle(int kernel, const double *origin, const double one_on_h, const double *rectangle, double *result, bool adaptive)
{
  // integrate over the rectangle
  double x0 = (origin[0] - rectangle[0]) * one_on_h;
  double x1 = (origin[0] - rectangle[2]) * one_on_h;
  double y0 = (origin[1] - rectangle[1]) * one_on_h;
  double y1 = (origin[1] - rectangle[3]) * one_on_h;

  if (kernel == 2) // gaussian kernel
    integrate_gaussian(x0, y0, x1, y1, adaptive, result);
  else if (kernel == 1) // biweight - note: reversal of x0,x1,y0,y1 to make the routine simpler (assumes x0 < x1, y0 < y1)
    integrate_biweight(x1, y1, x0, y0, adaptive, result);
}

void get_bc_rectangle(unsigned int num_points, const double *x, const int kernel, const double *h, unsigned int num_rectangles, const double *rectangles, const double *in_boundary, bool use_adaptive, int correction_type, double *correction)
{
  // simplest case:  No correction
  if (0 == correction_type)
  {
    for (unsigned int i = 0; i < num_points; i++)
    {
      correction[4*i] = 1;
      correction[4*i + 1] = 1;
      correction[4*i + 2] = 0;
      correction[4*i + 3] = 0;
    }
    return;
  }

  // run through the x grid and integrate
  for (unsigned int i = 0; i < num_points; i++)
  {
    if (in_boundary && !in_boundary[i])
    { // outside our boundary
      correction[4*i] = 1;
      correction[4*i + 1] = 1;
      correction[4*i + 2] = 0;
      correction[4*i + 3] = 0;
      continue;
    }
    // set origin + bandwidth
    const double origin[2] = { x[2*i], x[2*i + 1] };
    const double one_on_h = use_adaptive ? 1/h[i] : 1/h[0];

    // run through our rectangles and integrate
    double result[9];
    for (int j = 0; j < 9; j++)
      result[j] = 0;
    for (unsigned int r = 0; r < num_rectangles; r++)
      integrate_rectangle(kernel, origin, one_on_h, &rectangles[r*4], result, use_adaptive);

    // test for order 1 correction
    if (correction_type == 1)
    {
      correction[4*i] = result[0];
      correction[4*i+1] = 1;
      correction[4*i+2] = 0;
      correction[4*i+3] = 0;
    }
    else
    {
      // compute correction terms
      double b0 = result[4]*result[8] - result[5]*result[7];
      double b1 = result[2]*result[7] - result[1]*result[8];
      double b2 = result[1]*result[5] - result[2]*result[4];
      double d = b0*result[0] + b1*result[3] + b2*result[6];
      correction[4*i] = d;
      correction[4*i + 1] = b0;
      correction[4*i + 2] = b1;
      correction[4*i + 3] = b2;
    }
    printf("%.3f %.3f %1.4f %1.4f %1.4f %1.4f\n", origin[0], origin[1], correction[4*i+0], correction[4*i+1], correction[4*i+2], correction[4*i+3]);
  }
//  printf("Time taken for boundary calculation: %d ms\n", GetTickCount() - time);
}


#ifndef NMATLAB
// function bc = get_boundary_correction(data, kernel, h, boundary, in_boundary_index, correction_level)

void mexFunction( int nlhs,mxArray *plhs[],int nrhs, mxArray *prhs[])
{
  // Check for proper number of arguments.
  if(nrhs < 4)
  {
    mexErrMsgTxt("Four inputs required.");
  }
  else if(nlhs>1)
  {
    mexErrMsgTxt("Too many output arguments");
  }

  // calculate size of input data
  unsigned int num_points = mxGetN(prhs[0]);
  printf("Number of points: %i\n", num_points);

  // create boundary correction matrix
  plhs[0] = mxCreateDoubleMatrix(4, num_points, mxREAL);
  double *correction = mxGetPr(plhs[0]);

  // input variables
  double *points = mxGetPr(prhs[0]);
  int kernel = (int)*mxGetPr(prhs[1]);
  double *bandwidth = mxGetPr(prhs[2]);
  double *boundary = mxGetPr(prhs[3]);
  double *in_boundary = (nrhs < 5) ? NULL : mxGetPr(prhs[4]);
  int correction_level = (nrhs < 6) ? 2 : (int)*mxGetPr(prhs[5]);
  printf("Have In Boundary: %s\n", in_boundary ? "Yes" : "No");

  // check if it's adaptive correction
  bool use_adaptive = mxGetN(prhs[2]) == num_points;
  printf("Use adaptive: %s\n", use_adaptive ? "Yes" : "No");

  // check if it's triangular or rectangular boundary correction
  int boundary_type = mxGetM(prhs[3]);
  unsigned int boundary_size = mxGetN(prhs[3]);
  printf("Boundary Type: %i\n", boundary_type);
  printf("Boundary Size: %i\n", boundary_size);

  if (boundary_type == 2)
  {
/*  TODO: Triangular routines
    int num_triangles = mxGetN(prhs[2]) / 3;
    double *triangles = mxGetPr(prhs[2]);

    //  printf("triangles:\n");
    //  for (int i = 0; i < num_triangles * 6; i++)
    //    printf("%f\n", triangles[i]);

    // we assume for the dcutri code that the triangles are stored as follows:
    // x0
    //    y0
    // x1
    //    y1
    // x2
    //    y2
    use_adaptive_correction = nrhs > 3;
    // call our routine
    calc_boundary(grid, bandwidth, num_points, triangles, num_triangles, correction); */
  }
  else
  {
    get_bc_rectangle(num_points, points, kernel, bandwidth, boundary_size, boundary, in_boundary, use_adaptive, correction_level, correction);
  }
}
#endif
