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
#include <stdio.h>

#define NMATLAB

#ifndef NMATLAB
#include "mex.h"      // for matlab
#endif

double const quietNaN = std::numeric_limits<double>::quiet_NaN();

double eval_kernel(const int kernel, const double u, const double v)
{
  switch (kernel)
  {
  case 0:
    {
      // uniform
      static const double one_on_pi = 0.31830988618379067153776752674503;
      double abs_val = u*u + v*v;
      if (abs_val > 1)
        return 0;
      return one_on_pi;
    }
  case 1:
    {
      // biweight (rotationally symmetric)

      // z = (15/16)^2*(1-x.*x).^2.*(1-y.*y).^2;
      // z(find(abs(x) > 1)) = 0;
      // z(find(abs(y) > 1)) = 0;
      static const double three_on_pi = 0.95492965855137201461330258023509;
      double one_minus_abs_val = 1.0 - (u*u + v*v);
      if (one_minus_abs_val < 0)
        return 0;
      return three_on_pi * one_minus_abs_val * one_minus_abs_val;
    }
  case 2:
  default:
    {
      // gaussian
      static const double one_on_2pi = 0.15915494309189533576888376337251;
      return one_on_2pi*exp(-0.5*(u*u + v*v));
    }
  }
}

// evaluate_pdf(num_data, data, h_data, x, y, boundary_correction, ignore_entry = -1) - evaluates the pdf at the given point.
// num_data - number of data points.
// data - data points, of size num_data*2.
// h - bandwidth at each datapoint, of size num_data.
// x,y - point to evaluate the data at.
// boundary_correction - array of 4 boundary correction values (denominator and K, xK, yK terms)
// ignore_entry - if present, ignore this entries contribution to the pdf.
double evaluate_pdf(unsigned int num_data, const double *data, const int kernel, const double *h_data, double x, double y, const double *boundary_correction, unsigned int ignore_entry = (unsigned int)-1)
{
//  printf("evaluate_pdf called with first data point (%f,%f) and eval point (%f,%f), and bc = (%f,%f,%f,%f)\n", data[0], data[1], x, y, boundary_correction[0], boundary_correction[1], boundary_correction[2], boundary_correction[3]);
  double p = 0;
  // note: this can be vectorized - we don't care about order
  for (unsigned int i = 0; i < num_data; i++)
  {
    if (i == ignore_entry)
      continue;
    double u = (x - data[i*2]) / h_data[i];
    double v = (y - data[i*2 + 1]) / h_data[i];

    double z = eval_kernel(kernel, u, v);

    // TODO: Should this be + or - ???
    // TODO: we can add checking for accuracy here (where boundary_correction[0] is small)
    p += (boundary_correction[1] + boundary_correction[2]*u + boundary_correction[3]*v) * z / (boundary_correction[0] * h_data[i] * h_data[i]);
  }
  if (ignore_entry != (unsigned int)-1)
    num_data--;
  return (p > 0) ? p / num_data : 0;
}

// evaluate_density(num_data, data, h, num_points, points, boundary_correction, pdf) - evaluates the density at the given points.
// num_data - number of data points.
// data - data points, of size num_data*2.
// h - bandwidth at each datapoint, of size num_data.
// num_points - number of points to evaluate the density at.
// points - points to evaluate the data at, of size num_points*2.
// boundary_correction - array of 4 boundary correction values per evaluation point (denominator and K, xK, yK terms)
// in_boundary - integer array of size num_points.  If NULL, is assumed an array full of ones.  one == inside, 0 == outside.
// pdf - returned density, of size num_points
void evaluate_density(unsigned int num_data, const double *data, const int kernel, double *h, unsigned int num_points, const double *points, const double *boundary_correction, const double *in_boundary, double *pdf)
{
  if (in_boundary)
  {
    for (unsigned int i = 0; i < num_points; i++)
    {
      if (!in_boundary[i])
        pdf[i] = quietNaN;
      else
        pdf[i] = evaluate_pdf(num_data, data, kernel, h, points[i*2], points[i*2 + 1], &boundary_correction[i*4]);
    }
  }
  else
  { // faster with no boundary checking (saves bool check)
    for (unsigned int i = 0; i < num_points; i++)
      pdf[i] = evaluate_pdf(num_data, data, kernel, h, points[i*2], points[i*2 + 1], &boundary_correction[i*4]);
  }
}

#ifndef NMATLAB
// function pdf = evaluate_pdf(data, kernel, h_data, eval_points, bc_eval_points [, in_boundary_index])
void mexFunction( int nlhs,mxArray *plhs[],int nrhs, mxArray *prhs[])
{
  // Check for proper number of arguments.
  if(nrhs < 5)
  {
    mexErrMsgTxt("Five inputs required.");
  }
  else if(nlhs>1)
  {
    mexErrMsgTxt("Too many output arguments");
  }

  // calculate size of data points and grid points
  unsigned int num_data = mxGetN(prhs[0]);
  unsigned int num_points = mxGetN(prhs[3]);

  // create density matrix
  plhs[0] = mxCreateDoubleMatrix(1, num_points, mxREAL);
  double *pdf = mxGetPr(plhs[0]);

  // input variables
  double *data = mxGetPr(prhs[0]);
  int kernel = (int)*mxGetPr(prhs[1]);
  double *bandwidth = mxGetPr(prhs[2]);
  double *grid = mxGetPr(prhs[3]);
  double *bc_grid = mxGetPr(prhs[4]);
  double *in_boundary = (nrhs < 6) ? NULL : mxGetPr(prhs[5]);

  // and execute our function
  evaluate_density(num_data, data, kernel, bandwidth, num_points, grid, bc_grid, in_boundary, pdf);
}
#endif
