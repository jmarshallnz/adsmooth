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

#include "ConstructPDF.h"
#include "mathutils.h"
#include <math.h>
#include <limits>
#include <stdio.h>
//#include "R.h"

#define NMATLAB

#ifndef NMATLAB
#include "mex.h"      // for matlab
#endif

double const quietNaN = std::numeric_limits<double>::quiet_NaN();

double *alloc(unsigned int amount)
{
  double *mem = new double[amount];
  if (!mem)
    printf("unable to allocate %d doubles\n", amount);
  return mem;
}

void dealloc(double *mem)
{
  if (mem)
    delete[] mem;
}

void get_h_and_bc(unsigned int num_data, const double *data, int kernel, double bandwidth, unsigned int num_gridX, unsigned int num_gridY, const double *grid, unsigned int boundary_size, const double *boundary, const double *boundary_index, bool adaptive, int boundary_corrected, double **h_data, double **bc_data, double **bc_grid)
{
  if (adaptive)
  {
    // construct a pilot (TODO: This can probably go directly into get_adaptive_bandwidth, to save an allocation)
    double *pilot = alloc(num_gridX*num_gridY);
    construct_pdf(num_data, data, kernel, bandwidth, num_gridX, num_gridY, grid, boundary_size, boundary, boundary_index, false, boundary_corrected, pilot);
    // find our h_data values (adaptive bandwidth at each data point)
    double *h_grid = alloc(num_gridX*num_gridY);
    get_adaptive_bandwidth(num_data, data, bandwidth, num_gridX, num_gridY, grid, pilot, *h_data, h_grid);
    get_bc_rectangle(num_data, data, kernel, *h_data, boundary_size, boundary, NULL, adaptive, boundary_corrected, *bc_data);
    if (*bc_grid)
      get_bc_rectangle(num_gridX*num_gridY, grid, kernel, h_grid, boundary_size, boundary, boundary_index, adaptive, boundary_corrected, *bc_grid);
    dealloc(h_grid);
    dealloc(pilot);
  }
  else
  {
    for (unsigned int i = 0; i < num_data; i++)
      (*h_data)[i] = bandwidth;
    get_bc_rectangle(num_data, data, kernel, &bandwidth, boundary_size, boundary, NULL, false, boundary_corrected, *bc_data);
    if (*bc_grid)
      get_bc_rectangle(num_gridX*num_gridY, grid, kernel, &bandwidth, boundary_size, boundary, boundary_index, false, boundary_corrected, *bc_grid);
  }
}


double estimate_mise(unsigned int num_data, const double *data, int kernel, double bandwidth, unsigned int num_gridX, unsigned int num_gridY, const double *grid, unsigned int boundary_size, const double *boundary, const double *boundary_index, bool adaptive, int boundary_corrected)
{
  printf("estimate_mise called with num_data=%d, bandwidth=%f, num_gridX=%d, num_gridY=%d, boundary_size=%d, adaptive=%s, boundary_corrected=%s\n",
    num_data, bandwidth, num_gridX, num_gridY, boundary_size, adaptive ? "true" : "false", boundary_corrected ? ((boundary_corrected > 1) ? "full" : "half") : "false");
  double *h_data = alloc(num_data);
  double *bc_data = alloc(num_data * 4);
  double *bc_grid = alloc(num_gridX*num_gridY * 4);

  get_h_and_bc(num_data, data, kernel, bandwidth, num_gridX, num_gridY, grid, boundary_size, boundary, boundary_index, adaptive, boundary_corrected, &h_data, &bc_data, &bc_grid);

  // our cross-validation term
  double cv_sum = 0;
  for (unsigned int i = 0; i < num_data; i++)
    cv_sum += evaluate_pdf(num_data, data, kernel, h_data, data[2*i], data[2*i+1], &bc_data[4*i], i);

  // and the f_hat^2 term
  double *pdf = alloc(num_gridX*num_gridY);
  evaluate_density(num_data, data, kernel, h_data, num_gridX*num_gridY, grid, bc_grid, boundary_index, pdf);
  for (unsigned int i = 0; i < num_gridX*num_gridY; i++)
    if (!isnan(pdf[i]))
      pdf[i] *= pdf[i];
  double fs = trapezoid_rule(num_gridX, num_gridY, grid, pdf);
  dealloc(pdf);

  dealloc(bc_grid);
  dealloc(bc_data);
  dealloc(h_data);
//  printf("returning %f (cv_sum=%f,fs=%f)\n", fs - 2.0/num_data*cv_sum, cv_sum, fs);
  return fs - 2.0/num_data*cv_sum;
}

double estimate_rr_mise(unsigned int num_data_f, const double *data_f, unsigned int num_data_g, const double *data_g, int kernel, double bandwidth, double delta, unsigned int num_gridX, unsigned int num_gridY, const double *grid, unsigned int boundary_size, const double *boundary, const double *boundary_index, bool adaptive, int boundary_corrected, bool logarithmic)
{
//  Rprintf("estimate_rr_mise called with num_data=%d,%d, bandwidth=%f, num_gridX=%d, num_gridY=%d, boundary_size=%d, adaptive=%s, boundary_corrected=%i\n",
//          num_data_f, num_data_g, bandwidth, num_gridX, num_gridY, boundary_size, adaptive ? "true" : "false", boundary_corrected);

  // grab our h_data_f, bc_data_f, and bc_grid_f
  double *h_data_f = alloc(num_data_f);
  double *bc_data_f = alloc(num_data_f * 4);
  double *bc_grid = alloc(num_gridX*num_gridY * 4);
  get_h_and_bc(num_data_f, data_f, kernel, bandwidth, num_gridX, num_gridY, grid, boundary_size, boundary, boundary_index, adaptive, boundary_corrected, &h_data_f, &bc_data_f, &bc_grid);

  // and for g as well
  double *h_data_g = alloc(num_data_g);
  double *bc_data_g = alloc(num_data_g * 4);
  double *empty = NULL; // no need to recalculate bc_grid
  get_h_and_bc(num_data_g, data_g, kernel, bandwidth, num_gridX, num_gridY, grid, boundary_size, boundary, boundary_index, adaptive, boundary_corrected, &h_data_g, &bc_data_g, &empty);

  // and the rho_hat^2 term
  double *rho = alloc(num_gridX*num_gridY);
  double *pdf_f = alloc(num_gridX*num_gridY);
  evaluate_density(num_data_f, data_f, kernel, h_data_f, num_gridX*num_gridY, grid, bc_grid, boundary_index, pdf_f);

  double *pdf_g = alloc(num_gridX*num_gridY);
  evaluate_density(num_data_g, data_g, kernel, h_data_g, num_gridX*num_gridY, grid, bc_grid, boundary_index, pdf_g);

  if (delta == 0)
  { // we can't have delta == 0 for the routine below - we must have a non-zero minimum for f and g
    // in order for everything to work ok.  To do this we evaluate f_max / 100 as a reasonable starting point.
    // in general this will be too large
    double minFG = INT_MAX;  // probably should be float_max
    double maxFG = 0;
    for (unsigned int i = 0; i < num_gridX*num_gridY; i++)
    {
      if (!isnan(pdf_f[i]))
      {
        if (pdf_f[i] > 0 && pdf_f[i] < minFG)
          minFG = pdf_f[i];
        if (pdf_f[i] > maxFG)
          maxFG = pdf_f[i];
      }
      if (!isnan(pdf_g[i]))
      {
        if (pdf_g[i] > 0 && pdf_g[i] < minFG)
          minFG = pdf_g[i];
        if (pdf_g[i] > maxFG)
          maxFG = pdf_g[i];
      }
    }
//    Rprintf("setting delta=%g (minFG = %g)\n", maxFG * 0.01, minFG);
    delta = maxFG * 0.01;
  }
  // our cross-validation terms

  // first the log^2(f/g) term
  for (unsigned int i = 0; i < num_gridX*num_gridY; i++)
  {
    if (!isnan(pdf_f[i]) && !isnan(pdf_g[i]))
    {
      rho[i] = (pdf_f[i] + delta)/(pdf_g[i] + delta);
      if (logarithmic)
        rho[i] = log(rho[i]);
      rho[i] *= rho[i];
    }
    else
      rho[i] = quietNaN;
  }
  double fs = trapezoid_rule(num_gridX, num_gridY, grid, rho);
//  Rprintf("log^2(f/g) term evaluates to %g\n", fs);

  if (!logarithmic)
  { // the non-log case is quite a bit trickier - we need the following additional terms
    for (unsigned int i = 0; i < num_gridX*num_gridY; i++)
    {
      if (!isnan(pdf_f[i]) && !isnan(pdf_g[i]))
        rho[i] = (pdf_f[i] + delta)*pdf_f[i]/((pdf_g[i] + delta)*(pdf_g[i] + delta));
    }
    fs -= 2*trapezoid_rule(num_gridX, num_gridY, grid, rho);
    for (unsigned int i = 0; i < num_gridX*num_gridY; i++)
    {
      if (!isnan(pdf_f[i]) && !isnan(pdf_g[i]))
        rho[i] = (pdf_f[i] + delta)*(pdf_f[i]+delta)*pdf_g[i]/((pdf_g[i] + delta)*(pdf_g[i] + delta)*(pdf_g[i] + delta));
    }
    fs += 2*trapezoid_rule(num_gridX, num_gridY, grid, rho);
  }
  // now the leave-one-out bits
  double cv_sum_f = 0;
  for (unsigned int i = 0; i < num_data_f; i++)
  {
    double f_i = evaluate_pdf(num_data_f, data_f, kernel, h_data_f, data_f[2*i], data_f[2*i+1], &bc_data_f[4*i], i);
    double g = evaluate_pdf(num_data_g, data_g, kernel, h_data_g, data_f[2*i], data_f[2*i+1], &bc_data_f[4*i]);
//   Rprintf("term %i in cv_sum_f = %g, f_i=%g, g=%g, delta=%g\n", i, log((f_i + delta)/(g + delta))/(f_i+delta), f_i, g, delta);

    if (logarithmic)
      cv_sum_f += log((f_i + delta)/(g + delta))/(f_i+delta);
    else
      cv_sum_f += (f_i + delta)/((g + delta)*(g + delta));
  }

//  Rprintf("cv_sum_f term evaluates to %g\n", cv_sum_f);

  double cv_sum_g = 0;
  for (unsigned int i = 0; i < num_data_g; i++)
  {
    double f = evaluate_pdf(num_data_f, data_f, kernel, h_data_f, data_g[2*i], data_g[2*i+1], &bc_data_g[4*i]);
    double g_i = evaluate_pdf(num_data_g, data_g, kernel, h_data_g, data_g[2*i], data_g[2*i+1], &bc_data_g[4*i], i);
//    Rprintf("term %i in cv_sum_g = %g, f=%g, g_i=%g, delta=%g\n", i, log((f + delta)/(g_i + delta))/(g_i+delta), f, g_i, delta);

    if (logarithmic)
      cv_sum_g += log((f + delta)/(g_i + delta))/(g_i+delta);
    else
      cv_sum_g += (f + delta)*(f + delta) / ((g_i + delta) * (g_i + delta) * (g_i + delta));
  }

//  Rprintf("cv_sum_g term evaluates to %g\n", cv_sum_g);

  dealloc(pdf_f);
  dealloc(pdf_g);
  dealloc(rho);

  dealloc(bc_data_g);
  dealloc(h_data_g);

  dealloc(bc_grid);
  dealloc(bc_data_f);
  dealloc(h_data_f);
//  Rprintf("returning %g (cv_sum_f=%g,cv_sum_g=%g,fs=%g)\n", -fs - 2.0/num_data_f*cv_sum_f + 2.0/num_data_g*cv_sum_g, cv_sum_f, cv_sum_g, fs);
  return -fs - 2.0/num_data_f*cv_sum_f + 2.0/num_data_g*cv_sum_g;
}

double optimal_delta(double *p_f, double *p_g, unsigned int num_gridX, unsigned int num_gridY)
{
  double minFG = INT_MAX;
  double maxF = 0;
  double maxG = 0;
  for (unsigned int i = 0; i < num_gridX*num_gridY; i++)
  {
    if (!isnan(p_f[i]))
    {
      if (p_f[i] > 0 && p_f[i] < minFG)
        minFG = p_f[i];
      if (p_f[i] > maxF)
        maxF = p_f[i];
    }
    if (!isnan(p_g[i]))
    {
      if (p_g[i] > 0 && p_g[i] < minFG)
        minFG = p_g[i];
      if (p_g[i] > maxG)
        maxG = p_g[i];
    }
  }
  return maxG * 0.05;
}

double choose_delta(unsigned int num_data_f, const double *data_f, unsigned int num_data_g, const double *data_g, int kernel, double bandwidth, unsigned int num_gridX, unsigned int num_gridY, const double *grid, unsigned int boundary_size, const double *boundary, const double *boundary_index, bool adaptive, int boundary_corrected)
{
//  Rprintf("choose_delta called with num_data=%d,%d, bandwidth=%f, num_gridX=%d, num_gridY=%d, boundary_size=%d, adaptive=%s, boundary_corrected=%i\n",
//          num_data_f, num_data_g, bandwidth, num_gridX, num_gridY, boundary_size, adaptive ? "true" : "false", boundary_corrected);

  double *h_data_f = alloc(num_data_f);
  double *bc_data_f = alloc(num_data_f * 4);
  double *bc_grid = alloc(num_gridX*num_gridY * 4);
  get_h_and_bc(num_data_f, data_f, kernel, bandwidth, num_gridX, num_gridY, grid, boundary_size, boundary, boundary_index, adaptive, boundary_corrected, &h_data_f, &bc_data_f, &bc_grid);

  // and for g as well
  double *h_data_g = alloc(num_data_g);
  double *bc_data_g = alloc(num_data_g * 4);
  double *empty = NULL; // no need to recalculate bc_grid
  get_h_and_bc(num_data_g, data_g, kernel, bandwidth, num_gridX, num_gridY, grid, boundary_size, boundary, boundary_index, adaptive, boundary_corrected, &h_data_g, &bc_data_g, &empty);

  // and the rho_hat^2 term
  double *pdf_f = alloc(num_gridX*num_gridY);
  evaluate_density(num_data_f, data_f, kernel, h_data_f, num_gridX*num_gridY, grid, bc_grid, boundary_index, pdf_f);

  double *pdf_g = alloc(num_gridX*num_gridY);
  evaluate_density(num_data_g, data_g, kernel, h_data_g, num_gridX*num_gridY, grid, bc_grid, boundary_index, pdf_g);

  double delta = optimal_delta(pdf_f, pdf_g, num_gridX, num_gridY);

  dealloc(pdf_f);
  dealloc(pdf_g);

  dealloc(bc_data_g);
  dealloc(h_data_g);

  dealloc(bc_grid);
  dealloc(bc_data_f);
  dealloc(h_data_f);

  return delta;
}

#ifndef NMATLAB
// function = estimate_mise(data, kernel, h, grid, boundary, in_boundary_index, adaptive, boundary_corrected), h_start/5, h_start*5);
void mexFunction( int nlhs,mxArray *plhs[],int nrhs, mxArray *prhs[])
{
  // Check for proper number of arguments.
  if(nrhs < 6)
  {
    mexErrMsgTxt("At least six inputs required.");
  }
  else if(nlhs>1)
  {
    mexErrMsgTxt("Too many output arguments");
  }

  // calculate size of data points and grid points
  unsigned int num_data = mxGetN(prhs[0]);
  unsigned int num_grid = mxGetN(prhs[3]);
  unsigned int boundary_size = mxGetN(prhs[4]);
  unsigned int boundary_type = mxGetM(prhs[4]);

  // create mise output
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *mise = mxGetPr(plhs[0]);

  // input variables
  double *data = mxGetPr(prhs[0]);
  int kernel = (int)*mxGetPr(prhs[1]);
  double bandwidth = *mxGetPr(prhs[2]);
  double *grid = mxGetPr(prhs[3]);
  double *boundary = mxGetPr(prhs[4]);
  double *in_boundary = mxGetPr(prhs[5]);
  bool use_adaptive = (nrhs < 7) ? false : (*mxGetPr(prhs[6]) != 0);
  int boundary_corrected = (nrhs < 8) ? 0 : *mxGetPr(prhs[7]);

  unsigned int num_gridX = (unsigned int)sqrt((double)num_grid);
  unsigned int num_gridY = num_gridX;
  // and execute our function
  *mise = estimate_mise(num_data, data, kernel, bandwidth, num_gridX, num_gridY, grid, boundary_size, boundary, in_boundary, use_adaptive, boundary_corrected);
}
#endif

/* Relative risk starts here

#ifndef NMATLAB
void mexFunction( int nlhs,mxArray *plhs[],int nrhs, mxArray *prhs[])
{
  // Check for proper number of arguments.
  if(nrhs < 6)
  {
    mexErrMsgTxt("At least six inputs required.");
  }
  else if(nlhs>1)
  {
    mexErrMsgTxt("Too many output arguments");
  }

  // calculate size of data points and grid points
  unsigned int num_data_f = mxGetN(prhs[0]);
  unsigned int num_data_g = mxGetN(prhs[1]);
  unsigned int num_grid = mxGetN(prhs[3]);
  unsigned int boundary_size = mxGetN(prhs[4]);
  unsigned int boundary_type = mxGetM(prhs[4]);

  // create mise output
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *mise = mxGetPr(plhs[0]);

  // input variables
  double *data_f = mxGetPr(prhs[0]);
  double *data_g = mxGetPr(prhs[1]);
  double *h_delta = mxGetPr(prhs[2]);
  double *grid = mxGetPr(prhs[3]);
  double *boundary = mxGetPr(prhs[4]);
  double *in_boundary = mxGetPr(prhs[5]);
  bool use_adaptive = (nrhs < 7) ? false : (*mxGetPr(prhs[6]) != 0);
  bool boundary_corrected = (nrhs < 8) ? false : (*mxGetPr(prhs[7]) != 0);

  unsigned int num_gridX = (unsigned int)sqrt((double)num_grid);
  unsigned int num_gridY = num_gridX;
  // and execute our function
  *mise = estimate_rr_mise(num_data_f, data_f, num_data_g, data_g, h_delta[0], h_delta[1], num_gridX, num_gridY, grid, boundary_size, boundary, in_boundary, use_adaptive, boundary_corrected);
}
#endif*/