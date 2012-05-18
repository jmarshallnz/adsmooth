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
#include <math.h>
#include <limits>
#include <cstdio>

#define NMATLAB

#ifndef NMATLAB
#include "mex.h"      // for matlab
#endif

void construct_pdf(unsigned int num_data, const double *data, int kernel, double bandwidth, unsigned int num_gridX, unsigned int num_gridY, const double *grid, unsigned int boundary_size, const double *boundary, const double *boundary_index, bool adaptive, int boundary_corrected, double *pdf)
{
  double *index = (double *)boundary_index;
  if (!boundary_index)
  {
    index = new double[num_gridX*num_gridY];
    if (!index)
      printf("failed to allocate %i doubles for the boundary_index\n", num_gridX*num_gridY);
    get_boundary_index_rect(num_gridX*num_gridY, grid, boundary_size, boundary, index);
  }

  double *h_data = new double[num_data];
  if (!h_data)
    printf("failed to allocate %i doubles for the h_data\n", num_data);
  double *bc_grid = new double[num_gridX*num_gridY * 4];
  if (!bc_grid)
    printf("failed to allocate %i doubles for the bc_grid\n", num_gridX*num_gridY*4);
  if (adaptive)
  {
    // construct a pilot (TODO: This can probably go directly into get_adaptive_bandwidth, to save an allocation)
    double *pilot = new double[num_gridX*num_gridY];
    if (!pilot)
      printf("failed to allocate %i doubles for the pilot\n", num_gridX*num_gridY);
    construct_pdf(num_data, data, kernel, bandwidth, num_gridX, num_gridY, grid, boundary_size, boundary, index, false, boundary_corrected, pilot);
    // find our h_data values (adaptive bandwidth at each data point)
    double *h_grid = new double[num_gridX*num_gridY];
    if (!h_grid)
      printf("failed to allocate %i doubles for the h_grid\n", num_gridX*num_gridY);
    get_adaptive_bandwidth(num_data, data, bandwidth, num_gridX, num_gridY, grid, pilot, h_data, h_grid);
    get_bc_rectangle(num_gridX*num_gridY, grid, kernel, h_grid, boundary_size, boundary, index, adaptive, boundary_corrected, bc_grid);
    delete[] h_grid;
    delete[] pilot;
  }
  else
  {
    for (unsigned int i = 0; i < num_data; i++)
      h_data[i] = bandwidth;
    get_bc_rectangle(num_gridX*num_gridY, grid, kernel, &bandwidth, boundary_size, boundary, index, false, boundary_corrected, bc_grid);
//    Rprintf("Boundary correction major components are as follows:\n");
//    for (unsigned int i = 0; i < num_gridX*num_gridY; i++)
//      Rprintf("grid(%f,%f) has BC(%f,%f,%f,%f)\n", grid[2*i], grid[2*i+1], bc_grid[4*i], bc_grid[4*i+1], bc_grid[4*i+2], bc_grid[4*i+3]);
  }
  evaluate_density(num_data, data, kernel, h_data, num_gridX*num_gridY, grid, bc_grid, index, pdf);
  delete[] bc_grid;
  delete[] h_data;
  if (!boundary_index)
    delete[] index;
}

#ifndef NMATLAB
// function pdf = constructpdf(data, kernel, h, grid, boundary[, in_boundary_index, adaptive, boundary_corrected])
void mexFunction( int nlhs,mxArray *plhs[],int nrhs, mxArray *prhs[])
{
  // Check for proper number of arguments.
  if(nrhs < 5)
  {
    mexErrMsgTxt("At least five inputs required.");
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

  // create density matrix
  plhs[0] = mxCreateDoubleMatrix(1, num_grid, mxREAL);
  double *pdf = mxGetPr(plhs[0]);

  // input variables
  double *data = mxGetPr(prhs[0]);
  int kernel = (int)*mxGetPr(prhs[1]);
  double bandwidth = *mxGetPr(prhs[2]);
  double *grid = mxGetPr(prhs[3]);
  double *boundary = mxGetPr(prhs[4]);
  double *in_boundary = (nrhs < 6) ? NULL : mxGetPr(prhs[5]);
  bool use_adaptive = (nrhs < 7) ? false : (*mxGetPr(prhs[6]) != 0);
  int boundary_corrected = (nrhs < 8) ? 0 : (int)*mxGetPr(prhs[7]);

  unsigned int num_gridX = (unsigned int)sqrt((double)num_grid);
  unsigned int num_gridY = num_gridX;
  // and execute our function
  construct_pdf(num_data, data, kernel, bandwidth, num_gridX, num_gridY, grid, boundary_size, boundary, in_boundary, use_adaptive, boundary_corrected, pdf);
}
#endif
