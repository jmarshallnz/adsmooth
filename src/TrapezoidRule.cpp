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
#include "mathutils.h"

#define NMATLAB

#ifndef NMATLAB
#include "mex.h"      // for matlab
#endif

double eval_function(unsigned int offset, const double *function)
{
  double f = function[offset];
  if (isnan(f))
    return 0;
  return f;
}

extern "C" {

  void test_my_stuff(int *input, int *output)
  {
    *output = *input * 2;
  }
}

double trapezoid_rule(unsigned int gridSizeX, unsigned int gridSizeY, const double *grid, const double *function)
{
  double volume = 0;
  for (unsigned int x = 0; x < gridSizeX - 1; x++)
  {
    for (unsigned int y = 0; y < gridSizeY - 1; y++)
    {
      unsigned int offset = x * gridSizeY + y;
      double dx = grid[2*(offset + gridSizeY)] - grid[2*offset];
      double dy = grid[2*(offset + 1) + 1] - grid[2*offset + 1];
      //printf("calculating volume of trapezoid at (%f,%f)\n", grid[2*(x*gridSizeX + y)], grid[2*(x*gridSizeX + y) + 1]);
      double f1 = eval_function(offset, function);
      double f2 = eval_function(offset + 1, function);
      double f3 = eval_function(offset + gridSizeY, function);
      double f4 = eval_function(offset + gridSizeY + 1, function);
      volume += (f1 + f2 + f3 + f4) * dx * dy;
    }
  }
  return volume * 0.25;
}

#ifndef NMATLAB
// function a = trapezoid_rule(f, grid)

void mexFunction( int nlhs,mxArray *plhs[],int nrhs, mxArray *prhs[])
{
  // Check for proper number of arguments.
  if(nrhs < 2)
  {
    mexErrMsgTxt("Two inputs required.");
  }
  else if(nlhs>1)
  {
    mexErrMsgTxt("Too many output arguments");
  }

  // calculate size of input data
  unsigned int num_points = (unsigned int)sqrt((double)mxGetN(prhs[0]));
  //printf("Number of points: %i\n", num_points);

  // create area matrix
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *area = mxGetPr(plhs[0]);

  // input variables
  double *function = mxGetPr(prhs[0]);
  double *grid = mxGetPr(prhs[1]);

  *area = trapezoid_rule(num_points, num_points, grid, function);
}
#endif
