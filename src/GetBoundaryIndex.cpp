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
//#include "mathutils.h"

#define NMATLAB

#ifndef NMATLAB
#include "mex.h"      // for matlab
#endif

void get_boundary_index_rect(unsigned int num_points, const double *points, unsigned int boundary_size, const double *boundary, double *index)
{
  for (unsigned int i = 0; i < num_points; i++)
  {
    const double *x = &points[i*2];
    index[i] = 0;     // outside
    for (unsigned int j = 0; j < boundary_size; j++)
    { // check if this point is inside our rectangle
      const double *rect = &boundary[j*4];
      if (x[0] >= rect[0] && x[0] <= rect[2] && x[1] >= rect[1] && x[1] <= rect[3])
      {
        index[i] = 1; // inside
        break;
      }
    }
  }
}

unsigned int get_points_inside_rect(unsigned int num_points, const double *points, unsigned int boundary_size, const double *boundary, double *points_inside)
{
  unsigned int count = 0;
  for (unsigned int i = 0; i < num_points; i++)
  {
    const double *x = &points[i*2];
    for (unsigned int j = 0; j < boundary_size; j++)
    { // check if this point is inside our rectangle
      const double *rect = &boundary[j*4];
      if (x[0] >= rect[0] && x[0] <= rect[2] && x[1] >= rect[1] && x[1] <= rect[3])
      {
        points_inside[2*count] = points[2*i];
        points_inside[2*count+1] = points[2*i+1];
        count++;
        break;
      }
    }
  }
  return count;
}

double triangle_area(const double *a, const double *b, const double *c)
{
  return fabs(a[0]*(b[1] - c[1]) - a[1]*(b[0] - c[0]) + b[0]*c[1]-b[1]*c[0]);
}

void get_boundary_index_triangle(unsigned int num_points, const double *points, unsigned int boundary_size, const double *boundary, double *index)
{
  static const double tolerance = 1e-6;
  for (unsigned int i = 0; i < num_points; i++)
  {
    const double *x = &points[i*2];
    index[i] = 0;     // outside
    for (unsigned int j = 0; j < boundary_size; j++)
    { // check if this point is inside our triangle
      const double *p1 = &boundary[j*6];
      const double *p2 = &boundary[j*6+2];
      const double *p3 = &boundary[j*6+4];
      // compute the area of this triangle
      double a = triangle_area(p1,p2,p3);
      double a1 = triangle_area(x,p2,p3);
      double a2 = triangle_area(p1,x,p3);
      double a3 = triangle_area(p1,p2,x);
      if (a < tolerance) // basically zero area "triangle" -> ignore it.
        continue;
      if (a1 + a2 + a3 <= a * (1 + tolerance))
      { // inside
        index[i] = 1;
        break;
      }
    }
  }
}

#ifndef NMATLAB
// function index = get_boundary_index(grid, boundary)

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
  unsigned int num_points = mxGetN(prhs[0]);
  //printf("Number of points: %i\n", num_points);

  // create boundary index matrix
  plhs[0] = mxCreateDoubleMatrix(1, num_points, mxREAL);
  double *index = mxGetPr(plhs[0]);

  // input variables
  double *points = mxGetPr(prhs[0]);
  double *boundary = mxGetPr(prhs[1]);

  // check if it's triangular or rectangular boundary correction
  int boundary_type = mxGetM(prhs[1]);
  unsigned int boundary_size = mxGetN(prhs[1]);
  //printf("Boundary Type: %i\n", boundary_type);
  //printf("Boundary Size: %i\n", boundary_size);

  if (boundary_type == 2)
  {
    get_boundary_index_triangle(num_points, points, boundary_size / 3, boundary, index);
  }
  else
  {
    get_boundary_index_rect(num_points, points, boundary_size, boundary, index);
  }
}
#endif
