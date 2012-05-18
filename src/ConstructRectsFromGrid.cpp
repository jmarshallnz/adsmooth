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
#include <stdlib.h>
#include <math.h>
#include <limits>

#define NMATLAB

#ifndef NMATLAB
#include "mex.h"      // for matlab
#endif

unsigned int generate_x_rects_from_grid(unsigned int num_gridX, unsigned int num_gridY, const double *grid, const double *boundary_index, const double *offset, double *rects)
{
  unsigned int numRects = 0;
  for (unsigned int x = 0; x < num_gridX; x++)
  {
    unsigned int startY = (unsigned int)-1;
    for (unsigned int y = 0; y < num_gridY; y++)
    {
      const double *endPoint = &grid[2*(x*num_gridY + y)];
      if (startY != (unsigned int)-1 && boundary_index[x*num_gridY + y] == 0)
      { // finished a rect from in_rect to y
        const double *startPoint = &grid[2*(x*num_gridY + startY)];
        startY = (unsigned int)-1;
        if (rects)
        { // append a rectangle to our list
          rects[4*numRects] = startPoint[0] - offset[0];
          rects[4*numRects + 1] = startPoint[1] - offset[1];
          rects[4*numRects + 2] = endPoint[0] + offset[0];
          rects[4*numRects + 3] = endPoint[1] - offset[1];
        }
        numRects++;
      }
      else if (startY == (unsigned int)-1 && boundary_index[x*num_gridY + y] != 0)
      { // starting a rect
        startY = y;
      }
    }
    if (startY != (unsigned int)-1) // ending a rect
    {
      const double *startPoint = &grid[2*(x*num_gridY + startY)];
      const double *endPoint = &grid[2*(x*num_gridY + num_gridY - 1)];
      if (rects)
      { // append a rectangle to our list
        rects[4*numRects] = startPoint[0] - offset[0];
        rects[4*numRects + 1] = startPoint[1] - offset[1];
        rects[4*numRects + 2] = endPoint[0] + offset[0];
        rects[4*numRects + 3] = endPoint[1] + offset[1];
      }
      numRects++;
    }
  }
  return numRects;
}

unsigned int generate_y_rects_from_grid(unsigned int num_gridX, unsigned int num_gridY, const double *grid, const double *boundary_index, const double *offset, double *rects)
{
  unsigned int numRects = 0;
  for (unsigned int y = 0; y < num_gridY; y++)
  {
    unsigned int startX = (unsigned int)-1;
    for (unsigned int x = 0; x < num_gridX; x++)
    {
      const double *endPoint = &grid[2*(x*num_gridY + y)];
      if (startX != (unsigned int)-1 && boundary_index[x*num_gridY + y] == 0)
      { // finished a rect from startX to x
        const double *startPoint = &grid[2*(startX*num_gridY + y)];
        startX = (unsigned int)-1;
        if (rects)
        { // append a rectangle to our list
          rects[4*numRects] = startPoint[0] - offset[0];
          rects[4*numRects + 1] = startPoint[1] - offset[1];
          rects[4*numRects + 2] = endPoint[0] - offset[0];
          rects[4*numRects + 3] = endPoint[1] + offset[1];
        }
        numRects++;
      }
      else if (startX == (unsigned int)-1 && boundary_index[x*num_gridY + y] != 0)
      { // starting a rect
        startX = x;
      }
    }
    if (startX != (unsigned int)-1) // ending a rect
    {
      const double *startPoint = &grid[2*(startX*num_gridY + y)];
      const double *endPoint = &grid[2*((num_gridX - 1)*num_gridY + y)];
      if (rects)
      { // append a rectangle to our list
        rects[4*numRects] = startPoint[0] - offset[0];
        rects[4*numRects + 1] = startPoint[1] - offset[1];
        rects[4*numRects + 2] = endPoint[0] + offset[0];
        rects[4*numRects + 3] = endPoint[1] + offset[1];
      }
      numRects++;
    }
  }
  return numRects;
}

unsigned int construct_rects_from_grid(unsigned int num_gridX, unsigned int num_gridY, const double *grid, const double *boundary_index, double **rects)
{
  // run through the restricted grid and compute our rectangles - do it in
  // both directions to minimize the number of integrations required
  double offset[2];
  offset[0] = 0.5 * (grid[2*num_gridY] - grid[0]);
  offset[1] = 0.5 * (grid[2*1+1] - grid[1]);

  // first count the number of rectangles in each direction (we want to minimize the number)
  unsigned int numRectsX = generate_x_rects_from_grid(num_gridX, num_gridY, grid, boundary_index, offset, NULL);
  unsigned int numRectsY = generate_y_rects_from_grid(num_gridX, num_gridY, grid, boundary_index, offset, NULL);

  if (numRectsX < numRectsY)
  {
    *rects = new double[numRectsX * 4];
    return generate_x_rects_from_grid(num_gridX, num_gridY, grid, boundary_index, offset, *rects);
  }
  else
  {
    *rects = new double[numRectsY * 4];
    return generate_y_rects_from_grid(num_gridX, num_gridY, grid, boundary_index, offset, *rects);
  }
}

#ifndef NMATLAB
// function rect_boundary = construct_rects_from_grid(grid, triangle_boundary);
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

  // calculate size of data points and grid points
  unsigned int num_grid = mxGetN(prhs[0]);
  unsigned int boundary_size = mxGetN(prhs[1]) / 3;

  double *rects = NULL;  // will be created below

  // input variables
  double *grid = mxGetPr(prhs[0]);
  double *boundary = mxGetPr(prhs[1]);
  unsigned int num_gridX = (unsigned int)sqrt((double)num_grid);
  unsigned int num_gridY = num_gridX;

  // generate our boundary index
  double *boundary_index = new double[num_gridX * num_gridY];
  get_boundary_index_triangle(num_gridX*num_gridY, grid, boundary_size, boundary, boundary_index);

  // and execute our function
  unsigned int num_rects = construct_rects_from_grid(num_gridX, num_gridY, grid, boundary_index, &rects);

  delete[] boundary_index;

  // create output vector
  plhs[0] = mxCreateDoubleMatrix(4, num_rects, mxREAL);
  double *dest_rects = mxGetPr(plhs[0]);

  for (unsigned int i = 0; i < 4 * num_rects; i++)
    dest_rects[i] = rects[i];

  delete[] rects;
}
#endif