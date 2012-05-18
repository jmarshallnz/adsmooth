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
#include <cstdio>
#include "mathutils.h"

#define NMATLAB

#ifndef NMATLAB
#include "mex.h"      // for matlab
#endif

double interpolate(unsigned int num_gridX, unsigned int num_gridY, const double *grid, const double *function, const double *point)
{
  // first check whether it's outside our grid
  if (grid[0] > point[0] || grid[1] > point[1] ||
      grid[2*(num_gridX*num_gridY - 1)] < point[0] || grid[2*(num_gridX*num_gridY - 1) + 1] < point[1])
  { // outside - we can't extrapolate, so simply set to 0
    printf("ERROR: Point (%f,%f) is outside: (%f,%f)->(%f,%f)\n", point[0], point[1], grid[0], grid[1], grid[2*(num_gridX*num_gridY - 1)], grid[2*(num_gridX*num_gridY - 1) + 1]);
    return 0;
  }

  // ok, it's inside the grid we have covered, so find the grid location
  unsigned int x = num_gridX - 1;
  while (x && grid[2*x*num_gridY] > point[0])
    x--;

  unsigned int y = num_gridY - 1;
  while (y && grid[2*y+1] > point[1])
    y--;

//  printf("Point (%f,%f) is inside rect at (%i,%i): (%f,%f)->(%f,%f)\n", point[0], point[1], x, y, grid[2*(x*num_gridY + y)], grid[2*(x*num_gridY + y) + 1], grid[2*((x+1)*num_gridY + y+1)], grid[2*((x+1)*num_gridY + y+1) + 1]);

  // interpolate from here
  double w1 = (point[0] - grid[2*(x*num_gridY + y)]) * (point[1] - grid[2*(x*num_gridY + y) + 1]);
  double w2 = (grid[2*((x+1)*num_gridY + y+1)] - point[0]) * (point[1] - grid[2*(x*num_gridY + y) + 1]);
  double w3 = (point[0] - grid[2*(x*num_gridY + y)]) * (grid[2*((x+1)*num_gridY + y+1) + 1] - point[1]);
  double w4 = (grid[2*((x+1)*num_gridY + y+1)] - point[0]) * (grid[2*((x+1)*num_gridY + y+1) + 1] - point[1]);
  double w = (grid[2*((x+1)*num_gridY + y+1)] - grid[2*(x*num_gridY + y)]) * (grid[2*((x+1)*num_gridY + y+1) + 1] - grid[2*(x*num_gridY + y) + 1]);
//  printf("Weights are (%f,%f),(%f,%f)\n", w4/w, w3/w, w2/w, w1/w);

  return (w4 * function[x*num_gridY + y] + w3 * function[(x+1)*num_gridY + y] + w2 * function[x*num_gridY + y+1] + w1 * function[(x+1)*num_gridY + y+1])/w;
}

void get_adaptive_bandwidth(unsigned int num_points, const double *points, const double bandwidth, unsigned int num_gridX, unsigned int num_gridY, const double *grid, double *pilot, double *h_points, double *h_grid)
{
  // min and max adaptivity
  static const double min_h_scale = 0.3;
  static const double max_h_scale = 3.0;

  // first make sure our pilot is positive (set to the min), and dump any NaN's to zero
  double minPilot = std::numeric_limits<double>::max();
  for (unsigned int i = 0; i < num_gridX * num_gridY; i++)
    if (pilot[i] < minPilot && pilot[i] > 0) minPilot = pilot[i];
  //printf("minimum pilot = %f\n", minPilot);
  
  for (unsigned int i = 0; i < num_gridX * num_gridY; i++)
  {
    if (isnan(pilot[i]))
      pilot[i] = 0;
    else if (pilot[i] <= 0)
      pilot[i] = minPilot;
  }

  // now interpolate across the grid to the data points, and calculate the geometric mean
  double one_on_num_points = 1.0/num_points;
  double geom_mean = 1.0;
  for (unsigned int i = 0; i < num_points; i++)
  {
    h_points[i] = interpolate(num_gridX, num_gridY, grid, pilot, &points[2*i]);
    geom_mean *= pow(h_points[i], one_on_num_points);
  }
  //printf("geometric mean = %f\n", 1000000*geom_mean);
  //for (unsigned int i = 0; i < num_points; i++)
  //  printf("point (%f,%f) interpolates to %f\n", points[2*i], points[2*i+1], 1000000*h_points[i]);

  // scale by the geometric mean, bound the scale of adaptibility, and scale by the bandwidth
  for (unsigned int i = 0; i < num_points; i++)
  {
    h_points[i] = sqrt(geom_mean / h_points[i]);
    if (h_points[i] > max_h_scale) h_points[i] = max_h_scale;
    if (h_points[i] < min_h_scale) h_points[i] = min_h_scale;
    h_points[i] *= bandwidth;
  }
  //for (unsigned int i = 0; i < num_points; i++)
  //  printf("bandwidth at (%f,%f) is %f\n", points[2*i], points[2*i+1], 1000000*h_points[i]);

  // ok, now finally compute the bandwidth on the grid
  for (unsigned int i = 0; i < num_gridX * num_gridY; i++)
  {
    h_grid[i] = sqrt(geom_mean / pilot[i]);
    if (h_grid[i] > max_h_scale) h_grid[i] = max_h_scale;
    if (h_grid[i] < min_h_scale) h_grid[i] = min_h_scale;
    h_grid[i] *= bandwidth;
  }
}

#ifndef NMATLAB
// function [h_data h_grid] = get_adaptive_bandwidth(points, h, grid, pilot)

void mexFunction( int nlhs,mxArray *plhs[],int nrhs, mxArray *prhs[])
{
  // Check for proper number of arguments.
  if(nrhs < 4)
  {
    mexErrMsgTxt("Four inputs required.");
  }
  else if(nlhs>2)
  {
    mexErrMsgTxt("Too many output arguments");
  }

  // calculate size of input data
  unsigned int num_points = mxGetN(prhs[0]);
  //printf("Number of points: %i\n", num_points);

  unsigned int num_gridX = (unsigned int)sqrt((double)mxGetN(prhs[2]));
  unsigned int num_gridY = num_gridX;

  // create bandwidth matrix for points and grid bandwidths
  plhs[0] = mxCreateDoubleMatrix(1, num_points, mxREAL);
  double *h_points = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(1, num_gridX*num_gridY, mxREAL);
  double *h_grid = mxGetPr(plhs[1]);

  // input variables
  double *points = mxGetPr(prhs[0]);
  double bandwidth = *mxGetPr(prhs[1]);
  double *grid = mxGetPr(prhs[2]);
  double *pilot = mxGetPr(prhs[3]);

  // compute the adaptive bandwidth
  get_adaptive_bandwidth(num_points, points, bandwidth, num_gridX, num_gridY, grid, pilot, h_points, h_grid);
}
#endif
