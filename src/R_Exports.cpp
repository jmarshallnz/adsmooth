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

#define _DLL

#ifdef _DLL

#include <stdlib.h>
#include "ConstructPDF.h"
#include "mathutils.h"
//#include "R.h"

extern "C"
{
__declspec(dllexport) void construct_pdf_R(int *num_data, double *data, int *kernel, double *bandwidth, int *num_gridX, int *num_gridY, double *grid, int *boundary_size, double *boundary, int *type, double *pdf)
{
  bool adaptive = false;
  int boundary_corrected = 0;
  if (type)
  {
    adaptive = ((*type) & 4) == 4;
    boundary_corrected = (*type) & 3;
  }
//  Rprintf("calling construct_pdf_R adaptive=%s, boundary_corrected=%d\n", adaptive ? "true" : "false", boundary_corrected);

  // first compute the boundary_index
  double *boundary_index = new double[(*num_gridX)*(*num_gridY)];
  get_boundary_index_triangle((*num_gridX)*(*num_gridY), grid, *boundary_size, boundary, boundary_index);
//  Rprintf("successfully called get_boundary_index_triangle\n");

  // now construct our rectangular boundary
  double *rects = NULL;
  unsigned int num_rects = construct_rects_from_grid(*num_gridX, *num_gridY, grid, boundary_index, &rects);
//  Rprintf("successfully called construct_rects_from_grid, creating %d rectangles\n", num_rects);
//  for (unsigned int i = 0; i < num_rects; i++)
//    Rprintf("rect (i): (%f,%f)->(%f,%f)\n", rects[4*i], rects[4*i+1], rects[4*i+2], rects[4*i+3]);

  double *points = alloc(*num_data * 2);
  unsigned int num_points = get_points_inside_rect(*num_data, data, num_rects, rects, points);

  construct_pdf(num_points, points, *kernel, *bandwidth, (unsigned int)*num_gridX, (unsigned int)*num_gridY, grid, num_rects, rects, boundary_index, adaptive, boundary_corrected, pdf);

  dealloc(points);

  // and free up our created data
  delete[] rects;
  delete[] boundary_index;
}

__declspec(dllexport) void construct_boundary_index_R(int *num_gridX, int *num_gridY, double *grid, int *num_triangles, double *triangulation, double *index)
{
  get_boundary_index_triangle((*num_gridX)*(*num_gridY), grid, *num_triangles, triangulation, index);
}

__declspec(dllexport) void estimate_mise_R(int *num_data, double *data, int *kernel, double *bandwidth, int *num_gridX, int *num_gridY, double *grid, int *num_triangles, double *triangulation, int *type, double *mise)
{
//  Rprintf("calling construct_pdf_R with num_data = %d, num_gridX = %d, num_gridY = %d, boundary_size = %d, type = %d\n", *num_data, *num_gridX, *num_gridY, *boundary_size, *type);
  bool adaptive = false;
  int boundary_corrected = 0;
  if (type)
  {
    adaptive = ((*type) & 4) == 4;
    boundary_corrected = (*type) & 3;
  }
//  Rprintf("calling estimate_mise_R adaptive=%s, boundary_corrected=%d\n", adaptive ? "true" : "false", boundary_corrected);

  // first compute the boundary_index
  double *boundary_index = new double[(*num_gridX)*(*num_gridY)];
  get_boundary_index_triangle((*num_gridX)*(*num_gridY), grid, *num_triangles, triangulation, boundary_index);
//  Rprintf("successfully called get_boundary_index_triangle\n");

  // now construct our rectangular boundary
  double *rects = NULL;
  unsigned int num_rects = construct_rects_from_grid(*num_gridX, *num_gridY, grid, boundary_index, &rects);
//  Rprintf("successfully called construct_rects_from_grid, creating %d rectangles\n", num_rects);

  // throw away data outside our boundary
  double *points = alloc(*num_data * 2);
  unsigned int num_points = get_points_inside_rect(*num_data, data, num_rects, rects, points);

  *mise = estimate_mise(num_points, points, *kernel, *bandwidth, *num_gridX, *num_gridY, grid, num_rects, rects, boundary_index, adaptive, boundary_corrected);

  dealloc(points);

  // and free up our created data
  delete[] rects;
  delete[] boundary_index;
}

__declspec(dllexport) void estimate_rr_mise_R(int *num_data_f, double *data_f, int *num_data_g, double *data_g, int *kernel, double *bandwidth, double *delta, int *num_gridX, int *num_gridY, double *grid, int *num_triangles, double *triangulation, int *type, double *mise)
{
//  Rprintf("calling construct_pdf_R with num_data = %d, num_gridX = %d, num_gridY = %d, boundary_size = %d, type = %d\n", *num_data, *num_gridX, *num_gridY, *boundary_size, *type);
  bool adaptive = false;
  bool logarithmic = true;
  int boundary_corrected = 0;
  if (type)
  {
    logarithmic = ((*type) & 8) == 8;
    adaptive = ((*type) & 4) == 4;
    boundary_corrected = (*type) & 3;
  }
//  Rprintf("calling estimate_mise_R adaptive=%s, boundary_corrected=%d\n", adaptive ? "true" : "false", boundary_corrected);

  // first compute the boundary_index
  double *boundary_index = new double[(*num_gridX)*(*num_gridY)];
  get_boundary_index_triangle((*num_gridX)*(*num_gridY), grid, *num_triangles, triangulation, boundary_index);
//  Rprintf("successfully called get_boundary_index_triangle\n");

  // now construct our rectangular boundary
  double *rects = NULL;
  unsigned int num_rects = construct_rects_from_grid(*num_gridX, *num_gridY, grid, boundary_index, &rects);
//  Rprintf("successfully called construct_rects_from_grid, creating %d rectangles\n", num_rects);

  // throw away data outside our boundary
  double *points_f = alloc(*num_data_f * 2);
  unsigned int num_points_f = get_points_inside_rect(*num_data_f, data_f, num_rects, rects, points_f);
  double *points_g = alloc(*num_data_g * 2);
  unsigned int num_points_g = get_points_inside_rect(*num_data_g, data_g, num_rects, rects, points_g);

  *mise = estimate_rr_mise(num_points_f, points_f, num_points_g, points_g, *kernel, *bandwidth, *delta, *num_gridX, *num_gridY, grid, num_rects, rects, boundary_index, adaptive, boundary_corrected, logarithmic);

  dealloc(points_f);
  dealloc(points_g);

  // and free up our created data
  delete[] rects;
  delete[] boundary_index;
}

__declspec(dllexport) void choose_delta_R(int *num_data_f, double *data_f, int *num_data_g, double *data_g, int *kernel, double *bandwidth, int *num_gridX, int *num_gridY, double *grid, int *num_triangles, double *triangulation, int *type, double *delta)
{
//  Rprintf("calling construct_pdf_R with num_data = %d, num_gridX = %d, num_gridY = %d, boundary_size = %d, type = %d\n", *num_data, *num_gridX, *num_gridY, *boundary_size, *type);
  bool adaptive = false;
  int boundary_corrected = 0;
  if (type)
  {
    adaptive = ((*type) & 4) == 4;
    boundary_corrected = (*type) & 3;
  }
//  Rprintf("calling estimate_mise_R adaptive=%s, boundary_corrected=%d\n", adaptive ? "true" : "false", boundary_corrected);

  // first compute the boundary_index
  double *boundary_index = new double[(*num_gridX)*(*num_gridY)];
  get_boundary_index_triangle((*num_gridX)*(*num_gridY), grid, *num_triangles, triangulation, boundary_index);
//  Rprintf("successfully called get_boundary_index_triangle\n");

  // now construct our rectangular boundary
  double *rects = NULL;
  unsigned int num_rects = construct_rects_from_grid(*num_gridX, *num_gridY, grid, boundary_index, &rects);
//  Rprintf("successfully called construct_rects_from_grid, creating %d rectangles\n", num_rects);

  // throw away data outside our boundary
  double *points_f = alloc(*num_data_f * 2);
  unsigned int num_points_f = get_points_inside_rect(*num_data_f, data_f, num_rects, rects, points_f);
  double *points_g = alloc(*num_data_g * 2);
  unsigned int num_points_g = get_points_inside_rect(*num_data_g, data_g, num_rects, rects, points_g);
  
  // grab our h_data_f, bc_data_f, and bc_grid_f
  *delta = choose_delta(num_points_f, points_f, num_points_g, points_g, *kernel, *bandwidth, *num_gridX, *num_gridY, grid, num_rects, rects, boundary_index, adaptive, boundary_corrected);

  dealloc(points_f);
  dealloc(points_g);

  // and free up our created data
  delete[] rects;
  delete[] boundary_index;
}
}
#endif