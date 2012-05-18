#pragma once

double *alloc(unsigned int amount);
void dealloc(double *mem);
double evaluate_pdf(unsigned int num_data, const double *data, const int kernel, const double *h_data, double x, double y, const double *boundary_correction, unsigned int ignore_entry = -1);
void evaluate_density(unsigned int num_data, const double *data, const int kernel, double *h, unsigned int num_points, const double *points, const double *boundary_correction, const double *in_boundary, double *pdf);
void construct_pdf(unsigned int num_data, const double *data, const int kernel, double bandwidth, unsigned int num_gridX, unsigned int num_gridY, const double *grid, unsigned int boundary_size, const double *boundary, const double *boundary_index, bool adaptive, int boundary_corrected, double *pdf);
void get_boundary_index_rect(unsigned int num_points, const double *points, unsigned int boundary_size, const double *boundary, double *index);
unsigned int get_points_inside_rect(unsigned int num_points, const double *points, unsigned int boundary_size, const double *boundary, double *points_inside);
void get_boundary_index_triangle(unsigned int num_points, const double *points, unsigned int boundary_size, const double *boundary, double *index);
void get_adaptive_bandwidth(unsigned int num_points, const double *points, const double bandwidth, unsigned int num_gridX, unsigned int num_gridY, const double *grid, double *pilot, double *h_points, double *h_grid);
void get_bc_rectangle(unsigned int num_points, const double *x, const int kernel, const double *h, unsigned int num_rectangles, const double *rectangles, const double *in_boundary, bool use_adaptive, int correction_type, double *correction);
double trapezoid_rule(unsigned int gridSizeX, unsigned int gridSizeY, const double *grid, const double *function);
unsigned int construct_rects_from_grid(unsigned int num_gridX, unsigned int num_gridY, const double *grid, const double *boundary_index, double **rects);
double estimate_mise(unsigned int num_data, const double *data, const int kernel, double bandwidth, unsigned int num_gridX, unsigned int num_gridY, const double *grid, unsigned int boundary_size, const double *boundary, const double *boundary_index, bool adaptive, int boundary_corrected);
double estimate_rr_mise(unsigned int num_data_f, const double *data_f, unsigned int num_data_g, const double *data_g, int kernel, double bandwidth, double delta, unsigned int num_gridX, unsigned int num_gridY, const double *grid, unsigned int boundary_size, const double *boundary, const double *boundary_index, bool adaptive, int boundary_corrected, bool logarithmic);
double optimal_delta(double *p_f, double *p_g, unsigned int num_gridX, unsigned int num_gridY);
double choose_delta(unsigned int num_data_f, const double *data_f, unsigned int num_data_g, const double *data_g, int kernel, double bandwidth, unsigned int num_gridX, unsigned int num_gridY, const double *grid, unsigned int boundary_size, const double *boundary, const double *boundary_index, bool adaptive, int boundary_corrected);
