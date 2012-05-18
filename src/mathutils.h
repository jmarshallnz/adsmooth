#pragma once

#define ONE_ON_SQRT_2 0.70710678118654752440084436210485

static const double rel_error = 1E-12;

#ifndef HAS_R
double erfc(double x);
double erf(double x);
bool isnan(double x);
#endif
double phi(double x);

