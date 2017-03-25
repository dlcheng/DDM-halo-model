#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

typedef struct
{
  gsl_spline ** psp;  /* Point to spline arrays, size = nx */
  int nx;
  int ny;
  double * data_x;
  double * data_y;	
  double ** data_zx;
}spline_2d;	


spline_2d * spline_2d_alloc(int nx, int ny);	 
void spline_2d_free(spline_2d * p);	
void spline_2d_init(spline_2d *p, double *x, double *y, double *z);	
double spline_2d_eval(spline_2d * p, double x, double y);
