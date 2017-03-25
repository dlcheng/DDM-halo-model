#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "Spline_2d.h"

spline_2d * spline_2d_alloc(int nx, int ny)
{
   int i;
   
   spline_2d *p;
   p = (spline_2d *) malloc(sizeof(spline_2d));
   p->psp = (gsl_spline **) malloc(sizeof(gsl_spline *) * nx);
   p->nx = nx;
   p->ny = ny;
   p->data_x = (double *) malloc(sizeof(double) * nx);
   p->data_y = (double *) malloc(sizeof(double) * ny);
   p->data_zx = (double **) malloc(sizeof(double *) * nx);	
   
   for(i=0; i<nx; i++)
     {
	 p->data_zx[i] = (double *) malloc(sizeof(double) * ny); 	 
     p->psp[i] = gsl_spline_alloc(gsl_interp_cspline, ny);
     }
   return p;	
}	 

void spline_2d_free(spline_2d * p)
{
  int i;
  
  for(i=0; i<p->nx; i++)
    {
	gsl_spline_free(p->psp[i]);   /* release the all the splines first */
	free(p->data_zx[i]);
    }
    
  free(p->psp);
  free(p->data_x);
  free(p->data_y);
  free(p->data_zx);
  free(p);
}  	

void spline_2d_init(spline_2d *p, double *x, double *y, double *z)
{
   int i, j;
   
   for(i=0; i<p->nx; i++)
     p->data_x[i] = x[i];
     
   for(i=0; i<p->ny; i++)
     p->data_y[i] = y[i];
     
   for(i=0; i<p->nx; i++)
     for(j=0; j<p->ny; j++)
        p->data_zx[i][j] = z[i * (p->ny) + j];    
   
   for(i=0; i<p->nx; i++)
	  gsl_spline_init(p->psp[i], p->data_y, p->data_zx[i], p->ny); 
		    
}	

double spline_2d_eval(spline_2d * p, double x, double y)
{
  gsl_interp_accel * acc_globel = gsl_interp_accel_alloc();
  gsl_interp_accel * acc_local  = gsl_interp_accel_alloc();
  gsl_spline *local_sp = gsl_spline_alloc(gsl_interp_cspline, p->nx);
  
  int i;
  double *temp_data;
  double result;
  
  temp_data = (double *) malloc(sizeof(double) * (p->nx));
  
  for(i=0; i<p->nx; i++)
     temp_data[i] = gsl_spline_eval(p->psp[i], y, acc_globel);
     
  gsl_spline_init(local_sp, p->data_x, temp_data, p->nx);
  
  result = gsl_spline_eval(local_sp, x, acc_local);
  
  free(temp_data);
  gsl_spline_free(local_sp);
  gsl_interp_accel_free(acc_globel);
  gsl_interp_accel_free(acc_local);
  
  return result;		
}	

