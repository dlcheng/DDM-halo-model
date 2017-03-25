#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>

#include "allvars.h"
#include "proto.h"


void init_cl(double z)  /* the background galaxy at z */
{
  double a = 1.0 / (1.0 + z);
  int i;	
  int size = 100;
	
  a_cl = (double *) malloc(sizeof(double) * size);
  x_cl = (double *) malloc(sizeof(double) * size);
  
  double da;
  da = (1.0 - a) / ((double) size - 1.0);
  
  for(i=0; i<size; i++)
    {
    a_cl[i] = 1.0 - da * i;
	x_cl[i] = comoving_dis_ing(a_cl[i]);
    }
    
   xs = x_cl[size -1];
   
   sp_ax_cl = gsl_spline_alloc(gsl_interp_cspline, size);
   spacc_ax_cl  = gsl_interp_accel_alloc();
   
   gsl_spline_init(sp_ax_cl, x_cl, a_cl, size);  /* give x return a */
    
}  /* end init_cl */

void free_cl()
{
  free(a_cl);
  free(x_cl);
  gsl_spline_free(sp_ax_cl);
  gsl_interp_accel_free(spacc_ax_cl);  		
}  /* end free_cl */	

double a_from_x(double x)
{
  double result = gsl_spline_eval(sp_ax_cl, x, spacc_ax_cl);
  
  return result;	
	
}  /* end a_from_x */	

double comoving_dis_ing_kernel(double a, void * param)
{
  double result;
	
  result = 1.0 / a / a;
  result /= sqrt(Omega_m/pow(a, 3) + Omega_v);
  
  return result;
}   /* end comoving_dis_ing_kernel */

double comoving_dis_ing(double a)
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
             	
   double result, error;
   		
   gsl_function F;
   F.function = &comoving_dis_ing_kernel;
   
   gsl_integration_qag(&F, a, 1, 0, 1e-5, 10000, GSL_INTEG_GAUSS21, w, &result, &error);
   
   gsl_integration_workspace_free(w);
   
   return result * 3e3;			
		
}   /* end comoving_dis_ing */	

double cl_ing_kernel(double x, void * param)
{
  double *p = (double *) param;
  
  double l = p[0];
  int flag = (int) p[1];
  
  double a = a_from_x(x);
  double z = 1.0 / a - 1.0;
  double D1 = xs - x;
  double D2 = xs; 
  double Dx = x;
  double k;
  
  if(Dx < 1e-3)
    return 0;    
  else
    k = l / Dx;
  
  if(k > kmax_2d)
    return 0;
   
  double pk;
   
  if(flag == 0)
    pk = ps_nonlin_from_halofit(k, z);
  else 
    {
    if(flag == 1)
     pk = ps_ddm_from_spline(k, z);
    else
      {
      if(flag == 2);
        pk = ps_lin(k, z); 
      }
     }         
       
  return 1.0 / pow(a, 2) * pow(D1/D2, 2) * pk;    
   	
}  /* end cl_ing_kernel */	

void cl_cal(double l, double lm, double ln, int flag, double *p_result)
{

/*
 * flag == 0, CDM
 * flag == 1, DDM
 * flag == 2, Linear 
*/ 	
   double p[2];
   double cl;
   double delta_cl;
   double fsky = 0.5;
   double n_bar = 35 * (4.6656e8/ PI) / 4.0 / PI; /* n_bar original in unit of per square arcmin, translated to
                                                     per steradian.
                                                   */
   double gamma_r = 0.35;
   
   p[0] = l;
   p[1] = (double) flag;
   
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
             	
   double result, error;
   		
   gsl_function F;
   F.function = &cl_ing_kernel;
   F.params = p;
   
   gsl_integration_qag(&F, 0, xs, 0, 1e-5, 10000, GSL_INTEG_GAUSS21, w, &result, &error);
   
   gsl_integration_workspace_free(w);	
	
   cl = result / 36.0 * pow(Omega_m, 2) * 1e-12;
   
   double Nl = (ln + lm + 1) * (ln - lm + 1) / 2;
   
   delta_cl = sqrt(1.0 / Nl / fsky) * (cl + gamma_r * gamma_r / n_bar);
   
   p_result[0] = cl;
   // * l * (l + 1.0) / (2.0 * PI);
   p_result[1] = delta_cl / cl;
   
}   /* end cl_cal */	
