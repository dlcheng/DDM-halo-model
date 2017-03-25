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


/* the fraction of decayed dark matter */
double decay_fraction(double a)
{
  double time_at_a = time_ing(a) * 9.77792354298172 / H0;
  double result = -1.0 * dp_t_rev * log(2) * time_at_a;

  result = exp(result);
  result = 1.0 - result;
  
  result *= (1.0 - initial_dm_frac);
  
  return result; 
	
}  /* end decay_fraction */

double time_ing_kenerl(double a, void *param)
{
  double result = Omega_m / a + Omega_v * a * a;
  
  result = 1.0 / sqrt(result);
  
  return result;
}  /* end decay_fraction_ing_kenerl */

double time_ing(double a)
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
             	
  double result, error;
   		
  gsl_function F;
  F.function = &time_ing_kenerl;
   
  gsl_integration_qag(&F, 0, a, 0, 1e-8, 10000, GSL_INTEG_GAUSS21, w, &result, &error);
   
  gsl_integration_workspace_free(w);
  
  return result;			
	
}  /* end decay_fraction_ing */

double length_ing_kenerl(double a, void * param)
{
 double a0 = *(double *) param;
 double result;
 
 result = a0 / pow(a, 3) / sqrt(Omega_m/pow(a, 3) + Omega_v);
 
 return result;	
}  /* end length_ing_kenerl */
	
double length_ing(double a0, void * param)
{
  double a1 = *(double *) param;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
             	
  double result, error;
   		
  gsl_function F;
  F.function = &length_ing_kenerl;
  F.params = &a0;
   
  gsl_integration_qag(&F, a0, a1, 0, 1e-8, 10000, GSL_INTEG_GAUSS21, w, &result, &error);
   
  gsl_integration_workspace_free(w);  	
	
  return result;	
}  /* end length_ing */	

double length_ing_min(double a0, void * param)
{
  return -1.0 * length_ing(a0, param);	
}	

double max_length(double a)
{
  double a1 = a;
  double a_low = 0;
  double a_up = a1;
  double a_guess = (a_low + a_up) / 2.0;
  
  int status;
  int iter = 0, max_iter = 1000;
  const gsl_min_fminimizer_type * T;
  gsl_min_fminimizer *s;
  gsl_function F;
  
  F.function = &length_ing_min;
  F.params = &a1;
  
  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc(T);
  gsl_min_fminimizer_set(s, &F, a_guess, a_low, a_up);
  
  do
    {
	 iter++;
	 status = gsl_min_fminimizer_iterate(s);	
		
	 a_guess = gsl_min_fminimizer_x_minimum(s);
	 a_low   = gsl_min_fminimizer_x_lower(s);	
	 a_up    = gsl_min_fminimizer_x_upper(s);
	 
	 status  = gsl_min_test_interval(a_low, a_up, 0, 1e-4); 	        
	}	
  while(status == GSL_CONTINUE && iter < max_iter);
  
  gsl_min_fminimizer_free(s);	
	
  return dp_v / 100.0 * length_ing(a_guess, &a1);
}  /* end max_length */	

double log_max_mass(double a)
{
  double r_max = max_length(a);
  
  double m_max = mass_from_r(r_max);
  
  return log10(m_max);	
} /* end log_max_mass */
