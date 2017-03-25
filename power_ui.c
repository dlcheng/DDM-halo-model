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

double ps_lin(double k, double z)
{
  return Aps * pow(k, ns) * prim_tk(k) * prim_tk(k) * growth_factor(z) * growth_factor(z);	
}   /* end power_spectrum */

double ps_ddm_from_spline(double k, double z)
{
  double tk;	

  tk = tk_ddm_from_spline(k, z);
  
  return  ps_nonlin_from_halofit(k, z) * tk * tk;
}  /* end ps_ddm_from_spline */	

double tk_ddm_from_spline(double k, double z)
{
  double tk;
  double a = 1.0 / (1.0 + z);
  double log_a = log10(a);
  double log_k = log10(k);

  if(k < kmin_2d)
    {
	tk = 1;
	return tk;	
	}	
	
  tk = spline_2d_eval(sp_ddm_tk_2d, log_a, log_k);
  
  if(tk < 0)
    tk = 0;
  
  return  tk;
}  /* end ps_ddm_from_spline */	

/* Interface for halofit */
double ps_lin_from_halofit(double k, double z)
{
  return pow(2998, 3) * P_L(0.9999/(1.0+z), k*2998);
}   /* end ps_lin_from_halofit */	

double ps_nonlin_from_halofit(double k, double z)
{
  return pow(2998, 3) * P_NL(0.9999/(1.0+z), k*2998);	
}   /* end ps_nonlin_from_halofit */	

double ps_nonlin_from_mod_halofit(double k, double z)
{
  double y = k / 10.0;
  return ps_lin_from_halofit(k, z) + (ps_nonlin_from_halofit(k, z) - ps_lin_from_halofit(k, z)) * (1.0 + 2 * y * y) / (1.0 + 2 * y);
}  /* end Phalofit_plus_nonlin, with corrections at small scales */	
