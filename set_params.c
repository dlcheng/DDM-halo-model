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

void set_cosmological_parameters(double x1, double x2, double x3, double x4, double x5, double x6)
{
  Omega_m = x1;
  Omega_v = x2;
  Omega_b = x3;
  H0      = x4;
  ns      = x5;             
  sigma_8 = x6;
  delta_c = 1.686;	  
}    /* end set_cosmological_parameters */	

void set_2d_spline_parameters(double x1, double x2, int x3, double x4, double x5, int x6)
{ 
  double z_min, z_max;
  kmin_2d = x1;
  kmax_2d = x2;
  k_num = x3;  
  z_min = x4;
  z_max = x5;
  amin_2d = 1.0 / ( 1.0 + z_max);
  amax_2d = 1.0 / ( 1.0 + z_min);       
  a_num = x6;	  
}   /* end set_2d_spline_parameters */	

void set_decay_parameters(double x1, double x2, double x3)
{
  dp_t_rev = x1;        /* 1/tau in unit of Gyr^-1 */	
  dp_v     = x2;        /* kick velocity in unit of km/s */ 
  initial_dm_frac = x3; /* the initial fraction of the daughter particles */  
}    /* end set_decay_parameters */	
