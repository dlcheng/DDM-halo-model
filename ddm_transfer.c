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

/*  Contains function to describe the mass function and c-M relation 
    change from the fitting formula.
*/

double mf_ddm_t(double M)
{
#ifdef TURN_OFF_MFT
  return 1.0;
#endif	
	
  double Mx = pow(10, b_m);
  double f_eff = a_m * decay_f;
  double Mvk = pow(10, log_ddm_norm_mass);
  
  double mfac = pow(1.0 + Mx / (M / Mvk), f_eff);
	
  return  mfac;			
}         /* end mf_ddm_t */	

double c_ddm_t(double M)
{
#ifdef TURN_OFF_CMT
  return 1.0;
#endif	
	
  double Mx = pow(10, b_cm);
  double f_eff = a_cm * decay_f;
  double Mvk = pow(10, log_ddm_norm_mass);
  
  double cfac = pow(1.0 + Mx / (M / Mvk), f_eff);
	
  return cfac;	
		
}          /* end c_ddm_t */	


