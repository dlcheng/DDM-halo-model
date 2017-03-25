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

double factor_uk(double k, double M, double z, int flag)
{
/* 
   flag = 0      ------>  CDM	
   flag = 1      ------>  DDM
*/	
  double rvir = vir_radius(M, z);            /* the co-moving virial radius_from_mass */
  double c;
  double m_fac;
  
  if(flag == 0)
    c = c_cdm(M, z);
  else 
    c = c_ddm(M, z);
    
  m_fac = factor_mass(rvir, k, c);    
 	
  return m_fac;	
}	/* end factor_uk */

double factor_mass(double R, double k, double c)
{
  double rs = R / c;
  double part_1 = sin(k*rs) * si_intg(k*rs, (1+c)*k*rs);
  double part_2;
  double part_3;
  double part;
  double c_fac;
  double result;
  
  if(c*k*rs < 1e-8)
    part_2 = -1.0 * c / (1+c);
  else
    part_2 = -1.0 * sin(c*k*rs) / ((1+c)*k*rs);
    
  part_3 = cos(k*rs) * ci_intg(k*rs, (1+c)*k*rs);
  
  part = part_1 + part_2 + part_3;  
  c_fac = log(1+c) - (c/(1+c));
	
  result = part / c_fac;
  
  if(fabs(result - 1) < 1e-5)
     result = 1;
     
  return result;   	
}   /* end factor_mass */	

double ci_intg_kernel(double x, void * param)
{
  return cos(x) / x;		
}  /* end Ci_kernerl */

double ci_intg(double x1, double x2)
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
             	
   double result, error;
   		
   gsl_function F;
   F.function = &ci_intg_kernel;
   
   gsl_integration_qag(&F, x1, x2, 0, 1e-5, 10000, GSL_INTEG_GAUSS21, w, &result, &error);
   
   gsl_integration_workspace_free(w);
   return result;		
		
}   /* end ci_intg */		

double si_intg_kernel(double x, void * param)
{
   if(x < 1e-8)
     return 1;
   
   return sin(x) / x;  	
}   /* end si_intg_kernel */

double si_intg(double x1, double x2)
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
             	
   double result, error;
   		
   gsl_function F;
   F.function = &si_intg_kernel;
   
   gsl_integration_qag(&F, x1, x2, 0, 1e-5, 10000, GSL_INTEG_GAUSS21, w, &result, &error);
   
   gsl_integration_workspace_free(w);
   return result;		
	
}   /* end si_intg */		
