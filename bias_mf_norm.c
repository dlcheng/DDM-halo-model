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

double fst_v(double v)
{	
  double result = Ast * sqrt(2.0/PI) * sqrt(qst) * v;
  
  result *= 1.0 + pow(sqrt(qst) * v, -2.0 * pst);
  result *= exp(-qst / 2.0 * v * v);	
	
  return result;
}      /* end fst_v */	

double mf_norm_intg_kernel(double v, void *param)
{
   return fst_v(v) / v;	
}     /* end mf_norm_int_kernel */


void mf_normalization()
{
  double vmin;
  double result, error;
  
/* ns >=1, sigma(0) -> infinity  */  
  if(ns >=1)
    vmin = 0;
  else
    vmin = delta_z / sqrt(sigma_m_sq(0));  

   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);    
   gsl_function F;
   F.function = &mf_norm_intg_kernel;

   gsl_integration_qagiu(&F, vmin, 0, 1e-4, 10000, w, &result, &error);     	
   
   gsl_integration_workspace_free(w);
  
   mf_norm_fac = 1.0 / result;
}    /* end mf_normalization */

double bst_cdm(double v)
{ 
  double result = 1.0 + (qst * v * v - 1.0 ) / delta_z;
  result += 2.0 * pst / delta_z / (1.0 + pow(qst * v * v, pst));
  	
  return result;	
}       /* end bst_cdm */	

double bias_norm_intg_kernel(double v, void * param)
{
  return fst_v(v) / v * bst_cdm(v) * mf_norm_fac;	
}    /* end bias_norm_intg_kernel */	

void bias_normalization()
{
  double vmin, result, error;
  
  if(ns >=1)
    vmin  = 0;
  else
    vmin = delta_z / sqrt(sigma_m_sq(0));  	

   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);    
   gsl_function F;
   F.function = &bias_norm_intg_kernel;

   gsl_integration_qagiu(&F, vmin, 0, 1e-4, 10000, w, &result, &error);     	
   
   gsl_integration_workspace_free(w);
   
   bias_norm_fac = 1.0 / result;	
	
}    /* end bias_normalization */	
