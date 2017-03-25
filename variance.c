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

double sigma_m_int_kernel(double k, void *param)
{	
  double R = *(double *) param;		
  return 0.5 / PI / PI * pow(k, ns + 2.0) * prim_tk(k) * prim_tk(k) * win_1(k * R) * win_1(k * R);	
}      /* end sigma_m_before_int_kernel*/

double unnorm_sigma_m_sq(double M)
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
              
   double R = radius_from_mass(M);		
   double result, error;
   		
   gsl_function F;
   F.function = &sigma_m_int_kernel;
   F.params = &R;

   gsl_integration_qagiu(&F, 0, 0, 1e-6, 10000, w, &result, &error);   
   
   gsl_integration_workspace_free(w);   
   return result;
}     /* end sigma_m */	

void init_sigma_m()
{
   Aps = sigma_8 * sigma_8 / unnorm_sigma_m_sq(mass_from_r(8.0));	
}  /* end init_sigma_m */	

double sigma_m_sq(double M)
{
  return Aps * unnorm_sigma_m_sq(M);		
}    /* end sigma_m */	

