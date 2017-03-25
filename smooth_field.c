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


/**********************************     CDM  ***************************/
double cdm_halo_density_kernel(double M, void *param)
{
   return M * ncdm(M);	
}	/* end cdm_halo_density_kernel */

double cdm_halo_density()
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
             	
   double result, error;
   		
   gsl_function F;
   F.function = &cdm_halo_density_kernel;
   
   gsl_integration_qag(&F, Mmin, Mmax, 0, 1e-4, 10000, GSL_INTEG_GAUSS21, w, &result, &error);
   
   gsl_integration_workspace_free(w);   
   return result;		
}   /* end cdm_halo_density */	

double cdm_smooth_bias_kernel(double M, void *param)
{
   return M * ncdm(M) * bcdm(M);		
}   /* end ddm_smooth_bias_kernel */	

double cdm_smooth_bias()
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
             	
   double result, error;
   double beff;
   		
   gsl_function F;
   F.function = &cdm_smooth_bias_kernel;
   
   gsl_integration_qag(&F, Mmin, Mmax, 0, 1e-4, 10000, GSL_INTEG_GAUSS21, w, &result, &error);

   gsl_integration_workspace_free(w);
      
   beff = result / rho_halo_cdm;	  
   
   return (1.0 - halo_mass_frac_cdm * beff) / (1.0 - halo_mass_frac_cdm);	
}    /* end cdm_smooth_bias */	

/**********************************  DDM  ***************************/
double ddm_halo_density_kernel(double M, void *param)
{ 
   return M * nddm(M);	
}	/* end ddm_halo_density_kernel */

double ddm_halo_density()
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
             	
   double result, error;
   		
   gsl_function F;
   F.function = &ddm_halo_density_kernel;
   
   gsl_integration_qag(&F, Mmin, Mmax, 0, 1e-4, 10000, GSL_INTEG_GAUSS21, w, &result, &error);

   gsl_integration_workspace_free(w);
      
   return result;		
}   /* end ddm_halo_density */	

double ddm_smooth_bias_kernel(double M, void *param)
{
   return M * nddm(M) * bddm(M);		
}   /* end ddm_smooth_bias_kernel */	

double ddm_smooth_bias()
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
             	
   double result, error;
   double beff;
   		
   gsl_function F;
   F.function = &ddm_smooth_bias_kernel;
   
   gsl_integration_qag(&F, Mmin, Mmax, 0, 1e-4, 10000, GSL_INTEG_GAUSS21, w, &result, &error);

   gsl_integration_workspace_free(w);
      
   beff = result / rho_halo_ddm;	
		
   return (1.0 - halo_mass_frac_ddm * beff) / (1.0 - halo_mass_frac_ddm);	
}    /* end ddm_smooth_bias */	


/***************************  initialization function *******************/
void init_smooth_field()
{
  rho_halo_cdm = cdm_halo_density();
  halo_mass_frac_cdm = rho_halo_cdm / rho_m;
  bs_cdm = cdm_smooth_bias();
  
  rho_halo_ddm = ddm_halo_density();
  halo_mass_frac_ddm = rho_halo_ddm / rho_m;
  bs_ddm = ddm_smooth_bias();
  
}  /* end init_smooth_field */
