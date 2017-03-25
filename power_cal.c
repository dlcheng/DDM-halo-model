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


double tk_ddm_from_hm(double k, double z)
{
  double temp_cdm[6];
  double temp_ddm[6];
  
  calculate_pk(k, z, temp_cdm, 0);
  calculate_pk(k, z, temp_ddm, 1);
  
  return sqrt(temp_ddm[5]/temp_cdm[5]);	
}  	

void calculate_pk(double k, double z, double *temp, int flag)
{
/* 
   flag = 0   ----->   CDM	
   flag = 1   ----->   DDM
*/
  double f;
  double bs;
  double rho_h;
  double Pnl_kR;
  double HH_factor, H_factor;
  
  if(flag == 0)
    {
    f = halo_mass_frac_cdm;
    bs = bs_cdm;
    rho_h = rho_halo_cdm;
    }
  else 
    {
    f = halo_mass_frac_ddm;
    bs = bs_ddm;
    rho_h = rho_halo_ddm;
    } 

   Pnl_kR = ps_lin(k, z);
	
   HH_factor = halo_halo_factor(k, z, flag);
   H_factor  = halo_factor(k, z, flag);

   temp[0] = k;
   temp[1] = bs * bs * Pnl_kR;
   temp[1] *= (1.0 - f) * (1.0 - f);
   temp[2] = bs * Pnl_kR / rho_h * HH_factor;
   temp[2] *= 2.0 * f * (1.0 - f);
   temp[3] = Pnl_kR / rho_h / rho_h * HH_factor * HH_factor;
   temp[3] *= f * f;
   temp[4] = H_factor / rho_h / rho_h;
   temp[4] *= f * f;
   temp[5] = temp[1] + temp[2] + temp[3] + temp[4];
      
}    /* end calculate_pk */	


double halo_halo_factor_intg_kernel(double M, void *param)
{
   double *p = (double *) param;
   
   double k = p[0];
   double z = p[1];
   int flag = (int) p[2];
   
   double n, b, u;
   u = factor_uk(k, M, z, flag);
   
   if(flag < 1)
     {
	   n = ncdm(M);
	   b = bcdm(M);	 
	   }	 
   else
     {
	   n = nddm(M);
	   b = bddm(M);	 
	   }	   
	 
   return M * n * b * u;		
}     /* end halo_halo_factor_intg_kernel */


double halo_halo_factor(double k, double z, int flag)
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
   gsl_set_error_handler_off();
          	
   double result, error;
   double p[3];
   
   p[0] = k;
   p[1] = z;
   p[2] = flag;
   	
   gsl_function F;
   F.function = &halo_halo_factor_intg_kernel;
   F.params = p;
     
   gsl_integration_qag(&F, Mmin, Mmax, 0, 1e-4, 10000, GSL_INTEG_GAUSS21, w, &result, &error);

   gsl_integration_workspace_free(w);
   return result;			
		
}     /* end halo_halo_factor */

double halo_factor_intg_kernel(double M, void *param)
{
   double *p = (double *) param;
   
   double k = p[0];
   double z = p[1];
   int flag = (int) p[2];
   
   double n, u;
   
   u = factor_uk(k, M, z, flag);
   
   if(flag < 1)
     n = ncdm(M);
   else  
     n = nddm(M);

   return n * M * M * u * u;			
	
}     /* end halo_factor_intg_kernel */

double halo_factor(double k, double z, int flag)
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
   gsl_set_error_handler_off();
             	
   double result, error;
   double p[3];
   
   p[0] = k;
   p[1] = z;
   p[2] = flag;
   	
   gsl_function F;
   F.function = &halo_factor_intg_kernel;
   F.params = p;
     
   gsl_integration_qag(&F, Mmin, Mmax, 0, 1e-4, 10000, GSL_INTEG_GAUSS21, w, &result, &error);

   gsl_integration_workspace_free(w);   
   return result;			
}    /* end halo_factor */
