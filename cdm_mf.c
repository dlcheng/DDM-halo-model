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


/*********************************     Spline functions related to v(M) relation *************************/
void init_spline_mu()
{
   int i;	
   double sig_m;

/* initialize LogM */	
   for(i=0; i<NI; i++)
      {	
      SplogM[i] = log10(Mmin) + i * (log10(Mmax) - log10(Mmin)) / (double)(NI - 1);
      sig_m = sqrt(sigma_m_sq(pow(10, SplogM[i])));
      SpMu[i] = delta_z / sig_m;
      }
      
   spacc_mu = gsl_interp_accel_alloc();    
   sp_mu = gsl_spline_alloc(gsl_interp_cspline, NI);
   
   gsl_spline_init(sp_mu, SplogM, SpMu, NI);
		
}               /* end init_spline_mu */	

/* Interface */
double mu(double M)
{
  double temp_logm = log10(M);
  double result = gsl_spline_eval(sp_mu, temp_logm, spacc_mu);
  
  return result;	
}              /* end mu */
 
double dlogsigmasq_dlogm(double M)
{
  double temp_logm = log10(M);
  double result = -2.0 / mu(M) / log(10) * gsl_spline_eval_deriv(sp_mu, temp_logm, spacc_mu);
  
  return result;
}            /* end dlogsigmasq_dlogm */
 
/*********************************     Spline functions related to CDM mass function *************************/
void init_spline_ncdm()
{
  int i;
  double temp_mu;
  double temp_M;
  double temp_dsm_dm;
  
  for(i=0; i<NI; i++)
    {
	temp_M = pow(10, SplogM[i]);
	temp_mu = mu(temp_M);
	temp_dsm_dm = dlogsigmasq_dlogm(temp_M);
     SpNcdm[i] = fabs(-0.5 * rho_m / temp_M / temp_M * fst_v(temp_mu) * temp_dsm_dm * mf_norm_fac);
    }
    
   spacc_ncdm = gsl_interp_accel_alloc();    
   sp_ncdm = gsl_spline_alloc(gsl_interp_cspline, NI);
   
   gsl_spline_init(sp_ncdm, SplogM, SpNcdm, NI);    
}             /* end init_spline_ncdm */

/* Interface */
double ncdm(double M)
{
  double temp_logm = log10(M);	
  double result = gsl_spline_eval(sp_ncdm, temp_logm, spacc_ncdm);
  
  return result;	
}	         /* end ncdm */
