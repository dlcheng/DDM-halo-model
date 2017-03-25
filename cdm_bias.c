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

void init_spline_bcdm()
{
  int i;	
  double temp_M;
  double temp_mu;	
	
  for(i=0; i<NI; i++)	
    {
	temp_M = pow(10, SplogM[i]);
	temp_mu = mu(temp_M);
	SpBcdm[i] = bst_cdm(temp_mu) * bias_norm_fac;
    }
	
  spacc_bcdm = gsl_interp_accel_alloc();
  sp_bcdm = gsl_spline_alloc(gsl_interp_cspline, NI);
  
  gsl_spline_init(sp_bcdm, SplogM, SpBcdm, NI);
}                      /* end init_spline_bcdm */	

/* Interface */
double bcdm(double M)
{
  double temp_logm = log10(M);	
  double result = gsl_spline_eval(sp_bcdm, temp_logm, spacc_bcdm);
 	
  return result;
}      /* end bcdm */	
