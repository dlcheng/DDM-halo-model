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

void init_spline_nddm()
{
   int i;
   double temp_M;

   for(i=0; i<NI; i++)
     {
     temp_M = pow(10, SplogM[i]);
     SpNddm[i] = SpNcdm[i] * mf_ddm_t(temp_M);
     }

   spacc_nddm = gsl_interp_accel_alloc();    
   sp_nddm = gsl_spline_alloc(gsl_interp_cspline, NI);
   
   gsl_spline_init(sp_nddm, SplogM, SpNddm, NI);   
}                  /* end init_spline_nddm */

/* Interface */
double nddm(double M)
{
  double temp_logm = log10(M);	
  double result = gsl_spline_eval(sp_nddm, temp_logm, spacc_nddm);
  
  return result;		
}                  /* end nddm */	
