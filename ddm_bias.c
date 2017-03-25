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

void init_spline_bddm()
{
  int i;	
  double temp_Mf, temp_Mi;	
	
  for(i=0; i<NI; i++)
    {	
	  temp_Mf = pow(10, SplogM[i]);
	  temp_Mi = find_initial_mass(temp_Mf);
	  SpBddm[i] = bcdm(temp_Mi);
#ifdef DEBUG	
	  printf("i=%.3d, temp_Mf = %.6e, temp_Mi = %.6e, Acc_ddm = %.6e, Acc_cdm = %.6e, Bcdm = %.6e\n", 
	          i, temp_Mf, temp_Mi, acc_mf_ddm(temp_Mf), acc_mf_cdm(temp_Mi), bcdm(temp_Mi));
#endif	         
    }
  
  spacc_bddm = gsl_interp_accel_alloc();    
  sp_bddm = gsl_spline_alloc(gsl_interp_cspline, NI);
   
  gsl_spline_init(sp_bddm, SplogM, SpBddm, NI);   
}             /* end init_spline_bddm */	

/* Interface */
double bddm(double M)
{
  double temp_logm = log10(M);	
  double result = gsl_spline_eval(sp_bddm, temp_logm, spacc_bddm);
       
  return result;	
}             /* end bddm */	
