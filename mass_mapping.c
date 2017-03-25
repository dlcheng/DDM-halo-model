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

double ncdm_intg_kernel(double M, void *param)
{
   return ncdm(M);		
}       /* end ncdm_intg_kernel */	

double acc_mf_cdm(double Mi)
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
   
   double result, error;
   		
   gsl_function F;
   F.function = &ncdm_intg_kernel;
     
   gsl_integration_qag(&F, Mi, Mmax, 0, 1e-4, 10000, GSL_INTEG_GAUSS21, w, &result, &error);

   gsl_integration_workspace_free(w);
   return result;		
}       /* end acc_mf_cdm */	

double nddm_intg_kernel(double M, void *param)
{
   return nddm(M);		
}       /* end nddm_intg_kernel */	

double acc_mf_ddm(double Mf)
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);
   
   double result, error;
   		
   gsl_function F;
   F.function = &nddm_intg_kernel;
        
   gsl_integration_qag(&F, Mf, Mmax, 0, 1e-4, 10000, GSL_INTEG_GAUSS21, w, &result, &error);
   
   gsl_integration_workspace_free(w);
   return result;		
}       /* end acc_mf_ddm */	

/*************************************  local mass mapping ***************************/
double mass_mapping_my_func(double Mi, void *param)
{
  double right_hand = *(double *)param;	

  return acc_mf_cdm(Mi) - right_hand;	
	
}      /* end mass_mapping_my_func */	

double mass_mapping_my_diff_func(double Mi, void *param)
{  	
   return -1.0 * ncdm(Mi);		
}    /* end mass_mapping_my_diff_func */	

void mass_mapping_my_func_diff_func(double M, void * param, double *f, double *df)
{
  *f = mass_mapping_my_func(M, param);	
  *df = mass_mapping_my_diff_func(M, param);	
}    /* end mass_mapping_my_func_diff_func */	

double find_initial_mass(double Mf)
{
   double right_hand = acc_mf_ddm(Mf);	

/* If the mf difference is so small */
   if(fabs(mf_ddm_t(Mf) - 1.0) < 1e-4 || fabs(right_hand) < 1e-30)  
     return Mf;
 	   
   int status;
   int iter = 0, max_iter = 100;
   const gsl_root_fdfsolver_type *T;	
   gsl_root_fdfsolver *s;
   double Mt;   
   double Mg = Mf;
   
   gsl_function_fdf FDF;
   FDF.f = &mass_mapping_my_func;
   FDF.df = &mass_mapping_my_diff_func;
   FDF.fdf = &mass_mapping_my_func_diff_func;
   FDF.params = &right_hand;
   
   T = gsl_root_fdfsolver_newton;
   s = gsl_root_fdfsolver_alloc (T);
   gsl_root_fdfsolver_set(s, &FDF, Mg);
   
   do 
     {
	  iter++;
	  status =gsl_root_fdfsolver_iterate(s);	 
	  Mt = Mg;	 
	  Mg = gsl_root_fdfsolver_root(s);	 
	  status = gsl_root_test_delta(Mg, Mt, 0 , 1e-4);
	   }
   while(status == GSL_CONTINUE && iter < max_iter);
	
   gsl_root_fdfsolver_free(s);
        
   return Mg;
}          /* find_intitial_mass */	

