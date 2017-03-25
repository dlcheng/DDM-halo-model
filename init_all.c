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

int init_global()
{
 rho_crit = 3e10*MPCTOM/(8.0*PI*G0)/SUNTOKG;
 rho_m    = rho_crit*Omega_m;
 set_cosmological_parameters_to_default();                /* set parameters for halofit */	
	
 init_growth();	
 init_sigma_m();   /* Calculate the sigma_8 normalization */
 
 int nx = a_num;
 int ny = k_num;
 double ix, iy, temp;
 int i;
 
 ax = gsl_vector_alloc(nx);
 ky = gsl_vector_alloc(ny);
 tkz = gsl_vector_alloc(nx*ny);
 
 ix = log10(amax_2d/amin_2d) / ((double) nx - 1.0);
 iy = log10(kmax_2d/kmin_2d) / ((double) ny - 1.0);
 
 for(i=0; i<nx; i++)
   {
	temp = log10(amin_2d) + i * ix;   
//	temp = pow(10, temp);
	gsl_vector_set(ax, i, temp);
   }	      
   
 for(i=0; i<ny; i++)
   {
    temp = log10(kmin_2d) + i * iy;
//  temp = pow(10, temp);  
	gsl_vector_set(ky, i, temp);   
   }
   	     
 pre_2d_spline_flag = 0;  	     	   
 
 return 1; 	  	
}  /* end init_step_1 */	

int init_2d_spline_ddm()
{
  free_previous_2d_spline();
  
  int nx = a_num;
  int ny = k_num;
  double a_temp;
  double z_temp;
  double k_temp;
  int i, j;	
  double tk;
//  FILE *fp;  
//  fp = fopen("data.txt", "w+");
  
  sp_ddm_tk_2d = spline_2d_alloc(nx,ny);
	  	  
  for(i=0; i<nx; i++)
    {
#ifdef VERBOSE	
	printf("%d & %d\n", i+1, nx);
#endif	
	a_temp = gsl_vector_get(ax, i);
	a_temp = pow(10, a_temp);
	z_temp = 1.0 / a_temp - 1.0;
	init_same_time(a_temp);           /* initialize things at the same time */	
	
	for(j=0; j<ny; j++)
	  {
	  k_temp = gsl_vector_get(ky, j);
	  k_temp = pow(10, k_temp);
	  tk = tk_ddm_from_hm(k_temp, z_temp);
#ifdef VERBOSE	  
	  printf("%.6e	%.6e	%.6e\n", a_temp, k_temp, tk);
#endif	  
	  
      if(tk <= 1e-3)
         tk = 0;    
      
//      fprintf(fp, "%.6e	%.6e	%.6e\n", a_temp, k_temp, tk);       
	  gsl_vector_set(tkz, i*ny+j, tk);  
	  }	
	free_1d_spline();	
	} 		

//  fclose(fp); 
  spline_2d_init(sp_ddm_tk_2d, ax->data, ky->data, tkz->data);
  
  pre_2d_spline_flag = 1;
  
  return 1;
}     /* end init_step_2 */	

int init_same_time(double a)
{
 double z = 1.0 / a - 1.0;
 
 delta_z = delta_c / growth_factor(z);
 decay_f = decay_fraction(a);
 log_ddm_norm_mass = log_max_mass(a);
// printf("a=%.6e, decay_f = %.6e, log_ddm_norm_mass = %.6e\n", a, decay_f, log_ddm_norm_mass);
 
 mf_norm_fac = 1.0;
 bias_norm_fac = 1.0;
 
 init_spline_mu();
 init_spline_ncdm();
 init_spline_bcdm();
 init_spline_nddm();
 init_spline_bddm();
 init_smooth_field();	
 
 return 1;
}      /* end init_same_time */
