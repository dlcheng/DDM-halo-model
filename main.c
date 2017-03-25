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

int main(int argc, char **argv)
{
/* Step 1
 * 1) Cosmological parameters.
 * 2) 2D table details
 * 3) Global initialization 
 */ 	
 
  set_cosmological_parameters(0.3, 0.7, (0.024/0.7/0.7), 0.7, 0.96, 0.8);
  set_2d_spline_parameters(0.005, 150, 30, 0, 1, 10);
  init_global();
  
/* Explore the decay parameter space and compare with WDM Ly-alpha suppression */ 

#ifdef EXPLORE_PARAM_SPACE
  double ztemp = 3.0;   // the redshift the probe
  double wdmmass = 3.;   // the reference WDM mass


  char filename[100];
  double atemp = 1.0 / (1.0 + ztemp);
  FILE * fp;
  sprintf(filename,"%s%.1f%s","./Result/Ly-alpha-",wdmmass,"keV.dat"); 
  fp = fopen(filename, "w+");
  
  double tautemp;
  double vktemp;
  double tkddm30;
  double tkddm50;
  double lyflag;
  
  double taumax = 100;
  double taumin = 5;
  
  double vkmax = 500;   /* to vk = 500 km/s */
  double vkmin = 10;
  
  int imax = 20;
  int jmax = 100;
  int i, j;
  int j_pre = 0;
  
  double tkwdm30 = tk_wdm(30, wdmmass, ztemp);
  double tkwdm50 = tk_wdm(50, wdmmass, ztemp);
  
  for(i=0; i<imax; i++) 
  {
    tautemp = log10(taumin) + i * log10(taumax/taumin) / (double) (imax -1);
    tautemp = pow(10, tautemp);
    	  
    for(j=j_pre; j<jmax; j++)	  
     {
     vktemp = log10(vkmin) + j * log10(vkmax/vkmin) / (double) (jmax - 1);
     vktemp = pow(10, vktemp);
        
     set_decay_parameters(1.0/tautemp, vktemp, 0);     
     init_same_time(atemp);
     tkddm30 = tk_ddm_from_hm(30, ztemp);
     tkddm50 = tk_ddm_from_hm(50, ztemp);
  
     if(tkddm30 > tkwdm30 && tkddm50 > tkwdm50)
       lyflag = 1.0;     // good
     else
       lyflag = 0.0;     // rule out
 
     printf("%.6e\t%.6e\t%d\n", tautemp, vktemp, (int) lyflag);

     if(j == j_pre && lyflag == 0.0) //rule out at the first place
       break;

     if(j == (jmax - 1)) // at the last vk, no matter whether it is ruled out, write it down
      {
      fprintf(fp,"%.6e\t%.6e\n", vktemp, tautemp);
      break;
      }

     if((lyflag == 0.0 && j!= j_pre) && (j != (jmax - 1)) )
       {
       fprintf(fp,"%.6e\t%.6e\n", vktemp, tautemp);
       j_pre = j; 
       break;
       }
  }
  
  free_1d_spline();
  }
  
  fclose(fp);
#endif 

/* Output the transfer function for centain decay parameters */

#ifdef OUTPUT_TK  
  double tautemp = 13.79;
  double vktemp = 50.0;
  double ztemp = 0.0;

  if(argc == 4)
   {
    tautemp = atof(argv[1]);
    vktemp  = atof(argv[2]);
    ztemp   = atof(argv[3]);
   }
  else
   {
    printf("Need 3 arguments:\n");
    printf("e.g. ddm_tk 13.79 100 0\n");
    printf("where the first is tau, the second is vk and the last is redshift.\n");
    return 0;
   }

  double atemp = 1.0 / (1.0 + ztemp);    
  double kmintemp = 0.01;
  double kmaxtemp = 100.0;
  
  double ktemp;
  double tktemp;
  
  int i;
  int imax = 30;
  
  set_decay_parameters(1.0/tautemp, vktemp, 0); 
  init_same_time(atemp);
  
  printf("#k, tk (tau=%.2f, vk=%.2f)\n", tautemp, vktemp);

  for(i=0; i<imax; i++)
    {
    ktemp = log10(kmintemp) + i * log10(kmaxtemp/kmintemp) / (double) (imax - 1);
    ktemp = pow(10.0, ktemp);
				
    tktemp = tk_ddm_from_hm(ktemp, ztemp);
    printf("%.6e	%.6e\n", ktemp, tktemp);    
    }

#endif 

  return 0;
}   /* end main */	

double tk_wdm(double k, double mwdm, double z)
{
  double alpha = 0.0476 * pow(1.0 / mwdm, 1.85) * pow((1.0 + z) / 2.0, 1.3);
  double v = 3.0;
  double l = 0.6;
  double s = 0.4;
  
  double result = pow(alpha * k, v * l);
  result += 1.0;
  
  result = pow(result, (0.0 - s) / 2.0 / v);
  
  return result;
	
}  /* end tk_wdm */
