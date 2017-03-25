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

double win_1(double x)
{
  double result = 1;
  
  if(x < 1e-8)
    result = 1.0 - 0.1 * x * x;
  else  
	result = 3.0 / pow(x, 3) * (sin(x) - x*cos(x));
		
   return result;		
}           /* end win_1 */

	
double win_2(double x)
{
  double result = 0;	
	
  if(x < 1e-8)
	result = -0.2 * x * x;
  else
    result = 3.0 / pow(x, 3) * ((x * x - 3.0) * sin(x) + 3.0 * x * cos(x));	
	
  return result;	
}           /* end win_2 */
