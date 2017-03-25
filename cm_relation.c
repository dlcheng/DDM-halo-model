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


/* The Cold dark matter c-M relation */
double c_cdm(double M, double z)
{

  double alpha = -0.061094;
  double belta = 0.370894;
  double B = 4.710019;
  
  double result = B / pow(1.0 + z, belta) * pow(M/1e14, alpha);

  return result;	
 
}     /* end c_cdm */

/* The decay dark matter c-M relation */
double c_ddm(double M, double z)
{	
  return c_cdm(M, z) * c_ddm_t(M);
}     /* end c_ddm */
