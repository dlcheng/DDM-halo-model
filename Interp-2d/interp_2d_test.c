#include <gsl/gsl_vector.h>
#include <stdio.h>
#include "Spline_2d.h"

int main()

{   
	int nx = 20;
	int ny = 10;
	double kmin = 0.005;
	double kmax = 100;
	double amin = 0.5;
	double amax = 1;
	FILE *fp;         
    int i,j;
    
    fp = fopen("data.txt", "r");
        	
    spline_2d * i2d = spline_2d_alloc(nx,ny);  // step 1   
    gsl_vector * x = gsl_vector_alloc(nx);
    gsl_vector * y = gsl_vector_alloc(ny);
    gsl_vector * z = gsl_vector_alloc(nx*ny); 

    double ix, iy;
    double a_temp, k_temp, tk;
    ix = log10(amax/amin) / ((double) nx - 1.0);
    iy = log10(kmax/kmin) / ((double) ny - 1.0);
    
    for(i=0; i<nx; i++)
      {
	  a_temp = log10(amin) + i * ix;
//	  a_temp = pow(10, a_temp);
	  gsl_vector_set(x, i, a_temp);	  
	  }	  
    
    for(i=0; i<ny; i++)
      {
	  k_temp = log10(kmin) + i * iy;	  
//	  k_temp = pow(10, k_temp);
	  gsl_vector_set(y, i, k_temp);
	  }
	  	  
    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {
          fscanf(fp,"%lf	%lf	%lf\n", &a_temp, &k_temp, &tk);
          printf("a = %.6e, k = %.6e, tk = %.6e\n", a_temp, k_temp, tk);
          gsl_vector_set(z, i*ny+j, tk);
        }
    }
/*    
    for(i=0; i<nx*ny; i++)
      {
	  printf("Tk = %.6e\n", gsl_vector_get(z, i));		  
      }
*/      		  
    spline_2d_init(i2d,x->data,y->data,z->data);    //step 2
        
    double in_x, in_y;
    char flag;

    do
    {
    printf("In_a = ");
    scanf(" %lf", &in_x);
    printf("a = %lf\n", in_x);
    printf("In_k = ");
    scanf(" %lf", &in_y);   
    printf("k = %lf\n", in_y);
    printf("Est=%lf\n", spline_2d_eval(i2d, log10(in_x),log10(in_y)));
    printf("Continue(y/n)? ");
    if(scanf(" %c", &flag)); 
    }
    while(flag == 'y');
    
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_vector_free(z);
    

    spline_2d_free(i2d);
    
    return 0;
    
}
