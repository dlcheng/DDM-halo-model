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

/* Cosmological parameters */
 double Omega_m;
 double Omega_v;
 double Omega_b;
 double H0;
 double ns;         
 double sigma_8;
 double delta_c;

/* Decay parameters */
 double dp_t_rev;
 double dp_v;
 double initial_dm_frac;

/* 2D spline parameters */
 double kmin_2d;
 double kmax_2d;
 double amin_2d;
 double amax_2d;
 int k_num;    
 int a_num;

/* Decay related constants */
 double decay_f;                   /* fraction of decayed dark matter */
 double log_ddm_norm_mass;

/* Halo model  */
 double Aps;                       /* normalization factor for the PS */
 double Agr;                       /* normalization factor for the growth function */

 double halo_mass_frac_cdm;        /* The mass fraction of halos with the mass range of CDM */
 double halo_mass_frac_ddm;        /* The mass fraction of halos with the mass range of DDM */
 double rho_halo_cdm;              /* halo_mass_frac_cdm * rho_m */
 double rho_halo_ddm;              /* halo_mass_frac_ddm * rho_m */
 double bs_cdm;                    /* smooth field bias(small enough halos) for CDM */
 double bs_ddm;                    /* smooth field bias(small enough halos together with the real smooth field) for DDM*/
 
/* Density */
 double rho_crit;
 double rho_m;

/* others */
 double delta_z;
 double mf_norm_fac;
 double bias_norm_fac;

/* 1D Spline related */
 double SplogM[NI];
 double SpMu[NI];
 double SpNcdm[NI];
 double SpNddm[NI];
 double SpBcdm[NI];
 double SpBddm[NI];

 gsl_interp_accel * spacc_mu;
 gsl_interp_accel * spacc_ncdm;
 gsl_interp_accel * spacc_nddm;
 gsl_interp_accel * spacc_bcdm;
 gsl_interp_accel * spacc_bddm;

 gsl_spline *sp_mu;
 gsl_spline *sp_ncdm;
 gsl_spline *sp_nddm;
 gsl_spline *sp_bcdm;
 gsl_spline *sp_bddm;

/* 2D spline */
 gsl_vector * ax;               /* dimension of a */
 gsl_vector * ky;               /* dimension of k */
 gsl_vector * tkz;

 spline_2d *sp_ddm_tk_2d;
 int pre_2d_spline_flag;     

/* Cl calculation */
 double * a_cl;
 double * x_cl;
 gsl_spline *sp_ax_cl;           /* given x retrun a */
 gsl_interp_accel * spacc_ax_cl;
 double xs;
