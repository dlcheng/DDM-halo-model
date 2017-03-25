#ifndef ALLVAR_H
#define ALLVAR_H

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

#include "spline_2d.h"
#include "smith2.h"
#include "define.h"

/* Cosmological parameters */
extern double Omega_m;
extern double Omega_v;
extern double Omega_b;
extern double H0;
extern double ns;         
extern double sigma_8;
extern double delta_c;

/* Decay parameters */
extern double dp_t_rev;        /* in unit of Gyr^-1 */
extern double dp_v;            /* in unit of km/s */
extern double initial_dm_frac; /* the initial fraction of the daughter particles */

/* 2D spline parameters */
extern double kmin_2d;
extern double kmax_2d;
extern double amin_2d;
extern double amax_2d;
extern int k_num;    
extern int a_num;

/* Decay related constants */
extern double decay_f;                   /* fraction of decayed dark matter */
extern double log_ddm_norm_mass;

/* Halo model  */
extern double Aps;                       /* normalization factor for the PS */
extern double Agr;                       /* normalization factor for the growth function */

extern double halo_mass_frac_cdm;        /* The mass fraction of halos with the mass range of CDM */
extern double halo_mass_frac_ddm;        /* The mass fraction of halos with the mass range of DDM */
extern double rho_halo_cdm;              /* halo_mass_frac_cdm * rho_m */
extern double rho_halo_ddm;              /* halo_mass_frac_ddm * rho_m */
extern double bs_cdm;                    /* smooth field bias(small enough halos) for CDM */
extern double bs_ddm;                    /* smooth field bias(small enough halos together with the real smooth field)
                                              for DDM
                                           */
/* Density */
extern double rho_crit;
extern double rho_m;

/* others */
extern double delta_z;
extern double mf_norm_fac;
extern double bias_norm_fac;

/* 1D Spline related */
extern double SplogM[NI];
extern double SpMu[NI];
extern double SpNcdm[NI];
extern double SpNddm[NI];
extern double SpBcdm[NI];
extern double SpBddm[NI];

extern gsl_interp_accel * spacc_mu;
extern gsl_interp_accel * spacc_ncdm;
extern gsl_interp_accel * spacc_nddm;
extern gsl_interp_accel * spacc_bcdm;
extern gsl_interp_accel * spacc_bddm;

extern gsl_spline *sp_mu;
extern gsl_spline *sp_ncdm;
extern gsl_spline *sp_nddm;
extern gsl_spline *sp_bcdm;
extern gsl_spline *sp_bddm;

/* 2D spline */
extern gsl_vector * ax;               /* dimension of a */
extern gsl_vector * ky;               /* dimension of k */
extern gsl_vector * tkz;              /* dimension of z */

extern spline_2d *sp_ddm_tk_2d;
extern int pre_2d_spline_flag;         

/* Cl calculation */
extern double * a_cl;
extern double * x_cl;
extern gsl_spline *sp_ax_cl;           /* given x retrun a */
extern gsl_interp_accel * spacc_ax_cl;
extern double xs;                      /* maxium comoving distance */

#endif
