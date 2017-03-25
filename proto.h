#ifndef ALLVAR_H
 #include "allvars.h"
#endif

double fst_v(double v);
double mf_norm_intg_kernel(double v, void *param);
void mf_normalization();
double bst_cdm(double v);
double bias_norm_intg_kernel(double v, void * param);
void bias_normalization();

void init_spline_bcdm();
double bcdm(double M);

void init_spline_mu();
double mu(double M);
void init_spline_ncdm();
double ncdm(double M);

double c_cdm(double M, double z);
double c_ddm(double M, double z);

void init_spline_bddm();
double bddm(double M);

void init_spline_nddm();
double nddm(double M);

double mf_ddm_t(double M);
double c_ddm_t(double M);

double decay_fraction(double a);
double time_ing_kenerl(double a, void *param);
double time_ing(double a);
double length_ing_kenerl(double a, void * param);
double length_ing(double a0, void * param);
double length_ing_min(double a0, void * param);
double max_length(double a);
double log_max_mass(double a);

double g_factor(double z);
double unnorm_growth(double z);
void init_growth();
double growth_factor(double z);

int init_global();
int init_2d_spline_ddm();
int init_same_time(double a);

double ncdm_intg_kernel(double M, void *param);
double acc_mf_cdm(double Mi);
double nddm_intg_kernel(double M, void *param);
double acc_mf_ddm(double Mf);
double mass_mapping_my_func(double Mi, void *param);
double mass_mapping_my_diff_func(double Mi, void *param);
void mass_mapping_my_func_diff_func(double M, void * param, double *f, double *df);
double find_initial_mass(double Mf);

double factor_uk(double k, double M, double z, int flag);
double factor_mass(double R, double k, double c);
double ci_intg_kernel(double x, void * param);
double ci_intg(double x1, double x2);
double si_intg_kernel(double x, void * param);
double si_intg(double x1, double x2);

double tk_ddm_from_hm(double k, double z);
void calculate_pk(double k, double z, double *temp, int flag);
double halo_halo_factor_intg_kernel(double M, void *param);
double halo_halo_factor(double k, double z, int flag);
double halo_factor_intg_kernel(double M, void *param);
double halo_factor(double k, double z, int flag);

double ps_lin(double k, double z);
double ps_ddm_from_spline(double k, double z);
double tk_ddm_from_spline(double k, double z);
double ps_lin_from_halofit(double k, double z);
double ps_nonlin_from_halofit(double k, double z);
double ps_nonlin_from_mod_halofit(double k, double z);

void set_cosmological_parameters(double x1, double x2, double x3, double x4, double x5, double x6);
void set_2d_spline_parameters(double x1, double x2, int x3, double x4, double x5, int x6);
void set_decay_parameters(double x1, double x2, double x3);

double cdm_halo_density_kernel(double M, void *param);
double cdm_halo_density();
double cdm_smooth_bias_kernel(double M, void *param);
double cdm_smooth_bias();
double ddm_halo_density_kernel(double M, void *param);
double ddm_halo_density();
double ddm_smooth_bias_kernel(double M, void *param);
double ddm_smooth_bias();
void init_smooth_field();

double prim_tk(double k);

double prim_tk_bbks(double k);

double radius_from_mass(double M);
double vir_radius(double M, double z);
double mass_from_r(double R);
void free_1d_spline();
void free_previous_2d_spline();
void mark(char *s);

double sigma_m_int_kernel(double k, void *param);
double unnorm_sigma_m_sq(double M);
void init_sigma_m();
double sigma_m_sq(double M);

double win_1(double x);
double win_2(double x);

void init_cl(double a);
void free_cl();
double a_from_x(double x);
double comoving_dis_ing_kernel(double a, void * param);
double comoving_dis_ing(double a);
double cl_ing_kernel(double x, void * param);
void cl_cal(double l, double lm, double ln, int flag, double *p_result);


double tk_wdm(double k, double mwdm, double z);
