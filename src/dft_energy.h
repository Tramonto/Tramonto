/* This file was automatically generated.  Do not edit! */
void safe_free(void **ptr);
void safe_free(void **ptr);
void print_freen_profile_1D(double *freen_profile_1D,char *output_file);
double integrand_CMS_freen_bulk(int iunk,int inode_box,double **x);
double integrand_CMS_freen(int iunk,int inode_box,double **x);
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
#include <string.h>
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define CMS_SCFT     1
#define CMS          0
#define UNIFORM_INTERFACE  0
extern int Type_interface;
void print_to_file(FILE *fp,double val,char *var_label,int first);
#define SWITCH_BULK_OUTPUT_ALL 6
#define SWITCH_BULK_OUTPUT 5
extern int Print_rho_switch;
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
#define NO_SCREEN    4 
double integrand_WJDCcomp_freen_bulk(int iunk,int inode_box,double **x);
double integrand_WJDCcomp_freen(int iunk,int inode_box,double **x);
double integrand_WJDC_freen_bulk(int iunk,int inode_box,double **x);
double integrand_WJDC_freen(int iunk,int inode_box,double **x);
extern int Lseg_densities;
double integrand_WTC_freen_bulk(int iunk,int inode_box,double **x);
double integrand_WTC_freen(int iunk,int inode_box,double **x);
#define WTC          2
double integrand_elec_MSAcorr_freen_bulk(int iunk,int inode_box,double **x);
double integrand_elec_MSAcorr_freen(int iunk,int inode_box,double **x);
#define DELTAC_GENERAL 2
#define DELTAC_RPM     1 
double integrand_surface_charge(int iunk,int inode_box,int iwall,double **x);
double integrateOverSurface(double(*fp_integrand)(int,int,int,double **),int iunk,double **x,double *profile);
double integrand_elec_PB_freen(int iunk,int inode_box,double **x);
#define PI    3.141592653589793238462643383279502884197169399375
extern double Temp_elec;
#define NCOMP_MAX 5
extern double Rho_b[NCOMP_MAX];
extern double Charge_f[NCOMP_MAX];
extern int Ncomp;
extern int Type_coul;
double integrand_att_freen_bulk(int iunk,int inode_box,double **x);
double integrand_att_freen(int iunk,int inode_box,double **x);
extern int Type_attr;
double integrand_hs_freen_bulk(int iunk,int inode_box,double **x);
double integrand_hs_freen(int iunk,int inode_box,double **x);
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_func;
double integrand_vext_freen(int iunk,int inode_box,double **x);
extern int Nwall;
double integrand_mu_freen_bulk(int iunk,int inode_box,double **x);
double integrand_mu_freen(int iunk,int inode_box,double **x);
double integrateInSpace_SumInComp(double(*fp_integrand)(int,int,double **),int **nelhit,double **x,double *profile);
void print_to_screen(double val,char *var_label);
extern int **Nel_hit;
double integrand_ideal_gas_freen_bulk(int iunk,int inode_box,double **x);
extern int **Nel_hit2;
double integrand_ideal_gas_freen(int iunk,int inode_box,double **x);
double integrateInSpace(double(*fp_integrand)(int,int,double **),int iunk,int **nelhit,double **x,double *profile);
#define NEQ_TYPE       12 
extern int Phys2Unk_last[NEQ_TYPE];
#define DENSITY        0
extern int Phys2Unk_first[NEQ_TYPE];
#define WJDC3        5 
#define WJDC2        4 
#define WJDC         3
extern int Type_poly;
extern int L_HSperturbation;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define SCREEN_VERBOSE     3 
#define SCREEN_ERRORS_ONLY  0 
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Nnodes_per_proc;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern double *Integration_profile;
#define FILES_BASIC        0
extern int Iwrite_files;
extern int Ndim;
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
double calc_free_energy(FILE *fp,double **x);
double WJDCgraft_freen(double **x);
extern int Grafted_Logical;

