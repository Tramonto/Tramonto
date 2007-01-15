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
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define CMS_SCFT     2
#define CMS          0
void print_to_file(FILE *fp,double val,char *var_label,int first);
double integrand_WTC_freen_bulk(int iunk,int inode_box,double **x);
double integrand_WTC_freen(int iunk,int inode_box,double **x);
double integrand_surface_charge(int iunk,int inode_box,int iwall,double **x);
double integrateOverSurface(double(*fp_integrand)(int,int,int,double **),int iunk,double **x,double *profile);
double integrand_adsorption_bulk(int iunk,int inode_box,double **x);
double integrand_adsorption(int iunk,int inode_box,double **x);
double integrand_maxwell_stress_freen(int iunk,int inode_box,double **x);
#define PI    M_PI
extern double Temp_elec;
#define NCOMP_MAX 5
extern double Rho_b[NCOMP_MAX];
extern int Ncomp;
extern int Type_coul;
double integrand_att_freen_bulk(int iunk,int inode_box,double **x);
double integrand_att_freen(int iunk,int inode_box,double **x);
extern int Type_attr;
double integrand_hs_freen_bulk(int iunk,int inode_box,double **x);
double integrand_hs_freen(int iunk,int inode_box,double **x);
extern int Type_func;
double integrand_vext_freen(int iunk,int inode_box,double **x);
double integrand_mu_freen_bulk(int iunk,int inode_box,double **x);
double integrand_mu_freen(int iunk,int inode_box,double **x);
double integrateInSpace_SumInComp(double(*fp_integrand)(int,int,double **),int **nelhit,double **x,double *profile);
void print_to_screen(double val,char *var_label);
extern int **Nel_hit;
double integrand_ideal_gas_freen_bulk(int iunk,int inode_box,double **x);
extern int **Nel_hit2;
double integrand_ideal_gas_freen(int iunk,int inode_box,double **x);
double integrateInSpace(double(*fp_integrand)(int,int,double **),int iunk,int **nelhit,double **x,double *profile);
extern double Charge_f[NCOMP_MAX];
#define NMER_MAX     100
extern int Unk2Comp[NMER_MAX];
#define NEQ_TYPE       8
extern int Phys2Unk_last[NEQ_TYPE];
#define DENSITY        0
extern int Phys2Unk_first[NEQ_TYPE];
#define WTC          3
#define NONE       -1
#define NONE      -1
#define NONE -1
extern int Type_poly;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define NO_SCREEN    2 
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Nnodes_per_proc;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
extern double *Integration_profile;
#define VERBOSE      3 
extern int Iwrite;
extern int Ndim;
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
double calc_free_energy(FILE *fp,double **x);
