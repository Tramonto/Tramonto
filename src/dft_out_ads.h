/* This file was automatically generated.  Do not edit! */
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
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define NCOMP_MAX 5
extern double Charge_f[NCOMP_MAX];
extern double Rho_b[NCOMP_MAX];
#define NMER_MAX     100
extern double Rho_seg_b[NMER_MAX];
void print_to_file(FILE *fp,double val,char *var_label,int first);
void print_to_screen(double val,char *var_label);
double integrand_fluid_charge(int iunk,int inode_box,double **x);
double integrateInSpace_SumInComp(double(*fp_integrand)(int,int,double **),int **nelhit,double **x,double *profile);
void calc_fluid_charge(FILE *fp,double **x);
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int **Nel_hit;
double integrand_adsorption_bulk(int iunk,int inode_box,double **x);
void print_to_file_comp(FILE *fp,int icomp,double val,char *var_label,int first);
void print_to_screen_comp(int icomp,double val,char *var_label);
extern int **Nel_hit2;
double integrand_adsorption(int iunk,int inode_box,double **x);
double integrateInSpace(double(*fp_integrand)(int,int,double **),int iunk,int **nelhit,double **x,double *profile);
extern int Unk2Comp[NMER_MAX];
#define NEQ_TYPE       8
extern int Phys2Unk_last[NEQ_TYPE];
#define DENSITY        0
extern int Phys2Unk_first[NEQ_TYPE];
#define NO_SCREEN    2 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Ncomp;
extern int Nseg_tot;
extern int Lseg_densities;
extern double *Integration_profile;
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void calc_adsorption(FILE *fp,double **x);
