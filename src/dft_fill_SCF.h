/* This file was automatically generated.  Do not edit! */
double load_lambda_field(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
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
extern void *LinProbMgr_manager;
#define CALC_RESID_ONLY  3
#define INIT_GUESS_FLAG  2
extern double Rho_t;
#define NCOMP_MAX 5
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern int Ncomp;
#define SCF_CONSTR	   9
double load_SCF_field(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
extern double Rho_b[NCOMP_MAX];
#define NBLOCK_MAX   20 
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
double prefactor_rho_scft(int itype_mer);
#define G_CHAIN       11 
double resid_and_Jac_ChainDensity(int func_type,double **x,int iunk,int unk_B,int loc_inode,int inode_box,int resid_only_flag,double(*fp_prefactor)(int));
#define SCF_FIELD	  10
double fill_zero_value(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
extern double VEXT_MAX;
extern double **Vext;
extern int **Zero_density_TF;
#define DENSITY        0
#define NEQ_TYPE       13 
extern int Phys2Unk_first[NEQ_TYPE];
double load_SCF_density(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
extern int Phys2Nunk[NEQ_TYPE]; 