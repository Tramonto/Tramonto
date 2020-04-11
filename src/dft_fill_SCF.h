/* This file was automatically generated.  Do not edit! */
double load_Chain_Geqns_SCF(int func_type_field,int Njacobian_types,int Njacobian_sums,void(*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
double CMS_Resid_Bulk_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
double CMS_Resid_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
void CMS_Jacobian_GCHAIN_derivFIELD(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
void CMS_Jacobian_GCHAIN_derivG(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
double load_SCF_Geqns(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
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
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
extern void *LinProbMgr_manager;
#define CALC_RESID_ONLY  3
#define INIT_GUESS_FLAG  2
#define SCF_CONSTR	   9
#define NEQ_TYPE       12 
extern int Phys2Nunk[NEQ_TYPE];
extern double Rho_t;
#define NCOMP_MAX 6
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern int Ncomp;
double load_SCF_field(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
extern double Rho_b[NCOMP_MAX];
#define NBLOCK_MAX   20 
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
double prefactor_rho_scft(int itype_mer,int inode_box,double **x);
#define G_CHAIN       11 
double resid_and_Jac_ChainDensity(int func_type,double **x,int iunk,int unk_B,int loc_inode,int inode_box,int resid_only_flag,double(*fp_prefactor)(int,int,double **));
#define SCF_FIELD	  10
double fill_zero_value(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
extern double VEXT_MAX;
extern double **Vext;
extern int **Zero_density_TF;
#define DENSITY        0
extern int Phys2Unk_first[NEQ_TYPE];
double load_SCF_density(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
