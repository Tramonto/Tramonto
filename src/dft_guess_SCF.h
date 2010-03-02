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
#define INIT_GUESS_FLAG  2
double load_Chain_Geqns_SCF(int func_type_field,int Njacobian_types,int Njacobian_sums,void(*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#define G_CHAIN       11 
void node_box_to_ijk_box(int node_box,int *ijk_box);
extern void *LinProbMgr_manager;
#define NCOMP_MAX 5
extern int Geqn_start[NCOMP_MAX];
extern int ***Poly_to_Unk;
extern int **Nbond;
extern int ***Bonds;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int *Unk_to_Bond;
extern int *Unk_to_Seg;
extern int *Unk_to_Poly;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Nbonds;
double CMS_Resid_Bulk_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
double CMS_Resid_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
#define NBOND_MAX 4
#define NMER_MAX     200
void calc_init_polymer_G_SCF(double **xInBox,double **xOwned);
void calc_init_lambda(double **xInBox,double **xOwned);
extern double Rho_t;
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
#define NEQ_TYPE       13 
extern int Restart_field[NEQ_TYPE];
#define SCF_CONSTR	   9
extern int Phys2Nunk[NEQ_TYPE];
extern double **Vext;
extern int **Zero_density_TF;
void calc_init_SCFfield(double **xInBox,double **xOwned);
extern double Rho_b[NCOMP_MAX];
extern double VEXT_MAX;
#define SCF_FIELD	  10
#define DENSITY        0
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nmissing_densities;
#define RESTART_FEWERCOMP  4
extern int Restart;
extern int Ncomp;
extern int *L2B_node;
extern int Nnodes_per_proc;
void setup_polymer_SCF_field(double **xInBox,double **xOwned,int guess_type);
