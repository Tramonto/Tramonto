/* This file was automatically generated.  Do not edit! */
double load_Chain_Geqns(int func_type_field,int Njacobian_types,int Njacobian_sums,void(*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
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
extern int *Pol_Sym;
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
extern int *Unk_to_Poly;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
double WJDC_Resid_Bulk_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
double WJDC_Resid_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
#define NBOND_MAX 4
#define NMER_MAX     200
void calc_init_polymer_G_wjdc(double **xInBox,double **xOwned);
extern double G_WJDC_b[NMER_MAX *NBOND_MAX];
#define G_CHAIN       11 
extern int *Unk_to_Seg;
extern int Nbonds;
void setup_polymer_G_wjdc(double **xOwned);
void safe_free(void **ptr);
void safe_free(void **ptr);
#define INIT_GUESS_FLAG  2
typedef struct RB_Struct RB_Struct;
double load_euler_lagrange(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,struct RB_Struct *dphi_drb,int mesh_coarsen_flag_i,int resid_only_flag);
void node_box_to_ijk_box(int node_box,int *ijk_box);
void FMT1stDeriv_switch(double **x,struct RB_Struct *dphi_drb);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_func;
extern void *LinProbMgr_manager;
extern int *B2L_node;
extern int Nnodes_box;
extern int Nunk_per_node;
#define NDIM_MAX  3
struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};
void calc_init_WJDC_field(double **xInBox,double **xOwned);
extern double Field_WJDC_b[NMER_MAX];
extern int Unk2Comp[NMER_MAX];
extern int *L2B_node;
extern int **Zero_density_TF;
#define WJDC_FIELD     8
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nmissing_densities;
#define RESTART_FEWERCOMP  4
extern int Restart;
extern int Nnodes_per_proc;
extern int Ncomp;
#define WJDC3        5 
#define WJDC2        4 
extern int Nseg_tot;
#define WJDC         3
extern int Type_poly;
void setup_polymer_field_wjdc(double **xOwned);
