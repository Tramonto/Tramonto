/* This file was automatically generated.  Do not edit! */
double load_CMS_Geqns(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
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
#define CMS_G          2 
double load_CMS_field(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#define CMS_FIELD      1
double load_bond_wtc(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#define BONDWTC       7
double load_cavity_wtc(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#define CAVWTC     6
double load_nonlinear_transport_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x);
double load_linear_transport_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x);
extern int Linear_transport;
#define DIFFUSION      5
double load_poisson_control(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x);
#define POISSON        3
double load_rho_bar_v(double **x,int iunk,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag);
#define DELTA_FN       0
extern int Nrho_bar_s;
#define THETA_FN       1
double load_rho_bar_s(int sten_type,double **x,int iunk,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag);
#define NEQ_TYPE       8
extern int Phys2Unk_first[NEQ_TYPE];
#define HSRHOBAR       4
double load_CMS_density(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
typedef struct RB_Struct RB_Struct;
double load_euler_lagrange(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,struct RB_Struct *dphi_drb,int mesh_coarsen_flag_i,int resid_only_flag);
#define WTC          3
extern int Type_poly;
#define DENSITY        0
#define NCOMP_MAX 5
#define NMER_MAX     40
extern int Unk2Phys[3 *NCOMP_MAX+NMER_MAX+NMER_MAX *NMER_MAX+13];
void safe_free(void **ptr);
void safe_free(void **ptr);
double load_coarse_node_Ndim(int loc_inode,int inode_box,int iunk,double **x);
#define FLAG_PBELEC -777
#define FLAG_BULK   -888
void load_coarse_node_1dim(int loc_inode,int inode_box,int *ijk_box,int iunk,double **x);
#define FLAG_1DBC   -999
extern int *Mesh_coarsen_flag;
extern int L1D_bc;
extern int Nwall_type;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Mesh_coarsening;
void node_box_to_ijk_box(int node_box,int *ijk_box);
extern int *L2B_node;
extern int Nnodes_per_proc;
#define NODAL_FLAG -999
void print_profile_box(double **x,char *outfile);
#define VERBOSE      3 
void FMT3_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3,double *dphi_drb_loc);
#define FMT3       2
void FMT2_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3,double *dphi_drb_loc);
#define FMT2       1
void FMT1_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3,double *dphi_drb_loc);
void calc_FMT_derivatives(void(*fp_FMTderiv)(double *,double,double,double *,double *),int inode_box,double **x,struct RB_Struct *dphi_drb);
#define FMT1       0
extern int Nnodes_box;
#define NONE       -1
#define NONE      -1
#define NONE -1
extern int Type_func;
extern int Nunk_per_node;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
#define NO_SCREEN    2 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void fill_resid_and_matrix(double **x,int iter,int resid_only_flag,int unk_flag);
#define NDIM_MAX  3
struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};
void load_standard_node(int,int,int *,int,double **,struct RB_Struct *,double *,int,int);
void load_standard_node(int loc_inode,int inode_box,int *ijk_box,int iunk,double **x,struct RB_Struct *dphi_drb,double *resid_unk,int mesh_coarsen_flag_i,int resid_only_flag);
void print_residuals(int,int,double *);
void print_residuals(int loc_inode,int iunk,double *resid_unk);
