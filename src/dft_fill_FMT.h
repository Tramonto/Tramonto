/* This file was automatically generated.  Do not edit! */
typedef struct RB_Struct RB_Struct;
void calc_FMT_derivatives(void(*fp_FMTderiv)(double *,double,double,double *,double *),int inode_box,double **x,struct RB_Struct *dphi_drb);
double prefactor_rho_bar_v(int iunk,int jcomp,int *offset);
double load_rho_bar_v(double **x,int iunk,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag);
double constant_boundary(int iunk,int jnode_box);
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
extern int **Zero_density_TF;
#define DENSITY        0
#define NMER_MAX     100
extern int Unk2Comp[NMER_MAX];
extern int Lseg_densities;
double jac_rho_bar(int junk,int jnode_box,double **x);
double resid_rho_bar(int junk,int jnode_box,double **x);
double prefactor_rho_bar_s(int iunk,int jcomp,int *offset);
double resid_and_Jac_sten_fill_sum_Ncomp(int sten_type,double **x,int iunk,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag,int jzone_flag,double(*fp_prefactor)(int,int,int *),double(*fp_resid)(int,int,double **),double(*fp_jacobian)(int,int,double **));
extern int Nwall;
extern int Nlists_HW;
extern int Lhard_surf;
double load_rho_bar_s(int sten_type,double **x,int iunk,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag);
#define NDIM_MAX  3
struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};
struct RB_Struct FMT2ndDerivTheta_switch(double *n);
struct RB_Struct FMT2ndDerivDelta_switch(double *n,int *offset,double *sign,int icomp);
void solutionVec_to_nOrdering(double *rhoBar_SVOrdering,double *n);
extern int Nrho_bar_s;
#define HSRHOBAR       2
#define NEQ_TYPE       11 
extern int Phys2Unk_first[NEQ_TYPE];
extern void *LinProbMgr_manager;
#define INIT_GUESS_FLAG  2
extern double Dphi_Drhobar_RTF[10];
extern double Dphi_Drhobar_LBB[10];
extern double Dphi_Drhobar_b[10];
#define THETA_FN_R            1
#define NCOMP_MAX 5
extern double Inv_rad[NCOMP_MAX];
extern double Esize_x[NDIM_MAX];
extern double Inv_4pir[NCOMP_MAX];
extern double Inv_4pirsq[NCOMP_MAX];
extern double Fac_overlap_hs[NCOMP_MAX];
#define DELTA_FN_R            0
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
int find_jzone(int izone,int inode_box);
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Ndim;
struct Stencil_Struct {
  int        Length;      /* Number of nodes that interact with current 
                             node through this stencil                    */
  int      **Offset;      /* 2D array to be set to size [Length][Ndim] that
                             gives relative position of interacting node  */
  double    *Weight;      /* 1D array of size [Length] that gives weight
                             of interacting node to total stencil         */
  double   **HW_Weight;   /* 2D array of size [Length][Nnodes_per_el] that
                             holds the weights based on what element they
                             are being contributed from. Only used for Hard
                             Walls when stencil point is a boundary node  */
};
double load_nonlocal_hs_rosen_rb(int sten_type,int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,struct RB_Struct *dphi_drb,int resid_only_flag);
