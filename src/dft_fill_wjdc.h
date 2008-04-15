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
extern int **Poly_to_Unk_SegAll;
double d2y_dxi3_dxi2(double sigma_1,double sigma_2,double xi_2,double xi_3);
double d2y_dxi3_sq(double sigma_1,double sigma_2,double xi_2,double xi_3);
double d2y_dxi2_sq(double sigma_1,double sigma_2,double xi_2,double xi_3);
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define CALC_RESID_ONLY  3
#define INIT_GUESS_FLAG  2
#define NCOMP_MAX 5
extern double Fac_overlap[NCOMP_MAX][NCOMP_MAX];
#define PI    M_PI
extern int **Bonds_SegAll;
extern int *Nbonds_SegAll;
double calc_dens_seg(int iseg,int inode_box,double **x);
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
extern int **Nodes_2_boundary_wall;
extern int Lhard_surf;
extern int Nlists_HW;
extern int *Pol_Sym_Seg;
extern int Nseg_tot;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
#define THETA_FN_SIG          5
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
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
double load_polyWJDC_cavityEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
double y_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
double constant_boundary(int iunk,int jnode_box);
double dy_dxi3_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
double dy_dxi2_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
#define CAVWTC         4
extern void *LinProbMgr_manager;
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
double yterm_wjdc(int icomp,int jcomp,int jnode_box,double **x);
extern int **Nbond;
double load_Chain_Geqns(int func_type_field,int Njacobian_types,int Njacobian_sums,void(*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
double WJDC_Resid_Bulk_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
double WJDC_Resid_GCHAIN(int iunk,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
void WJDC_Jacobian_GCHAIN_derivCAVITY(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
void WJDC_Jacobian_GCHAIN_derivFIELD(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
void WJDC_Jacobian_GCHAIN_derivG(int iunk,int loc_inode,int pol_num,int jseg,int unk_B,int inode_box,int jnode_box,int nunk,int *unk,double weight,double **x);
double load_WJDC_Geqns(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
#define NMER_MAX     100
extern double Betamu_chain[NMER_MAX];
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int Nmer[NCOMP_MAX];
extern int Npol_comp;
double prefactor_rho_wjdc(int iseg);
#define G_CHAIN        9 
double resid_and_Jac_ChainDensity(int func_type,double **x,int iunk,int unk_B,int loc_inode,int inode_box,int resid_only_flag,double(*fp_prefactor)(int));
#define WJDC_FIELD     8
double fill_zero_value(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
extern double VEXT_MAX;
extern double **Vext;
extern int **Zero_density_TF;
extern int Unk2Comp[NMER_MAX];
#define DENSITY        0
#define NEQ_TYPE       11 
extern int Phys2Unk_first[NEQ_TYPE];
double load_WJDC_density(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
typedef struct RB_Struct RB_Struct;
double load_euler_lagrange(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,struct RB_Struct *dphi_drb,int mesh_coarsen_flag_i,int resid_only_flag);
#define NDIM_MAX  3
struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};
double load_WJDC_field(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,struct RB_Struct *dphi_drb,int mesh_coarsen_flag_i,int resid_only_flag);
