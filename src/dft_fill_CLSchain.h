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
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define NCOMP_MAX 5
extern double Rho_b_RTF[NCOMP_MAX];
extern double Rho_b_LBB[NCOMP_MAX];
#define NMER_MAX     200
extern double Rho_seg_RTF[NMER_MAX];
extern double Rho_seg_LBB[NMER_MAX];
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
extern int Lhard_surf;
extern int **Nbond;
extern int Nlists_HW;
double load_Chain_Geqns_SCF(int func_type_field,int Njacobian_types,int Njacobian_sums,void(*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
double load_polymer_recursion(int sten_type,int func_type_field,int Njacobian_types,int Njacobian_sums,void(*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),int iunk,int loc_inode,int inode_box,int unk_B,int itype_mer,int izone,int *ijk_box,double **x,int resid_only_flag);
extern double Bond_ff[NCOMP_MAX][NCOMP_MAX];
extern int *L2G_node;
double fill_zero_value(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
extern double VEXT_MAX;
extern double **Vext;
extern int *Unk_to_Bond;
extern int *Unk_to_Seg;
extern int *Unk_to_Poly;
#define G_CHAIN       11 
double load_Chain_Geqns(int func_type_field,int Njacobian_types,int Njacobian_sums,void(*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **),double(*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **),int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,int resid_only_flag);
extern double Field_WJDC_b[NMER_MAX];
extern double Rho_seg_b[NMER_MAX];
#define NBOND_MAX 4
extern double G_WJDC_b[NMER_MAX *NBOND_MAX];
double resid_and_Jac_ChainDensity_WJDC2(int func_type,double **x,int iunk,int unk_B,int loc_inode,int inode_box,int resid_only_flag,double(*fp_prefactor)(int));
double dy_dxi3_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
double dy_dxi2_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
extern int **Zero_density_TF;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
extern int Ncomp;
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
#define NDIM_MAX  3
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
double y_cav(double sigma_1,double sigma_2,double xi_2,double xi_3);
#define CAVWTC         4
#define DELTA_FN_BOND         6
void grafted_jac(int sten_type,int itype_mer,int iunk,int loc_inode,int *ijk_box,int izone,int unk_G,int inode_boxl,double prefac,double ysqrt,double **x);
extern int *Mesh_coarsen_flag;
extern int L1D_bc;
extern int Nwall_type;
extern int Mesh_coarsening;
void node_box_to_ijk_box(int node_box,int *ijk_box);
int node_to_node_box(int inode);
int position_to_node(double *NodePos);
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern double *Poly_graft_dist;
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int WallType[NWALL_MAX];
extern int Graft_wall[NCOMP_MAX];
void node_to_position(int inode,double *NodePos);
extern int **Nodes_2_boundary_wall;
extern int ***Bonds;
extern int Grafted[NCOMP_MAX];
extern int *L2B_node;
extern int Nnodes_per_proc;
extern int ***Poly_to_Unk;
extern int Geqn_start[NCOMP_MAX];
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern int *Pol_Sym;
extern int **Poly_to_Unk_SegAll;
extern double Gsum[NCOMP_MAX];
extern int *Nbonds_SegAll;
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int Lseg_densities;
extern int Nmer[NCOMP_MAX];
#define NBLOCK_MAX   20 
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
#define WJDC3        5 
#define CMS_SCFT     1
#define CMS          0
extern int Unk2Comp[NMER_MAX];
#define WJDC2        4 
#define WJDC         3
extern int Type_poly;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern void *LinProbMgr_manager;
#define CALC_RESID_ONLY  3
#define INIT_GUESS_FLAG  2
double fill_constant_density_chain(int iunk,int icomp,int iseg,double fac_FIELD,int loc_inode,int inode_box,double **x,int resid_only_flag);
#define DENSITY        0
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern double Esize_x[NDIM_MAX];
extern int Grad_dim;
extern double Size_x[NDIM_MAX];
extern int *B2G_node;
#define PHASE_INTERFACE 2
extern int Type_interface;
extern int Lconstrain_interface;
double resid_and_Jac_ChainDensity(int func_type,double **x,int iunk,int unk_B,int loc_inode,int inode_box,int resid_only_flag,double(*fp_prefactor)(int));
